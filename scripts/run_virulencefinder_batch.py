#!/usr/bin/env python3
"""Manifest-bound, resumable VirulenceFinder 3.2.1 batch runner.

The authoritative input is the selected Longcycler assembly manifest.  One
VirulenceFinder process is launched per assembly because VirulenceFinder 3.2.1
silently uses only the first value passed to ``-ifa``.

User-changeable settings live in config/virulencefinder_sensitivity.toml.  The
effective analytical settings are included in every cache context, so changing
a threshold, database, overlap, tool build, or input FASTA invalidates reuse.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import datetime as dt
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import tomllib
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


RUNNER_SCHEMA_VERSION = 1
EXPECTED_PROFILE_NAMES = (
    "web_default_id90_cov60",
    "matched_id80_cov80",
)
FASTA_ALLOWED = frozenset("ACGTURYSWKMBDHVN.-*")
METADATA_COLUMNS = (
    "Participant_id",
    "tp_lab",
    "episode_key",
    "Assembly_ID",
    "fasta_sha256",
    "UTI_Status",
    "Event_type",
    "full_path",
)


class ContractError(RuntimeError):
    """Raised when a scientific, provenance, or output contract is violated."""


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_text(value: str) -> str:
    return sha256_bytes(value.encode("utf-8"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
        os.replace(tmp_name, path)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def atomic_write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fields})
        os.replace(tmp_name, path)
    finally:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "t", "1", "yes", "y"}


def resolve_path(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = root / path
    # Keep virtual-environment executable symlinks intact: resolving
    # venv/bin/python to the Homebrew target bypasses the venv site-packages.
    return Path(os.path.abspath(path))


def run_capture(command: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def load_config(config_path: Path, root: Path) -> dict[str, Any]:
    with config_path.open("rb") as handle:
        config = tomllib.load(handle)
    if int(config.get("schema_version", -1)) != RUNNER_SCHEMA_VERSION:
        raise ContractError(
            f"Config schema {config.get('schema_version')} does not match runner schema "
            f"{RUNNER_SCHEMA_VERSION}."
        )
    for key, value in list(config["paths"].items()):
        config["paths"][key] = str(resolve_path(root, str(value)))
    profiles = config.get("profiles", {})
    if tuple(profiles.keys()) != EXPECTED_PROFILE_NAMES:
        raise ContractError(
            "Profiles must be declared in the approved order: " + ", ".join(EXPECTED_PROFILE_NAMES)
        )
    for name, profile in profiles.items():
        identity = float(profile["minimum_identity"])
        coverage = float(profile["minimum_coverage"])
        if not (0 <= identity <= 1 and 0 <= coverage <= 1):
            raise ContractError(f"Invalid identity/coverage fraction for {name}.")
    return config


def validate_fasta_and_hash(path: Path) -> dict[str, Any]:
    digest = hashlib.sha256()
    contig_ids: set[str] = set()
    n_records = 0
    total_bases = 0
    current_id: str | None = None
    current_bases = 0
    with path.open("rb") as handle:
        for line_number, raw in enumerate(handle, start=1):
            digest.update(raw)
            stripped = raw.strip()
            if not stripped:
                continue
            try:
                line = stripped.decode("ascii")
            except UnicodeDecodeError as exc:
                raise ContractError(f"Non-ASCII FASTA content in {path}:{line_number}") from exc
            if line.startswith(">"):
                if current_id is not None and current_bases == 0:
                    raise ContractError(f"FASTA record {current_id!r} has no sequence in {path}.")
                header = line[1:].strip()
                if not header:
                    raise ContractError(f"Empty FASTA header in {path}:{line_number}")
                current_id = header.split()[0]
                if current_id in contig_ids:
                    raise ContractError(f"Duplicate FASTA contig identifier {current_id!r} in {path}.")
                contig_ids.add(current_id)
                n_records += 1
                current_bases = 0
                continue
            if current_id is None:
                raise ContractError(f"Sequence data precede the first FASTA header in {path}:{line_number}")
            sequence = re.sub(r"\s+", "", line).upper()
            invalid = set(sequence) - FASTA_ALLOWED
            if invalid:
                raise ContractError(
                    f"Invalid nucleotide character(s) {sorted(invalid)!r} in {path}:{line_number}"
                )
            current_bases += len(sequence)
            total_bases += len(sequence)
    if n_records == 0:
        raise ContractError(f"No FASTA records found in {path}.")
    if current_id is not None and current_bases == 0:
        raise ContractError(f"FASTA record {current_id!r} has no sequence in {path}.")
    return {
        "computed_sha256": digest.hexdigest(),
        "validated_contigs": n_records,
        "validated_total_bases": total_bases,
        "bytes": path.stat().st_size,
    }


def build_database_catalog(config: dict[str, Any]) -> tuple[list[dict[str, Any]], list[str], dict[str, Any]]:
    database = Path(config["paths"]["database"])
    descriptions: dict[str, str] = {}
    with (database / "config").open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ContractError(f"Malformed VirulenceFinder database config line: {line}")
            descriptions[parts[0].strip()] = parts[2].strip()

    functions: dict[str, str] = {}
    with (database / "notes.txt").open("r", encoding="utf-8") as handle:
        for raw in handle:
            if ":" not in raw:
                continue
            gene, description = raw.split(":", 1)
            functions.setdefault(gene.strip(), description.strip().strip(":"))

    rows: list[dict[str, Any]] = []
    families_by_db: dict[str, set[str]] = defaultdict(set)
    duplicate_headers: dict[str, Counter[str]] = {}
    for db_name in config["software"]["databases"]:
        fasta_path = database / f"{db_name}.fsa"
        headers: list[str] = []
        with fasta_path.open("r", encoding="utf-8") as handle:
            for raw in handle:
                if raw.startswith(">"):
                    headers.append(raw[1:].strip())
        counts = Counter(headers)
        duplicate_headers[db_name] = counts
        for raw_ref_id in sorted(counts):
            parts = raw_ref_id.split(":")
            family = parts[0].strip()
            allele = parts[1].strip() if len(parts) > 1 else ""
            accession = ":".join(parts[2:]).strip() if len(parts) > 2 else ""
            families_by_db[db_name].add(family)
            rows.append(
                {
                    "database": db_name,
                    "database_description": descriptions.get(db_name, ""),
                    "raw_ref_id": raw_ref_id,
                    "gene_family": family,
                    "allele": allele,
                    "accession": accession,
                    "reference_header_occurrences": counts[raw_ref_id],
                    "function": functions.get(family, ""),
                }
            )

    expected_ecoli = int(config["cohort"]["expected_ecoli_gene_families"])
    expected_stx = int(config["cohort"]["expected_stx_gene_families"])
    if len(families_by_db["virulence_ecoli"]) != expected_ecoli:
        raise ContractError(
            f"virulence_ecoli family count changed: {len(families_by_db['virulence_ecoli'])} != {expected_ecoli}"
        )
    if len(families_by_db["stx"]) != expected_stx:
        raise ContractError(f"stx family count changed: {len(families_by_db['stx'])} != {expected_stx}")
    family_universe = sorted(set().union(*families_by_db.values()))
    observed_overlap = sorted(
        families_by_db["virulence_ecoli"] & families_by_db["stx"]
    )
    expected_overlap = sorted(config["cohort"]["expected_cross_database_family_overlap"])
    if observed_overlap != expected_overlap:
        raise ContractError(
            f"Cross-database gene-family overlap changed: {observed_overlap!r}; "
            f"expected {expected_overlap!r}."
        )
    expected_combined = int(config["cohort"]["expected_combined_gene_families"])
    if len(family_universe) != expected_combined:
        raise ContractError(
            f"Combined gene-family universe is {len(family_universe)}, expected {expected_combined}."
        )
    expected_duplicate = config["cohort"]["known_duplicate_reference_header"]
    observed_duplicates = sorted(
        header
        for counts in duplicate_headers.values()
        for header, n in counts.items()
        if n > 1
    )
    if observed_duplicates != [expected_duplicate]:
        raise ContractError(
            f"Official database duplicate-header set changed: {observed_duplicates!r}; "
            f"expected {[expected_duplicate]!r}."
        )
    diagnostics = {
        "families_by_database": {key: len(value) for key, value in families_by_db.items()},
        "combined_family_count": len(family_universe),
        "cross_database_family_overlap": observed_overlap,
        "duplicate_reference_headers": observed_duplicates,
    }
    return rows, family_universe, diagnostics


def collect_environment(config: dict[str, Any], root: Path) -> dict[str, Any]:
    paths = config["paths"]
    software = config["software"]
    python = Path(paths["python"])
    blastn = Path(paths["blastn"])
    database = Path(paths["database"])
    required = [python, blastn, database / "config", database / "VERSION"]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise ContractError("Missing runtime component(s): " + ", ".join(missing))

    tool_version = run_capture([str(python), "-m", "virulencefinder", "--version"], root)
    python_version = run_capture([str(python), "--version"], root)
    blast_version = run_capture([str(blastn), "-version"], root)
    pip_freeze = run_capture([str(python), "-m", "pip", "freeze"], root)
    db_commit = run_capture(["git", "-C", str(database), "rev-parse", "HEAD"], root)
    commands = [tool_version, python_version, blast_version, pip_freeze, db_commit]
    if any(command.returncode != 0 for command in commands):
        labels = ("virulencefinder", "python", "blastn", "pip_freeze", "database_git")
        detail = {
            label: {
                "returncode": command.returncode,
                "stdout": command.stdout.strip(),
                "stderr": command.stderr.strip(),
            }
            for label, command in zip(labels, commands)
            if command.returncode != 0
        }
        raise ContractError("One or more runtime provenance commands failed: " + stable_json(detail))

    observed_tool = tool_version.stdout.strip()
    observed_python = python_version.stdout.strip().replace("Python ", "")
    observed_blast = blast_version.stdout.strip().splitlines()[0]
    observed_db_version = (database / "VERSION").read_text(encoding="utf-8").strip()
    observed_commit = db_commit.stdout.strip()
    checks = {
        "virulencefinder_version": (observed_tool, software["virulencefinder_version"]),
        "python_version": (observed_python, software["python_version"]),
        "blast_version_prefix": (observed_blast, software["blast_version_prefix"]),
        "database_version": (observed_db_version, software["database_version"]),
        "database_commit": (observed_commit, software["database_commit"]),
        "database_config_sha256": (sha256_file(database / "config"), software["database_config_sha256"]),
        "virulence_ecoli_sha256": (
            sha256_file(database / "virulence_ecoli.fsa"),
            software["virulence_ecoli_sha256"],
        ),
        "stx_sha256": (sha256_file(database / "stx.fsa"), software["stx_sha256"]),
    }
    main_path = next((python.parent.parent / "lib").glob("python*/site-packages/virulencefinder/__main__.py"), None)
    if main_path is None:
        raise ContractError("Could not locate the installed VirulenceFinder __main__.py.")
    checks["virulencefinder_main_sha256"] = (
        sha256_file(main_path),
        software["virulencefinder_main_sha256"],
    )
    failed = {
        key: {"observed": observed, "expected": expected}
        for key, (observed, expected) in checks.items()
        if (not str(observed).startswith(str(expected)) if key == "blast_version_prefix" else str(observed) != str(expected))
    }
    if failed:
        raise ContractError("Pinned runtime contract failed: " + stable_json(failed))

    free_gib = shutil.disk_usage(root).free / (1024**3)
    minimum = float(config["execution"]["min_free_gib"])
    if free_gib < minimum:
        raise ContractError(f"Only {free_gib:.1f} GiB free; at least {minimum:.1f} GiB is required.")

    patterns = [str(x) for x in config["execution"].get("heavy_process_patterns", [])]
    heavy_lines: list[str] = []
    if patterns:
        proc = run_capture(["pgrep", "-fl", "|".join(re.escape(x) for x in patterns)], root)
        if proc.returncode == 0:
            heavy_lines = [line for line in proc.stdout.splitlines() if line.strip()]
    return {
        "virulencefinder_version": observed_tool,
        "python_version": observed_python,
        "blast_version": blast_version.stdout.strip(),
        "database_version": observed_db_version,
        "database_commit": observed_commit,
        "database_config_sha256": checks["database_config_sha256"][0],
        "virulence_ecoli_sha256": checks["virulence_ecoli_sha256"][0],
        "stx_sha256": checks["stx_sha256"][0],
        "virulencefinder_main_sha256": checks["virulencefinder_main_sha256"][0],
        "pip_freeze": sorted(line for line in pip_freeze.stdout.splitlines() if line.strip()),
        "free_gib_at_preflight": round(free_gib, 3),
        "heavy_process_matches": heavy_lines,
    }


def load_and_validate_manifest(
    config: dict[str, Any], root: Path, workers: int
) -> list[dict[str, Any]]:
    manifest_path = Path(config["paths"]["manifest"])
    status_path = Path(config["paths"]["clinical_status"])
    if not manifest_path.exists() or not status_path.exists():
        raise ContractError("Manifest or clinical status file is missing.")
    manifest = read_csv(manifest_path)
    status = read_csv(status_path)
    cohort = config["cohort"]
    expected_n = int(cohort["expected_assemblies"])
    if len(manifest) != expected_n:
        raise ContractError(f"Manifest has {len(manifest)} rows; expected {expected_n}.")

    status_by_key: dict[str, dict[str, str]] = {}
    for row in status:
        key = f"{row.get('Participant_id', '')}||{row.get('tp_lab', '')}"
        if key in status_by_key:
            raise ContractError(f"Duplicate clinical status key: {key}")
        status_by_key[key] = row

    seen_paths: set[str] = set()
    seen_assemblies: set[str] = set()
    seen_keys: set[str] = set()
    prepared: list[dict[str, Any]] = []
    required_assembler = str(cohort["required_assembler"]).lower()
    required_organism = str(cohort["required_organism"])
    forbidden = [str(x) for x in cohort.get("forbidden_path_fragments", [])]
    for index, raw in enumerate(manifest):
        row: dict[str, Any] = dict(raw)
        assembler = (raw.get("assembler") or raw.get("Assembler") or "").lower()
        if assembler != required_assembler:
            raise ContractError(f"Non-{required_assembler} assembly in manifest row {index + 1}.")
        for field in ("QC_PASS", "selected_canonical", "file_exists", "usable_fasta"):
            if not parse_bool(raw.get(field)):
                raise ContractError(f"Manifest row {index + 1} has {field} != TRUE.")
        if raw.get("Clinical_Organism") != required_organism:
            raise ContractError(
                f"Manifest row {index + 1} organism {raw.get('Clinical_Organism')!r} != {required_organism!r}."
            )
        fasta = Path(raw.get("full_path") or raw.get("fasta_path") or "").resolve()
        if not fasta.is_file():
            raise ContractError(f"Selected FASTA does not exist: {fasta}")
        if any(fragment in str(fasta) for fragment in forbidden):
            raise ContractError(f"Forbidden input path detected: {fasta}")
        assembly = raw.get("Assembly_ID", "")
        key = f"{raw.get('Participant_id', '')}||{raw.get('tp_lab', '')}"
        if not assembly or key.startswith("||"):
            raise ContractError(f"Missing Assembly_ID or episode key in manifest row {index + 1}.")
        if str(fasta) in seen_paths or assembly in seen_assemblies or key in seen_keys:
            raise ContractError(f"Duplicate selected path, assembly, or episode key at row {index + 1}.")
        seen_paths.add(str(fasta))
        seen_assemblies.add(assembly)
        seen_keys.add(key)
        clinical = status_by_key.get(key)
        if clinical is None:
            raise ContractError(f"Selected episode is absent from clinical status map: {key}")
        row.update(
            {
                "manifest_index": index + 1,
                "full_path": str(fasta),
                "episode_key": key,
                "UTI_Status": clinical.get("UTI_Status", ""),
                "Event_type": clinical.get("Event_type", raw.get("Event_type", "")),
                "assembler": assembler,
            }
        )
        prepared.append(row)

    if len({row["Participant_id"] for row in prepared}) != int(cohort["expected_residents"]):
        raise ContractError("Selected resident denominator changed.")
    uti_counts = Counter(row["UTI_Status"] for row in prepared)
    if uti_counts["UTI"] != int(cohort["expected_uti"]) or uti_counts["Not_UTI"] != int(
        cohort["expected_not_uti"]
    ):
        raise ContractError(f"Operational UTI denominator changed: {dict(uti_counts)}")

    print(f"Validating and hashing {len(prepared)} selected FASTAs ...", flush=True)
    validations: dict[str, dict[str, Any]] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, min(workers, 4))) as pool:
        future_map = {
            pool.submit(validate_fasta_and_hash, Path(row["full_path"])): row for row in prepared
        }
        for completed, future in enumerate(concurrent.futures.as_completed(future_map), start=1):
            row = future_map[future]
            result = future.result()
            stored = str(row.get("fasta_sha256", "")).lower()
            if result["computed_sha256"] != stored:
                raise ContractError(
                    f"FASTA SHA-256 mismatch for {row['Assembly_ID']}: "
                    f"{result['computed_sha256']} != {stored}"
                )
            validations[row["Assembly_ID"]] = result
            if completed % 50 == 0 or completed == len(prepared):
                print(f"  FASTA validation: {completed}/{len(prepared)}", flush=True)
    for row in prepared:
        row.update(validations[row["Assembly_ID"]])
    return prepared


def build_context(
    config: dict[str, Any], environment: dict[str, Any], row: dict[str, Any], profile_name: str
) -> tuple[dict[str, Any], str]:
    profile = config["profiles"][profile_name]
    context = {
        "runner_schema_version": RUNNER_SCHEMA_VERSION,
        "assembly_id": row["Assembly_ID"],
        "episode_key": row["episode_key"],
        "full_path": row["full_path"],
        "fasta_sha256": row["fasta_sha256"],
        "profile": profile_name,
        "minimum_identity": float(profile["minimum_identity"]),
        "minimum_coverage": float(profile["minimum_coverage"]),
        "overlap_nt": int(config["execution"]["overlap_nt"]),
        "databases": list(config["software"]["databases"]),
        "virulencefinder_version": environment["virulencefinder_version"],
        "python_version": environment["python_version"],
        "blast_version": environment["blast_version"],
        "database_version": environment["database_version"],
        "database_commit": environment["database_commit"],
        "database_config_sha256": environment["database_config_sha256"],
        "virulence_ecoli_sha256": environment["virulence_ecoli_sha256"],
        "stx_sha256": environment["stx_sha256"],
        "virulencefinder_main_sha256": environment["virulencefinder_main_sha256"],
    }
    return context, sha256_text(stable_json(context))


def slugify(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._") or "sample"


def build_command(
    config: dict[str, Any], row: dict[str, Any], profile_name: str, output_dir: Path, json_path: Path
) -> list[str]:
    profile = config["profiles"][profile_name]
    command = [
        config["paths"]["python"],
        "-m",
        "virulencefinder",
        "-ifa",
        row["full_path"],
        "-o",
        str(output_dir),
        "-j",
        str(json_path),
        "-p",
        config["paths"]["database"],
        "-d",
        ",".join(config["software"]["databases"]),
        "-b",
        config["paths"]["blastn"],
        "-t",
        f"{float(profile['minimum_identity']):.6g}",
        "-l",
        f"{float(profile['minimum_coverage']):.6g}",
        "--overlap",
        str(int(config["execution"]["overlap_nt"])),
        "-q",
    ]
    if command.count("-ifa") != 1 or command[command.index("-ifa") + 1] != row["full_path"]:
        raise ContractError("Command builder violated the one-FASTA-per-process contract.")
    return command


def canonical_hits_from_data(
    data: dict[str, Any], config: dict[str, Any], row: dict[str, Any], profile_name: str
) -> list[dict[str, Any]]:
    software = config["software"]
    profile = config["profiles"][profile_name]
    if str(data.get("software_version")) != str(software["virulencefinder_version"]):
        raise ContractError("VirulenceFinder output software version mismatch.")
    databases = data.get("databases")
    if not isinstance(databases, dict) or len(databases) != 1:
        raise ContractError("VirulenceFinder output database metadata are missing or ambiguous.")
    database_meta = next(iter(databases.values()))
    if str(database_meta.get("database_version")) != str(software["database_version"]):
        raise ContractError("VirulenceFinder output database version mismatch.")
    if str(database_meta.get("database_commit")) != str(software["database_commit"]):
        raise ContractError("VirulenceFinder output database commit mismatch.")
    executions = data.get("software_executions")
    if not isinstance(executions, dict) or len(executions) != 1:
        raise ContractError("VirulenceFinder output must contain exactly one software execution.")
    parameters = next(iter(executions.values())).get("parameters", {})
    if str(Path(parameters.get("inputfasta", "")).resolve()) != row["full_path"]:
        raise ContractError("VirulenceFinder output input path does not equal the selected FASTA.")
    if parameters.get("inputfastq_1") is not None or parameters.get("aligner") != "blastn":
        raise ContractError("VirulenceFinder output is not an assembled-FASTA BLAST run.")
    if set(parameters.get("databases", [])) != set(software["databases"]):
        raise ContractError("VirulenceFinder output selected database set changed.")
    numeric_checks = {
        "min_id": float(profile["minimum_identity"]),
        "min_cov": float(profile["minimum_coverage"]),
        "overlap": float(config["execution"]["overlap_nt"]),
    }
    for key, expected in numeric_checks.items():
        if not math.isclose(float(parameters.get(key, math.nan)), expected, rel_tol=0, abs_tol=1e-9):
            raise ContractError(f"VirulenceFinder output parameter {key} changed.")

    regions = data.get("seq_regions", {})
    phenotypes = data.get("phenotypes", {})
    if not isinstance(regions, dict) or not isinstance(phenotypes, dict):
        raise ContractError("VirulenceFinder output result collections are malformed.")
    min_identity_pct = 100 * float(profile["minimum_identity"])
    min_coverage_pct = 100 * float(profile["minimum_coverage"])
    deduplicated: dict[tuple[Any, ...], dict[str, Any]] = {}
    result_keys: dict[tuple[Any, ...], list[str]] = defaultdict(list)
    for raw_result_key, region in regions.items():
        if not isinstance(region, dict):
            raise ContractError("VirulenceFinder seq_region entry is not an object.")
        identity = float(region.get("identity", math.nan))
        fragment_coverage = float(region.get("coverage", math.nan))
        if not math.isfinite(identity) or identity + 1e-7 < min_identity_pct:
            raise ContractError(f"Hit below identity threshold in {row['Assembly_ID']}: {identity}")
        if not math.isfinite(fragment_coverage):
            raise ContractError(f"Hit has non-finite coverage in {row['Assembly_ID']}")
        ref_id = str(region.get("ref_id", ""))
        family = str(region.get("name", "")) or ref_id.split(":", 1)[0]
        if not ref_id or not family:
            raise ContractError("VirulenceFinder hit lacks ref_id or gene-family name.")
        ref_databases = sorted(
            {str(value).split(":")[-1] for value in region.get("ref_database", [])}
        )
        if not ref_databases or not set(ref_databases).issubset(set(software["databases"])):
            raise ContractError(f"Unexpected hit database(s): {ref_databases}")
        parts = ref_id.split(":")
        allele = parts[1] if len(parts) > 1 else ""
        accession = ":".join(parts[2:]) if len(parts) > 2 else str(region.get("ref_acc", ""))
        functions = sorted(
            {
                str(phenotypes[key].get("function", ""))
                for key in region.get("phenotypes", [])
                if key in phenotypes and phenotypes[key].get("function")
            }
        )
        hit = {
            "profile": profile_name,
            "Participant_id": row["Participant_id"],
            "tp_lab": row["tp_lab"],
            "episode_key": row["episode_key"],
            "Assembly_ID": row["Assembly_ID"],
            "full_path": row["full_path"],
            "fasta_sha256": row["fasta_sha256"],
            "database": ";".join(ref_databases),
            "gene_family": family,
            "allele": allele,
            "accession": accession,
            "raw_ref_id": ref_id,
            "identity_pct": identity,
            # VirulenceFinder may emit multiple low-coverage fragments for a
            # putative split gene. The aggregate threshold is applied below;
            # both the reported fragment value and aggregate value are kept.
            "reported_fragment_coverage_pct": fragment_coverage,
            "coverage_pct": fragment_coverage,
            "grade": region.get("grade", ""),
            "query_id": region.get("query_id", ""),
            "query_start": region.get("query_start_pos", ""),
            "query_end": region.get("query_end_pos", ""),
            "reference_start": region.get("ref_start_pos", ""),
            "reference_end": region.get("ref_end_pos", ""),
            "alignment_length": region.get("alignment_length", ""),
            "reference_length": region.get("ref_seq_length", ""),
            "phenotype_function": "; ".join(functions),
        }
        biological_key = (
            hit["database"],
            family,
            ref_id,
            hit["query_id"],
            hit["query_start"],
            hit["query_end"],
            hit["reference_start"],
            hit["reference_end"],
        )
        result_keys[biological_key].append(str(raw_result_key))
        existing = deduplicated.get(biological_key)
        if existing is None or (identity, fragment_coverage) > (
            float(existing["identity_pct"]),
            float(existing["coverage_pct"]),
        ):
            deduplicated[biological_key] = hit

    # Reapply the coverage contract at the raw reference level. This handles
    # legitimate split genes while excluding a documented cgecore edge case
    # where low-coverage fragments can survive the tool's internal filter.
    grouped: dict[tuple[str, str], list[tuple[tuple[Any, ...], dict[str, Any]]]] = defaultdict(list)
    for key, hit in deduplicated.items():
        grouped[(hit["database"], hit["raw_ref_id"])].append((key, hit))

    retained: dict[tuple[Any, ...], dict[str, Any]] = {}
    for group_hits in grouped.values():
        ref_lengths = {int(item[1]["reference_length"]) for item in group_hits}
        if len(ref_lengths) != 1:
            raise ContractError("Split-reference hits disagree on reference length.")
        ref_length = next(iter(ref_lengths))
        intervals = sorted(
            (
                min(int(item[1]["reference_start"]), int(item[1]["reference_end"])),
                max(int(item[1]["reference_start"]), int(item[1]["reference_end"])),
            )
            for item in group_hits
        )
        covered = 0
        current_start: int | None = None
        current_end: int | None = None
        for start, end in intervals:
            if current_start is None:
                current_start, current_end = start, end
            elif start <= int(current_end) + 1:
                current_end = max(int(current_end), end)
            else:
                covered += int(current_end) - int(current_start) + 1
                current_start, current_end = start, end
        if current_start is not None:
            covered += int(current_end) - int(current_start) + 1
        aggregate_coverage = 100 * covered / ref_length if ref_length > 0 else math.nan
        if not math.isfinite(aggregate_coverage) or aggregate_coverage + 1e-7 < min_coverage_pct:
            continue
        for key, hit in group_hits:
            hit["aggregate_reference_coverage_pct"] = aggregate_coverage
            hit["coverage_pct"] = aggregate_coverage
            hit["coverage_filter_basis"] = (
                "individual_hit"
                if float(hit["reported_fragment_coverage_pct"]) + 1e-7 >= min_coverage_pct
                else "split_reference_union"
            )
            retained[key] = hit

    hits = []
    for key in sorted(retained, key=lambda item: tuple(str(x) for x in item)):
        hit = retained[key]
        hit["raw_result_keys"] = ";".join(sorted(result_keys[key]))
        hit["collapsed_duplicate_result_keys"] = max(0, len(result_keys[key]) - 1)
        hits.append(hit)
    return hits


def validate_result_file(
    path: Path, config: dict[str, Any], row: dict[str, Any], profile_name: str
) -> tuple[dict[str, Any], list[dict[str, Any]], str]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"Could not parse VirulenceFinder JSON: {path}") from exc
    hits = canonical_hits_from_data(data, config, row, profile_name)
    canonical_hash = sha256_text(stable_json(hits))
    return data, hits, canonical_hash


def job_paths(
    config: dict[str, Any], row: dict[str, Any], profile_name: str, context_hash: str
) -> tuple[Path, Path, Path]:
    output = Path(config["paths"]["output"])
    stem = f"{slugify(row['Assembly_ID'])}.{row['fasta_sha256'][:16]}.{context_hash[:12]}"
    directory = output / "cache" / profile_name / stem
    return directory, directory / "result.json", directory / "meta.json"


def validate_cached_job(
    config: dict[str, Any], row: dict[str, Any], profile_name: str, context_hash: str
) -> tuple[bool, dict[str, Any] | None, list[dict[str, Any]] | None, str]:
    directory, result_path, meta_path = job_paths(config, row, profile_name, context_hash)
    if not result_path.exists() or not meta_path.exists():
        return False, None, None, "missing_cache_files"
    try:
        with meta_path.open("r", encoding="utf-8") as handle:
            meta = json.load(handle)
        if meta.get("status") != "success" or meta.get("context_hash") != context_hash:
            return False, meta, None, "cache_metadata_context_mismatch"
        _, hits, canonical_hash = validate_result_file(result_path, config, row, profile_name)
        if meta.get("canonical_hits_sha256") != canonical_hash:
            return False, meta, None, "cache_canonical_hit_hash_mismatch"
        if meta.get("result_json_sha256") != sha256_file(result_path):
            return False, meta, None, "cache_json_hash_mismatch"
        return True, meta, hits, "valid"
    except (OSError, json.JSONDecodeError, ContractError) as exc:
        return False, None, None, f"cache_validation_error:{exc}"


def write_failure_diagnostic(
    config: dict[str, Any], row: dict[str, Any], profile_name: str, context_hash: str,
    attempt: int, command: list[str], completed: subprocess.CompletedProcess[str] | None,
    error: str, result_path: Path | None,
) -> None:
    base = Path(config["paths"]["output"]) / "failures" / profile_name
    stem = f"{slugify(row['Assembly_ID'])}.{context_hash[:12]}.attempt{attempt}"
    diagnostic = {
        "generated_at": utc_now(),
        "Assembly_ID": row["Assembly_ID"],
        "episode_key": row["episode_key"],
        "profile": profile_name,
        "context_hash": context_hash,
        "attempt": attempt,
        "command": command,
        "returncode": completed.returncode if completed is not None else None,
        "error": error,
    }
    atomic_write_json(base / f"{stem}.json", diagnostic)
    if completed is not None:
        atomic_write_text(base / f"{stem}.stdout.txt", completed.stdout or "")
        atomic_write_text(base / f"{stem}.stderr.txt", completed.stderr or "")
    if result_path is not None and result_path.exists():
        shutil.copy2(result_path, base / f"{stem}.result.json")


def run_one_job(
    config: dict[str, Any], environment: dict[str, Any], row: dict[str, Any],
    profile_name: str, resume: bool, attempt: int,
) -> dict[str, Any]:
    context, context_hash = build_context(config, environment, row, profile_name)
    directory, result_path, meta_path = job_paths(config, row, profile_name, context_hash)
    if resume:
        valid, meta, hits, reason = validate_cached_job(
            config, row, profile_name, context_hash
        )
        if valid and meta is not None and hits is not None:
            return {
                "manifest_index": row["manifest_index"],
                "profile": profile_name,
                "Assembly_ID": row["Assembly_ID"],
                "episode_key": row["episode_key"],
                "fasta_sha256": row["fasta_sha256"],
                "context_hash": context_hash,
                "cache_dir": str(directory),
                "status": "success",
                "cache_reused": True,
                "attempt": attempt,
                "elapsed_seconds": meta.get("elapsed_seconds", ""),
                "n_hits": len(hits),
                "failure_reason": "",
            }

    if directory.exists():
        shutil.rmtree(directory)
    output = Path(config["paths"]["output"])
    tmp_root = output / ".tmp"
    tmp_root.mkdir(parents=True, exist_ok=True)
    tmp_dir = Path(tempfile.mkdtemp(prefix=f"{slugify(row['Assembly_ID'])}.{profile_name}.", dir=tmp_root))
    tmp_json = tmp_dir / "result.json"
    command = build_command(config, row, profile_name, tmp_dir, tmp_json)
    started = time.monotonic()
    completed: subprocess.CompletedProcess[str] | None = None
    error = ""
    try:
        completed = run_capture(command, Path.cwd())
        elapsed = time.monotonic() - started
        if completed.returncode != 0:
            raise ContractError(f"VirulenceFinder exit status {completed.returncode}")
        _, hits, canonical_hash = validate_result_file(tmp_json, config, row, profile_name)
        directory.mkdir(parents=True, exist_ok=False)
        shutil.copy2(tmp_json, result_path)
        atomic_write_text(directory / "stderr.txt", completed.stderr or "")
        meta = {
            "status": "success",
            "generated_at": utc_now(),
            "context": context,
            "context_hash": context_hash,
            "command": command,
            "attempt": attempt,
            "elapsed_seconds": round(elapsed, 6),
            "n_hits": len(hits),
            "zero_hit_success": len(hits) == 0,
            "canonical_hits_sha256": canonical_hash,
            "result_json_sha256": sha256_file(result_path),
        }
        atomic_write_json(meta_path, meta)
        return {
            "manifest_index": row["manifest_index"],
            "profile": profile_name,
            "Assembly_ID": row["Assembly_ID"],
            "episode_key": row["episode_key"],
            "fasta_sha256": row["fasta_sha256"],
            "context_hash": context_hash,
            "cache_dir": str(directory),
            "status": "success",
            "cache_reused": False,
            "attempt": attempt,
            "elapsed_seconds": round(elapsed, 6),
            "n_hits": len(hits),
            "failure_reason": "",
        }
    except Exception as exc:  # noqa: BLE001 - diagnostic boundary around external process
        error = str(exc)
        elapsed = time.monotonic() - started
        write_failure_diagnostic(
            config, row, profile_name, context_hash, attempt, command, completed, error, tmp_json
        )
        return {
            "manifest_index": row["manifest_index"],
            "profile": profile_name,
            "Assembly_ID": row["Assembly_ID"],
            "episode_key": row["episode_key"],
            "fasta_sha256": row["fasta_sha256"],
            "context_hash": context_hash,
            "cache_dir": str(directory),
            "status": "error",
            "cache_reused": False,
            "attempt": attempt,
            "elapsed_seconds": round(elapsed, 6),
            "n_hits": "",
            "failure_reason": error,
        }
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


def task_key(result: dict[str, Any]) -> tuple[int, str]:
    return int(result["manifest_index"]), str(result["profile"])


def run_jobs(
    config: dict[str, Any], environment: dict[str, Any], rows: list[dict[str, Any]],
    workers: int, resume: bool,
) -> list[dict[str, Any]]:
    tasks = [(row, profile) for row in rows for profile in EXPECTED_PROFILE_NAMES]
    total = len(tasks)
    results: dict[tuple[int, str], dict[str, Any]] = {}
    last_report = time.monotonic()
    report_jobs = int(config["execution"]["progress_every_jobs"])
    report_seconds = int(config["execution"]["progress_every_seconds"])
    started = time.monotonic()
    print(f"Launching {total} jobs with {workers} worker(s) ...", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        future_map = {
            pool.submit(run_one_job, config, environment, row, profile, resume, 1): (row, profile)
            for row, profile in tasks
        }
        for completed_n, future in enumerate(concurrent.futures.as_completed(future_map), start=1):
            result = future.result()
            results[task_key(result)] = result
            now = time.monotonic()
            if completed_n % report_jobs == 0 or now - last_report >= report_seconds or completed_n == total:
                elapsed = max(now - started, 1e-9)
                rate = completed_n / elapsed
                eta = (total - completed_n) / rate if rate > 0 else math.nan
                failures = sum(value["status"] != "success" for value in results.values())
                reused = sum(parse_bool(value["cache_reused"]) for value in results.values())
                print(
                    f"  progress {completed_n}/{total}; failures={failures}; reused={reused}; "
                    f"elapsed={elapsed / 60:.1f} min; ETA={eta / 60:.1f} min",
                    flush=True,
                )
                last_report = now

    retry_count = int(config["execution"]["retry_failures"])
    failures = [value for value in results.values() if value["status"] != "success"]
    for retry_index in range(1, retry_count + 1):
        if not failures:
            break
        print(f"Retrying {len(failures)} failed job(s) sequentially (retry {retry_index}) ...", flush=True)
        new_failures = []
        for old in failures:
            row_by_manifest_index = {int(item["manifest_index"]): item for item in rows}
            row = row_by_manifest_index[int(old["manifest_index"])]
            result = run_one_job(
                config, environment, row, str(old["profile"]), False, retry_index + 1
            )
            results[task_key(result)] = result
            if result["status"] != "success":
                new_failures.append(result)
        failures = new_failures
    return sorted(results.values(), key=task_key)


def select_pilot_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ranked = sorted(rows, key=lambda row: (int(float(row.get("n_contigs") or 0)), row["Assembly_ID"]))
    candidates = [ranked[0], ranked[(len(ranked) - 1) // 2], ranked[-1]]
    unique: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in candidates + ranked:
        if row["Assembly_ID"] not in seen:
            unique.append(row)
            seen.add(row["Assembly_ID"])
        if len(unique) == 3:
            break
    return unique


def run_pilot(
    config: dict[str, Any], environment: dict[str, Any], rows: list[dict[str, Any]], workers: int
) -> None:
    pilot_rows = select_pilot_rows(rows)
    first_pass = run_jobs(config, environment, pilot_rows, min(workers, 3), True)
    if any(result["status"] != "success" for result in first_pass):
        raise ContractError("Pilot execution failed; full batch is blocked.")
    second_pass = run_jobs(config, environment, pilot_rows, min(workers, 3), True)
    if not all(parse_bool(result["cache_reused"]) for result in second_pass):
        raise ContractError("Pilot resume check did not reuse all six valid cache entries.")
    first_hits = {(r["manifest_index"], r["profile"]): int(r["n_hits"]) for r in first_pass}
    second_hits = {(r["manifest_index"], r["profile"]): int(r["n_hits"]) for r in second_pass}
    if first_hits != second_hits:
        raise ContractError("Pilot canonical hit counts changed across forced reparse/cache validation.")
    pilot_report = []
    for result in second_pass:
        row = rows[int(result["manifest_index"]) - 1]
        pilot_report.append(
            {
                **result,
                "n_contigs": row.get("n_contigs", ""),
                "pilot_cache_reparse_pass": True,
            }
        )
    fields = list(pilot_report[0].keys())
    atomic_write_csv(Path(config["paths"]["output"]) / "pilot_status.csv", pilot_report, fields)
    atomic_write_text(
        Path(config["paths"]["output"]) / "PILOT_COMPLETE.txt",
        f"VirulenceFinder pilot: PASS\nGenerated: {utc_now()}\nJobs: 6\nCache reparse: PASS\n",
    )


def marker_value(present: set[str], genes: list[str]) -> int:
    return int(any(gene in present for gene in genes))


def build_primary_concordance(
    config: dict[str, Any], rows: list[dict[str, Any]], family_universe: list[str],
    matched_presence: dict[str, set[str]], concordance_dir: Path,
) -> dict[str, Any]:
    primary_path = Path(config["paths"]["primary_vf_presence_absence"])
    primary_rows = read_csv(primary_path)
    if len(primary_rows) != int(config["cohort"]["expected_assemblies"]):
        raise ContractError("Primary VFDB presence/absence row count changed.")
    primary_meta = {"Participant_id", "tp_lab", "episode_key", "Assembly_ID", "fasta_sha256"}
    primary_genes = [field for field in primary_rows[0] if field not in primary_meta]
    primary_casefold: dict[str, list[str]] = defaultdict(list)
    cge_casefold: dict[str, list[str]] = defaultdict(list)
    for gene in primary_genes:
        primary_casefold[gene.casefold()].append(gene)
    for gene in family_universe:
        cge_casefold[gene.casefold()].append(gene)
    crosswalk: list[dict[str, Any]] = []
    mapping: dict[str, str] = {}
    for cge in family_universe:
        if cge in primary_genes:
            primary = cge
            kind = "exact"
        elif len(cge_casefold[cge.casefold()]) == 1 and len(primary_casefold.get(cge.casefold(), [])) == 1:
            primary = primary_casefold[cge.casefold()][0]
            kind = "case_normalized_one_to_one"
        else:
            primary = ""
            kind = "unmapped"
        if primary:
            mapping[cge] = primary
        crosswalk.append(
            {
                "cge_gene_family": cge,
                "primary_vfdb_gene": primary,
                "mapping_type": kind,
                "eligible_for_concordance": bool(primary),
            }
        )
    atomic_write_csv(
        concordance_dir / "cge_to_primary_vfdb_gene_crosswalk.csv",
        crosswalk,
        list(crosswalk[0].keys()),
    )
    primary_by_key = {row["episode_key"]: row for row in primary_rows}
    expected_keys = {row["episode_key"] for row in rows}
    if set(primary_by_key) != expected_keys:
        raise ContractError("Primary VFDB matrix episode keys do not equal the selected manifest keys.")
    long_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    per_episode: list[dict[str, Any]] = []
    for cge, primary in sorted(mapping.items()):
        counts = Counter()
        for row in rows:
            key = row["episode_key"]
            cge_value = int(cge in matched_presence[key])
            primary_value = int(float(primary_by_key[key].get(primary, 0) or 0) > 0)
            counts[(primary_value, cge_value)] += 1
            long_rows.append(
                {
                    "episode_key": key,
                    "Assembly_ID": row["Assembly_ID"],
                    "fasta_sha256": row["fasta_sha256"],
                    "cge_gene_family": cge,
                    "primary_vfdb_gene": primary,
                    "primary_present": primary_value,
                    "cge_matched_80_80_present": cge_value,
                    "agreement": int(primary_value == cge_value),
                }
            )
        tn, fp, fn, tp = counts[(0, 0)], counts[(0, 1)], counts[(1, 0)], counts[(1, 1)]
        n = tn + fp + fn + tp
        observed = (tn + tp) / n
        primary_rate = (fn + tp) / n
        cge_rate = (fp + tp) / n
        expected = primary_rate * cge_rate + (1 - primary_rate) * (1 - cge_rate)
        kappa = (observed - expected) / (1 - expected) if expected < 1 else math.nan
        summaries.append(
            {
                "cge_gene_family": cge,
                "primary_vfdb_gene": primary,
                "n": n,
                "true_positive": tp,
                "true_negative": tn,
                "primary_only": fn,
                "cge_only": fp,
                "agreement": observed,
                "cohen_kappa": kappa,
                "primary_prevalence": primary_rate,
                "cge_prevalence": cge_rate,
            }
        )
    shared_cge = sorted(mapping)
    shared_primary = [mapping[gene] for gene in shared_cge]
    for row in rows:
        key = row["episode_key"]
        cge_set = {gene for gene in shared_cge if gene in matched_presence[key]}
        primary_set = {
            cge for cge, primary in zip(shared_cge, shared_primary)
            if int(float(primary_by_key[key].get(primary, 0) or 0) > 0)
        }
        union = cge_set | primary_set
        per_episode.append(
            {
                "episode_key": key,
                "Assembly_ID": row["Assembly_ID"],
                "fasta_sha256": row["fasta_sha256"],
                "shared_gene_universe_n": len(shared_cge),
                "primary_shared_burden": len(primary_set),
                "cge_shared_burden": len(cge_set),
                "intersection_n": len(cge_set & primary_set),
                "union_n": len(union),
                "jaccard": len(cge_set & primary_set) / len(union) if union else 1.0,
                "exact_shared_profile_match": int(cge_set == primary_set),
            }
        )
    if long_rows:
        atomic_write_csv(
            concordance_dir / "primary_shared_gene_episode_concordance.csv",
            long_rows,
            list(long_rows[0].keys()),
        )
        atomic_write_csv(
            concordance_dir / "primary_shared_gene_summary.csv",
            summaries,
            list(summaries[0].keys()),
        )
        atomic_write_csv(
            concordance_dir / "primary_shared_profile_episode_summary.csv",
            per_episode,
            list(per_episode[0].keys()),
        )
    return {"mapped_shared_gene_families": len(mapping)}


def combine_outputs(
    config: dict[str, Any], environment: dict[str, Any], rows: list[dict[str, Any]],
    catalog: list[dict[str, Any]], family_universe: list[str], run_results: list[dict[str, Any]],
) -> dict[str, Any]:
    output = Path(config["paths"]["output"])
    output.mkdir(parents=True, exist_ok=True)
    atomic_write_text(
        output / "RESTRICTED_INTERNAL_DATA.txt",
        "RESTRICTED INTERNAL DATA\n"
        "\n"
        "The batch-level manifests, hit tables, presence/absence matrices and episode metrics "
        "contain internal episode keys and/or participant identifiers needed for validated joins.\n"
        "Research-facing RQ06-RQ08 pair and case tables are separately deidentified below "
        "rq_sensitivity/ using RES, EPI and PAIR labels. Do not publish the internal tables "
        "without the study's approved disclosure controls.\n",
    )
    expected_jobs = int(config["cohort"]["expected_assemblies"]) * len(EXPECTED_PROFILE_NAMES)
    if len(run_results) != expected_jobs or any(result["status"] != "success" for result in run_results):
        raise ContractError("Cannot combine: the 1,064-job run manifest is incomplete or contains failures.")
    result_by_key = {task_key(result): result for result in run_results}
    if len(result_by_key) != expected_jobs:
        raise ContractError("Duplicate or missing assembly/profile jobs in run manifest.")

    all_hits: list[dict[str, Any]] = []
    presence: dict[str, dict[str, set[str]]] = {
        profile: {} for profile in EXPECTED_PROFILE_NAMES
    }
    for row in rows:
        for profile in EXPECTED_PROFILE_NAMES:
            context, context_hash = build_context(config, environment, row, profile)
            del context
            valid, meta, hits, reason = validate_cached_job(config, row, profile, context_hash)
            if not valid or meta is None or hits is None:
                raise ContractError(f"Cache became invalid during combination: {row['Assembly_ID']} {profile}: {reason}")
            all_hits.extend(hits)
            presence[profile][row["episode_key"]] = {hit["gene_family"] for hit in hits}

    hit_fields = [
        "profile", "Participant_id", "tp_lab", "episode_key", "Assembly_ID", "full_path",
        "fasta_sha256", "database", "gene_family", "allele", "accession", "raw_ref_id",
        "identity_pct", "reported_fragment_coverage_pct", "aggregate_reference_coverage_pct",
        "coverage_pct", "coverage_filter_basis", "grade", "query_id", "query_start", "query_end",
        "reference_start", "reference_end", "alignment_length", "reference_length",
        "phenotype_function", "raw_result_keys", "collapsed_duplicate_result_keys",
    ]
    all_hits.sort(
        key=lambda hit: (
            EXPECTED_PROFILE_NAMES.index(hit["profile"]),
            next(row["manifest_index"] for row in rows if row["episode_key"] == hit["episode_key"]),
            hit["gene_family"], hit["raw_ref_id"], str(hit["query_id"]), str(hit["query_start"]),
        )
    )
    atomic_write_csv(output / "hits_long.csv", all_hits, hit_fields)
    atomic_write_csv(
        output / "gene_reference_catalog.csv",
        catalog,
        [
            "database", "database_description", "raw_ref_id", "gene_family", "allele",
            "accession", "reference_header_occurrences", "function",
        ],
    )

    metrics: list[dict[str, Any]] = []
    markers = {name: list(genes) for name, genes in config["markers"].items()}
    for profile in EXPECTED_PROFILE_NAMES:
        pa_rows: list[dict[str, Any]] = []
        for row in rows:
            key = row["episode_key"]
            present = presence[profile][key]
            base = {column: row.get(column, "") for column in METADATA_COLUMNS}
            pa_rows.append({**base, **{gene: int(gene in present) for gene in family_universe}})
            group_values = {name: marker_value(present, genes) for name, genes in markers.items()}
            metrics.append(
                {
                    "profile": profile,
                    **base,
                    "raw_hit_count": sum(
                        1 for hit in all_hits if hit["profile"] == profile and hit["episode_key"] == key
                    ),
                    "distinct_cge_gene_family_count": len(present),
                    **{f"cge_{name}_group": value for name, value in group_values.items()},
                    "CGE_available_ExPEC_group_count": sum(group_values.values()),
                    "stx1_present": int("stx1" in present),
                    "stx2_present": int("stx2" in present),
                }
            )
        atomic_write_csv(
            output / f"presence_absence_{profile}.csv",
            pa_rows,
            list(METADATA_COLUMNS) + family_universe,
        )
    metric_fields = list(metrics[0].keys())
    atomic_write_csv(output / "episode_metrics.csv", metrics, metric_fields)

    concordance_dir = output / "concordance"
    concordance_dir.mkdir(parents=True, exist_ok=True)
    within_rows: list[dict[str, Any]] = []
    for row in rows:
        key = row["episode_key"]
        web = presence["web_default_id90_cov60"][key]
        matched = presence["matched_id80_cov80"][key]
        union = web | matched
        within_rows.append(
            {
                "episode_key": key,
                "Assembly_ID": row["Assembly_ID"],
                "fasta_sha256": row["fasta_sha256"],
                "web_default_burden": len(web),
                "matched_80_80_burden": len(matched),
                "matched_minus_web_burden": len(matched) - len(web),
                "intersection_n": len(web & matched),
                "union_n": len(union),
                "jaccard": len(web & matched) / len(union) if union else 1.0,
                "exact_profile_match": int(web == matched),
            }
        )
    atomic_write_csv(
        concordance_dir / "within_cge_profile_episode_concordance.csv",
        within_rows,
        list(within_rows[0].keys()),
    )
    prevalence_rows: list[dict[str, Any]] = []
    for gene in family_universe:
        web_n = sum(gene in presence["web_default_id90_cov60"][row["episode_key"]] for row in rows)
        matched_n = sum(gene in presence["matched_id80_cov80"][row["episode_key"]] for row in rows)
        prevalence_rows.append(
            {
                "gene_family": gene,
                "web_default_present_n": web_n,
                "web_default_prevalence": web_n / len(rows),
                "matched_80_80_present_n": matched_n,
                "matched_80_80_prevalence": matched_n / len(rows),
                "matched_minus_web_prevalence": (matched_n - web_n) / len(rows),
            }
        )
    atomic_write_csv(
        concordance_dir / "within_cge_profile_gene_prevalence.csv",
        prevalence_rows,
        list(prevalence_rows[0].keys()),
    )
    primary_diag = build_primary_concordance(
        config,
        rows,
        family_universe,
        presence["matched_id80_cov80"],
        concordance_dir,
    )

    run_fields = [
        "manifest_index", "profile", "Assembly_ID", "episode_key", "fasta_sha256",
        "context_hash", "cache_dir", "status", "cache_reused", "attempt",
        "elapsed_seconds", "n_hits", "failure_reason",
    ]
    atomic_write_csv(output / "run_manifest.csv", run_results, run_fields)
    atomic_write_text(
        output / "BATCH_COMPLETE.txt",
        f"VirulenceFinder batch: PASS\nGenerated: {utc_now()}\n"
        f"Assemblies: {len(rows)}\nProfiles: {len(EXPECTED_PROFILE_NAMES)}\nJobs: {len(run_results)}\n",
    )
    return {
        "assemblies": len(rows),
        "profiles": len(EXPECTED_PROFILE_NAMES),
        "jobs": len(run_results),
        "canonical_hits": len(all_hits),
        "gene_families": len(family_universe),
        **primary_diag,
    }


def write_preflight_outputs(
    config_path: Path, config: dict[str, Any], environment: dict[str, Any],
    rows: list[dict[str, Any]], database_diagnostics: dict[str, Any], root: Path,
) -> None:
    output = Path(config["paths"]["output"])
    output.mkdir(parents=True, exist_ok=True)
    preflight_rows = [
        {
            "manifest_index": row["manifest_index"],
            "Participant_id": row["Participant_id"],
            "tp_lab": row["tp_lab"],
            "episode_key": row["episode_key"],
            "Assembly_ID": row["Assembly_ID"],
            "full_path": row["full_path"],
            "fasta_sha256": row["fasta_sha256"],
            "validated_contigs": row["validated_contigs"],
            "validated_total_bases": row["validated_total_bases"],
            "bytes": row["bytes"],
            "UTI_Status": row["UTI_Status"],
            "Event_type": row["Event_type"],
        }
        for row in rows
    ]
    atomic_write_csv(output / "preflight_manifest.csv", preflight_rows, list(preflight_rows[0].keys()))
    effective = {
        "generated_at": utc_now(),
        "project_root": str(root),
        "source_config": str(config_path),
        "source_config_sha256": sha256_file(config_path),
        "runner_schema_version": RUNNER_SCHEMA_VERSION,
        "config": config,
        "environment": environment,
        "database_diagnostics": database_diagnostics,
        "manifest_sha256": sha256_file(Path(config["paths"]["manifest"])),
        "preflight_manifest_sha256": sha256_file(output / "preflight_manifest.csv"),
    }
    atomic_write_json(output / "effective_config.json", effective)
    atomic_write_json(output / "provenance.json", effective)
    atomic_write_text(
        output / "PREFLIGHT_COMPLETE.txt",
        f"VirulenceFinder preflight: PASS\nGenerated: {utc_now()}\nAssemblies: {len(rows)}\n",
    )


def determine_workers(config: dict[str, Any], environment: dict[str, Any], override: int | None) -> int:
    workers = int(override if override is not None else config["execution"]["workers"])
    if environment.get("heavy_process_matches"):
        workers = min(workers, int(config["execution"]["workers_when_heavy_job_detected"]))
    if workers < 1:
        raise ContractError("Worker count must be at least one.")
    return workers


def verify_primary_release(config: dict[str, Any]) -> None:
    release = Path(config["paths"]["primary_release_dir"])
    marker = release / "RUN_COMPLETE.txt"
    checks = release / "final_contract_checks.csv"
    statuses = release / "final_question_status.csv"
    if not marker.exists() or "PASS" not in marker.read_text(encoding="utf-8"):
        raise ContractError("Primary research-question release marker is missing or not PASS.")
    check_rows = read_csv(checks)
    if len(check_rows) != 18 or not all(parse_bool(row.get("pass")) for row in check_rows):
        raise ContractError("Primary release no longer has 18 passing final contract checks.")
    status_rows = read_csv(statuses)
    expected = {f"RQ{i:02d}" for i in range(1, 11)}
    observed = {row.get("research_question") for row in status_rows if row.get("status") == "complete"}
    if observed != expected:
        raise ContractError("Primary RQ01-RQ10 completion status changed.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        default="config/virulencefinder_sensitivity.toml",
        help="Commented TOML file containing all user-changeable and locked settings.",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--preflight-only", action="store_true", help="Validate inputs/runtime and stop.")
    mode.add_argument("--pilot-only", action="store_true", help="Run and reparse the six-job pilot only.")
    mode.add_argument("--combine-only", action="store_true", help="Rebuild combined tables from valid caches.")
    mode.add_argument("--verify-only", action="store_true", help="Validate existing caches and combined outputs.")
    parser.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Reuse only fully validated context-matched caches (default comes from TOML).",
    )
    parser.add_argument("--workers", type=int, help="Override execution.workers from the TOML file.")
    parser.add_argument(
        "--retry-failures",
        action="store_true",
        help="Resume the full task set and retry any job lacking a valid success cache.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    os.chdir(root)
    config_path = resolve_path(root, args.config)
    config = load_config(config_path, root)
    output = Path(config["paths"]["output"])
    output.mkdir(parents=True, exist_ok=True)
    verify_primary_release(config)
    environment = collect_environment(config, root)
    workers = determine_workers(config, environment, args.workers)
    rows = load_and_validate_manifest(config, root, workers)
    catalog, family_universe, db_diagnostics = build_database_catalog(config)
    write_preflight_outputs(config_path, config, environment, rows, db_diagnostics, root)
    print(
        f"Preflight PASS: {len(rows)} Longcycler assemblies; {len(family_universe)} gene families; "
        f"workers={workers}.",
        flush=True,
    )
    if args.preflight_only:
        return 0
    if args.pilot_only:
        run_pilot(config, environment, rows, workers)
        print("Pilot PASS.", flush=True)
        return 0

    resume = bool(config["execution"]["resume"]) if args.resume is None else bool(args.resume)
    if args.retry_failures:
        resume = True
    if args.combine_only or args.verify_only:
        run_results: list[dict[str, Any]] = []
        for row in rows:
            for profile in EXPECTED_PROFILE_NAMES:
                _, context_hash = build_context(config, environment, row, profile)
                valid, meta, hits, reason = validate_cached_job(config, row, profile, context_hash)
                run_results.append(
                    {
                        "manifest_index": row["manifest_index"],
                        "profile": profile,
                        "Assembly_ID": row["Assembly_ID"],
                        "episode_key": row["episode_key"],
                        "fasta_sha256": row["fasta_sha256"],
                        "context_hash": context_hash,
                        "cache_dir": str(job_paths(config, row, profile, context_hash)[0]),
                        "status": "success" if valid else "error",
                        "cache_reused": valid,
                        "attempt": meta.get("attempt", "") if meta else "",
                        "elapsed_seconds": meta.get("elapsed_seconds", "") if meta else "",
                        "n_hits": len(hits) if hits is not None else "",
                        "failure_reason": "" if valid else reason,
                    }
                )
    else:
        run_results = run_jobs(config, environment, rows, workers, resume)

    run_fields = [
        "manifest_index", "profile", "Assembly_ID", "episode_key", "fasta_sha256",
        "context_hash", "cache_dir", "status", "cache_reused", "attempt",
        "elapsed_seconds", "n_hits", "failure_reason",
    ]
    atomic_write_csv(output / "run_manifest.csv", run_results, run_fields)
    failures = [result for result in run_results if result["status"] != "success"]
    if failures:
        raise ContractError(
            f"{len(failures)} VirulenceFinder job(s) remain failed; combined/RQ outputs are blocked."
        )
    summary = combine_outputs(config, environment, rows, catalog, family_universe, run_results)
    provenance_path = output / "provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["batch_completed_at"] = utc_now()
    provenance["batch_summary"] = summary
    provenance["run_manifest_sha256"] = sha256_file(output / "run_manifest.csv")
    provenance["hits_long_sha256"] = sha256_file(output / "hits_long.csv")
    atomic_write_json(provenance_path, provenance)
    print("VirulenceFinder batch and combined outputs: PASS", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ContractError as exc:
        print(f"BLOCKED: {exc}", file=sys.stderr, flush=True)
        raise SystemExit(2)
