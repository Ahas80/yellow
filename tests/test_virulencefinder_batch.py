"""Regression tests for the manifest-bound VirulenceFinder batch runner."""

from __future__ import annotations

import copy
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "virulencefinder_batch", ROOT / "scripts" / "run_virulencefinder_batch.py"
)
assert SPEC and SPEC.loader
VF = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VF)


def minimal_config(tmp: Path) -> dict:
    return {
        "paths": {
            "python": "/pinned/python",
            "blastn": "/pinned/blastn",
            "database": "/pinned/database",
            "output": str(tmp / "out"),
            "manifest": str(tmp / "manifest.csv"),
            "clinical_status": str(tmp / "status.csv"),
        },
        "software": {
            "virulencefinder_version": "3.2.1",
            "database_version": "2.1.0",
            "database_commit": "commit",
            "databases": ["virulence_ecoli", "stx"],
        },
        "cohort": {
            "expected_assemblies": 1,
            "expected_residents": 1,
            "expected_uti": 0,
            "expected_not_uti": 1,
            "required_assembler": "longcycler",
            "required_organism": "Escherichia coli",
            "forbidden_path_fragments": ["Rowenas analysis"],
        },
        "profiles": {
            "web_default_id90_cov60": {"minimum_identity": 0.90, "minimum_coverage": 0.60},
            "matched_id80_cov80": {"minimum_identity": 0.80, "minimum_coverage": 0.80},
        },
        "execution": {"overlap_nt": 30},
    }


def minimal_row(fasta: Path) -> dict:
    return {
        "manifest_index": 1,
        "Participant_id": "R1",
        "tp_lab": "T0",
        "episode_key": "R1||T0",
        "Assembly_ID": "A1",
        "full_path": str(fasta.resolve()),
        "fasta_sha256": VF.sha256_file(fasta),
    }


def result_document(row: dict, profile: str, regions: dict | None = None) -> dict:
    thresholds = {
        "web_default_id90_cov60": (0.90, 0.60),
        "matched_id80_cov80": (0.80, 0.80),
    }
    identity, coverage = thresholds[profile]
    return {
        "software_version": "3.2.1",
        "databases": {
            "db": {"database_version": "2.1.0", "database_commit": "commit"}
        },
        "software_executions": {
            "random-execution-key": {
                "parameters": {
                    "inputfasta": row["full_path"],
                    "inputfastq_1": None,
                    "aligner": "blastn",
                    "databases": ["virulence_ecoli", "stx"],
                    "min_id": identity,
                    "min_cov": coverage,
                    "overlap": 30,
                }
            }
        },
        "seq_regions": regions or {},
        "phenotypes": {"geneA": {"function": "test function"}},
    }


def region(
    *, key: str = "geneA:1:ACC", identity: float = 99.0,
    coverage: float = 100.0, ref_start: int = 1, ref_end: int = 100,
    query_start: int = 1, query_end: int = 100, query: str = "contig1",
) -> dict:
    return {
        "phenotypes": ["geneA"],
        "ref_database": ["VirulenceFinder-2.1.0:virulence_ecoli"],
        "ref_id": key,
        "name": key.split(":", 1)[0],
        "ref_acc": "ACC",
        "identity": identity,
        "alignment_length": abs(ref_end - ref_start) + 1,
        "ref_seq_length": 100,
        "ref_start_pos": ref_start,
        "ref_end_pos": ref_end,
        "query_id": query,
        "query_start_pos": query_start,
        "query_end_pos": query_end,
        "coverage": coverage,
        "grade": 2,
    }


class RunnerUnitTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.tmp = Path(self.temp.name)
        self.fasta = self.tmp / "assembly.fasta"
        self.fasta.write_text(">contig1\nACGTACGT\n", encoding="ascii")
        self.config = minimal_config(self.tmp)
        self.row = minimal_row(self.fasta)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_command_contains_exactly_one_fasta_argument(self) -> None:
        command = VF.build_command(
            self.config, self.row, "web_default_id90_cov60", self.tmp / "job", self.tmp / "job.json"
        )
        self.assertEqual(command.count("-ifa"), 1)
        self.assertEqual(command[command.index("-ifa") + 1], self.row["full_path"])

    def test_zero_hit_json_is_a_successful_parse(self) -> None:
        data = result_document(self.row, "web_default_id90_cov60")
        self.assertEqual(VF.canonical_hits_from_data(data, self.config, self.row, "web_default_id90_cov60"), [])

    def test_random_keys_duplicate_hits_and_multiple_loci_are_stable(self) -> None:
        a = region(query_start=10, query_end=109)
        regions = {
            "geneA-random-one": a,
            "geneA-random-two": copy.deepcopy(a),
            "geneA-second-locus": region(query_start=500, query_end=599),
        }
        data = result_document(self.row, "web_default_id90_cov60", regions)
        hits = VF.canonical_hits_from_data(data, self.config, self.row, "web_default_id90_cov60")
        self.assertEqual(len(hits), 2)
        collapsed = [hit for hit in hits if int(hit["query_start"]) == 10][0]
        self.assertEqual(collapsed["collapsed_duplicate_result_keys"], 1)
        self.assertEqual(collapsed["raw_result_keys"], "geneA-random-one;geneA-random-two")

    def test_split_reference_union_enforces_coverage_threshold(self) -> None:
        regions = {
            "fragment-z": region(coverage=40.0, ref_start=61, ref_end=100, query_start=300, query_end=339),
            "fragment-a": region(coverage=40.0, ref_start=1, ref_end=40, query_start=1, query_end=40),
        }
        data = result_document(self.row, "matched_id80_cov80", regions)
        hits = VF.canonical_hits_from_data(data, self.config, self.row, "matched_id80_cov80")
        self.assertEqual(len(hits), 2)
        self.assertTrue(all(hit["aggregate_reference_coverage_pct"] == 80.0 for hit in hits))
        self.assertTrue(all(hit["coverage_filter_basis"] == "split_reference_union" for hit in hits))

        data["seq_regions"]["fragment-z"].update(
            {"ref_start_pos": 81, "ref_end_pos": 100, "coverage": 20.0}
        )
        self.assertEqual(
            VF.canonical_hits_from_data(data, self.config, self.row, "matched_id80_cov80"), []
        )

    def test_below_identity_and_malformed_json_are_rejected(self) -> None:
        data = result_document(
            self.row, "web_default_id90_cov60", {"bad": region(identity=89.99)}
        )
        with self.assertRaises(VF.ContractError):
            VF.canonical_hits_from_data(data, self.config, self.row, "web_default_id90_cov60")
        broken = self.tmp / "broken.json"
        broken.write_text('{"truncated":', encoding="utf-8")
        with self.assertRaises(VF.ContractError):
            VF.validate_result_file(broken, self.config, self.row, "web_default_id90_cov60")

    def test_cache_context_changes_for_every_scientific_input(self) -> None:
        environment = {
            "virulencefinder_version": "3.2.1", "python_version": "3.14.0",
            "blast_version": "2.17.0", "database_version": "2.1.0",
            "database_commit": "commit", "database_config_sha256": "c",
            "virulence_ecoli_sha256": "e", "stx_sha256": "s",
            "virulencefinder_main_sha256": "v",
        }
        _, baseline = VF.build_context(
            self.config, environment, self.row, "web_default_id90_cov60"
        )
        mutations = []
        for path, value in [
            (("profiles", "web_default_id90_cov60", "minimum_identity"), 0.91),
            (("profiles", "web_default_id90_cov60", "minimum_coverage"), 0.61),
            (("execution", "overlap_nt"), 31),
            (("software", "databases"), ["virulence_ecoli"]),
        ]:
            cfg = copy.deepcopy(self.config)
            target = cfg
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            mutations.append(VF.build_context(cfg, environment, self.row, "web_default_id90_cov60")[1])
        changed_row = dict(self.row, fasta_sha256="f" * 64)
        mutations.append(VF.build_context(self.config, environment, changed_row, "web_default_id90_cov60")[1])
        changed_env = dict(environment, database_commit="new")
        mutations.append(VF.build_context(self.config, changed_env, self.row, "web_default_id90_cov60")[1])
        self.assertTrue(all(value != baseline for value in mutations))

    def test_fasta_parser_rejects_empty_records_duplicate_ids_and_bad_sequence(self) -> None:
        for contents in (
            "",
            ">one\nACGT\n>one\nACGT\n",
            ">one\n",
            ">one\nACGTZ\n",
            "ACGT\n",
        ):
            self.fasta.write_text(contents, encoding="ascii")
            with self.subTest(contents=contents), self.assertRaises(VF.ContractError):
                VF.validate_fasta_and_hash(self.fasta)


class ManifestContractTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.tmp = Path(self.temp.name)
        self.config = minimal_config(self.tmp)
        self.fasta = self.tmp / "assembly.fasta"
        self.fasta.write_text(">c1\nACGT\n", encoding="ascii")
        self.base = {
            "Assembly_ID": "A1", "Participant_id": "R1", "tp_lab": "T0",
            "assembler": "longcycler", "QC_PASS": "TRUE", "selected_canonical": "TRUE",
            "file_exists": "TRUE", "usable_fasta": "TRUE",
            "Clinical_Organism": "Escherichia coli", "full_path": str(self.fasta),
            "fasta_sha256": VF.sha256_file(self.fasta), "Event_type": "Routine",
        }
        with Path(self.config["paths"]["clinical_status"]).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["Participant_id", "tp_lab", "UTI_Status", "Event_type"])
            writer.writeheader()
            writer.writerow({"Participant_id": "R1", "tp_lab": "T0", "UTI_Status": "Not_UTI", "Event_type": "Routine"})

    def tearDown(self) -> None:
        self.temp.cleanup()

    def write_manifest(self, row: dict) -> None:
        with Path(self.config["paths"]["manifest"]).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)

    def test_valid_manifest_passes(self) -> None:
        self.write_manifest(self.base)
        result = VF.load_and_validate_manifest(self.config, self.tmp, 1)
        self.assertEqual(len(result), 1)

    def test_manifest_rejects_nonlongcycler_qc_noncanonical_missing_and_stale(self) -> None:
        variants = [
            {"assembler": "flye"},
            {"QC_PASS": "FALSE"},
            {"selected_canonical": "FALSE"},
            {"file_exists": "FALSE"},
            {"usable_fasta": "FALSE"},
            {"fasta_sha256": "0" * 64},
        ]
        for variant in variants:
            row = dict(self.base, **variant)
            self.write_manifest(row)
            with self.subTest(variant=variant), self.assertRaises(VF.ContractError):
                VF.load_and_validate_manifest(self.config, self.tmp, 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
