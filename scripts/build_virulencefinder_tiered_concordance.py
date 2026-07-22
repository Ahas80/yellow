#!/usr/bin/env python3
"""Build a conservative, tiered CGE/VFDB comparison without forcing mappings.

The primary ABRicate/VFDB release is read-only.  This script writes a separate
exploratory comparison in which individual-gene statistics are limited to
sequence/annotation-supported direct pairs, broader relationships are reduced
to explicitly configured binary system endpoints, and unavailable mappings
remain NA rather than being converted to absence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import statistics
import tempfile
import tomllib
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "virulencefinder_tiered_crosswalk.toml"
DEFAULT_OUTPUT = (
    ROOT
    / "results"
    / "virulencefinder_cge_3_2_1"
    / "concordance"
    / "tiered_cross_database_comparison"
)
PRIMARY_MATRIX = ROOT / "results" / "research_questions" / "_inputs" / "vf_presence_absence_532.csv"
CGE_WEB_MATRIX = ROOT / "results" / "virulencefinder_cge_3_2_1" / "presence_absence_web_default_id90_cov60.csv"
CGE_MATCHED_MATRIX = ROOT / "results" / "virulencefinder_cge_3_2_1" / "presence_absence_matched_id80_cov80.csv"

PRIMARY_META = {"Participant_id", "tp_lab", "episode_key", "Assembly_ID", "fasta_sha256"}
CGE_META = PRIMARY_META | {"UTI_Status", "Event_type", "full_path"}
PROFILE_PATHS = {
    "web_default_id90_cov60": CGE_WEB_MATRIX,
    "matched_id80_cov80": CGE_MATCHED_MATRIX,
}


class ContractError(RuntimeError):
    """Raised when an input or mapping contract is not satisfied."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ContractError(f"CSV has no header: {path}")
        rows = list(reader)
        return rows, list(reader.fieldnames)


def binary(value: Any, label: str) -> int:
    try:
        numeric = float(value or 0)
    except (TypeError, ValueError) as exc:
        raise ContractError(f"Non-numeric binary value for {label}: {value!r}") from exc
    if numeric not in (0.0, 1.0):
        raise ContractError(f"Expected 0/1 for {label}, found {value!r}")
    return int(numeric)


def safe_ratio(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator else None


def binary_metrics(primary: Sequence[int], cge: Sequence[int]) -> dict[str, Any]:
    if len(primary) != len(cge) or not primary:
        raise ContractError("Binary comparison vectors must have equal non-zero length.")
    n11 = sum(a == 1 and b == 1 for a, b in zip(primary, cge))
    n00 = sum(a == 0 and b == 0 for a, b in zip(primary, cge))
    primary_only = sum(a == 1 and b == 0 for a, b in zip(primary, cge))
    cge_only = sum(a == 0 and b == 1 for a, b in zip(primary, cge))
    n = len(primary)
    agreement = (n11 + n00) / n
    primary_prevalence = (n11 + primary_only) / n
    cge_prevalence = (n11 + cge_only) / n
    expected = (
        primary_prevalence * cge_prevalence
        + (1 - primary_prevalence) * (1 - cge_prevalence)
    )
    kappa = safe_ratio(agreement - expected, 1 - expected)
    return {
        "n": n,
        "both_present_n": n11,
        "both_absent_n": n00,
        "primary_only_n": primary_only,
        "cge_only_n": cge_only,
        "primary_prevalence": primary_prevalence,
        "cge_prevalence": cge_prevalence,
        "cge_minus_primary_prevalence": cge_prevalence - primary_prevalence,
        "positive_jaccard": safe_ratio(n11, n11 + primary_only + cge_only),
        "positive_agreement_dice": safe_ratio(2 * n11, 2 * n11 + primary_only + cge_only),
        "negative_agreement": safe_ratio(2 * n00, 2 * n00 + primary_only + cge_only),
        "overall_agreement": agreement,
        "cohen_kappa": kappa,
    }


def rank(values: Sequence[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda pair: pair[1])
    result = [0.0] * len(values)
    start = 0
    while start < len(indexed):
        end = start + 1
        while end < len(indexed) and indexed[end][1] == indexed[start][1]:
            end += 1
        mean_rank = (start + 1 + end) / 2
        for position in range(start, end):
            result[indexed[position][0]] = mean_rank
        start = end
    return result


def pearson(first: Sequence[float], second: Sequence[float]) -> float | None:
    if len(first) != len(second) or not first:
        return None
    mean_first = statistics.fmean(first)
    mean_second = statistics.fmean(second)
    numerator = sum((a - mean_first) * (b - mean_second) for a, b in zip(first, second))
    denom_first = sum((a - mean_first) ** 2 for a in first)
    denom_second = sum((b - mean_second) ** 2 for b in second)
    return safe_ratio(numerator, math.sqrt(denom_first * denom_second))


def spearman(first: Sequence[float], second: Sequence[float]) -> float | None:
    return pearson(rank(first), rank(second))


def atomic_write_csv(path: Path, rows: list[dict[str, Any]], fields: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with tempfile.NamedTemporaryFile(
        mode="w", newline="", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        handle.write(text)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def load_inputs(config_path: Path) -> dict[str, Any]:
    with config_path.open("rb") as handle:
        config = tomllib.load(handle)
    primary_rows, primary_fields = read_csv(PRIMARY_MATRIX)
    primary_features = [field for field in primary_fields if field not in PRIMARY_META]
    profiles: dict[str, dict[str, Any]] = {}
    for profile, path in PROFILE_PATHS.items():
        rows, fields = read_csv(path)
        profiles[profile] = {
            "path": path,
            "rows": rows,
            "fields": fields,
            "families": [field for field in fields if field not in CGE_META],
        }
    expected_samples = int(config["contract"]["expected_samples"])
    if len(primary_rows) != expected_samples:
        raise ContractError(f"Primary matrix has {len(primary_rows)} rows; expected {expected_samples}.")
    if len(primary_features) != int(config["contract"]["expected_primary_features"]):
        raise ContractError("Primary feature universe changed.")
    primary_keys = [row["episode_key"] for row in primary_rows]
    if len(set(primary_keys)) != expected_samples:
        raise ContractError("Primary episode keys are missing or duplicated.")
    for profile, payload in profiles.items():
        if len(payload["rows"]) != expected_samples:
            raise ContractError(f"{profile} row count changed.")
        if len(payload["families"]) != int(config["contract"]["expected_cge_families"]):
            raise ContractError(f"{profile} CGE family universe changed.")
        if payload["families"] != profiles["web_default_id90_cov60"]["families"]:
            raise ContractError("CGE profile family columns differ or changed order.")
        by_key = {row["episode_key"]: row for row in payload["rows"]}
        if set(by_key) != set(primary_keys):
            raise ContractError(f"{profile} episode-key universe differs from primary.")
        payload["rows"] = [by_key[key] for key in primary_keys]
        for primary_row, cge_row in zip(primary_rows, payload["rows"]):
            for field in ("Assembly_ID", "fasta_sha256"):
                if primary_row[field] != cge_row[field]:
                    raise ContractError(f"{profile} {field} does not match primary for {primary_row['episode_key']}.")
    return {
        "config": config,
        "primary_rows": primary_rows,
        "primary_features": primary_features,
        "profiles": profiles,
        "cge_families": profiles["web_default_id90_cov60"]["families"],
    }


def expand_system_rules(config: dict[str, Any], primary: list[str], cge: list[str]) -> list[dict[str, Any]]:
    primary_set, cge_set = set(primary), set(cge)
    expanded: list[dict[str, Any]] = []
    for raw in config.get("system_rule", []):
        rule = dict(raw)
        cge_members = list(rule.get("cge", []))
        prefix = rule.get("cge_prefix")
        if prefix:
            cge_members = [family for family in cge if family.startswith(prefix)]
            expected = int(rule["expected_cge_members"])
            if len(cge_members) != expected:
                raise ContractError(
                    f"System {rule['id']} prefix {prefix!r} expanded to {len(cge_members)}, expected {expected}."
                )
        primary_members = list(rule["primary"])
        missing_primary = sorted(set(primary_members) - primary_set)
        missing_cge = sorted(set(cge_members) - cge_set)
        if missing_primary or missing_cge:
            raise ContractError(
                f"System {rule['id']} has unavailable members: primary={missing_primary}, CGE={missing_cge}."
            )
        if rule.get("mode") != "any" or not primary_members or not cge_members:
            raise ContractError(f"System {rule['id']} must use non-empty mode='any' members.")
        rule["primary_members"] = primary_members
        rule["cge_members"] = cge_members
        expanded.append(rule)
    return expanded


def mapping_definitions(data: dict[str, Any]) -> dict[str, Any]:
    config = data["config"]
    primary_features = data["primary_features"]
    cge_families = data["cge_families"]
    direct_exact = [dict(item) for item in config["direct_exact"]["mappings"]]
    direct_alias = [dict(item) for item in config.get("direct_alias", [])]
    high_risk = [dict(item) for item in config.get("high_risk_exact", [])]
    direct = [dict(item, mapping_type="direct_exact") for item in direct_exact] + [
        dict(item, mapping_type="direct_alias") for item in direct_alias
    ]
    primary_direct = [item["primary"] for item in direct]
    cge_direct = [item["cge"] for item in direct]
    primary_high = [item["primary"] for item in high_risk]
    cge_high = [item["cge"] for item in high_risk]
    if len(set(primary_direct)) != len(primary_direct) or len(set(cge_direct)) != len(cge_direct):
        raise ContractError("Strict direct mappings must be one-to-one and collision-free.")
    if set(primary_direct) & set(primary_high) or set(cge_direct) & set(cge_high):
        raise ContractError("Strict-direct and high-risk exact tiers overlap.")
    missing_primary = sorted((set(primary_direct) | set(primary_high)) - set(primary_features))
    missing_cge = sorted((set(cge_direct) | set(cge_high)) - set(cge_families))
    if missing_primary or missing_cge:
        raise ContractError(f"Configured direct mappings are unavailable: {missing_primary}; {missing_cge}")
    systems = expand_system_rules(config, primary_features, cge_families)
    system_primary = set().union(*(set(rule["primary_members"]) for rule in systems))
    system_cge = set().union(*(set(rule["cge_members"]) for rule in systems))
    any_primary = set(primary_direct) | set(primary_high) | system_primary
    any_cge = set(cge_direct) | set(cge_high) | system_cge
    contract = config["contract"]
    expected = {
        "direct_exact": len(direct_exact),
        "direct_alias": len(direct_alias),
        "high_risk_exact": len(high_risk),
        "system_rules": len(systems),
        "system_primary_memberships": sum(len(rule["primary_members"]) for rule in systems),
        "system_cge_memberships": sum(len(rule["cge_members"]) for rule in systems),
        "any_tier_primary": len(any_primary),
        "any_tier_cge": len(any_cge),
        "non_comparable_primary": len(set(primary_features) - any_primary),
        "non_comparable_cge": len(set(cge_families) - any_cge),
    }
    for name, actual in expected.items():
        configured = int(contract[f"expected_{name}"])
        if actual != configured:
            raise ContractError(f"{name} reconciled to {actual}; contract requires {configured}.")
    return {
        "direct": direct,
        "high_risk": high_risk,
        "systems": systems,
        "primary_direct": set(primary_direct),
        "cge_direct": set(cge_direct),
        "primary_high": set(primary_high),
        "cge_high": set(cge_high),
        "system_primary": system_primary,
        "system_cge": system_cge,
        "any_primary": any_primary,
        "any_cge": any_cge,
        "counts": expected,
    }


def make_catalogs(data: dict[str, Any], defs: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    config = data["config"]
    technical = config["direct_exact"].get("technical_flags", {})
    direct_by_primary = {item["primary"]: item for item in defs["direct"]}
    direct_by_cge = {item["cge"]: item for item in defs["direct"]}
    high_by_primary = {item["primary"]: item for item in defs["high_risk"]}
    high_by_cge = {item["cge"]: item for item in defs["high_risk"]}
    primary_systems: dict[str, list[str]] = defaultdict(list)
    cge_systems: dict[str, list[str]] = defaultdict(list)
    for rule in defs["systems"]:
        for member in rule["primary_members"]:
            primary_systems[member].append(rule["id"])
        for member in rule["cge_members"]:
            cge_systems[member].append(rule["id"])

    primary_catalog: list[dict[str, Any]] = []
    for feature in data["primary_features"]:
        if feature in direct_by_primary:
            item = direct_by_primary[feature]
            tier = item["mapping_type"]
            state = "strict_direct"
            counterpart = item["cge"]
            rationale = item.get("rationale", config["direct_exact"]["default_rationale"])
            flag = technical.get(counterpart, "")
        elif feature in high_by_primary:
            item = high_by_primary[feature]
            tier, state, counterpart = "high_risk_exact", "ambiguous_not_direct", item["cge"]
            rationale, flag = item["rationale"], ""
        elif feature in defs["system_primary"]:
            tier, state, counterpart = "system_only", "aggregate_only", "NA"
            rationale, flag = "Comparable only through the listed system-level OR rule(s).", ""
        else:
            tier, state, counterpart = "not_comparable", "not_represented_or_unresolved", "NA"
            rationale, flag = "No defensible configured CGE counterpart; retain as method-specific.", ""
        primary_catalog.append(
            {
                "primary_feature": feature,
                "exclusive_tier": tier,
                "comparison_state": state,
                "cge_counterpart": counterpart,
                "eligible_for_direct_concordance": state == "strict_direct",
                "system_rules": "|".join(sorted(primary_systems.get(feature, []))),
                "technical_flag": flag,
                "rationale": rationale,
            }
        )

    cge_catalog: list[dict[str, Any]] = []
    for family in data["cge_families"]:
        if family in direct_by_cge:
            item = direct_by_cge[family]
            tier = item["mapping_type"]
            state = "strict_direct"
            counterpart = item["primary"]
            rationale = item.get("rationale", config["direct_exact"]["default_rationale"])
            flag = technical.get(family, "")
        elif family in high_by_cge:
            item = high_by_cge[family]
            tier, state, counterpart = "high_risk_exact", "ambiguous_not_direct", item["primary"]
            rationale, flag = item["rationale"], ""
        elif family in defs["system_cge"]:
            tier, state, counterpart = "system_only", "aggregate_only", "NA"
            rationale, flag = "Comparable only through the listed system-level OR rule(s).", ""
        else:
            tier, state, counterpart = "not_comparable", "not_represented_or_unresolved", "NA"
            rationale, flag = "No defensible configured primary counterpart; retain as method-specific.", ""
        cge_catalog.append(
            {
                "cge_family": family,
                "exclusive_tier": tier,
                "comparison_state": state,
                "primary_counterpart": counterpart,
                "eligible_for_direct_concordance": state == "strict_direct",
                "system_rules": "|".join(sorted(cge_systems.get(family, []))),
                "technical_flag": flag,
                "rationale": rationale,
            }
        )
    return primary_catalog, cge_catalog


def calls(rows: Sequence[dict[str, str]], field: str) -> list[int]:
    return [binary(row[field], field) for row in rows]


def direct_outputs(data: dict[str, Any], defs: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    primary_rows = data["primary_rows"]
    technical = data["config"]["direct_exact"].get("technical_flags", {})
    summaries: list[dict[str, Any]] = []
    episode_rows: list[dict[str, Any]] = []
    profile_summaries: list[dict[str, Any]] = []
    for profile, payload in data["profiles"].items():
        cge_rows = payload["rows"]
        primary_matrix = [calls(primary_rows, item["primary"]) for item in defs["direct"]]
        cge_matrix = [calls(cge_rows, item["cge"]) for item in defs["direct"]]
        for item, primary_values, cge_values in zip(defs["direct"], primary_matrix, cge_matrix):
            summaries.append(
                {
                    "profile": profile,
                    "mapping_type": item["mapping_type"],
                    "primary_feature": item["primary"],
                    "cge_family": item["cge"],
                    "technical_flag": technical.get(item["cge"], ""),
                    **binary_metrics(primary_values, cge_values),
                }
            )
        primary_burdens: list[int] = []
        cge_burdens: list[int] = []
        episode_jaccards: list[float] = []
        exact_profiles = 0
        for index, primary_row in enumerate(primary_rows):
            primary_vector = [matrix[index] for matrix in primary_matrix]
            cge_vector = [matrix[index] for matrix in cge_matrix]
            intersection = sum(a == 1 and b == 1 for a, b in zip(primary_vector, cge_vector))
            union = sum(a == 1 or b == 1 for a, b in zip(primary_vector, cge_vector))
            jaccard = intersection / union if union else 1.0
            differences = sum(a != b for a, b in zip(primary_vector, cge_vector))
            exact_profiles += differences == 0
            primary_burden, cge_burden = sum(primary_vector), sum(cge_vector)
            primary_burdens.append(primary_burden)
            cge_burdens.append(cge_burden)
            episode_jaccards.append(jaccard)
            episode_rows.append(
                {
                    "profile": profile,
                    "episode_key": primary_row["episode_key"],
                    "Assembly_ID": primary_row["Assembly_ID"],
                    "fasta_sha256": primary_row["fasta_sha256"],
                    "direct_target_universe_n": len(defs["direct"]),
                    "primary_present_n": primary_burden,
                    "cge_present_n": cge_burden,
                    "profile_difference_n": differences,
                    "positive_jaccard": jaccard,
                    "exact_profile": differences == 0,
                }
            )
        profile_summaries.append(
            {
                "profile": profile,
                "n_assemblies": len(primary_rows),
                "direct_target_universe_n": len(defs["direct"]),
                "median_episode_jaccard": statistics.median(episode_jaccards),
                "mean_episode_jaccard": statistics.fmean(episode_jaccards),
                "exact_profile_n": exact_profiles,
                "exact_profile_proportion": exact_profiles / len(primary_rows),
                "restricted_burden_spearman": spearman(primary_burdens, cge_burdens),
                "median_primary_present_n": statistics.median(primary_burdens),
                "median_cge_present_n": statistics.median(cge_burdens),
            }
        )
    return summaries, episode_rows, profile_summaries


def high_risk_outputs(data: dict[str, Any], defs: dict[str, Any]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for profile, payload in data["profiles"].items():
        for item in defs["high_risk"]:
            result.append(
                {
                    "profile": profile,
                    "comparison_role": "diagnostic_only_not_direct_concordance",
                    "primary_feature": item["primary"],
                    "cge_family": item["cge"],
                    "recommended_system_rule": item.get("recommended_system_rule", ""),
                    "rationale": item["rationale"],
                    **binary_metrics(
                        calls(data["primary_rows"], item["primary"]),
                        calls(payload["rows"], item["cge"]),
                    ),
                }
            )
    return result


def any_call(row: dict[str, str], members: Iterable[str]) -> int:
    return int(any(binary(row[member], member) for member in members))


def system_outputs(data: dict[str, Any], defs: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summaries: list[dict[str, Any]] = []
    call_matrix: list[dict[str, Any]] = []
    primary_calls: dict[str, list[int]] = {}
    for rule in defs["systems"]:
        primary_calls[rule["id"]] = [any_call(row, rule["primary_members"]) for row in data["primary_rows"]]
    for index, primary_row in enumerate(data["primary_rows"]):
        out: dict[str, Any] = {
            "episode_key": primary_row["episode_key"],
            "Assembly_ID": primary_row["Assembly_ID"],
            "fasta_sha256": primary_row["fasta_sha256"],
        }
        for rule in defs["systems"]:
            out[f"primary__{rule['id']}"] = primary_calls[rule["id"]][index]
            for profile, payload in data["profiles"].items():
                out[f"{profile}__{rule['id']}"] = any_call(payload["rows"][index], rule["cge_members"])
        call_matrix.append(out)
    for rule in defs["systems"]:
        for profile, payload in data["profiles"].items():
            cge_values = [any_call(row, rule["cge_members"]) for row in payload["rows"]]
            summaries.append(
                {
                    "profile": profile,
                    "system_id": rule["id"],
                    "system_label": rule["label"],
                    "support": rule["support"],
                    "aggregation": "any_member_present",
                    "primary_member_n": len(rule["primary_members"]),
                    "cge_member_n": len(rule["cge_members"]),
                    "primary_members": "|".join(rule["primary_members"]),
                    "cge_members": "|".join(rule["cge_members"]),
                    "rationale": rule["rationale"],
                    **binary_metrics(primary_calls[rule["id"]], cge_values),
                }
            )
    return summaries, call_matrix


def scope_summary(data: dict[str, Any], primary_catalog: list[dict[str, Any]], cge_catalog: list[dict[str, Any]]) -> list[dict[str, Any]]:
    labels = {
        "direct_exact": "strict_direct",
        "direct_alias": "strict_direct",
        "high_risk_exact": "ambiguous_exact_name",
        "system_only": "aggregate_only",
        "not_comparable": "not_comparable",
    }
    rows: list[dict[str, Any]] = []
    for method, catalog, field in (
        ("primary_ABRicate_VFDB", primary_catalog, "primary_feature"),
        ("CGE_VirulenceFinder", cge_catalog, "cge_family"),
    ):
        counts: dict[str, int] = defaultdict(int)
        for item in catalog:
            counts[labels[item["exclusive_tier"]]] += 1
        total = len(catalog)
        for category in ("strict_direct", "ambiguous_exact_name", "aggregate_only", "not_comparable"):
            rows.append(
                {
                    "method": method,
                    "unit": "feature" if field == "primary_feature" else "gene_family",
                    "mapping_category": category,
                    "n": counts[category],
                    "proportion": counts[category] / total,
                    "total_universe_n": total,
                }
            )
    return rows


def non_comparable_prevalence(data: dict[str, Any], defs: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    n = len(data["primary_rows"])
    for feature in sorted(set(data["primary_features"]) - defs["any_primary"]):
        values = calls(data["primary_rows"], feature)
        rows.append(
            {
                "method": "primary_ABRicate_VFDB",
                "profile": "primary_id80_cov80",
                "target": feature,
                "target_type": "feature",
                "counterpart_state": "NA_not_comparable",
                "present_n": sum(values),
                "prevalence": sum(values) / n,
            }
        )
    for profile, payload in data["profiles"].items():
        for family in sorted(set(data["cge_families"]) - defs["any_cge"]):
            values = calls(payload["rows"], family)
            rows.append(
                {
                    "method": "CGE_VirulenceFinder",
                    "profile": profile,
                    "target": family,
                    "target_type": "gene_family",
                    "counterpart_state": "NA_not_comparable",
                    "present_n": sum(values),
                    "prevalence": sum(values) / n,
                }
            )
    return rows


def build(config_path: Path, output: Path) -> None:
    data = load_inputs(config_path)
    defs = mapping_definitions(data)
    primary_catalog, cge_catalog = make_catalogs(data, defs)
    direct_summary, direct_episode, direct_profile = direct_outputs(data, defs)
    high_risk = high_risk_outputs(data, defs)
    systems, system_calls = system_outputs(data, defs)
    scope = scope_summary(data, primary_catalog, cge_catalog)
    unavailable = non_comparable_prevalence(data, defs)

    checks = [
        ("primary_catalog_rows", len(primary_catalog), 227),
        ("cge_catalog_rows", len(cge_catalog), 681),
        ("strict_direct_pairs", len(defs["direct"]), 48),
        ("direct_summary_rows", len(direct_summary), 96),
        ("direct_episode_rows", len(direct_episode), 1064),
        ("direct_profile_rows", len(direct_profile), 2),
        ("high_risk_summary_rows", len(high_risk), 18),
        ("system_summary_rows", len(systems), 14),
        ("system_call_rows", len(system_calls), 532),
        ("scope_rows", len(scope), 8),
        ("not_comparable_rows", len(unavailable), 138 + 2 * 515),
    ]
    if any(actual != expected for _, actual, expected in checks):
        raise ContractError(f"Output reconciliation failed: {checks}")

    output.mkdir(parents=True, exist_ok=True)
    atomic_write_csv(output / "primary_feature_mapping_catalog.csv", primary_catalog)
    atomic_write_csv(output / "cge_family_mapping_catalog.csv", cge_catalog)
    atomic_write_csv(output / "strict_direct_gene_concordance.csv", direct_summary)
    atomic_write_csv(output / "strict_direct_episode_concordance.csv", direct_episode)
    atomic_write_csv(output / "strict_direct_profile_summary.csv", direct_profile)
    atomic_write_csv(output / "ambiguous_exact_name_diagnostics.csv", high_risk)
    atomic_write_csv(output / "system_level_concordance.csv", systems)
    atomic_write_csv(output / "system_level_call_matrix.csv", system_calls)
    atomic_write_csv(output / "mapping_scope_summary.csv", scope)
    atomic_write_csv(output / "not_comparable_prevalence.csv", unavailable)
    atomic_write_csv(
        output / "validation_checks.csv",
        [
            {"check": name, "expected": expected, "observed": actual, "status": "PASS"}
            for name, actual, expected in checks
        ],
    )

    provenance = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_role": "exploratory tiered cross-database concordance; primary VFDB release unchanged",
        "schema_version": data["config"]["contract"]["schema_version"],
        "config": str(config_path),
        "config_sha256": sha256(config_path),
        "inputs": {
            str(path): sha256(path)
            for path in (PRIMARY_MATRIX, CGE_WEB_MATRIX, CGE_MATCHED_MATRIX)
        },
        "counts": defs["counts"],
        "important_limit": "A complete feature-by-family mapping is not possible; unavailable targets remain NA.",
    }
    atomic_write_text(output / "provenance.json", json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    atomic_write_text(
        output / "README.md",
        "\n".join(
            [
                "# Tiered CGE/VFDB comparison",
                "",
                "This is an exploratory harmonization layer; it does not alter the primary ABRicate/VFDB release or the completed VirulenceFinder sensitivity analysis.",
                "",
                "- 48 sequence/annotation-supported direct pairs may be compared at the individual-gene level. They remain provisionally curated rather than a universal ontology.",
                "- Nine same-name pairs are explicitly ambiguous and are shown only as diagnostics.",
                "- Seven configured system rules compare binary any-member presence; component counts and database-wide burden are not compared.",
                "- 138/227 primary features and 515/681 CGE families have no defensible configured counterpart and remain NA, not absent.",
                "- Neither method is treated as a gold standard; `primary_only` and `cge_only` are neutral discordance labels.",
                "",
                "Edit the commented mapping contract in `config/virulencefinder_tiered_crosswalk.toml`, then rerun `python3 scripts/build_virulencefinder_tiered_concordance.py`.",
                "Generate the publication figure with `Rscript scripts/visualize_tiered_virulence_comparison.R`.",
                "",
                "Main figure outputs: `tiered_virulence_comparison.png` and `tiered_virulence_comparison.pdf`.",
                "",
            ]
        ),
    )
    print(f"Tiered CGE/VFDB concordance: PASS ({output})")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    build(arguments.config.resolve(), arguments.output.resolve())
