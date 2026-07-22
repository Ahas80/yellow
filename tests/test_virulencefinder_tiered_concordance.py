#!/usr/bin/env python3

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_virulencefinder_tiered_concordance.py"
SPEC = importlib.util.spec_from_file_location("tiered_concordance", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


class TieredConcordanceTests(unittest.TestCase):
    def test_binary_metrics_are_symmetric_and_reconciled(self):
        result = MODULE.binary_metrics([1, 1, 0, 0], [1, 0, 1, 0])
        self.assertEqual(result["both_present_n"], 1)
        self.assertEqual(result["both_absent_n"], 1)
        self.assertEqual(result["primary_only_n"], 1)
        self.assertEqual(result["cge_only_n"], 1)
        self.assertAlmostEqual(result["positive_jaccard"], 1 / 3)
        self.assertAlmostEqual(result["positive_agreement_dice"], 0.5)
        self.assertAlmostEqual(result["overall_agreement"], 0.5)
        self.assertAlmostEqual(result["cohen_kappa"], 0.0)

    def test_full_build_preserves_na_and_contract_counts(self):
        primary_hash_before = MODULE.sha256(MODULE.PRIMARY_MATRIX)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "tiered"
            MODULE.build(MODULE.DEFAULT_CONFIG, output)
            with (output / "validation_checks.csv").open(newline="", encoding="utf-8") as handle:
                checks = list(csv.DictReader(handle))
            self.assertTrue(checks)
            self.assertTrue(all(row["status"] == "PASS" for row in checks))

            with (output / "primary_feature_mapping_catalog.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                primary = list(csv.DictReader(handle))
            with (output / "cge_family_mapping_catalog.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                cge = list(csv.DictReader(handle))
            self.assertEqual(len(primary), 227)
            self.assertEqual(len(cge), 681)
            self.assertEqual(sum(row["comparison_state"] == "strict_direct" for row in primary), 48)
            self.assertEqual(sum(row["comparison_state"] == "strict_direct" for row in cge), 48)
            self.assertTrue(
                all(
                    row["cge_counterpart"] == "NA"
                    for row in primary
                    if row["comparison_state"] == "not_represented_or_unresolved"
                )
            )
            self.assertTrue(
                all(
                    row["primary_counterpart"] == "NA"
                    for row in cge
                    if row["comparison_state"] == "not_represented_or_unresolved"
                )
            )
        self.assertEqual(primary_hash_before, MODULE.sha256(MODULE.PRIMARY_MATRIX))


if __name__ == "__main__":
    unittest.main()
