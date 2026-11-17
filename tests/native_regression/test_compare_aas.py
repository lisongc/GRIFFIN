#!/usr/bin/env python3
"""Purpose: test rejection of missing mixed corrections and malformed comparison data.

Merged-version authors: Jackson Wallace and Lisong Chen.
Last updated: 2026-11-16.
"""
from copy import deepcopy
import unittest

from compare_aas import FIELDS, verify


class ExtendedComparisonTests(unittest.TestCase):
    def setUp(self):
        old = {key: 1.0 for key in FIELDS}
        old.update(id="case", sigma_pb=10.0, sm_baseline_pb=9.0, delta_sigma_pb=1.0,
                   amplitudes=[[1e-4, 0.0] for _ in range(4)],
                   contact_amplitudes=[[1e-8, 0.0] for _ in range(4)])
        self.old = [old]
        self.aas = deepcopy(self.old)
        self.aas[0].update(sigma_pb=10.01, sm_baseline_pb=9.005, delta_sigma_pb=1.005)
        self.aas[0]["amplitudes"][0][1] = 1e-6
        self.public = deepcopy(self.aas)
        self.cases = [{"id": "case", "group": "unit", "expected": {
            key: old[key] for key in ("sigma_pb", "sm_baseline_pb", "delta_sigma_pb",
                                    "amplitudes", "contact_amplitudes", "MW_pole_GeV", "FA_e")}}]

    def report(self):
        return verify(self.cases, self.public, self.old, self.aas)

    def test_valid_extension_passes(self):
        report = self.report()
        self.assertTrue(report["passed"])
        self.assertEqual(report["cases_with_resolved_extension_signal"], 1)

    def test_old_public_amplitude_fails(self):
        self.public = deepcopy(self.old)
        self.assertFalse(self.report()["passed"])

    def test_unchanged_sm_result_fails(self):
        self.public[0]["sm_baseline_pb"] = self.old[0]["sm_baseline_pb"]
        self.assertFalse(self.report()["passed"])

    def test_wrong_finite_shift_fails(self):
        self.public[0]["delta_sigma_pb"] = self.old[0]["delta_sigma_pb"]
        self.assertFalse(self.report()["passed"])

    def test_corrupt_old_baseline_fails(self):
        self.old[0]["sigma_pb"] += 0.1
        self.assertFalse(self.report()["passed"])

    def test_common_mass_change_does_not_evade_archive(self):
        for rows in (self.public, self.old, self.aas):
            rows[0]["MW_pole_GeV"] += 0.1
        self.assertFalse(self.report()["passed"])

    def test_contact_change_fails(self):
        self.public[0]["contact_amplitudes"][0][0] *= 2.0
        self.assertFalse(self.report()["passed"])

    def test_missing_record_rejected(self):
        self.public = []
        with self.assertRaises(ValueError): self.report()

    def test_duplicate_reference_rejected(self):
        self.aas *= 2
        with self.assertRaises(ValueError): self.report()

    def test_missing_reference_field_rejected(self):
        del self.aas[0]["FV_e"]
        with self.assertRaises(ValueError): self.report()

    def test_nonfinite_reference_rejected(self):
        self.aas[0]["amplitudes"][0][0] = float("nan")
        with self.assertRaises(ValueError): self.report()

    def test_no_extension_in_either_path_fails(self):
        self.aas = deepcopy(self.old)
        self.public = deepcopy(self.old)
        self.assertFalse(self.report()["passed"])


if __name__ == "__main__":
    unittest.main()
