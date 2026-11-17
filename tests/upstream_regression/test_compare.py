#!/usr/bin/env python3
"""Negative checks that keep the independent regression sensitive to drift."""
import copy
import unittest

from common import RESOLVED_INDICES, finite_json_records
from compare import verify


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        inputs = [1.0] * len(RESOLVED_INDICES)
        amplitudes = [[1e-5, -2e-6] for _ in range(4)]
        self.fixture = dict(id="case", group="strict", s_GeV2=10000.0,
                            resolved_inputs=inputs,
                            expected=dict(amplitudes=amplitudes, sigma_pb=10.0))
        self.actual = {"case": dict(id="case", s_GeV2=10000.0,
                                   resolved_inputs=inputs[:], sm_inputs=inputs[:],
                                   total_amplitudes=copy.deepcopy(amplitudes),
                                   sm_amplitudes=copy.deepcopy(amplitudes),
                                   total_sigma_pb=10.0, sm_sigma_pb=10.0,
                                   model_correction_pb=0.0)}
        self.manifest = dict(strict_cases=1, diagnostic_cases=0,
                             near_pole_policy="diagnostic", independence="external")

    def test_identical_reference_passes(self):
        result = verify([self.fixture], self.actual, self.manifest)
        self.assertTrue(result["passed"])
        self.assertEqual(result["strict_reference_comparisons"], 18)

    def test_common_amplitude_drift_fails(self):
        # Both public views can agree with one another while disagreeing with
        # the independently frozen kernel; this must still fail.
        for field in ("total_amplitudes", "sm_amplitudes"):
            self.actual["case"][field][0][0] += 1e-10
        self.assertFalse(verify([self.fixture], self.actual, self.manifest)["passed"])

    def test_cross_section_and_mass_drift_fail(self):
        self.actual["case"]["sm_sigma_pb"] += 1e-6
        self.actual["case"]["resolved_inputs"][0] += 1e-6
        self.assertFalse(verify([self.fixture], self.actual, self.manifest)["passed"])

    def test_near_pole_difference_is_reported_separately(self):
        self.fixture["group"] = "near_pole_diagnostic"
        self.manifest.update(strict_cases=0, diagnostic_cases=1)
        for field in ("total_amplitudes", "sm_amplitudes"):
            self.actual["case"][field][0][0] += 1e-10
        result = verify([self.fixture], self.actual, self.manifest)
        self.assertTrue(result["passed"])
        self.assertGreater(result["field_statistics"]["near_pole_diagnostic_total_amplitudes"]["max_absolute_error"], 0)

    def test_missing_duplicate_and_nonfinite_rows_rejected(self):
        for text in ("", '{"id":"case","x":NaN}',
                     '{"id":"case"}\n{"id":"case"}'):
            with self.assertRaises(ValueError):
                finite_json_records(text, {"case"})

    def test_missing_amplitude_is_rejected(self):
        self.actual["case"]["total_amplitudes"].pop()
        with self.assertRaises(ValueError):
            verify([self.fixture], self.actual, self.manifest)


if __name__ == "__main__":
    unittest.main()
