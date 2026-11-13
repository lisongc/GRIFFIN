#!/usr/bin/env python3
"""Unit-test comparison failure modes without a physics executable."""
import unittest
from compare import compare, input_csv


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.cases = [{"id": "sample", "group": "unit", "expected": {
            "sigma_pb": 10., "delta_sigma_pb": .01,
            "amplitudes": [[1e-5, 0.], [0., 0.], [0., 0.], [0., 0.]]}}]
        self.actual = [{"id": "sample", **self.cases[0]["expected"]}]

    def test_identical_passes(self):
        report = compare(self.cases, self.actual)
        self.assertTrue(report["passed"])
        self.assertEqual(report["numeric_comparisons"], 10)

    def test_wrong_shift_fails_despite_correct_full_cross_section(self):
        self.actual[0]["delta_sigma_pb"] += 1e-6
        self.assertFalse(compare(self.cases, self.actual)["passed"])

    def test_missing_result_fails(self):
        self.assertFalse(compare(self.cases, [])["passed"])

    def test_missing_observable_fails(self):
        del self.actual[0]["delta_sigma_pb"]
        self.assertFalse(compare(self.cases, self.actual)["passed"])

    def test_duplicate_id_fails(self):
        self.assertFalse(compare(self.cases, self.actual * 2)["passed"])

    def test_unexpected_id_fails(self):
        self.assertFalse(compare(self.cases, self.actual + [{"id": "unexpected"}])["passed"])

    def test_nonfinite_fails(self):
        self.actual[0]["sigma_pb"] = float("nan")
        self.assertFalse(compare(self.cases, self.actual)["passed"])

    def test_shape_mismatch_fails(self):
        self.actual[0]["amplitudes"] = [[1e-5, 0.]]
        self.assertFalse(compare(self.cases, self.actual)["passed"])

    def test_csv_zero_fills_unused_coefficients(self):
        case = {"id": "a", "inputs": {"scheme": "gmu", "energy_GeV": 50,
                "cos_theta": .5, "Lambda_GeV": 246, "coefficients": {"CphiD": .001}}}
        value = input_csv([case], ["CphiD", "CphiWB"])
        self.assertEqual(value.splitlines()[1], "a,gmu,50,0.5,246,0.001,0.0")


if __name__ == "__main__":
    unittest.main()
