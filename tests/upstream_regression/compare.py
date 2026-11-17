#!/usr/bin/env python3
"""Offline comparison with independently compiled upstream v1.1.1 fixtures."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

from common import HERE, RESOLVED_INDICES, UPSTREAM_COMMIT, UPSTREAM_TREE, run_public


def verify(fixtures, actual, manifest):
    failures = []
    statistics = {}
    counts = {"strict_cases": 0, "diagnostic_cases": 0, "strict_reference_comparisons": 0}

    def compare(label, identity, observed, expected, absolute, relative, assertion=True):
        if isinstance(observed, bool) or not isinstance(observed, (int, float)) or not math.isfinite(observed):
            raise ValueError(identity + ": nonnumeric or nonfinite " + label)
        difference = abs(observed-expected)
        fraction = difference/(absolute+relative*abs(expected))
        stat = statistics.setdefault(label, {"comparisons": 0, "max_absolute_error": 0.0,
                                            "max_tolerance_fraction": 0.0})
        stat["comparisons"] += 1
        if difference > stat["max_absolute_error"]:
            stat["max_absolute_error"] = difference
            stat["worst_absolute_id"] = identity
        stat["max_tolerance_fraction"] = max(stat["max_tolerance_fraction"], fraction)
        if assertion and difference > absolute+relative*abs(expected):
            failures.append(dict(id=identity, field=label, actual=observed, expected=expected,
                                 absolute_error=difference, tolerance_fraction=fraction))

    required = {"id", "s_GeV2", "resolved_inputs", "sm_inputs", "total_amplitudes",
                "sm_amplitudes", "total_sigma_pb", "sm_sigma_pb", "model_correction_pb"}
    for fixture in fixtures:
        identity = fixture["id"]
        row = actual[identity]
        if set(row) != required:
            raise ValueError(identity + ": unexpected public output fields")
        strict = fixture["group"] == "strict"
        counts["strict_cases" if strict else "diagnostic_cases"] += 1
        prefix = "strict" if strict else "near_pole_diagnostic"
        compare("s_GeV2", identity, row["s_GeV2"], fixture["s_GeV2"], 1e-9, 2e-12)
        for field in ("resolved_inputs", "sm_inputs"):
            if len(row[field]) != len(RESOLVED_INDICES):
                raise ValueError(identity + ": wrong resolved input count")
            for index, observed, expected in zip(RESOLVED_INDICES, row[field], fixture["resolved_inputs"]):
                compare("resolved_input_"+str(index), identity, observed, expected,
                        5e-10 if index in (0,20) else 5e-12, 2e-12)
        for field in ("total_amplitudes", "sm_amplitudes"):
            if len(row[field]) != 4 or any(len(pair) != 2 for pair in row[field]):
                raise ValueError(identity + ": wrong amplitude shape")
            for k in range(4):
                for part in range(2):
                    compare(prefix+"_"+field, identity, row[field][k][part],
                            fixture["expected"]["amplitudes"][k][part],
                            2e-15, 2e-11, assertion=strict)
                    if strict:
                        counts["strict_reference_comparisons"] += 1
        for field in ("total_sigma_pb", "sm_sigma_pb"):
            compare(prefix+"_"+field, identity, row[field], fixture["expected"]["sigma_pb"],
                    1e-9, 2e-11, assertion=strict)
            if strict:
                counts["strict_reference_comparisons"] += 1
        # These zero-Wilson identities remain assertions also in the diagnostic region.
        compare("zero_wilson_model_correction_pb", identity, row["model_correction_pb"], 0.0, 1e-15, 0.0)
        compare("zero_wilson_total_vs_sm_sigma_pb", identity, row["total_sigma_pb"], row["sm_sigma_pb"], 1e-9, 2e-11)
        for k in range(4):
            for part in range(2):
                compare("zero_wilson_total_vs_sm_amplitudes", identity,
                        row["total_amplitudes"][k][part], row["sm_amplitudes"][k][part], 2e-15, 2e-11)
    if counts["strict_cases"] != manifest["strict_cases"] or counts["diagnostic_cases"] != manifest["diagnostic_cases"]:
        raise ValueError("unexpected scan coverage")
    return dict(passed=not failures, **counts, failure_count=len(failures), failures=failures[:30],
                field_statistics=statistics, near_pole_policy=manifest["near_pole_policy"],
                independence=manifest["independence"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    try:
        manifest = json.loads((HERE / "manifest.json").read_text())
        data = (HERE / "fixtures.jsonl").read_bytes()
        if hashlib.sha256(data).hexdigest() != manifest["fixtures_sha256"]:
            raise ValueError("frozen fixture hash mismatch")
        if manifest["upstream_commit"] != UPSTREAM_COMMIT or manifest["upstream_tree"] != UPSTREAM_TREE:
            raise ValueError("unexpected reference revision")
        if hashlib.sha256((HERE / "reference_driver.cc").read_bytes()).hexdigest() != manifest["reference_driver_sha256"]:
            raise ValueError("reference driver changed since fixture generation")
        fixtures = [json.loads(line) for line in data.decode().splitlines() if line.strip()]
        if len({row["id"] for row in fixtures}) != len(fixtures):
            raise ValueError("duplicate fixture ids")
        actual = run_public(args.runner, fixtures, manifest["raw_inputs"])
        report = verify(fixtures, actual, manifest)
        report["upstream_commit"] = manifest["upstream_commit"]
        report["fixtures_sha256"] = manifest["fixtures_sha256"]
    except (OSError, ValueError, RuntimeError, KeyError, TypeError, subprocess.TimeoutExpired) as error:
        report = {"passed": False, "error": str(error)}
    encoded = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.summary:
        args.summary.write_text(encoded)
    sys.stdout.write(encoded)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
