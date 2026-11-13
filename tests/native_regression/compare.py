#!/usr/bin/env python3
"""Execute a new native evaluator and compare against archived Jackson results."""
import argparse
import csv
import hashlib
import io
import json
import math
import subprocess
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
# Absolute + relative tolerances. Complex arrays are compared componentwise.
# The absolute shift floor also covers cancellation in sigma(C)-sigma(0).
TOLERANCES = {
    "sigma_pb": (1e-9, 2e-11),
    "sm_baseline_pb": (1e-9, 2e-11),
    "delta_sigma_pb": (1e-9, 2e-11),
    "amplitudes": (2e-15, 2e-11),
    "contact_amplitudes": (2e-20, 2e-12),
    "MW_pole_GeV": (5e-10, 2e-12),
    "MW_running_GeV": (5e-10, 2e-12),
    "MZ_pole_GeV": (5e-12, 2e-12),
    "GW_pole_GeV": (5e-11, 2e-11),
    "GW_running_GeV": (5e-11, 2e-11),
    "FA_e": (1e-12, 2e-11),
    "FA_mu": (1e-12, 2e-11),
    "SW_e": (1e-12, 2e-11),
    "SW_mu": (1e-12, 2e-11),
}


def input_csv(cases, coefficients):
    buf = io.StringIO()
    names = ["id", "scheme", "energy_GeV", "cos_theta", "Lambda_GeV"] + coefficients
    writer = csv.DictWriter(buf, fieldnames=names, lineterminator="\n")
    writer.writeheader()
    for case in cases:
        inputs = case["inputs"]
        row = {"id": case["id"], **{k: inputs[k] for k in names[1:5]}}
        row.update({k: inputs["coefficients"].get(k, 0.0) for k in coefficients})
        writer.writerow(row)
    return buf.getvalue()


def compare(cases, actual_rows, tolerances=None):
    """Independent comparison; fail on omissions, duplicates and nonfinite data."""
    tolerances = TOLERANCES if tolerances is None else tolerances
    errors, statistics = [], {}
    actual = {}
    for row in actual_rows:
        if not isinstance(row, dict) or not isinstance(row.get("id"), str):
            errors.append({"error": "Result must be an object with a string id"})
            continue
        if row["id"] in actual:
            errors.append({"id": row["id"], "error": "Duplicate result id"})
        actual[row["id"]] = row
    expected_ids = {c["id"] for c in cases}
    for extra in sorted(set(actual) - expected_ids):
        errors.append({"id": extra, "error": "Unexpected result id"})
    count = 0

    def scalar(case_id, field, path, expected, found):
        nonlocal count
        count += 1
        if isinstance(found, bool) or not isinstance(found, (int, float)) or not math.isfinite(found):
            errors.append({"id": case_id, "field": path, "error": "Nonfinite or nonnumeric result"})
            return
        assert math.isfinite(expected)
        atol, rtol = tolerances[field]
        error = abs(found - expected)
        limit = atol + rtol * abs(expected)
        relative = error / abs(expected) if expected != 0 else None
        stat = statistics.setdefault(field, {"count": 0, "max_absolute_error": 0.0,
                                             "max_relative_error_nonzero": 0.0,
                                             "max_tolerance_fraction": 0.0})
        stat["count"] += 1
        if error >= stat["max_absolute_error"]:
            stat.update(max_absolute_error=error, worst_absolute_id=case_id)
        if relative is not None and relative >= stat["max_relative_error_nonzero"]:
            stat.update(max_relative_error_nonzero=relative, worst_relative_id=case_id)
        if error / limit >= stat["max_tolerance_fraction"]:
            stat.update(max_tolerance_fraction=error / limit, worst_tolerance_id=case_id)
        if error > limit:
            errors.append({"id": case_id, "field": path, "expected": expected, "actual": found,
                           "absolute_error": error, "allowed_error": limit})

    def visit(case_id, field, path, expected, found):
        if isinstance(expected, list):
            if not isinstance(found, list) or len(found) != len(expected):
                errors.append({"id": case_id, "field": path, "error": "Array shape mismatch"})
                return
            for i, (x, y) in enumerate(zip(expected, found)):
                visit(case_id, field, f"{path}[{i}]", x, y)
        else:
            scalar(case_id, field, path, expected, found)

    for case in cases:
        row = actual.get(case["id"])
        if row is None:
            errors.append({"id": case["id"], "error": "Missing result"})
            continue
        for field, expected in case["expected"].items():
            if field not in row:
                errors.append({"id": case["id"], "field": field, "error": "Missing observable"})
            else:
                visit(case["id"], field, field, expected, row[field])
    return {
        "passed": not errors, "expected_cases": len(cases), "actual_cases": len(actual),
        "numeric_comparisons": count, "groups": dict(Counter(c["group"] for c in cases)),
        "tolerances": {k: {"absolute": a, "relative": r} for k, (a, r) in tolerances.items()},
        "field_statistics": statistics, "failure_count": len(errors), "failures": errors[:100],
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--runner", type=Path, help="New evaluator executable; stdin CSV, stdout JSONL")
    mode.add_argument("--actual", type=Path, help="Previously captured new-evaluator JSONL")
    mode.add_argument("--emit-input", action="store_true", help="Emit runner CSV only")
    ap.add_argument("--summary", type=Path, help="Also save machine-readable comparison report")
    ap.add_argument("--save-actual", type=Path, help="Save new-evaluator stdout for diagnostics")
    ap.add_argument("--group", action="append", help="Select fixture group(s); default all")
    ap.add_argument("--timeout", type=float, default=300, help="Evaluator time limit in seconds")
    args = ap.parse_args()
    manifest = json.loads((HERE / "manifest.json").read_text())
    raw = (HERE / "fixtures.jsonl").read_bytes()
    if hashlib.sha256(raw).hexdigest() != manifest["fixtures_sha256"]:
        ap.error("Fixture integrity hash mismatch; regenerate intentionally from archived sources")
    cases = [json.loads(line) for line in raw.decode().splitlines()]
    if args.group:
        unknown = set(args.group) - set(manifest["groups"])
        if unknown:
            ap.error("Unknown fixture group(s): " + ", ".join(sorted(unknown)))
        cases = [c for c in cases if c["group"] in args.group]
    csv_input = input_csv(cases, manifest["csv_coefficients"])
    if args.emit_input:
        sys.stdout.write(csv_input)
        return 0
    try:
        if args.runner:
            result = subprocess.run([str(args.runner.resolve())], input=csv_input,
                                    text=True, capture_output=True, timeout=args.timeout)
            if result.returncode:
                raise RuntimeError(f"Evaluator exit {result.returncode}: {result.stderr[-8000:]}")
            output = result.stdout
        else:
            output = args.actual.read_text()
        if args.save_actual:
            args.save_actual.write_text(output)
        actual = [json.loads(line) for line in output.splitlines() if line.strip()]
        report = compare(cases, actual)
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
        report = {"passed": False, "error": str(exc), "expected_cases": len(cases)}
    report["fixture_sha256"] = manifest["fixtures_sha256"]
    encoded = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.summary:
        args.summary.write_text(encoded)
    sys.stdout.write(encoded)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
