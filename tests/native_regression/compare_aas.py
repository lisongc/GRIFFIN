#!/usr/bin/env python3
"""Check the extended public API without rewriting the historical Jackson archive.

The two reference modes independently assemble low-level SM amplitudes, but
share the public finite input solution and existing physics kernels. This is
an assembly/regression check, not a new independent Jackson calculation of
the alpha_f*alpha_s extension.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

from compare import HERE, TOLERANCES, compare, input_csv

# Existing acceptance criteria are unchanged. FV was not a historical Jackson
# target; the new assembly comparison applies the existing FA/SW tolerance to it.
EXTENDED_TOLERANCES = dict(TOLERANCES, FV_e=(1e-12, 2e-11), FV_mu=(1e-12, 2e-11))
FIELDS = (
    "sigma_pb", "sm_baseline_pb", "delta_sigma_pb", "MW_pole_GeV",
    "MZ_pole_GeV", "MW_running_GeV", "GW_pole_GeV", "GW_running_GeV",
    "FA_e", "FA_mu", "SW_e", "SW_mu", "FV_e", "FV_mu",
    "amplitudes", "contact_amplitudes",
)
CHANGED_FIELDS = {"sigma_pb", "sm_baseline_pb", "delta_sigma_pb", "amplitudes"}


def indexed_complete(rows, ids, label):
    result = {}
    for row in rows:
        if not isinstance(row, dict) or set(row) != {"id", *FIELDS}:
            raise ValueError(f"{label}: incomplete or unexpected output fields")
        if row["id"] not in ids or row["id"] in result:
            raise ValueError(f"{label}: unexpected or duplicate result id")
        for field in FIELDS:
            value = row[field]
            if field in ("amplitudes", "contact_amplitudes"):
                if not isinstance(value, list) or len(value) != 4 or any(
                    not isinstance(pair, list) or len(pair) != 2 for pair in value
                ):
                    raise ValueError(f"{label}: malformed {field}")
                values = [part for pair in value for part in pair]
            else:
                values = [value]
            if any(isinstance(v, bool) or not isinstance(v, (int, float))
                   or not math.isfinite(v) for v in values):
                raise ValueError(f"{label}: nonfinite or nonnumeric {field}")
        result[row["id"]] = row
    if set(result) != ids:
        raise ValueError(f"{label}: missing output records")
    return result


def verify(cases, public_rows, old_rows, aas_rows):
    ids = {case["id"] for case in cases}
    public = indexed_complete(public_rows, ids, "public API")
    old = indexed_complete(old_rows, ids, "SMNNLO reference")
    aas = indexed_complete(aas_rows, ids, "SMaas reference")
    extended_cases = [dict(case, expected={key: aas[case["id"]][key] for key in FIELDS})
                      for case in cases]
    invariant_cases = [dict(case, expected={key: old[case["id"]][key]
                                           for key in FIELDS if key not in CHANGED_FIELDS})
                       for case in cases]
    archive_invariants = [dict(case, expected={key: value for key, value in case["expected"].items()
                                              if key not in CHANGED_FIELDS}) for case in cases]
    stages = {
        "historical_smnnlo_vs_jackson_archive": compare(cases, old_rows),
        "public_aas_vs_low_level_assembly": compare(extended_cases, public_rows, EXTENDED_TOLERANCES),
        "unchanged_mass_ewpo_contact_outputs": compare(invariant_cases, public_rows, EXTENDED_TOLERANCES),
        "unchanged_fields_vs_jackson_archive": compare(archive_invariants, public_rows),
    }
    changed = sum(any(abs(aas[identity]["amplitudes"][k][part]
                          - old[identity]["amplitudes"][k][part]) > 2e-15
                      for k in range(4) for part in range(2)) for identity in ids)
    return {
        "passed": all(stage["passed"] for stage in stages.values()) and changed > 0,
        "cases": len(cases),
        "historical_numeric_comparisons": stages["historical_smnnlo_vs_jackson_archive"]["numeric_comparisons"],
        "extended_numeric_comparisons": stages["public_aas_vs_low_level_assembly"]["numeric_comparisons"],
        "cases_with_resolved_extension_signal": changed,
        "stages": stages,
        "reference_scope": "Independent amplitude assembly using the existing SM kernels and the public finite mass solution; not new Jackson extended golden data.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--save-actual", type=Path, help="Save unmodified public-API JSONL")
    parser.add_argument("--save-reference-smnnlo", type=Path)
    parser.add_argument("--save-reference-aas", type=Path)
    parser.add_argument("--timeout", type=float, default=300, help="Seconds per runner invocation")
    args = parser.parse_args()
    manifest = json.loads((HERE / "manifest.json").read_text())
    fixture_bytes = (HERE / "fixtures.jsonl").read_bytes()
    if hashlib.sha256(fixture_bytes).hexdigest() != manifest["fixtures_sha256"]:
        parser.error("Historical fixture integrity hash mismatch")
    cases = [json.loads(line) for line in fixture_bytes.decode().splitlines()]
    csv_input = input_csv(cases, manifest["csv_coefficients"])

    def run(option, destination):
        command = [str(args.runner.resolve())] + ([option] if option else [])
        result = subprocess.run(command, input=csv_input, text=True, capture_output=True, timeout=args.timeout)
        if result.returncode:
            raise RuntimeError(f"Runner {option or 'public'} exit {result.returncode}: {result.stderr[-8000:]}")
        if destination:
            destination.write_text(result.stdout)
        return [json.loads(line) for line in result.stdout.splitlines() if line.strip()]

    try:
        public = run(None, args.save_actual)
        old = run("--reference-smnnlo", args.save_reference_smnnlo)
        aas = run("--reference-aas", args.save_reference_aas)
        report = verify(cases, public, old, aas)
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        report = {"passed": False, "error": str(error), "cases": len(cases)}
    report["fixture_sha256"] = manifest["fixtures_sha256"]
    encoded = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.summary:
        args.summary.write_text(encoded)
    sys.stdout.write(encoded)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
