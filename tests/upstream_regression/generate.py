#!/usr/bin/env python3
"""Regenerate fixtures with a separate, unmodified upstream v1.1.1 build."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile

from common import (HERE, RAW_INPUTS, RESOLVED_INDICES, UPSTREAM_COMMIT,
                    UPSTREAM_TREE, cases, finite_json_records, run_public)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--upstream-repo", type=Path, required=True,
                        help="Local Git clone containing frozen v1.1.1; no network is used")
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--work-dir", type=Path,
                        help="New or empty build directory (defaults to a retained temporary directory)")
    parser.add_argument("--cxx", default="c++")
    args = parser.parse_args()
    work = args.work_dir or Path(tempfile.mkdtemp(prefix="griffin-v111-reference-"))
    work.mkdir(parents=True, exist_ok=True)
    if any(work.iterdir()):
        parser.error("--work-dir must be empty")
    source, build = work / "source", work / "build"
    source.mkdir()
    tree = subprocess.check_output(["git", "-C", str(args.upstream_repo), "rev-parse",
                                    UPSTREAM_COMMIT + "^{tree}"], text=True).strip()
    if tree != UPSTREAM_TREE:
        raise RuntimeError("unexpected upstream source tree")
    archive = subprocess.check_output(["git", "-C", str(args.upstream_repo), "archive",
                                       "--format=tar", UPSTREAM_COMMIT])
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
        for member in tar.getmembers():
            path = Path(member.name)
            if path.is_absolute() or ".." in path.parts or not (member.isfile() or member.isdir()):
                raise RuntimeError("unsupported archive member: " + member.name)
        tar.extractall(source)
    source_hashes = {p.relative_to(source).as_posix(): sha(p.read_bytes())
                     for p in sorted(source.rglob("*")) if p.is_file()}
    source_digest = sha(json.dumps(source_hashes, sort_keys=True, separators=(",", ":")).encode())
    compiler = shutil.which(args.cxx)
    if compiler is None:
        raise RuntimeError("C++ compiler not found: " + args.cxx)
    driver = work / "upstream_reference_driver"
    with (work / "build.log").open("w") as log:
        subprocess.run(["cmake", "-S", str(source), "-B", str(build),
                        "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_CXX_COMPILER=" + compiler],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run(["cmake", "--build", str(build), "-j", "4"],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
        subprocess.run([compiler, "-std=c++11", "-O2", "-I", str(source),
                        str(HERE / "reference_driver.cc"), str(build / "libgriffin.a"),
                        "-o", str(driver)], stdout=log, stderr=subprocess.STDOUT, check=True)
    scan = cases()
    actual = run_public(args.runner, scan, RAW_INPUTS)
    reference_input = ""
    for case in scan:
        row = actual[case["id"]]
        reference_input += ",".join(str(v) for v in [case["id"], case["final"],
            row["s_GeV2"], case["cos_theta"], *row["sm_inputs"]]) + "\n"
    (work / "reference_input.csv").write_text(reference_input)
    result = subprocess.run([str(driver)], input=reference_input, text=True,
                            capture_output=True, check=True, timeout=180)
    (work / "reference_output.jsonl").write_text(result.stdout)
    expected = finite_json_records(result.stdout, {case["id"] for case in scan})
    fixtures = []
    for case in scan:
        row = actual[case["id"]]
        reference = expected[case["id"]]
        fixtures.append(dict(case, s_GeV2=row["s_GeV2"], resolved_inputs=row["sm_inputs"],
                             expected={"amplitudes": reference["amplitudes"],
                                       "sigma_pb": reference["sigma_pb"]}))
    fixture_bytes = ("\n".join(json.dumps(row, separators=(",", ":"), allow_nan=False)
                               for row in fixtures) + "\n").encode()
    manifest = {
        "schema_version": 1,
        "upstream_repository": "https://github.com/lisongc/GRIFFIN",
        "upstream_tag": "v1.1.1", "upstream_commit": UPSTREAM_COMMIT,
        "upstream_tree": tree, "upstream_archive_sha256": sha(archive),
        "upstream_source_files_sha256": source_digest,
        "fixtures_sha256": sha(fixture_bytes),
        "reference_driver_sha256": sha((HERE / "reference_driver.cc").read_bytes()),
        "generator_sha256": sha(Path(__file__).read_bytes()),
        "public_runner_source_sha256": sha((HERE.parent / "upstream_regression_runner.cc").read_bytes()),
        "compiler": subprocess.check_output([compiler, "--version"], text=True).splitlines()[0],
        "reference_library_build": "CMake Release; upstream sources unmodified",
        "reference_driver_flags": "-std=c++11 -O2",
        "raw_inputs": RAW_INPUTS, "resolved_input_indices": RESOLVED_INDICES,
        "amplitude_order": ["VV", "AV", "VA", "AA"],
        "cross_section_unit": "pb", "cross_section_conversion": .38937966e9,
        "strict_cases": sum(c["group"] == "strict" for c in scan),
        "diagnostic_cases": sum(c["group"] != "strict" for c in scan),
        "tolerances": {"amplitude": {"absolute": 2e-15, "relative": 2e-11},
                       "sigma_pb": {"absolute": 1e-9, "relative": 2e-11}},
        "independence": "Amplitudes and angular squares come only from a separately linked frozen upstream library. Resolved zero-Wilson inputs, including the Gmu W mass, are supplied by the public calculator and frozen; its mass solver is not independently validated here.",
        "near_pole_policy": "Offsets +/-0.01 and +/-0.03 GeV^2 lie inside v1.2's changed one-loop interpolation. Their external-reference differences are diagnostic, not equality assertions or physics validation. Strict tolerances are unchanged for the 560 off-resonance cases.",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "fixtures.jsonl").write_bytes(fixture_bytes)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"work_directory": str(work), "strict_cases": manifest["strict_cases"],
                      "diagnostic_cases": manifest["diagnostic_cases"],
                      "fixtures_sha256": manifest["fixtures_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
