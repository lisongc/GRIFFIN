"""Input protocol and scan definition for the independent v1.1.1 reference."""
import csv
import io
import json
import math
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
UPSTREAM_COMMIT = "a18f53830283a78155750a23a3594ada8de4c074"
UPSTREAM_TREE = "55144aa54da8005033f87aacba95b60d54ce542e"
RESOLVED_INDICES = list(range(15)) + [16, 17, 18, 20, 21, 22, 23]
HEADER = "id,scheme,final,mode,coordinate,cos_theta,alpha,MZ_running,MW_running,GammaZ_running,GammaW_running,MH,MT,MB,MC,MS,MU,MD,ML,MM,ME,DeltaAlpha,alphaS,Gmu"
RAW_INPUTS = dict(zip(HEADER.split(",")[6:], [
    1/137.03599976, 91.1876, 80.377, 2.4952, 2.085, 125.1, 172.5,
    2.87, 1.27, .342, .246, .246, 1.777, .1057, .000511, .059, .1179, 1.166379e-5,
]))
FINALS = {"mu": 13, "tau": 15, "u": 2, "d": 1, "s": 3, "c": 4, "b": 5}


def cases():
    result = []
    for scheme in ("alpha", "gmu"):
        for flavor, final in FINALS.items():
            for energy in (30, 50, 70, 90, 95, 120, 160, 190):
                for angle_index, angle in enumerate((-.8, -.3, 0, .5, .8)):
                    result.append(dict(id=f"{scheme}_{flavor}_E{energy}_A{angle_index}",
                                       scheme=scheme, final=final, mode="energy",
                                       coordinate=energy, cos_theta=angle, group="strict"))
            for offset_index, offset in enumerate((-.03, -.01, .01, .03)):
                result.append(dict(id=f"{scheme}_{flavor}_pole_D{offset_index}",
                                   scheme=scheme, final=final, mode="delta_s",
                                   coordinate=offset, cos_theta=.5, group="near_pole_diagnostic"))
    return result


def csv_input(scan, raw_inputs):
    output = io.StringIO()
    output.write(HEADER + "\n")
    writer = csv.writer(output, lineterminator="\n")
    for case in scan:
        writer.writerow([case[k] for k in ("id", "scheme", "final", "mode", "coordinate", "cos_theta")]
                        + [raw_inputs[name] for name in HEADER.split(",")[6:]])
    return output.getvalue()


def finite_json_records(text, expected_ids):
    def finite(value):
        if isinstance(value, dict):
            return all(finite(v) for v in value.values())
        if isinstance(value, list):
            return all(finite(v) for v in value)
        return not isinstance(value, (int, float)) or (not isinstance(value, bool) and math.isfinite(value))
    rows = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        identity = row.get("id")
        if identity not in expected_ids or identity in rows or not finite(row):
            raise ValueError("unexpected, duplicate, or nonfinite record: " + str(identity))
        rows[identity] = row
    if set(rows) != set(expected_ids):
        raise ValueError("missing output records")
    return rows


def run_public(runner, scan, raw_inputs):
    result = subprocess.run([str(Path(runner).resolve())], input=csv_input(scan, raw_inputs),
                            text=True, capture_output=True, timeout=180)
    if result.returncode:
        raise RuntimeError(f"public runner exit {result.returncode}: {result.stderr[-4000:]}")
    return finite_json_records(result.stdout, {case["id"] for case in scan})
