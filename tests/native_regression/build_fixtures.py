#!/usr/bin/env python3
"""Mechanically copy archived native Jackson predictions; never re-evaluate them."""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
COEFFICIENTS = [
    "CphiD", "CphiWB", "Cphil1_11", "Cphil1_22", "Cphil3_11",
    "Cphil3_22", "Cphie_11", "Cphie_22", "Cll1122", "Cll1221",
    "Cee1122", "Cle1122", "Cle2211",
]
RAW_EXTRA = ["Cll2211", "Cll2112", "Cee1221", "Cee2211", "Cee2112"]
SOURCES = {
    "individual": "comparison_individual_20261026/jackson/raw.jsonl",
    "contact_mapping": "comparison_individual_20261026/jackson/mapping.jsonl",
    "native_tables67": "comparison_native_note_20261017/jackson/raw.jsonl",
    "earlier_tables67": "comparison_tables4_7_highsm_20261015/jackson/raw.jsonl",
    "native_table4": "comparison_native_note_20261017/ewpo/raw_jackson.jsonl",
    "earlier_table4": "comparison_tables4_7_highsm_20261015/ewpo/raw_jackson.csv",
}
COMMON_INPUTS = {
    "alpha_inverse": 137.03599976, "Gmu_GeVminus2": 1.166379e-5,
    "MZ_GeV": 91.1876, "GammaZ_GeV": 2.4952,
    "MW_GeV_alpha_input": 80.377, "GammaW_GeV_alpha_input": 2.085,
    "MH_GeV": 125.1, "mt_GeV": 172.5, "mb_GeV": 2.87,
    "mc_GeV": 0.652, "mtau_GeV": 1.777, "other_fermion_masses_GeV": 0,
    "DeltaAlpha": 0.059, "alpha_s": 0.1179,
    "process": "unpolarized massless e+ e- -> mu+ mu-",
    "cross_section": "d sigma / d cos(theta), pb",
    "GeVminus2_to_pb": 0.38937966e9,
}


def canonical(raw):
    """Translate physical contact strengths; not a modification of golden data."""
    out = {k: raw.get(k, 0.0) for k in COEFFICIENTS}
    out["Cll1122"] += raw.get("Cll2211", 0.0)
    out["Cll1221"] += raw.get("Cll2112", 0.0)
    out["Cee1122"] += sum(raw.get(k, 0.0) for k in RAW_EXTRA[2:])
    return {k: v for k, v in out.items() if v != 0.0}


def scenario_coefficients(scenario, strength):
    c = 0.060516 * strength
    out = {}
    if scenario in ("vertex", "combined"):
        out.update({k: c for k in COEFFICIENTS[:8]})
    if scenario in ("contact", "combined"):
        out.update({"Cll" + s: c for s in ("1122", "1221", "2211", "2112")})
        out.update({"Cee" + s: c / 3 for s in ("1122", "1221", "2211", "2112")})
        out.update(Cle1122=c / 2, Cle2211=c / 4)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archive-root", type=Path, required=True)
    args = ap.parse_args()
    cases, source_info = [], {}

    def load(group):
        p = args.archive_root / SOURCES[group]
        data = p.read_bytes()
        source_info[group] = {"path": SOURCES[group], "sha256": hashlib.sha256(data).hexdigest()}
        if p.suffix == ".csv":
            return list(csv.DictReader(data.decode().splitlines()))
        return [json.loads(line) for line in data.decode().splitlines() if line.strip()]

    def add(group, index, raw, expected, scheme="gmu", energy=90, angle=0.5, lam=246, **tags):
        cases.append({
            "id": f"{group}_{index:04d}", "group": group,
            "inputs": {"scheme": scheme, "energy_GeV": energy, "cos_theta": angle,
                       "Lambda_GeV": lam, "coefficients": canonical(raw)},
            "original_raw_coefficients": {k: v for k, v in raw.items() if v != 0.0},
            "expected": expected, "source": tags,
        })

    for group in ("individual", "contact_mapping", "native_tables67", "earlier_tables67"):
        for i, row in enumerate(load(group)):
            raw = row.get("raw_coefficients")
            if raw is None:
                raw = scenario_coefficients(row["scenario"], row["t"])
            expected = {}
            renames = {"mw_pole": "MW_pole_GeV", "mw_running": "MW_running_GeV",
                       "gw_pole": "GW_pole_GeV", "gw_running": "GW_running_GeV",
                       "mz_pole": "MZ_pole_GeV"}
            for key in ("sigma_pb", "sm_baseline_pb", "delta_sigma_pb", "MW_pole_GeV", "MW_running_GeV", "GW_pole_GeV",
                        "GW_running_GeV", "MZ_pole_GeV", "FA_e", "FA_mu",
                        "amplitudes", "contact_amplitudes"):
                if key in row:
                    expected[key] = row[key]
            for old, new in renames.items():
                if old in row:
                    expected[new] = row[old]
            if "sigma_nb" in row:
                expected["sigma_pb"] = 1000 * row["sigma_nb"]
                expected["sm_baseline_pb"] = 1000 * row["sm_baseline_nb"]
                expected["delta_sigma_pb"] = 1000 * row["delta_sigma_nb"]
            for key, val in row.get("form_factors", {}).items():
                assert val[1] == 0.0, "Nonreal FA/SW needs an explicit protocol extension"
                expected[key] = val[0]
            add(group, i, raw, expected, row["scheme"], row["energy_GeV"],
                row.get("cos_theta", .5), row.get("Lambda_GeV", 246),
                row=i + 1, coefficient=row.get("coefficient"),
                scenario=row.get("scenario"), population=row.get("population"), t=row["t"])

    # Native Table 4: combine the separately archived SW and FA records.
    paired = defaultdict(dict)
    for line, row in enumerate(load("native_table4"), 1):
        if row["observable"] != "FV":
            paired[(row["coefficient"], row["t"])][row["observable"]] = (line, row)
    for i, ((name, t), obs) in enumerate(paired.items()):
        assert set(obs) == {"SW", "FA"}
        row = obs["FA"][1]
        key = name.replace("Cll_", "Cll")
        add("native_table4", i, {key: row["C"]},
            {"FA_e": obs["FA"][1]["value"], "SW_e": obs["SW"][1]["value"],
             "MW_pole_GeV": row["MW_native_GeV"]},
            coefficient=name, t=t, rows=[obs[o][0] for o in ("SW", "FA")])

    # Earlier Table 4 stored finite +/- C=.001 directly, alongside derivatives.
    # Copy only those actual finite evaluations; never use the derivative columns.
    paired = defaultdict(dict)
    for line, row in enumerate(load("earlier_table4"), 2):
        if row["observable"] != "FV" and float(row["mw_step"]) == .001:
            paired[row["coefficient"]][row["observable"]] = (line, row)
    index = 0
    for name, obs in paired.items():
        assert set(obs) == {"SW", "FA"}
        for c, field in ((0., None), (.001, "native_plus_shift"), (-.001, "native_minus_shift")):
            expected = {}
            for observable, (_, row) in obs.items():
                expected[observable + "_e"] = float(row["native_zero"]) + (float(row[field]) if field else 0.)
            if field is None:
                expected["MW_pole_GeV"] = float(obs["FA"][1]["mw_reference"])
            add("earlier_table4", index, {name: c}, expected, coefficient=name,
                coefficient_dimensionless=c, rows=[obs[o][0] for o in ("SW", "FA")])
            index += 1

    # Archived mapping controls must match in amplitudes, masses and cross section.
    mapping = [c for c in cases if c["group"] == "contact_mapping"]
    for old, new in zip(mapping[::2], mapping[1::2]):
        assert old["inputs"] == new["inputs"], "Physical contact mapping mismatch"
        assert old["expected"] == new["expected"], "Native mapping control mismatch"
    counts = dict(Counter(c["group"] for c in cases))
    assert counts == {"individual": 520, "contact_mapping": 80, "native_tables67": 120,
                      "earlier_tables67": 168, "native_table4": 40, "earlier_table4": 15}, counts
    fixture = HERE / "fixtures.jsonl"
    fixture.write_text("".join(json.dumps(c, separators=(",", ":"), allow_nan=False) + "\n" for c in cases))
    manifest = {"format_version": 1, "cases": len(cases), "groups": counts,
                "common_inputs": COMMON_INPUTS, "csv_coefficients": COEFFICIENTS,
                "amplitude_order": ["VV", "AV", "VA", "AA"], "sources": source_info,
                "fixtures_sha256": hashlib.sha256(fixture.read_bytes()).hexdigest(),
                "excluded": "Historical standalone FV and all derivative/projection values.",
                "mapping_controls_exact_pairs": len(mapping) // 2}
    (HERE / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"cases": len(cases), "groups": counts, "mapping_controls_exact_pairs": 40}))


if __name__ == "__main__":
    main()
