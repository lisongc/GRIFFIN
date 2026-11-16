<!--
File: README.md
Purpose: Describe archived regression inputs, comparison assemblies and rerun commands.
Authors (merged-version documentation): Jackson Wallace and Lisong Chen
Last updated: 2026-11-16
-->

# Finite-result regression: historical baseline and extended amplitude

The current `mat_SMaas` calculation includes selected fermionic mixed
electroweak–QCD terms in both the pole residue and off-resonance background,
for both full and zero-coefficient predictions. The fixtures preserve the
earlier `mat_SMNNLO` results and their existing canonical-strength mapping;
they are bundled in this local directory. The runner's
default mode uses `SMEFTCalculator` directly; no generic engine or
perturbative-series projection is involved. Current test-run results are recorded in
[REGRESSION_TEST_REPORT.md](../../docs/REGRESSION_TEST_REPORT.md).

These fixtures preserve **Jackson's full, finite predictions**, not their
linearization, a fitted derivative, or an independently reconstructed amplitude.
The golden sources and manifest are never modified to accommodate the added
SM terms. They describe the earlier `mat_SMNNLO` scattering setup, so the
current public cross sections must not be compared to them as though the
SM content had stayed fixed.

The extended regression separates three checks:

1. A test-only assembly using `mat_SMNNLO` is checked against all original
   archived assertions, retaining the historical Jackson comparison.
2. Public full/SM amplitudes, cross sections and finite shifts are
   checked against a separately assembled `mat_SMaas` calculation.
3. Quantities unaffected by the scattering-class change—masses, widths,
   FA/SW/FV and contact amplitudes—are checked for consistency between the
   public calculation and comparison assemblies. The historical standalone Jackson FV
   remains excluded; FV comparisons here use the merged FA/SW construction.

The test-only assemblies construct the pole inputs and scattering amplitude
separately from the public amplitude method, but use the same resolved
mass inputs and underlying SM kernels. Thus these tests validate assembly,
matching and the intended extension. They are not a new independent Jackson
benchmark or an independent derivation of the SM corrections. No test-only
mode is exposed by the library's public calculator.

The current mixed background includes the upstream v1.1.1 correction to
`xscaas.cc` and `aasoffz.grid`. These historical fixtures and their tolerances
were not changed for that correction. The separate
[`test_near_pole.cc`](../sm/test_near_pole.cc) checks fixed numbers from an
independent upstream v1.1.1 evaluation; it does not use the shared-kernel
assembly as its expected answer. Its near-pole coverage and the remaining
exact-point limitation are documented in the regression report above.

## Historical archive coverage

| Fixture group | Cases | Content |
| --- | ---: | --- |
| `individual` | 520 | 13 coefficients, 2 schemes, 4 energies, 5 strengths |
| `contact_mapping` | 80 | 40 original/relabelled contact-strength control pairs |
| `native_tables67` | 120 | Vertex/contact/combined benchmarks, full finite evaluation |
| `earlier_tables67` | 168 | Earlier finite benchmark values, including smaller strengths |
| `native_table4` | 40 | Paired finite electron FA/SW and W mass |
| `earlier_table4` | 15 | Original finite electron FA/SW at dimensionless C=0,+0.001,-0.001 |
| Total | 943 | |

The standalone historical FV formula is intentionally **not** a golden target:
the merged implementation reconstructs FV from FA/SW. No archived derivative,
odd projection, or alternative fixed-mass result is included. Outputs absent
from an old record are not invented or inferred. The original FA/SW complex
records are real; their real values are retained.

## Common inputs

All tests use massless unpolarized `e+ e- -> mu+ mu-`, with cross sections
`d sigma / d cos(theta)` in pb, `cos(theta)=0.5`. Cross-section cases use
`sqrt(s)=50,90,95,120 GeV`; EWPO-only cases use a dummy 90 GeV evaluation point.

| Input | Value |
| --- | --- |
| alpha inverse; Gmu | 137.03599976; 1.166379e-5 GeV^-2 |
| MZ; GammaZ | 91.1876; 2.4952 GeV |
| MW; GammaW (alpha inputs) | 80.377; 2.085 GeV |
| MH; mt; mb; mc; mtau | 125.1; 172.5; 2.87; 0.652; 1.777 GeV |
| Other fermion masses | 0 |
| DeltaAlpha; alpha_s | 0.059; 0.1179 |
| Lambda | 246 GeV |
| GeV^-2 to pb | 0.38937966e9 |

Mass/width numbers above are the original input numbers; archived outputs
separately label pole and running quantities. The Gmu finite solver predicts
the W mass and width. All unlisted Wilson coefficients are zero. Usual
individual strengths obey dimensionless `C=0.060516*t`, where `t` is the
coefficient in TeV^-2: `0,+0.001,-0.001,+0.01,-0.01`. Earlier cross-section
fixtures also contain `t=+/-0.0001,+/-0.00001`.

## Explicit contact convention mapping

The merged evaluator keeps the named input API. Cll current-exchange aliases
and Cee current-exchange/Fierz aliases share one total unique-operator
strength. The fixtures already translate Jackson's independently populated
raw slots to these canonical strengths:

```
L.Cll1122 = J.Cll1122 + J.Cll2211
L.Cll1221 = J.Cll1221 + J.Cll2112
L.Cee1122 = J.Cee1122 + J.Cee1221 + J.Cee2211 + J.Cee2112
L.Cle1122 = J.Cle1122
L.Cle2211 = J.Cle2211
```

This preserves both LL/RR contact strengths and the muon-decay contribution:

$$
\Delta_\mu=C_{\phi l,11}^{(3)}+C_{\phi l,22}^{(3)}
-\frac12C_{ll,L}^{1221},\qquad
-\frac12C_{ll,L}^{1221}
=-\frac12(C_{ll,J}^{1221}+C_{ll,J}^{2112}).
$$

Vertex coefficients are unchanged. Each fixture keeps
`original_raw_coefficients` alongside the translated `inputs.coefficients`.
It is an explicit test-input translation, **not** automatic addition by input
setters. The runner assigns each canonical value once. Alias setters use
latest-write-wins semantics, so applying raw pieces sequentially through
aliases, or zeroing an alias after setting the total, would corrupt the
import. Do not set all 81 raw tuples including zeros; initialize with a new
`Input` or `clearWilsonCoefficients()` and assign one aggregated value per
unique operator. There is no new public importer API.

For example, old independent `Cll1221=a`, `Cll2112=b` must be imported with
one `setCll(first,second,second,first,a+b)` call, not two alias assignments.
Cll `1122`/`2211` and `1221`/`2112` form separate groups; Cee
`1122`/`2211`/`1221`/`2112` form one group. More generally the supported
orbits are `{prst,stpr}` for Cll and `{prst,stpr,ptsr,srpt}` for Cee, for all
generation labels. No extra Hermiticity or other-coefficient relations are
introduced. The pre-existing reordered Ceq/Cqe alias is unaffected.

The historical multi-coefficient contact benchmark has four `Cll=c` slots,
four `Cee=c/3` slots, `Cle1122=c/2` and `Cle2211=c/4`. Its translated input is
`Cll1122=Cll1221=2*c`, `Cee1122=4*c/3`, with Cle unchanged. All 40 archived
original-vs-canonical native controls agree exactly in both translated input
and expected output; regeneration checks this before writing fixtures.

## Runner protocol

The executable reads CSV from stdin and emits one JSON object per line on
stdout. Diagnostics belong on stderr. Input header:

```
id,scheme,energy_GeV,cos_theta,Lambda_GeV,CphiD,CphiWB,Cphil1_11,Cphil1_22,Cphil3_11,Cphil3_22,Cphie_11,Cphie_22,Cll1122,Cll1221,Cee1122,Cle1122,Cle2211
```

`scheme` is `alpha` or `gmu`. Coefficients are dimensionless. Output must echo
the string `id` and provide the relevant numeric fields:

```
sigma_pb, sm_baseline_pb, delta_sigma_pb,
MW_pole_GeV, MW_running_GeV, GW_pole_GeV, GW_running_GeV, MZ_pole_GeV,
FA_e, FA_mu, SW_e, SW_mu, amplitudes, contact_amplitudes
```

Each amplitude array contains four `[real,imaginary]` pairs ordered
**VV, AV, VA, AA** (initial current changes first). Extra output fields such as
FV are not compared to the deliberately excluded old Jackson FV; they can
be checked between the extended public calculation and comparison assembly. Only fields
actually archived for each case are required by the historical comparison.
Missing cases, missing observables, duplicates, unexpected ids, nonfinite
numbers and malformed arrays fail the test.

From the project root, run the current extended regression:

```sh
python3 tests/native_regression/compare_aas.py --runner build/release/tests/native_fixture_runner --summary build/release/native_regression_summary.json
ctest --test-dir build/release --output-on-failure -R regression
```

`compare_aas.py` runs the default public evaluator and the test-only
`--reference-smnnlo` and `--reference-aas` runner modes. Optional
`--save-actual path.jsonl`, `--save-reference-smnnlo path.jsonl` and
`--save-reference-aas path.jsonl` retain those three outputs separately.
Here "reference" means the test-only comparison assembly, which evaluates
both nonzero and zero Wilson coefficients; it is not the public zero-coefficient
SM result. The latter is read through `smResult()`, `smInputValues()` and
`smAmplitudes(point)`.
The summary separates the original archived assertions, the extended
public-versus-low-level assembly comparison, and unchanged-output checks.
The mass/width/FA/SW/contact subset is also checked directly against the
historical archive. FV is checked between merged evaluations only.

The older archive verifier remains available for historical SMNNLO output:

```sh
python3 tests/native_regression/compare.py --emit-input
python3 tests/native_regression/compare.py --actual saved_reference_smnnlo.jsonl --summary historical_summary.json
python3 tests/native_regression/test_compare.py
python3 tests/native_regression/test_compare_aas.py
```

Do not run `compare.py --runner` on the default public evaluator expecting
old scattering numbers: the extra SM terms intentionally change them.
The old verifier's `--group individual` option (repeatable) restricts its
archive checks; its default is all 943 cases. Reports give counts, failure
details and maximum absolute/relative differences per observable, including
the worst tolerance fraction. The fixture SHA-256 is checked before use.
The drivers have only Python standard-library dependencies.

Numerical comparison is `abs(actual-expected) <= atol + rtol*abs(expected)`.
Cross sections and shifts have a 1e-9 pb absolute floor; amplitudes 2e-15;
contacts 2e-20; FA/SW 1e-12; MW 5e-10 GeV. Relative tolerances are 2e-11,
except masses and contacts where they are 2e-12. Full exact settings appear
in `compare.py`, the extended-field settings in `compare_aas.py`, and the
machine-readable reports. The **shift is checked
separately**, so agreement in a large full cross section cannot mask an
incorrect small coefficient-induced change.

## Rebuilding fixtures intentionally

```sh
python3 tests/native_regression/build_fixtures.py --archive-root /path/to/griffin_jakcson
```

This mechanically copies native stored values, with only nb-to-pb conversion,
pairing of FA/SW records, the explicit contact-input mapping, and addition of
the earlier stored native shift to its stored native zero. The manifest
contains every original relative source path and SHA-256; fixture creation
does not invoke either physics implementation. Normal regression runs do not
need those external archive directories.
