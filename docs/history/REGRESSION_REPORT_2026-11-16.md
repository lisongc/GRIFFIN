<!--
Purpose: archive the 16 November 2026 verification, scope and reproduction.
Merged-version authors: Jackson Wallace and Lisong Chen.
Last updated: 2026-11-16.
-->

# Regression audit — 16 November 2026

**Recorded verification: 16 November 2026 — PASS, 15/15 CTests.** The fresh
`build/docs-cleanup` build passed in 23.67 seconds with both example options
disabled, including the scattering test moved into `tests/`. This layout
check supplements the physics-import verification below; it adds no new
physics claim. Sanitizers were not rerun for the upstream import or cleanup.

These are software/assembly and independent upstream-reference checks, not
validation of precision near-pole scattering. The calculator uses NNLO+ SM
form factors, `mat_SMaas`, Jackson's finite Gmu mass prescription and the full
amplitude square. Explicit dimension-six pole-input shifts are linear at
fixed mass; subsequent generated nonlinear dependence is retained. See the
[exact prescription](../USER_GUIDE.md#exact-prescription) and
[supported scope](../KNOWN_ISSUES.md#1-corrected-mixed-background-and-remaining-pole-restriction).

## Upstream import and recorded results

The imported upstream v1.2 commit is
`ca2d5234d7712063b192e079bd8e3a4346fe9a7e`; the local package version is 1.2.
It includes the v1.1.1 mixed self-energy correction, v1.2 one-loop near-pole
interpolation, separate `mat_SMNNLO2` class and grids, box-selection flags,
and `psobsfix::setval()`. Include paths/guards and trailing whitespace were
adapted, and the new example's uninitialized energy was fixed. Numerical
expressions and upstream disabled sectors were preserved.

`SMEFTCalculator` still uses `mat_SMaas`. Wilson matching, the finite Gmu mass
solution, FA/SW/FV insertion and contact normalization were not changed. The
updated SM background enters both full and zero-Wilson scattering results.

The original import Release build passed **15/15 CTests in 20.92 seconds**,
including **25 public headers** in installed-package checks. All eight
examples built at that stage exited successfully, with no reported nonfinite
numerical output. The existing 943-case regression passed all four stages;
historical fixture values and tolerances were unchanged.

| Comparison | Recorded coverage | Result |
|---|---|---|
| Historical `mat_SMNNLO` assembly versus Jackson archive | 943 cases; 22,475 numerical comparisons | PASS |
| Public calculator versus separate `mat_SMaas` assembly | 943 cases; 28,290 numerical comparisons | PASS |
| Independent upstream v1.1.1 reference | 560 strict cases; 10,080 comparisons | PASS |
| v1.1.1/v1.2 interpolation-region comparison | 56 diagnostic cases; equality not required | Finite/input/identity checks pass |
| Independent upstream v1.2 `mat_SMNNLO2` reference | 36 current-channel configurations; five complex outputs each | PASS |

Evidence: [import results](validation/results_2026-11-16_upstream_v12.json)
and [near-pole diagnostic output](validation/results_2026-11-16_v12.json).
Dated pre-import audits, hashes, timings and the 12 November amplitude-impact
study are retained in the [pre-import chronology](REGRESSION_LOG_2026-11-16_PRE_IMPORT.md).

## Independent upstream references

The [v1.1.1 regression](../../tests/upstream_regression/README.md) separately
compiles unmodified commit `a18f53830283a78155750a23a3594ada8de4c074`.
Its reference executable links only that library and evaluates `mat_SMaas`
plus an independent angular square. Checked-in fixtures make normal CTest
runs independent of Git, a reference checkout and network access.

The strict scan covers both electroweak schemes, electron initial states,
muon/tau/u/d/s/c/b final states, energies 30, 50, 70, 90, 95, 120, 160 and
190 GeV, and cosines -0.8, -0.3, 0, 0.5 and 0.8. All four complex amplitudes
in VV, AV, VA, AA order are checked through both public views, together with
both total/SM cross sections. There were **zero failures**; the maximum
amplitude-component difference was **4.07e-20** and the maximum cross-section
difference was **4.55e-13 pb**. Six harness checks detect drift or invalid output.

The reference receives the calculator's separately solved zero-Wilson SM
input values, frozen and checked on subsequent runs. Thus the scattering
kernels and angular square are independent, but this test does **not**
independently validate the Gmu mass solver. Existing mass/Jackson tests
retain that separate role. All Wilson coefficients are zero and Lambda is
1000 GeV; the manifest records full inputs, hashes, compiler and units.

The 56 diagnostics sample `s-MZ_pole^2 = +/-0.01, +/-0.03 GeV^2`, at cosine
0.5 for both schemes and all seven final flavors. These points lie inside
v1.2's changed one-loop interpolation, so cross-version equality is not an
assertion. Maximum differences were **2.67e-8** in an amplitude component
and **0.04658 pb** in a differential cross section. Finite outputs, frozen
resolved inputs and zero-Wilson total/SM identities remain assertions;
strict tolerances were not widened to hide the diagnostic differences.

The [v1.2 reference](../../tests/upstream_v12_reference/README.md) was generated
by the unmodified imported commit, using full-precision output. Its 36
configurations span muon/up/down final states, three energy/angle pairs and
four currents. Each checks the new total, a separate `mat_SMaas` total,
`resoffZ2f()`, `resoffZfb()` and `resvpol()`. `test_sm_nnlo2` also checks
total/background decomposition, inherited residue, typed constructors,
default box flags and borrowed fixed-observable updates. This verifies the
port of active upstream sectors, not complete off-resonance NNLO accuracy.

## Coverage and interpretation

- **Inputs/contacts:** coefficient and scale validation; Cll/Cee aliases,
  overwrite/clear and raw-to-canonical import; no double counting; Gmu factor.
- **Finite prescription:** direct vertices, Gmu equation/residual, independent
  zero-coefficient mass solution, FA/SW-based FV, full angular square,
  contact-quadratic terms and scale rescaling.
- **Extended amplitude:** 980 points across two schemes, seven final flavors,
  five energies, two angles and seven coefficient configurations; all four
  channels, full/SM results, residue/background additions and vanishing added
  mixed terms at zero strong coupling.
- **Public API/install:** snapshot ownership, copy/move and temporary lifetimes,
  rejected domains/processes, 25 standalone/double-included public headers,
  isolated installed clients, private-helper compile probes and absence of
  retired interfaces. Clients cover the example, public API, conversion,
  shared setup and separate EWPO/scattering headers.
- **Harness:** nine original and twelve extended self-tests, plus the six
  independent-reference checks, detect omissions, corrupted values, incorrect
  small differences and malformed/missing/nonfinite output.

The 943-case archive covers coefficient scans, contact controls and earlier
Tables 4/6/7 comparisons. Its common process is massless, unpolarized
$e^+e^-\to\mu^+\mu^-$, $\cos\theta=0.5$, at $\sqrt{s}=50,90,95,120$ GeV,
with differential cross sections in pb. See the [fixture README](../../tests/native_regression/README.md)
and [manifest](../../tests/native_regression/manifest.json) for complete inputs.

The historical Jackson cross sections use `mat_SMNNLO`; the added
`mat_SMaas` terms intentionally change them. Separate test-only assemblies
share the SM kernels and resolved mass inputs: they test assembly and
preserved matching, not an independent kernel derivation. Historical Jackson
standalone FV is not a golden target; merged FV is constructed from FA/SW
and the existing SM imaginary-ratio term. All comparisons use finite results,
not a linear extraction.

The corrected default background gives Im Q(0.01) about **-5.45e-14**, versus
about **6.74e-4** in the historical v1.1 diagnostic. This changes the old
observed behavior; it does not establish interpolation-boundary derivatives
or precision Z-lineshape accuracy. The exact real pole remains unsupported.

The separate `near_pole_nnlo2` diagnostic shows growth in the new box sector:
for electron annihilation to muons at cosine 0.5, offsets -0.1 and -0.01 GeV^2
give 0.96130 and 1.21544 nb with `mat_SMNNLO2`, versus 0.95498 and 0.95497 nb
with corrected `mat_SMaas`. This class remains separate from the SMEFT
calculator. Tests do not require the problematic behavior to persist.
Disabled sectors, unsupported flavors, positive-mass requirements and grid
extrapolation limitations are documented in [known issues](../KNOWN_ISSUES.md).

## Numerical tolerances

Each scalar satisfies `abs(actual-expected) <= atol + rtol*abs(expected)`.
The historical/assembly and v1.1.1 comparisons use:

| Quantity | Absolute tolerance | Relative tolerance |
|---|---:|---:|
| Full/SM/difference cross sections | $10^{-9}$ pb | $2\times10^{-11}$ |
| Amplitude real/imaginary components | $2\times10^{-15}$ | $2\times10^{-11}$ |
| Contact components | $2\times10^{-20}$ | $2\times10^{-12}$ |
| Pole/running W mass | $5\times10^{-10}$ GeV | $2\times10^{-12}$ |
| Pole Z mass | $5\times10^{-12}$ GeV | $2\times10^{-12}$ |
| W widths | $5\times10^{-11}$ GeV | $2\times10^{-11}$ |
| FA/SW; FV in extended assembly | $10^{-12}$ | $2\times10^{-11}$ |

The v1.2 frozen amplitudes use `atol=2e-17`, `rtol=2e-10`, with tighter
same-build decomposition checks. Missing/duplicate identifiers, missing
observables, malformed arrays and nonfinite values are rejected. Do not
regenerate fixtures merely to pass a failing test; review source provenance
and mismatches first. Reference READMEs above retain generation commands
and hashes; generation is an explicit maintenance action outside CTest.

## Reproduce the checks

Recorded environment: macOS arm64, AppleClang 16.0.0.16000026, C++11.
Python 3's standard library suffices. From the repository root:

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DGRIFFIN_BUILD_EXAMPLES=ON
cmake --build build/release -j 4
ctest --test-dir build/release --output-on-failure
./build/release/examples/smeft_example
```

The main examples are `smeft_example` and `testmatel2`. To build the five
older SM examples, additionally set `-DGRIFFIN_BUILD_LEGACY_EXAMPLES=ON`;
their executables are under `build/release/examples/legacy/`. The scattering
equivalence test runs with `BUILD_TESTING`, independently of example options.

Standalone summaries and optional diagnostics:

```bash
python3 tests/native_regression/compare_aas.py --runner build/release/tests/native_fixture_runner --summary build/release/regression_summary.json
python3 tests/upstream_regression/compare.py --runner build/release/tests/upstream_regression_runner --summary build/release/upstream_v111_regression_summary.json
cmake --build build/release --target near_pole_background near_pole_nnlo2 -j 4
./build/release/tests/near_pole_background
./build/release/tests/near_pole_nnlo2
```

The last recorded sanitizer run predates the import: **11/11 passed on
13 November**, with leak detection disabled for legacy allocating classes.
Its commands, timing and scope remain in the [historical record](REGRESSION_LOG_2026-11-16_PRE_IMPORT.md#reproduce-the-checks).
No leak-freedom or current sanitizer claim follows from that run.

These checks establish behavior within the documented prescription and
sampled domain, not complete higher-order SMEFT accuracy or validity at
arbitrary coefficients. Earlier revisions and the fresh Jackson rebuild
remain in the [integration log through 12 November](REGRESSION_LOG_2026-11-12.md).
