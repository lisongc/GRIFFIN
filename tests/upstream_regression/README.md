# Independent upstream v1.1.1 regression

This test compares the public zero-Wilson `SMEFTCalculator`, which uses
`mat_SMaas`, with frozen predictions from **a separately compiled upstream
library**. The reference is the immutable `v1.1.1` commit
`a18f53830283a78155750a23a3594ada8de4c074` (tree
`55144aa54da8005033f87aacba95b60d54ce542e`). Its source is extracted with
`git archive` and built without edits. `reference_driver.cc` links only to
that upstream archive; it does not include merged headers or call the
merged scattering kernels.

The reference includes the upstream v1.1.1 mixed electroweak–QCD stability
correction. It is appropriate for the retained `mat_SMaas` calculation;
it is not a reference for the new `mat_SMNNLO2` additions.

## Input matching and independence

Both `{alpha, MW, MZ}` and `{Gmu, alpha, MZ}` schemes are covered. All Wilson
coefficients are zero and Lambda is 1000 GeV. `manifest.json` records the
original running-width inputs, the legacy indices used for resolved inputs,
the compiler, source hashes, units, and comparison tolerances.

Generation asks the public calculator for its **separately resolved
zero-Wilson SM inputs** (`smInputValues()`). These exact numeric masses,
widths, couplings, and fermion masses are passed to a plain upstream `inval`.
The upstream reference therefore performs no second W-mass solution or
pole-mass conversion. In particular, the Gmu comparison isolates the
scattering calculation from a mismatch between old and new W-mass solvers.
The resolved inputs are also frozen in every fixture and checked during
each offline run, for both `resolvedInput()` and `smInputValues()`.

**The amplitude kernels and angular-square reference are independent; the
Gmu mass solver is not independently validated by this test.** The existing
mass-solution and archived Jackson tests cover that separate responsibility.
The cross section is the massless external-fermion unpolarized
`d sigma / d cos(theta)` in pb, using `0.38937966e9` for GeV^-2 to pb.

## Coverage and acceptance

The strict scan has 560 cases:

- Both electroweak input schemes.
- Electron initial states and muon, tau, u, d, s, c, b final states.
- Energies 30, 50, 70, 90, 95, 120, 160, 190 GeV.
- Cosines -0.8, -0.3, 0, 0.5, 0.8.
- All four complex amplitudes in VV, AV, VA, AA order, checked through both
  `amplitudes()` and `smAmplitudes()`, plus both cross-section views.

These give 10,080 direct numerical comparisons with the independent
reference, in addition to resolved-input and zero-Wilson identity checks.
The acceptance rule is `abs(actual-expected) <= atol + rtol*abs(expected)`.
Amplitude components use `atol=2e-15`, `rtol=2e-11`; cross sections use
`atol=1e-9 pb`, `rtol=2e-11`. These match the existing regression tolerances;
they were not widened to accommodate this port.

A separate 56-case diagnostic samples `s-MZ_pole^2 = +/-0.01, +/-0.03 GeV^2`
at cosine 0.5 for both schemes and all seven final states. These points lie
inside v1.2's deliberately changed one-loop near-pole interpolation. The
summary records their differences from v1.1.1, but **does not require
cross-version equality there**. Finite output, fixed inputs, and equality
of the two zero-Wilson public views are still required. Neither diagnostic
values nor a passing regression validate precision Z-lineshape predictions.
The exact real pole is excluded.

## Running offline

Build the `upstream_regression_runner` target, then run:

```sh
python3 tests/upstream_regression/compare.py \
  --runner build/release/tests/upstream_regression_runner \
  --summary build/release/upstream_v111_regression_summary.json
python3 tests/upstream_regression/test_compare.py
```

CTest uses the checked-in fixtures; it needs Python's standard library,
but neither Git access, a reference checkout, nor network access. Fixture
bytes and the reference driver's source are checked against SHA-256 values
in the manifest. The harness negative checks verify that common drift in
both public amplitude views, mass drift, cross-section drift, malformed
records, and nonfinite values are detected.

## Reproducing the reference

Use a local Git clone containing the frozen commit and a built public runner:

```sh
python3 tests/upstream_regression/generate.py \
  --runner build/release/tests/upstream_regression_runner \
  --upstream-repo upstream/GRIFFIN-v1.2 \
  --work-dir /private/tmp/griffin-v111-reference-reproduction \
  --output-dir /private/tmp/griffin-v111-reference-fixtures
```

The work directory must be empty. The script extracts and builds the frozen
source, compiles the independent driver, feeds it the public runner's
resolved inputs, and writes fixtures plus provenance. Its retained work
directory contains the source, build log, and exact reference-driver input
and output. `--cxx` selects another compiler. Different compilers may produce
different final floating-point digits or fixture hashes; assess numerical
differences with the stated tolerances.

Fixture regeneration is an explicit maintenance action, never part of
CTest. Review mismatches and source provenance before replacing fixtures;
do not regenerate merely to make a failing regression pass.
