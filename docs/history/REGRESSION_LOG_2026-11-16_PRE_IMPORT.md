<!--
Purpose: preserve the dated pre-upstream-import regression chronology.
Merged-version authors: Jackson Wallace and Lisong Chen.
Archived: 2026-11-16.
-->

# Regression chronology — 13–16 November 2026, before the upstream v1.2 import

These sections were moved from the current regression report on 16 November.
They record earlier snapshots, in reverse chronological order; numerical
results, source hashes and verification dates retain their original meaning.
References to “current”, unchanged v1.1 kernels, local commit state or header
counts apply only to the dated stage. The later upstream v1.1.1/v1.2 import
changes the near-pole behavior and brings the installed header count to 25.
Use the [16 November regression audit](REGRESSION_REPORT_2026-11-16.md) for current evidence.

The earlier [integration log](REGRESSION_LOG_2026-11-12.md) covers the preceding
merge and layout work. The 12 November impact subsection below is retained
because it explains the distinct comparison baselines in these later checks.
Migration links lead to the consolidated current guide; the obsolete
three-row diagnostic reference has been replaced with its recorded values.

## EWPO/scattering implementation separation — 16 November

`ScatteringPoint` and `CurrentChannel` now live in `xsc_smeft.h`, with
unchanged definitions and namespaces. FA/SW/FV evaluation is implemented
in `src/smeft/EWPO_smeft.cc`. The calculator shares the same private FA/SW
helper between observable and scattering evaluations; no second public
evaluator or extra helper header was introduced. Source review confirms
byte-identical moved function bodies, unchanged validation order and
unchanged FV construction. The finite Gmu solver and SM kernels are untouched.

A fresh release build passed **11/11 suites** in **19.74 seconds**. All
**24 public headers** compile independently and when included twice;
**five installed clients** pass. A new client tests EWPO selections and
scattering requests through their respective headers separately. Compile-only
probes confirm that the extracted calculator helpers remain private.

All **943 complete numerical records** are byte-identical before and after,
with SHA-256
`e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`.
All 22,475 historical and 28,290 extended comparisons pass. Fixture values,
tolerances and the example's numerical output are unchanged. See the
[machine-readable record](validation/results_2026-11-16_ewpo_split.json)
and [include/build migration note](../MIGRATION.md#header-path-mapping).
Sanitizers were not rerun for this source reorganization. These checks do
not validate or repair the unchanged `mat_SMaas` near-pole limitation.
Changes remain local, uncommitted and unpushed.

## SMEFT observable header rename — 16 November

`smeft/Observables.h` is now `smeft/EWPO_smeft.h`, without an old-path
wrapper. The filename, include guard and calculator include were updated;
all enum, struct and function definitions and namespaces are unchanged.
Current guides and migration instructions use the new path. Source and
installed-package checks require it and reject the retired header.

A fresh release build passed **11/11 suites** in **15.91 seconds**, including
independent/double inclusion of **24 public headers** and **four installed
clients**. All **943 complete numerical records** are byte-identical before
and after, with SHA-256
`e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`.
The 22,475 historical and 28,290 extended comparisons pass. Fixture values,
tolerances and example results are unchanged. Source review against the
pre-change snapshot confirms identical type/function definitions.

See the [machine-readable record](validation/results_2026-11-16_ewpo_header.json)
and [migration instructions](../MIGRATION.md#header-path-mapping).
No physics or `mat_SMaas` changes were made. Sanitizers were not rerun for
this header rename. Changes remain local, uncommitted and unpushed.

## Contact-amplitude API cleanup — 16 November

The local contact helper is now `fourFermionContactAmplitude`, with the same
typed and integer overloads. The unused, identically zero derivative helper
was removed, without aliases. Contact matching, normalization, validation,
finite Gmu solving and the full amplitude square are unchanged. See the
[migration note](../MIGRATION.md#api-name-changes).

A fresh release build passed **11/11 suites** in **15.81 seconds**, including
all **24 public headers** and **four installed clients**. Tests cover the new
name through both overloads, coefficient aliases, all four vector/axial
channels, scale normalization and rejected processes/currents. The installed
public-interface client also checks a nonzero contact amplitude. The API
cleanliness test rejects either retired function name in production code.

All **943 complete numerical records** are byte-identical before and after,
with SHA-256
`e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`.
The regression passes 22,475 historical and 28,290 extended comparisons;
fixtures, tolerances and the example's numerical output are unchanged.
Independent source review against the pre-change working-tree snapshot
confirms unchanged matching formulas and addition of contacts exactly once.

See the [machine-readable record](validation/results_2026-11-16_contact_api.json).
Changes remain local, uncommitted and unpushed. Sanitizers were not rerun
for this interface-only cleanup; the `mat_SMaas` calculation is unchanged.

## Explicit SM result and input names — 16 November

The zero-coefficient outputs now use `smResult()`, `smInputValues()` and
`smAmplitudes(point)`. The internal stored SM result and mass solution use
matching names. No old-name aliases remain; `total()` and `modelCorrection()`
are unchanged. The SM input values still come from a separate zero-coefficient
Gmu solve, and the SM amplitudes use the requested energy and angle.
See the [accessor migration table](../MIGRATION.md#api-name-changes).

A fresh release build passed **11/11 suites** in 16.71 seconds. After the
final test-name cleanup, the complete suite passed again in **15.32 seconds**.
All **24 public headers** and **four installed clients** pass. New compile-time
checks verify the accessor return types, and runtime checks compare the SM
input values and amplitudes with the zero-coefficient calculation. The API
cleanliness test rejects reintroduced old accessor names.

All **943 complete numerical records** are byte-identical before and after,
with SHA-256
`e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`.
The regression passes 22,475 historical and 28,290 extended comparisons.
The example's numerical output, fixture data, JSON output fields/order and
tolerances are unchanged. Source comparison after only the intended identifier
substitutions and comment removal confirms unchanged calculator code.

Current and historical documentation use SM terminology for zero coefficients;
historical numbers and verification dates are preserved. Genuine benchmark
references and test-only comparison assemblies remain distinct: those routines
evaluate nonzero as well as zero coefficients and are not the SM result.
File-purpose/authorship/update metadata is retained or added to touched files.

See the [machine-readable record](validation/results_2026-11-16_sm_names.json).
Changes remain local, uncommitted and unpushed. Sanitizers were not rerun for
this naming-only change; no physics formula or `mat_SMaas` calculation changed.

## Shared electroweak setup and file descriptions — 16 November

`EWInputScheme` and its name function now share the single
`sm/ElectroweakSetup.h` header with `ElectroweakSetup`. The old two SMEFT
headers are removed without forwarding wrappers. The class body, enum and
name-function body are byte-identical to their preceding definitions; all
namespaces, borrowing semantics and compatibility checks are unchanged.
The calculator and mass-solution header include the combined SM header.

A fresh release build passed **11/11 suites** in 16.22 seconds. All **24
remaining public headers** compile independently and when included twice.
The installed package passes four downstream tests, including a new client
that includes only the shared setup header, checks both schemes and their
names, and verifies both input-type mismatch errors without using SMEFT.
Packaging checks require the new header and reject the obsolete paths.

All **943 complete numerical records** are byte-identical before and after
this consolidation, with SHA-256
`e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`.
The regression passes 22,475 historical and 28,290 extended comparisons;
the SMEFT example, fixtures and tolerances are unchanged. The SM kernels,
grids, finite Gmu mass-solving prescription and contact matching are unchanged.

All 30 surviving changed/new files in the pending cleanup have top-of-file
purpose, attribution and last-update metadata (16 November 2026). This covers
source, headers, build/test files, current guides and JSON verification records;
JSON uses a leading metadata object rather than invalid comment syntax.
Merged files without clear prior authorship credit Jackson Wallace and Lisong
Chen, as requested. Inherited GRIFFIN credits are preserved, including Ayres
Freitas and Lisong Chen's original utility-code record, separately from the
integration contributors. Historical verification results are not rewritten.

See the [machine-readable record](validation/results_2026-11-16_electroweak_setup.json)
and [migration note](../MIGRATION.md#header-path-mapping).
These are local uncommitted changes; nothing was pushed. Sanitizers were not
rerun for this header/comment-only consolidation, and the near-pole
`mat_SMaas` limitation remains outside this verification's scope.

## Cross-section conversion function — 16 November

The separate `FermionPairCrossSection` class/header/source were removed.
The same angular square and normalization now live in
`differentialCrossSectionFromAmplitudes()` in `sm/tools.h` and `sm/tools.cc`.
The calculator's public `differentialCrossSection(point)` method and
`Amplitudes` array type are unchanged. No old class alias or forwarding header
remains; see the [migration note](../MIGRATION.md#api-name-changes).

A fresh release build passed **11/11 suites** in 15.23 seconds, including
independent compilation of all **25 remaining public headers**. The isolated
installation runs three downstream tests: the SMEFT example, public-interface
tests, and a new client that includes only `sm/tools.h` and uses the conversion
function. The removed header is confirmed absent from the installation.

Tests independently check the angular formula, all four amplitude channels,
quadratic terms, all seven supported final states, color factors, units,
angular endpoints and rejected processes/kinematics. The standalone function
now explicitly rejects nonfinite kinematic arguments. Its inclusive angular
endpoints do not change the calculator's stricter endpoint restrictions.

All **943 complete numerical records** are byte-identical before and after
this replacement, with SHA-256
`e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`.
The regression passes 22,475 historical and 28,290 extended comparisons;
fixtures and tolerances are unchanged. The SMEFT example is also identical.
Independent review confirms the same channel order, angular expression and
normalization. Existing width utilities, coefficient matching, mass solving,
SM amplitude kernels and grids are unchanged; this does not repair or validate
the separate near-pole `mat_SMaas` issue.

See the [machine-readable record](validation/results_2026-11-16_cross_section_function.json).
Sanitizers were not rerun for this refactor. All edits remain local,
uncommitted and unpushed.

## Private particle lookup — 16 November

The separate `SMEFTFlavor.h` helper was removed. The lookup and its metadata
types are now private members of `smeft::Input`, declared in `SMEFTInput.h`
and implemented in `SMEFTInput.cc`. Only the existing source-defined vertex
and contact coefficient-selection helpers have access. No public method,
coefficient convention or input storage was added or changed.

A fresh release build passed **11/11 suites** in 15.65 seconds, including
the isolated installed-package test and all **26 standalone public headers**.
Compile-only installation probes confirm that the public input accessor is
usable while the lookup and its two metadata types are rejected as private
members. Expanded tests reject negative, out-of-range and unsupported
in-range particle identifiers in both direct shifts, FA/FV/SW and both
contact functions, before any legacy array lookup. Existing valid-generation
and contact-alias tests also pass.

All **943 complete numerical records** are byte-identical to the preceding
header-layout version, with SHA-256
`e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`.
The regression passes 22,475 historical and 28,290 extended comparisons.
The SMEFT example is unchanged. Review confirms that all 11 particle mappings,
vertex/contact formulas and normalization factors are preserved. The SM
kernels, grids and finite Gmu solver are untouched. Fixtures and tolerances
were not changed; sanitizers were not rerun for this refactor.

See the [machine-readable record](validation/results_2026-11-16_private_lookup.json)
and [migration note](../MIGRATION.md#header-path-mapping).
The edits remain local, uncommitted and unpushed. This does not address the
separate near-pole `mat_SMaas` issue.

## SMEFT header grouping — 16 November

At this earlier stage, the calculator, setup, input-scheme and observable
headers were grouped under `include/griffin/smeft/`. The setup and scheme
were subsequently consolidated under `sm/`, as recorded above. Only the
convenience entry `include/griffin/SMEFT.h` remains at the include root.
No forwarding headers or aliases were added.
The umbrella's include order and all class/namespace names are unchanged;
see the [include-path migration table](../MIGRATION.md#header-path-mapping).

A fresh release build of the updated working tree passed **11/11 suites**
in 14.48 seconds. A separate fresh clone of the pre-change commit
`d3e0fdb692176d7f492632dffbc99bdf040f09c4` was also built and passed
11/11 suites. The installed-package test used a fresh prefix, compiled
each of the **26 public headers** independently and twice, and ran the
downstream example and API tests. New checks enforce the module paths
and reject leftover root headers in both the source and installed package.

All **943 complete numerical records** are byte-identical before and after
the move, with SHA-256
`e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`.
Both runs pass 22,475 historical and 28,290 extended comparisons; fixtures
and tolerances are unchanged. The SMEFT example also gives identical results.

An exact source comparison allowing only the include-path substitutions
confirms unchanged C++ definitions and calculations. The SM kernels, grids,
finite Gmu solver, direct shifts and contacts are untouched. In particular,
this layout change does **not** apply the upstream `mat_SMaas` correction
or validate the near-pole background. Sanitizers were not rerun for this
path-only change. Independent review found no stale active include paths
or broken links in the updated guides.

Machine-readable record:
[results_2026-11-16_header_layout.json](validation/results_2026-11-16_header_layout.json).
The layout edits were tested locally without committing or pushing them.

## W-mass solution rename — 14 November

The input-solution class/header/source are now named `WMassSolution`, and
the solver settings are `WMassSolverOptions`. No compatibility alias or
old-name forwarding header remains. The file comments explain the Gmu
solution, alpha-scheme input copy, ownership and automatic construction by
the calculator; the guides explain public, private and test-only headers.

A fresh release build passed **11/11 suites** in 15.01 seconds, including
the installed-package test and independent compilation of all 26 public
headers. The installed package contains `WMassSolution.h`, not the retired
header. The SMEFT example passed. All **943 complete numerical records**
are byte-identical to the pre-rename commit
`177205835f9c238dd45be6e9952667190c209d05`, with SHA-256
`e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`.
The regression passed 22,475 historical and 28,290 extended comparisons;
fixtures and tolerances were unchanged.

Source comparison, allowing only the intended names, comments and formatting,
confirms unchanged solver/calculator numerical code and defaults. The
`mat_SMaas` calculation and warning are unchanged. Sanitizers were not rerun
for this naming-only change. See the
[machine-readable record](validation/results_2026-11-14_wmass_rename.json)
and [migration note](../MIGRATION.md#header-path-mapping).

## Committed collaborator checkout — 14 November

The merged implementation and updated warning were committed on the local
branch `codex/merged-smeft`, at
`66b0edd9d9e8ad599a152b857d2b1b4f77ea2eb4`. A separate fresh local clone
(`git clone --no-hardlinks --branch codex/merged-smeft`) was built using only
committed files. No existing build products or external comparison archives
were needed. Nothing was pushed.

Running `./scripts/run_checks.sh` in that clone passed **11/11 suites** in
14.96 seconds, including installed-package and 26 standalone-header checks;
the SMEFT example also passed. The 943-case regression passed **22,475
historical** and **28,290 extended** numerical comparisons. Its complete
output SHA-256 is
`e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`,
identical to the pre-warning output. Fixtures and tolerances were unchanged.

The optional `near_pole_background` target built and reproduced the historical
v1.1 diagnostic: Im Q(h) was approximately 6.738e-4 at h = 0.03, 0.01 and
0.003 GeV^2. The exact three values are retained in the
[machine-readable record](validation/results_2026-11-14_fresh_checkout.json).
This confirmed the limitation of that snapshot; it did not validate near-pole
predictions or fix the matching. Hash checks confirm unchanged `xscaas.cc`,
`imaas.in`, `aasoffz.grid`, `xscnnlo.cc` and `delrho.cc` relative to the start
of this handoff. The exact-point rejection condition is unchanged; only its
message, source comments and documentation were clarified.

Both the source repository and fresh checkout were clean after testing.
This report and its [machine-readable record](validation/results_2026-11-14_fresh_checkout.json)
are a subsequent documentation-only commit. Sanitizers were not rerun for
this handoff; their most recent run is recorded below.

## Cross-section helper rename — 14 November

On 14 November, the public helper was renamed to `FermionPairCrossSection`,
with matching `.h` and `.cc` files under `core/`. The class, includes, build
references, tests and guides were updated directly, without old-name aliases
or forwarding headers. This intermediate class was replaced by the conversion
function on 16 November, as recorded above.
See [MIGRATION.md](../MIGRATION.md#api-name-changes).

This is a naming-only change: the angular formula, channel order, supported
processes, massless external-fermion approximation, normalization, finite
Gmu solution and full-square evaluation are unchanged. Source comparison
confirms that the core implementation and calculator code differ only by
the intended identifier/include substitutions and one local variable name.

A fresh macOS arm64/AppleClang 16 release build passed **11/11 test suites**
in 14.82 seconds, including the isolated installed-package/public-header
checks. The SMEFT example also passed. All **943 complete output records**
are byte-for-byte identical before and after the rename; their common
SHA-256 is `e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`.
The archived and extended regressions passed 22,475 and 28,290 numerical
comparisons respectively. Fixtures and tolerances were unchanged.

Machine-readable verification:
[results_2026-11-14_helper_rename.json](validation/results_2026-11-14_helper_rename.json).
Sanitizers were not rerun for this name-only change; the sanitizer result
below belongs to the 13 November run. These checks do not resolve or
independently validate the near-pole SM background issues.

## Previous documentation cleanup — 13 November

The 13 November documentation/code-comment cleanup changes no numerical
formula or public API. All **943 complete output records are byte-for-byte
identical** before and after this cleanup. Release and sanitizer builds each
pass **11/11 test suites**.

## What the regression establishes

The older Jackson cross-section fixtures used `mat_SMNNLO`. The current
`mat_SMaas` calculation intentionally adds SM terms, so it must not be
claimed to reproduce those historical cross sections unchanged.

| Check | Coverage | Result |
|---|---|---|
| Historical `mat_SMNNLO` assembly versus unchanged Jackson archive | 943 cases; 22,475 numerical comparisons | PASS |
| Current public calculation versus separate `mat_SMaas` assembly | 943 cases; 28,290 numerical comparisons | PASS |
| Unchanged masses, widths, EWPOs and contacts | Public/assembly checks and available archived values | PASS |
| Complete public output before/after documentation cleanup | All 943 records, including amplitudes and cross sections | Byte-identical |

The test-only assemblies construct the pole inputs, amplitude and angular
square separately from the public calculation. They share the underlying SM
kernels and the public resolved mass inputs. This checks assembly and
preserved matching; it is **not** an independent derivation of the SM kernels
or a new set of Jackson-supplied extended benchmarks. All comparisons use
finite results, not a linear extraction. The historical standalone Jackson
FV is not a golden target: the merged FV is constructed from its FA/SW and
the existing SM imaginary-ratio term.

Machine-readable results for the 13 November cleanup:
[results_2026-11-13_docs.json](validation/results_2026-11-13_docs.json).
Fixture values and tolerances were not changed.

## Documentation and code cleanup — 13 November

All current guides use the same NNLO+/finite-prescription terminology.
They distinguish fixed-mass explicit shifts from generated nonlinear terms,
and `modelCorrection()` from a derivative or linear interference.
A [documentation index](../README.md) supplies the reading order; the older
development narrative is preserved under [history/](README.md).

The source cleanup touches comments and exception-message wording only.
A token audit confirms unchanged C++ code after excluding comments and
reviewed string literals; the only changed literals are clearer SMEFT error
messages. All analytic expression fragments are unchanged. Public names,
include order and numerical formulas are unchanged.

Additional checks:

- All local Markdown file/heading links in the guides, current report,
  archived records, root README and fixture README resolve.
- The user-guide C++ snippet compiles and runs after supplying exactly the
  SM input initialization it references in the working example.
- Independent review found no blocking mismatch between the current
  documentation, code comments and implementation.
- The before/after public-output SHA-256 is
  `e03a1f2521d06eedfeb95068f0d9ca00c712308b015064736c0ac0fb4f81ec21`.
- The fixture SHA-256 remains
  `956234c20636f83187518926da5aa0e8b097297d7a269f012eb520f3bad919f1`.

The snapshot, audit script, snippet check and complete before/after outputs
are preserved in `../diagnostics/documentation_cleanup_20261113_CYx3Lk/`,
relative to the merged repository root. At that 13 November stage nothing
was committed or pushed, and neither original source repository was modified.
The later local integration commit is recorded above.

## Numerical effect of the preceding amplitude extension

This subsection records the **12 November physics-content update**, not a
change caused by the 13 November cleanup. The extension affected scattering
amplitudes; masses, widths, FA, SW, FV and contact amplitudes were unchanged.

For the fixture SM inputs in the Gmu scheme, with the process/angle above:

| Energy [GeV] | Previous `mat_SMNNLO` [pb] | Current `mat_SMaas` [pb] | Change [%] |
|---:|---:|---:|---:|
| 50 | 14.10596034 | 14.10210341 | -0.027343 |
| 90 | 448.26992780 | 447.59089166 | -0.151479 |
| 95 | 132.59554209 | 132.64568459 | +0.037816 |
| 120 | 9.77048132 | 9.77193638 | +0.014892 |

Here the percentage is $100(\sigma_{\rm new}-\sigma_{\rm old})/\sigma_{\rm old}$.
With only $C_{\phi WB}/\Lambda^2=+0.01\,\mathrm{TeV}^{-2}$ at 90 GeV,
the finite change $\sigma(C)-\sigma(0)$ went from -0.5688478552 pb to
-0.5638142269 pb: a 0.8849% reduction in magnitude.

Across the 943 scanned inputs, the largest relative change in a **full cross
section** was 0.15343%. The largest relative change in a **finite SMEFT
difference** was 6.2165%, measured against the absolute old difference:
-0.1258462224 to -0.1180230383 pb in a combined vertex/contact case. That
percentage uses a small cancellation-sensitive difference as its denominator;
it is not a 6% full-cross-section change or a bound for arbitrary coefficients.

See [the saved impact study](validation/aas_impact_2026-11-12.json)
and [extension test results](validation/results_2026-11-12_aas.json).

## Reproduce the checks

Environment used: macOS arm64, AppleClang 16.0.0.16000026, C++11. Run from
the merged repository root; Python 3's standard library is sufficient.

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DGRIFFIN_BUILD_EXAMPLES=ON
cmake --build build/release -j 4
ctest --test-dir build/release --output-on-failure
./build/release/examples/smeft_example
```

13 November release: **11/11 pass**, 15.57 seconds.
CTest writes `build/release/native_regression_summary.json`. To retain a
standalone numerical summary:

```bash
python3 tests/native_regression/compare_aas.py --runner build/release/tests/native_fixture_runner --summary build/release/regression_summary.json
```

Sanitizer configuration:

```bash
cmake -S . -B build/sanitize -DBUILD_TESTING=ON -DGRIFFIN_BUILD_EXAMPLES=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined'
cmake --build build/sanitize -j 4
ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=halt_on_error=1 ctest --test-dir build/sanitize --output-on-failure
```

13 November sanitizer: **11/11 pass**, 84.95 seconds; no reported address or
undefined-behavior errors. Leak detection was disabled for known legacy
allocating classes. This is not a leak-freedom claim.
