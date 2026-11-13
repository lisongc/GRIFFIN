# Migrating to the cleaned layout

Updated 14 November 2026. Rebuild the library and client programs after
updating the local merged version. The source-layout/interface cleanup does
not change the finite Gmu equation, coefficient normalization, FA/SW insertion,
FV construction or full-square treatment. The separate 12 November
scattering-content update changes `mat_SMNNLO` to `mat_SMaas`, as described
below. The 13 November documentation/comment cleanup does not change the
calculation or public class names.

## Cross-section helper rename — 14 November

The low-level helper is now `griffin::FermionPairCrossSection`, declared in
`<griffin/core/FermionPairCrossSection.h>` and implemented in
`src/core/FermionPairCrossSection.cc`. Replace
`MasslessFermionCrossSectionBuilder` and either old `xscmassless.h` include
spelling with these names; no compatibility alias or forwarding header is
retained.

This is a name-only change. The helper still computes the unpolarized
cross section for supported charged-fermion-pair annihilation with massless
external-fermion kinematics. No additional processes or corrections are
introduced. Calls to `SMEFTCalculator::differentialCrossSection()` and the
calculator's public `Amplitudes` type are unchanged.

## Scattering-content update

Existing `SMEFTCalculator` calls now include the additional selected fermionic
mixed electroweak–QCD pole-residue and off-resonance terms supplied by `mat_SMaas`.
No new caller option is required: both the full result and its independently
computed zero-coefficient reference use the extension at their resolved
inputs. Contact amplitudes still occur once, and the full square is retained.

Cross sections and their finite coefficient-induced shifts can therefore
change. The W-mass matching and standalone FA/SW/FV formulas are unchanged.
The original Jackson fixture numbers remain unchanged: the test suite
distinguishes the archived `mat_SMNNLO` baseline from the extended-amplitude
checks. See the
[regression report](REGRESSION_TEST_REPORT.md) for numerical effects and
verification. The runtime exact-real-pole restriction and process coverage
remain unchanged. The 14 November handoff adds a warning and optional
diagnostic for unresolved behaviour throughout the nearby region; it does
not repair or change `mat_SMaas`. See
[known issues](KNOWN_ISSUES.md#1-unresolved-near-pole-scattering-background).

Use **NNLO+ SM form factors with a finite, unlinearized dimension-six
prescription** to describe the current calculation. The FA/SW classes include
selected terms beyond NNLO. Explicit coefficient shifts at a fixed mass
remain linear, while the mass solution, pole construction and full square
retain their generated nonlinear dependence. No explicit $(\delta a)^2$ term
was added to the FA input. This is not a claim of complete off-resonance NNLO
or of every possible higher-power SMEFT contribution; see the
[exact prescription](USER_GUIDE.md#exact-prescription).

## Recommended includes and build target

For ordinary finite-SMEFT use:

```cpp
#include <griffin/SMEFT.h>
```

This supplies the calculator, setup, observable/kinematic types, coefficient
input and SM input classes. The fully qualified class names remain
`griffin::prediction::ElectroweakSetup` and
`griffin::prediction::SMEFTCalculator`; `griffin::smeft::Input` is unchanged.
No namespace aliases or class-name wrappers were introduced.

Include this convenience header first when also using advanced SM headers.
It loads the relevant inline SM formulas before the legacy mass-macro
redefinitions. Prefer typed `InputPar` setters; arbitrary legacy input-header
ordering remains a separate [known issue](KNOWN_ISSUES.md#5-ownership-and-legacy-error-handling).

With the repository as a CMake subdirectory, link to `GRIFFIN::griffin`.
For an installation, use:

```cmake
find_package(GRIFFIN CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE GRIFFIN::griffin)
```

Set `CMAKE_PREFIX_PATH` to the install prefix when configuring your client.
For manual builds, use `-I/path/to/installation/include` and link the installed
library. Do not add the repository root, `src/`, or `data/` as public include
directories. C++11 remains the minimum language level; the cleaned build
requires CMake 3.20 or newer.

## Header path mapping

Each new path below is the spelling used in `#include <...>`. On disk it
is under `include/`. Old include locations are not retained as wrappers.

| Old include | New include |
|---|---|
| `prediction/SMEFTCalculator.h` | `griffin/SMEFTCalculator.h` |
| `prediction/ElectroweakSetup.h` | `griffin/ElectroweakSetup.h` |
| `prediction/EWInputScheme.h` | `griffin/EWInputScheme.h` |
| `prediction/Observables.h` | `griffin/Observables.h` |
| `classes.h` | `griffin/core/classes.h` |
| `Cplx.h` | `griffin/core/Cplx.h` |
| `xscmassless.h` or `griffin/core/xscmassless.h` | `griffin/core/FermionPairCrossSection.h` |
| `SMval.h` | `griffin/sm/SMval.h` |
| `SMvalG.h` | `griffin/sm/SMvalG.h` |
| `EWPOZ.h` | `griffin/sm/EWPOZ.h` |
| `EWPOZ2.h` | `griffin/sm/EWPOZ2.h` |
| `deltar.h` | `griffin/sm/deltar.h` |
| `delrho.h` | `griffin/sm/delrho.h` |
| `ff.h` | `griffin/sm/ff.h` |
| `ff0.h` | `griffin/sm/ff0.h` |
| `xscnnlo.h` | `griffin/sm/xscnnlo.h` |
| `xscaas.h` | `griffin/sm/xscaas.h` |
| `tools.h` | `griffin/sm/tools.h` |
| `oneloop.h` | `griffin/numerics/oneloop.h` |
| `li.h` | `griffin/numerics/li.h` |
| `linex.h` | `griffin/numerics/linex.h` |
| `smeft/SMEFTInput.h` | `griffin/smeft/SMEFTInput.h` |
| `smeft/NativeInputPoint.h` | `griffin/smeft/NativeInputPoint.h` |
| `smeft/ff0_smeft.h` | `griffin/smeft/ff0_smeft.h` |
| `smeft/xsc_smeft.h` | `griffin/smeft/xsc_smeft.h` |

The old `smeft/SMEFTFlavor.h` was a helper rather than a supported user
interface. It is now private at `src/smeft/detail/SMEFTFlavor.h` and is not
installed. Client code should use the public `Fermion` and `Generation`
types instead of depending on that helper.

## Removed infrastructure

The following generic components are removed, not optional or deprecated:

- `PredictionEngine`, `TheoryProvider` and `SMProvider`.
- `PerturbativeSeries`, `OrderTag`, `LoopOrder`, `Accuracy` and their sector metadata.
- `PredictionResult`, `CombinationRule` and `CombinationKind`.
- The provider-only `ElectroweakSetup::inputAddress()` and `compatibleWith()` helpers.
- The `testprediction` example and synthetic generic-provider tests.

For SMEFT, call `SMEFTCalculator` directly. `FinitePrediction` is its result
type; `modelCorrection()` means the finite total minus the independently
computed zero-coefficient reference, not a coefficient of a truncated series.
There is no `Accuracy` replacement because this calculator intentionally
uses fixed SM content, now including the `mat_SMaas` scattering extension.

For SM-only work, use the retained SM classes directly from
`griffin/sm/`. The existing SM examples demonstrate those interfaces. A
zero-Wilson `SMEFTCalculator` also supplies its fixed-prescription SM limit,
but it is not an order-selection wrapper.

Earlier cleanups removed `linearEFT()`, the interference-only cross-section
helpers, the old linear SMEFT providers/observables and the unimplemented
`dr_SMNLO::compmw()` declaration. No compatibility layer restores them.
Use `SMvalGmu`/`invalGmu` for the SM mass solution and the finite calculator
for the merged SMEFT solution.

If migrating from older public names, replace `CalculationContext` with
`ElectroweakSetup` and `NativeSMEFTPrediction` with `SMEFTCalculator`.
The direct vertex functions are `directZAxialShift` and `directZVectorShift`;
they contain no separate linear W-mass response.

## Contact inputs still require canonical strengths

The 9 November Cll/Cee alias rules are unchanged by the folder move. Setters
overwrite one canonical value; they do not add equivalent labels. When
importing independent redundant raw entries, sum them first and assign
each unique coefficient once. See the
[contact-input migration examples](SMEFT_IMPLEMENTATION.md#contact-aliases-and-canonical-strengths).

## Examples, documents and build products

| Previous location/name | Current location/name |
|---|---|
| `examples/testnativeprediction.cc`, target `testnativeprediction` | `examples/smeft_example.cc`, target `smeft_example` |
| `scripts/run_prediction_architecture_checks.sh` | `scripts/run_checks.sh` |
| `NATIVE_SMEFT.md` | `docs/USER_GUIDE.md` |
| `SMEFT_IMPLEMENTATION.md` | `docs/SMEFT_IMPLEMENTATION.md` |
| `KNOWN_ISSUES.md` | `docs/KNOWN_ISSUES.md` |
| `REGRESSION_TEST_REPORT.md` | `docs/REGRESSION_TEST_REPORT.md` |
| `SMEFT_JACKSON_COMPARISON.md` | `docs/history/SMEFT_JACKSON_COMPARISON.md` |
| `PREDICTION_ENGINE.md` | Replaced by `docs/ARCHITECTURE.md` |

The recommended fresh build directories are `build/release` and
`build/sanitize`. Do not reuse a relocated CMake cache: its absolute paths
refer to the old layout. The pre-cleanup source snapshot and archived old
build trees are recoverable outside the source tree under
`../diagnostics/layout_cleanup_20261112_ObbFpc/` (old build trees in its
`old_builds/` subdirectory).

The retained `.in` expression fragments and `.grid` tables are compiled
into the library. Their relocation does not require runtime file searches
or a data-directory setting in external programs.
