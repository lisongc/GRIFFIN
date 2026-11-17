<!--
Purpose: list practical changes for external users of earlier merged GRIFFIN versions.
Merged-version authors: Jackson Wallace and Lisong Chen.
Last updated: 2026-11-17.
-->

# Updating an existing program

External users should rebuild the library and their programs in a fresh
CMake build directory and install prefix. Installing over an older prefix
leaves obsolete headers behind;
compatibility aliases and forwarding headers are not supplied. Include
`<griffin/SMEFT.h>` and link to `GRIFFIN::griffin` for ordinary SMEFT use.

The FV and `matel` ownership repairs change class layouts; `psobs` also has
virtual destruction. Rebuild external programs with the updated library.
FV constructors take `Fermion`; copies have independent FA/SW objects.
After moving an FV object, reassign its source before evaluation or setters.
Custom FV subclasses must treat protected `fa`/`sw` as owning pointers;
custom FA/SW child types must override `clone()` to preserve their type on copy.

Both typed `matel` constructors remain available. Numeric-value construction
owns its four fixed form factors; object-based construction borrows them.
Copies duplicate owned objects or retain borrowed references; moves transfer
ownership. Reassign moved-from matrix elements before evaluation. Keep SM
inputs and borrowed form factors alive for as long as they are used.

## Header path mapping

Umbrella-header users need no include changes. For individual headers:

| Previous include | Current include |
|---|---|
| `prediction/SMEFTCalculator.h`, `griffin/SMEFTCalculator.h` | `griffin/smeft/SMEFTCalculator.h` |
| Separate `ElectroweakSetup.h` / `EWInputScheme.h` under `prediction/`, `griffin/` or `griffin/smeft/` | `griffin/sm/ElectroweakSetup.h` |
| `Observables.h` under those same directories | `griffin/smeft/EWPO_smeft.h` for observables; `griffin/smeft/xsc_smeft.h` for scattering types |
| `smeft/NativeInputPoint.h`, `griffin/smeft/NativeInputPoint.h` | `griffin/smeft/WMassSolution.h` |
| Root SM headers such as `SMval.h`, `EWPOZ2.h`, `xscaas.h` | Same filename under `griffin/sm/` |
| `classes.h`, `Cplx.h` | Same filename under `griffin/core/` |
| `oneloop.h`, `li.h`, `linex.h` | Same filename under `griffin/numerics/` |
| `xscmassless.h`, `griffin/core/FermionPairCrossSection.h` | `griffin/sm/tools.h` |

## Typed labels

All input access, particle labels and current labels now use `InputPar`,
`Fermion` and `Current`. Integer-label overloads and the old macros are removed.
For example, use `FA_SMNNLO(Fermion::electron, input)` and
`input.get(InputPar::WMassComplex)`; do not pass `index(...)` into these APIs.
The `inval(int)` constructor still takes an allocation size.

Choose mass keys by convention, not by header order:

| Value | Pole key | Running-width key |
|---|---|---|
| W mass | `WMassComplex` (0) | `WMassExperiment` (20) |
| Z mass | `ZMassComplex` (1) | `ZMassExperiment` (21) |
| W width | `WWidthComplex` (17) | `WWidthExperiment` (22) |
| Z width | `ZWidthComplex` (18) | `ZWidthExperiment` (23) |

These are `InputPar::` members. Fill the running-width keys in `SMval`;
it derives the pole values. SM formulas read the pole keys explicitly.
`SMvalGmu` derives the W mass and width from its Gmu inputs.
The former mass-macro header-order restriction is removed.

## API name changes

| Previous API | Replacement |
|---|---|
| `result.reference()` | `result.smResult()` |
| `calculation.referenceInput()` | `calculation.smInputValues()` |
| `calculation.referenceAmplitudes(point)` | `calculation.smAmplitudes(point)` |
| `smeft::coeffS4fD6Tree(...)` | `smeft::fourFermionContactAmplitude(...)` |
| Cross-section builder classes and `squaredPrediction()` | `griffin::differentialCrossSectionFromAmplitudes(...)` |

The calculator's `differentialCrossSection(point)`, `total()` and
`modelCorrection()` remain available. For the standalone conversion's
arguments and units, see [the user guide](USER_GUIDE.md#converting-supplied-amplitudes).

## Coefficient inputs

Cll/Cee aliases identify canonical operator strengths. Setters overwrite;
sum independently populated raw aliases first, then assign each unique
coefficient once. See the [normalization and import rules](SMEFT_IMPLEMENTATION.md#contact-aliases-and-canonical-strengths).

## Examples

`GRIFFIN_BUILD_EXAMPLES=ON` builds `smeft_example` and `testmatel2`.
The five older SM examples moved to `examples/legacy/` and require
`GRIFFIN_BUILD_LEGACY_EXAMPLES=ON`. Scattering equivalence is now a test,
independent of example options. Follow the [current build/test instructions](../README.md#testing).
