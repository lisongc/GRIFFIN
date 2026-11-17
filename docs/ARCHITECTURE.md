<!--
File: ARCHITECTURE.md
Purpose: Explain calculation flow, ownership and source/build boundaries.
Authors (merged-version documentation): Jackson Wallace and Lisong Chen
Last updated: 2026-11-16
-->

# Architecture and source map

The finite SMEFT calculator uses NNLO+ SM form factors and `mat_SMaas`.
Its SM kernels incorporate upstream `v1.2` snapshot `ca2d523`; the additional
`mat_SMNNLO2` class is exposed separately for direct SM calculations.

## The public calculation

Start with `<griffin/SMEFT.h>` and [the example](../examples/smeft_example.cc).

| Object | Responsibility |
|---|---|
| `SMval` / `SMvalGmu` | Ordinary SM inputs, including pole-mass conversion |
| `smeft::Input` | Wilson coefficients and Lambda |
| `prediction::ElectroweakSetup` | Pair an electroweak scheme with compatible SM inputs |
| `prediction::SMEFTCalculator` | Resolve inputs and evaluate the finite prescription |
| `ScatteringPoint` / `CurrentChannel` | Species, kinematics, units and V/A channel |
| `PseudoObservable` | Select WMass, FA, SW or FV |
| `FinitePrediction` | Return `smResult()`, `total()` and `modelCorrection()` |

`griffin::prediction` is the public calculator's namespace. There is no
generic prediction engine or public SM-order selector.

## How a finite prediction is made

```text
SM inputs + Wilson coefficients + electroweak scheme
                         |
                  SMEFTCalculator
                         |
        two independently resolved input sets:
        requested coefficients / zero coefficients
                         |
        NNLO+ FA/SW + direct SMEFT vertex shifts
                         |
        mat_SMaas (includes mat_SMNNLO)
          pole/background + contacts once
                         |
                complete angular square
                         |
          total / SM result / their difference
```

`WMassSolution` solves the finite Gmu mass equation with SM Delta r
reevaluated at each iterate; the alpha scheme keeps the supplied masses.
FA/SW feed both observable requests and the scattering pole construction.
FV uses the same FA/SW and SM imaginary-ratio term. The formula definitions
are maintained in the [user guide](USER_GUIDE.md#exact-prescription) and
[SMEFT implementation](SMEFT_IMPLEMENTATION.md), rather than duplicated here.

The scattering step follows the implementation in
[`SMEFTCalculator.cc`](../src/smeft/SMEFTCalculator.cc). The classes are cumulative:

| Class | Scattering content |
|---|---|
| `mat_SMNNLO` | NNLO pole construction with the NLO off-resonance remainder |
| `mat_SMaas : mat_SMNNLO` | Adds selected fermionic mixed electroweak–QCD residue and background terms |
| `mat_SMNNLO2 : mat_SMaas` | Adds the active v1.2 fermionic two-loop off-resonance and vacuum-polarization terms |

Using `mat_SMNNLO` alone would omit the mixed scattering extension already
included in the calculator. The v1.2 import exposes `mat_SMNNLO2` separately;
it does not replace the calculator's `mat_SMaas` assembly. Selecting it there
would change the calculation and require validation of its added terms with
the finite SMEFT pole inputs and supported flavors. Its disabled sectors and
near-pole box behavior are listed in [known issues](KNOWN_ISSUES.md#upstream-v12-sm-class).
Local SMEFT contacts are absent from all three SM classes and are added once.

`differentialCrossSectionFromAmplitudes()` in `sm/tools.h` converts complete
VV, AV, VA, AA amplitudes into an unpolarized cross section, including color,
angular factors, normalization and supported-process checks. It adds no
radiative corrections. `SMEFTCalculator::differentialCrossSection()` combines
amplitude evaluation with this conversion.

## Input ownership

`ElectroweakSetup` borrows its SM input, which must remain alive while the
descriptor is used. `SMEFTCalculator` copies both SM and coefficient values
at construction; later caller mutations cannot alter it. Construct another
calculator for a changed input point and reuse it for energies/angles at
the same point.

The calculator owns two `WMassSolution` values and stack-owned pole inputs.
`resolvedInput()` and `smInputValues()` expose its requested and zero-Wilson
solutions; `massSolution()` exposes convergence diagnostics.

The direct SM `FV_SMLO`, `FV_SMNLO` and `FV_SMNNLO` classes own their FA/SW
pair with `std::unique_ptr`. Copies duplicate the owned form factors while
continuing to borrow the same SM input; moves transfer the pair. Inputs must
outlive these FV objects.

`matel` and its derived SM classes own the four `psobsfix` objects created
from numeric values; the other constructor borrows supplied form factors.
Copies duplicate owned objects or retain borrowed references; moves transfer
ownership. Reassign moved-from objects before evaluation. SM inputs and
borrowed form factors must outlive the matrix element.

## Folder responsibilities

| Location | Contents |
|---|---|
| `include/griffin/SMEFT.h` | Convenience entry header for SMEFT calculations |
| `include/griffin/smeft/` | Calculator, coefficients, EWPO/scattering requests, solver and shift interfaces |
| `include/griffin/sm/` | SM inputs, form factors, amplitudes and shared `ElectroweakSetup.h` |
| `include/griffin/core/`, `numerics/` | Shared types, base classes and numerical interfaces |
| `src/smeft/EWPO_smeft.cc` | Shared private FA/SW evaluation and FV construction |
| `src/smeft/SMEFTCalculator.cc` | Input coordination and complete scattering assembly |
| `src/smeft/xsc_smeft.cc` | Local four-fermion contacts; other SMEFT sources handle inputs, masses and shifts |
| `src/sm/`, `src/core/`, `src/numerics/` | SM physics, shared implementation and numerical algorithms |
| `src/sm/expressions/`, `src/sm/grids/` | Formula fragments and tables included during compilation |
| `examples/` | Main SMEFT and v1.2 SM examples; older SM demonstrations in `legacy/` |
| `tests/` | API, assembly, scattering-equivalence, independent-reference and diagnostic checks |
| `docs/history/` | Historical evidence, separate from current API instructions |

Particle-to-generation lookup is private bookkeeping within `smeft::Input`.
The EWPO helpers remain private calculator members, shared by observables
and scattering. Splitting these files creates no additional calculation path.

## Build boundary

`GRIFFIN::griffin` exposes only `include/` to external users. SM sources include
their expressions and grids locally. Installation supplies the archive, public headers
and CMake package metadata. Expressions and grids are compiled into the
archive, so external users need no runtime copies. Generated files stay under `build/`.

`GRIFFIN_BUILD_EXAMPLES` enables the two main examples.
`GRIFFIN_BUILD_LEGACY_EXAMPLES` independently enables `examples/legacy/`;
both default to OFF. `BUILD_TESTING` controls the test suite, including
scattering equivalence. Test reference headers are neither library inputs
nor installed interfaces.

Input access and particle/current interfaces use `InputPar`, `Fermion` and
`Current` throughout. SM kernels read pole mass/width keys; `SMval` converts
the separate running-width inputs. Header order cannot redefine these keys.

Tests distinguish historical Jackson fixtures, assembly checks sharing the
production kernels, and independent upstream references at matched inputs.
See [testing](../README.md#testing) for commands and reference provenance.
