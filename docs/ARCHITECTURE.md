<!--
File: ARCHITECTURE.md
Purpose: Explain calculation flow, input ownership and the public/private source layout.
Authors (merged-version documentation): Jackson Wallace and Lisong Chen
Last updated: 2026-11-17
-->

# Architecture and source map

The library has one finite SMEFT calculation path. Start with
`<griffin/SMEFT.h>` and [the working example](../examples/smeft_example.cc).
For formulas and conventions, use [SMEFT_IMPLEMENTATION.md](SMEFT_IMPLEMENTATION.md);
this page describes which code is responsible for each step.

## Calculation flow

```text
SM inputs + Wilson coefficients + ElectroweakSetup
                         |
                  SMEFTCalculator
                         |
          two independent input solutions:
          requested coefficients / all coefficients zero
                         |
          NNLO+ SM FA/SW at each resolved MW
              + direct SMEFT vertex shifts
                         |
          mat_SMaas pole and background amplitude
              + contact amplitude exactly once
                         |
                  full angular square
                         |
          total / SM result / total minus SM result
```

In Gmu, `WMassSolution` solves Jackson's finite mass equation and
re-evaluates SM Delta r at each iterate. In alpha, it retains supplied masses.
There is no extra linear W-mass derivative. The same FA/SW calculation supplies
both observable requests and scattering pole inputs; FV is constructed from
FA/SW and the existing SM imaginary-ratio term.

Both input solutions use `mat_SMaas`, including its selected mixed
electroweak–QCD residue/background terms and the upstream v1.1.1 code/grid
correction. Explicit dimension-six FA/SW shifts are linear at fixed mass;
the later finite mass solution, pole construction and full square retain
their generated nonlinear dependence. This is not every higher-power SMEFT
term or complete NNLO away from the resonance.

`differentialCrossSectionFromAmplitudes()` in `sm/tools.h` and
`src/sm/tools.cc` converts the complete `VV, AV, VA, AA` amplitudes to
an unpolarized differential cross section, supplying angular factors, color
multiplicity and normalization. It does not compute amplitudes or add
corrections: `mat_SMNNLO` and `mat_SMaas` perform that separate job.

## Public objects and ownership

| Object | Responsibility |
|---|---|
| `SMval` / `SMvalGmu` | Ordinary SM inputs |
| `InputPar` / `Fermion` / `Current` | The only SM parameter, particle and current labels; no parallel integer-label interface |
| `smeft::Input` | Wilson coefficients, aliases and Lambda |
| `prediction::ElectroweakSetup` | Pair SM inputs with the chosen electroweak scheme |
| `prediction::SMEFTCalculator` | Resolve inputs and evaluate finite full/SM predictions |
| `ScatteringPoint` | Initial/final species, energy squared, angle cosine and unit factor |
| `PseudoObservable` | Select WMass, FA, SW or FV |
| `FinitePrediction` | Return `total()`, `smResult()` and `modelCorrection()` |

`ElectroweakSetup` borrows its SM input, which must remain alive while
the setup is used. `SMEFTCalculator` copies inputs and coefficients at
construction; later caller changes do not affect it. Reuse a calculator for
energy/angle scans, and construct a new calculator for another input or
coefficient point. The namespace `griffin::prediction` is not a separate
engine or directory.

The low-level SM ownership repair follows `merged_griffin_v1.2`, without
importing its scattering formulas. `FV_SMLO`, `FV_SMNLO` and `FV_SMNNLO`
own their selected FA/SW pair through `std::unique_ptr`; construction no
longer creates lower-order pairs and then abandons them. Copies duplicate
the pair while preserving its LO/NLO/NNLO type. The SM input is still borrowed.

`matel` owns the four fixed form factors created by its numerical-value
constructor. Its object-based constructor instead borrows the supplied
form factors. Copies duplicate owned values or retain borrowed references;
moves transfer ownership. `setinput()` updates owned fixed form factors,
but does not modify external borrowed form factors. These rules are inherited
by `mat_SMNNLO` and `mat_SMaas`. Reassign moved-from FV or matrix-element
objects before using them again. Inputs and borrowed form factors must
outlive the objects that use them.

`psobs` has a virtual destructor so deleting an observable through a base
pointer runs the correct cleanup. This does not transfer ownership of its
borrowed SM input. Remaining input-class precautions are in
[known limitations](KNOWN_ISSUES.md#5-ownership-and-legacy-error-handling).

## Source responsibilities

| Location | Contents |
|---|---|
| `include/griffin/SMEFT.h` | Convenience entry header collecting the usual SMEFT/SM interfaces |
| `include/griffin/sm/ElectroweakSetup.h` | Both `EWInputScheme` and `ElectroweakSetup`; no Wilson coefficients |
| `include/griffin/smeft/SMEFTCalculator.h` | Calculator and finite-result types |
| `include/griffin/smeft/EWPO_smeft.h` | Observable selections |
| `include/griffin/smeft/xsc_smeft.h` | Scattering request types and local contact interface |
| Other `include/griffin/smeft/` headers | Coefficient storage, mass-solver diagnostics and direct shifts |
| `include/griffin/sm/`, `core/`, `numerics/` | SM, shared and numerical interfaces |
| `src/smeft/EWPO_smeft.cc` | Shared FA/SW evaluation and FA/SW-based FV |
| `src/smeft/SMEFTCalculator.cc` | Input coordination and complete scattering assembly |
| `src/smeft/xsc_smeft.cc` | Local four-fermion matching, not the full amplitude |
| Other `src/smeft/` files | Coefficient storage, finite MW solution and direct vertex shifts |
| `src/sm/`, `src/core/`, `src/numerics/` | SM physics, common implementation and numerical algorithms |
| `src/sm/expressions/`, `data/grids/` | Analytic fragments and retained interpolation tables |
| `examples/smeft_example.cc` | Public SMEFT workflow |
| `examples/sm/` | Retained SM demonstration programs |
| `tests/sm/` | Near-pole and manual/helper/calculator scattering checks |
| Other `tests/` files | API, SMEFT, installed-client and archived-result checks |

The FA/SW helpers are private calculator members. Particle-to-generation
lookup is private to `smeft::Input`; it selects coefficients, not a flavor
mixing model. `WMassSolution` owns the resolved inputs and convergence
diagnostics. Its header is installed because the calculator stores it by
value and exposes `massSolution()`; ordinary users do not construct it.

`tests/reference_scattering.h` is the only source-controlled header outside
`include/`. It is a test-only amplitude assembly, not installed or used by
the library. Its shared SM kernels make it an assembly check, not independent
validation of those kernels.

## Build and validation boundary

Clients link to `GRIFFIN::griffin` and need only the public `include/`
path. Source/data include paths are private. The selected expressions and
grids are compiled into the library; installed clients do not load them at
runtime. All grid files are retained, including those not selected by the
current build. Generated build/install products are not additional source
interfaces.

Typed labels are used throughout the SM and SMEFT implementations, not just
at the outer calculator interface. Experimental and complex-pole masses have
different `InputPar` labels, and including an input header cannot redefine
their meaning. Integer array indices remain an internal storage detail,
not a second public input convention.

The public cross section supports massless electron annihilation to muons,
taus and $u,d,s,c,b$. Exactly $s=M_Z^2$ is still rejected; the upstream
near-pole correction is not an implementation of every exact-point limit.
See [known limitations](KNOWN_ISSUES.md) and
[the regression report](REGRESSION_TEST_REPORT.md). Historical Jackson
fixtures check the earlier SMNNLO setup, shared-kernel checks test the current
assembly, and frozen upstream values test the corrected near-pole behaviour.
