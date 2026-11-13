# Architecture and source map

Updated 14 November 2026. Public headers, implementation and numerical data
are separated below. The calculator uses GRIFFIN's NNLO+ SM form factors
and `mat_SMaas` scattering amplitude, with dimension-six corrections
inserted through the finite, unlinearized prescription.

## The public calculation

Start with `<griffin/SMEFT.h>` and [the complete example](../examples/smeft_example.cc).
Only four concepts are needed for an ordinary calculation:

| Object | Meaning |
|---|---|
| `SMval` or `SMvalGmu` | Ordinary SM inputs |
| `smeft::Input` | Wilson coefficients and Lambda |
| `prediction::ElectroweakSetup` | Which electroweak inputs define the calculation |
| `prediction::SMEFTCalculator` | The finite calculation at those inputs |

`ScatteringPoint` supplies the initial/final species, energy squared, angle
cosine and unit conversion when requesting a cross section. `PseudoObservable`
selects WMass, FA, SW or FV. `FinitePrediction` returns `reference()`,
`total()` and their difference, `modelCorrection()`.

The C++ namespace `griffin::prediction` is retained to avoid renaming the
working public classes. It does **not** mean that a `PredictionEngine` is
present. The generic engine, model providers, tagged series, `Accuracy`,
`LoopOrder` and combination rules have been removed.

## How a finite prediction is made

```text
SM inputs + Wilson coefficients + electroweak input scheme
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
          total / reference / total minus reference
```

In Gmu, `NativeInputPoint` solves Jackson's finite mass equation, evaluating
SM Delta r again at each mass iterate. In the alpha scheme it keeps the
supplied masses. No separate linear W-mass derivative is appended. FV is
constructed from the same FA/SW values and the existing SM imaginary-ratio
term. The exact equations are in the [user guide](USER_GUIDE.md#exact-prescription).

The explicit dimension-six FA/SW shifts are linear in the coefficients
when the resolved mass is held fixed. For example, the FA input contains
$2a_0\,\delta a$, not an additional $(\delta a)^2$. Nonlinear coefficient
dependence is retained when solving the finite mass equation, constructing
the pole amplitude and squaring the full amplitude. This does not supply
every possible higher-power SMEFT term. `modelCorrection()` is the finite
total minus reference, not a first derivative or a linear-only contribution.

`mat_SMaas` extends `mat_SMNNLO` with selected fermionic mixed
electroweak–QCD terms in both the complete pole residue and the off-resonance
background. Both input solutions above use this extension; it is not applied
only to the nonzero-coefficient result. Contacts are added once, after the
pole/background amplitude and before the full square. The SM kernels are
evaluated at each solution's mass; their resulting coefficient dependence
is part of the finite prescription, not a new set of SMEFT loop amplitudes.

The shared `FermionPairCrossSection` helper converts the four complete
`VV, AV, VA, AA` amplitudes into an unpolarized differential cross section.
It supplies the angular square, color multiplicity, normalization and
supported-process checks; it does not calculate the SM or SMEFT amplitudes.
Its current scope is electron annihilation to the supported charged fermion
pairs, with external-fermion masses neglected in the kinematics. The name
does not imply support for massive external particles, polarization, Bhabha
scattering or neutrino final states. Ordinary users call
`SMEFTCalculator::differentialCrossSection()` rather than constructing this
helper. Its 14 November rename changes neither the formula nor its scope.

`SMEFTCalculator` copies its input values at construction. Later changes to
the caller's objects cannot change an existing calculation. Construct a new
calculator for another coefficient/input point; reuse it for different
scattering energies and angles at the same input point.

`ElectroweakSetup` itself is a lightweight borrowing descriptor: its SM input
must remain alive while that descriptor is used. Once the calculator has
been constructed, it no longer borrows either caller input or setup.

## Folder responsibilities

| Location | Contents |
|---|---|
| `include/griffin/SMEFT.h` | Convenience entry header for the finite calculation |
| `include/griffin/SMEFTCalculator.h` | Calculator and finite-result types |
| `include/griffin/ElectroweakSetup.h`, `EWInputScheme.h`, `Observables.h` | Setup and shared request types |
| `include/griffin/smeft/` | Coefficients, input-solver diagnostics, direct-shift and contact interfaces |
| `include/griffin/sm/` | Existing SM input, form-factor and scattering interfaces |
| `include/griffin/core/` | Common input/fermion types, complex numbers and `FermionPairCrossSection` helper |
| `include/griffin/numerics/` | Advanced loop-integral, interpolation and special-function interfaces |
| `src/smeft/` | Finite calculator, finite input solution and SMEFT matching implementations |
| `src/smeft/detail/` | Private flavor helpers; not part of the installed API |
| `src/sm/` | SM physics implementation |
| `src/sm/expressions/` | Analytic C++ formula fragments included during compilation |
| `src/core/`, `src/numerics/` | Common implementation and numerical algorithms |
| `data/grids/` | Numerical tables included during compilation, not runtime input files |
| `examples/` | `smeft_example.cc` first; retained SM examples for lower-level usage |
| `tests/` | Public API, matching, finite-physics and archived-result checks |
| `docs/history/` | Historical numerical records, explicitly separated from current instructions |

The low-level `NativeInputPoint` header remains installed because the public
calculator exposes solver diagnostics through it. In contrast,
`SMEFTFlavor.h` is an implementation-only helper. These classifications are
about access, not different physical prescriptions.

The old short SM filenames are retained within the `sm/` directories so
that links to the original GRIFFIN code remain recognizable. Moving a file
does not rename its C++ classes or change its numerical formulas.

## Build boundary

Client programs see `include/` through the `GRIFFIN::griffin` target. They
do not need the repository root, `src/` or `data/` on their include path.
Implementation/data include paths are private to the library build.

The installed package contains the library, public headers and CMake package
metadata. The `.in` expressions and `.grid` tables are already compiled into
the library; client programs do not load them at runtime. Generated files
belong under `build/`, and are not source material to edit or distribute as
part of the public API.

Include `<griffin/SMEFT.h>` first when mixing it with advanced SM headers.
It loads the relevant inline SM definitions before the legacy input headers
redefine mass macros. Use typed `InputPar` keys in new code, but note that
these keys alone do not fix arbitrary legacy header ordering; see
[known issues](KNOWN_ISSUES.md#5-ownership-and-legacy-error-handling).

## SM content and scope

The fixed SM ingredients are `FA_SMNNLO`, `SW_SMNNLO`, `mat_SMaas` and
`dr_SMNNLO`. The FA/SW classes contain NNLO contributions plus selected
three- and four-loop terms: **NNLO+**, rather than uniformly NNLO everywhere.
The mixed electroweak–QCD amplitude extension does not make the entire
off-resonance calculation complete at NNLO. There is no SM-order switch;
the original SM classes remain available directly for SM-only comparisons.

The supported cross sections are massless electron annihilation to muons,
taus or $u,d,s,c,b$. Bhabha and neutrino final states are not supplied.
The calculator rejects exactly $s=M_Z^2$, but its unchanged `mat_SMaas`
background also has unresolved single-/double-real-pole behaviour nearby.
The finite complex resonant denominator does not regularize those background
terms. No automatic exclusion window or fallback is installed; this branch
is not validated for precision Z-lineshape use. See the
[full-background warning and diagnostic](KNOWN_ISSUES.md#1-unresolved-near-pole-scattering-background).

Historical Jackson fixtures use `mat_SMNNLO` and remain unchanged. Current
tests distinguish that archived baseline from the extended public amplitude;
see [the regression report](REGRESSION_TEST_REPORT.md).
