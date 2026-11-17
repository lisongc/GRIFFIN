<!--
File: KNOWN_ISSUES.md
Purpose: Record supported-domain limits and unresolved issues.
Authors (merged-version documentation): Jackson Wallace and Lisong Chen
Last updated: 2026-11-17
-->

# Known issues and current scope

The supported SMEFT path is [`SMEFTCalculator`](USER_GUIDE.md), with fixed
SM content and a finite prescription. The retained nonlinear terms are
deliberate. See [migration](MIGRATION.md) for retired interfaces and input changes.

## 1. Corrected mixed background and remaining pole restriction

**The previously observed v1.1 mixed-background residual is resolved in the
imported implementation.** The updated `aasoffz.grid` is present, together
with the associated analytic changes. It should not be listed as an open
`mat_SMaas` defect. The changes are:

- v1.1.1 updates `aasoffz.grid`, removes the linear imaginary parts of the
  massless `zah0`/`zvh0` functions, and prescribes
  `zzhataas(s) = zzhataas(MZ²+10.01) * ((s-MZ²)/10.01)²`
  for `|s-MZ²| < 10 GeV²`.
- v1.2 linearly interpolates the complete one-loop pole-subtracted s-channel
  remainders between `s=MZ²(1±10⁻⁵)` for `|1-s/MZ²| < 10⁻⁵`, including
  the bosonic box terms in that subtraction.

The recorded residual diagnostic gives Im Q(0.01) about $-5.45\times10^{-14}$,
and the corrected `mat_SMaas` cross section is stable at the sampled
near-pole points. See the [dated audit](history/REGRESSION_REPORT_2026-11-16.md).

Exactly $s=M_Z^2$ still requires an explicit analytic limit: `resoffZaas()`
contains real-pole divisions, and `SMEFTCalculator` rejects that point.
Precision accuracy and smooth derivatives at interpolation boundaries have
not been established by the regression scans. These are remaining domain
and validation restrictions, not evidence that the old grid defect persists.
Pole pseudo-observables remain available.

To inspect the corrected background and the separate new-class diagnostic:

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/release --target near_pole_background near_pole_nnlo2 -j 4
./build/release/tests/near_pole_background
./build/release/tests/near_pole_nnlo2
```

## 2. Finite SMEFT scope and conventions

Vertices and contacts are real tree-level dimension-six corrections. The
explicit FA/SW shifts are linear at fixed resolved mass; the finite mass
solution, pole construction and full square retain their generated higher
powers. This is not a complete higher-order SMEFT expansion or a dimension-eight
calculation. See the [exact prescription](USER_GUIDE.md#exact-prescription).

Not implemented: complex/CP-odd coefficients, coefficient running or scale
matching, WCxf import/export, further tensor/Hermiticity constraints, a full
set of Higgs/dipole/charged-current observables, or uncertainty/covariance output.

Coefficient aliases identify unique operators; setters overwrite rather than
accumulate. Importers must aggregate independent raw contributions before
assigning a canonical strength. Follow the authoritative
[alias and normalization rules](SMEFT_IMPLEMENTATION.md#contact-aliases-and-canonical-strengths).

## 3. Process and observable coverage

The finite calculator supports massless electron annihilation to muons,
taus and $u,d,s,c,b$, with vector/axial channels. It rejects Bhabha scattering,
neutrino final states, non-electron initial states and collinear endpoints.
SM Bhabha support and local neutrino contacts do not supply the missing
crossed diagrams or complete SMEFT processes. SW for neutral fermions and a
finite DeltaR observable are also unavailable.

The cross section is a partonic two-body quantity without beam polarization,
real-radiation integration, ISR convolution, fragmentation or detector cuts.
Color multiplicity three is not an inclusive QCD correction. The caller
supplies the desired unit conversion.

## 4. SM accuracy is fixed in the SMEFT calculator

The calculator uses `FA_SMNNLO`, `SW_SMNNLO`, `mat_SMaas` and `dr_SMNNLO`,
with no SM-order selector. FA/SW contain NNLO plus selected three-/four-loop
terms (**NNLO+**) and their own mixed electroweak–QCD corrections.
`mat_SMaas` adds selected fermionic mixed terms in both residue and background,
for the full and zero-Wilson results at their separately resolved masses.
This is not complete off-resonance NNLO, nor a set of additional SMEFT loops.

Historical Jackson fixtures check the older amplitude assembly. Current
assembly tests share the production SM kernels. The additional
[upstream v1.1.1 regression](../tests/upstream_regression/README.md) uses
independently compiled scattering kernels at matched resolved inputs;
it does not independently validate the Gmu mass solver. Near-pole differences
are reported separately. Direct SM classes remain available under `sm/`.

### Upstream v1.2 SM class

`mat_SMNNLO2` is available from `<griffin/sm/xscnnlo2.h>`. Its `result()`
adds `resoffZ2f()`, `resoffZfb()` and `resvpol()` to the inherited construction.
The imported snapshot has these material limits:

- `resoffZfb()` returns `compboxfb(s)`; `compoffZfb(s)` is disabled.
- `resvpol()` returns `resaaz`; its computed `resaaa` term is disabled.
- The new fermion–boson box sector reports identical initial/final flavors,
  bottom quarks and neutrinos as unimplemented, then returns zero for that
  sector. A finite total can therefore have incomplete correction coverage.
- `resvpol()` takes logarithms of every charged-fermion mass from electron
  through top, requiring strictly positive masses. Interpolation routines can
  extrapolate outside their tabulated axes without an accuracy guarantee.

The new `resoffZfb()` contribution has observed near-pole growth. In this
snapshot it consists of the active `compboxfb(s)` box contribution; the
separate `compoffZfb(s)` term is disabled. For $e^+e^-\to\mu^+\mu^-$ at
$\cos\theta=0.5$, the recorded **differential cross section**
$d\sigma/d\cos\theta$ is:

| $s-M_Z^2$ (GeV$^2$, pole mass) | Corrected `mat_SMaas` (nb) | `mat_SMNNLO2` (nb) |
|---|---:|---:|
| $-0.1$ | 0.95498 | 0.96130 |
| $-0.01$ | 0.95497 | 1.21544 |

The second point differs by about **27%**. The
[component diagnostic](history/validation/results_2026-11-16_v12.json)
localizes the growing contribution to `resoffZfb()`, but does not establish
its underlying cause or a validity/exclusion window for all channels.
These are differential values at one angle, not integrated cross sections.

**This unresolved near-pole behavior is an important reason for keeping
`mat_SMNNLO2` separate from `SMEFTCalculator`.** The disabled terms and missing
bottom-quark box contribution are additional concerns: the calculator
supports bottom final states. It retains its existing `mat_SMaas`
prescription. Using `mat_SMNNLO2` with SMEFT would require resolving or
establishing the validity range of the added terms and validating their
combination with the finite SMEFT pole inputs and supported flavors.
The class name alone does not establish complete off-resonance NNLO coverage.

## 5. Fixed memory leaks and remaining low-level input errors

**The FV and `matel` memory leaks were fixed in `smeft_v1.1.2`.**

`SMEFTCalculator` has no identified ownership leak. It owns copies of its
SM/coefficient inputs and uses temporary pole form factors whose lifetimes
cover their evaluation. It checks required inputs, kinematics and nonfinite
results, reporting detected errors with exceptions that external users can catch.
These validation errors are separate from memory leaks.

`ElectroweakSetup` borrows its input until the calculator copies it at
construction. Keep that input alive while using the setup; references
returned by the calculator require the calculator itself to remain alive.
These are lifetime requirements, not unresolved allocation leaks.

The remaining issues below concern the **low-level SM input classes**:

| Location | Remaining behavior | Proposed repair |
|---|---|---|
| `inval::get()` / `set()` | Invalid accesses can call `exit(1)`, terminating the program. | Throw `std::out_of_range` for invalid keys and `std::domain_error` for invalid reads. |
| `inval` destructor | Deleting a derived input through an `inval*` does not have safe polymorphic destruction. | Add a virtual destructor; rebuild external programs. |
| Custom-capacity `SMval` / `SMvalGmu` / `invalGmu` inputs | Their `compute()` routines access storage directly; undersized containers can access memory out of bounds before calculator validation. | Require sufficient storage in these constructors, including construction from another `inval`, or remove custom-capacity construction for these fixed-layout inputs. |

Use the default-size input constructors. The calculator's checked input
copies and derived-W storage check protect its evaluation path; they cannot
repair an input object that was already constructed or used incorrectly.

Some older numerical routines also print an error and return zero for
unsupported cases, for example `C0` and `D0`. Such a finite value cannot be
detected by a nonfinite-output check. Their error paths need explicit failure
reporting as a separate change from memory ownership.

`SMval.h` no longer redefines mass macros. All input accesses now use explicit
`InputPar` keys, so the former pole/running-width header-order issue is fixed.
The distinct slots and external-user changes are listed in
[the migration guide](MIGRATION.md#typed-labels).
