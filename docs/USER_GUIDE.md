<!--
File: USER_GUIDE.md
Purpose: Explain inputs, supported observables and usage of the merged GRIFFIN calculation.
Authors (merged-version documentation): Jackson Wallace and Lisong Chen
Last updated: 2026-11-16
-->

# User guide: finite SMEFT calculation

The public entry point is `SMEFTCalculator` in `griffin::prediction`.
It combines Lisong's typed inputs with Jackson's finite prescription:
NNLO+ SM form factors, `mat_SMaas`, direct dimension-six shifts and the full
amplitude square. **Finite** means evaluation at the requested coefficients;
generated nonlinear terms are retained without claiming every higher-power
SMEFT contribution. See [known issues](KNOWN_ISSUES.md) for scattering scope.

## Inputs and mass treatment

Use `SMval`/`SMvalGmu` with typed `InputPar` setters and a separate
`smeft::Input`. Real dimensionless coefficients default to zero; Lambda in
GeV must be supplied. In the Gmu scheme, **the W mass is solved iteratively
at the requested finite coefficients**, with SM Delta r reevaluated at each
iterate. The zero-coefficient SM mass is solved separately. In the alpha
scheme the supplied W mass is retained. The calculator uses corrected
`mat_SMaas`; [the class comparison](ARCHITECTURE.md#how-a-finite-prediction-is-made)
explains its relation to `mat_SMNNLO` and `mat_SMNNLO2`.

Replacing Lambda by `k*Lambda` and coefficients by `k*k*C` leaves predictions
unchanged. Contact labels share canonical strengths as explained
[below](#contact-inputs).

## Exact prescription

Let $M$ be the pole W mass and $M_Z$ the supplied pole Z mass. Define

$$
k=\frac{\pi\alpha}{\sqrt{2}G_\mu M_Z^2},\qquad
c_0^2=\frac12+\sqrt{\frac14-k},\qquad s_0^2=1-c_0^2,
$$

$$
\Delta_\mu=C_{\phi l,11}^{(3)}+C_{\phi l,22}^{(3)}
-\frac12C_{ll}^{1221},
$$

$C_{ll}^{1221}$ is the canonical total strength; $2112$ is its alias,
not another contribution to sum.

$$
\delta x_6=-\frac{k}{2\sqrt{1/4-k}}
\frac{1}{\sqrt{2}G_\mu\Lambda^2}
\left[2\frac{c_0}{s_0}C_{\phi WB}
+\frac{c_0^2}{2s_0^2}C_{\phi D}+\Delta_\mu\right].
$$

Here $\delta x_6$ is the dimension-six input-matching term in the mass
equation, not a linearized prediction for $M_W(C)-M_W(0)$. The finite Gmu
equation is iterated as

$$
M_{n+1}^2=M_Z^2\left[\frac12+
\sqrt{\frac14-k\bigl(1+\operatorname{Re}\Delta r_{\rm SM}(M_n)\bigr)}
+\delta x_6\right].
$$

The additive shift is **outside** the square root. The seed is 80 GeV;
the default stopping criterion is
$|M_{n+1}^2-M_n^2|\leq10^{-4}\,{\rm GeV}^2$. The returned mass is the last
evaluated iterate $M_n$, not the next candidate. `massSolution()` exposes
the iteration count and residual; options control tolerance and iteration
limit. Invalid domains are rejected. The Gmu solve also updates the W width
and running mass/width with Jackson's formula. In the alpha scheme the
supplied masses and widths are unchanged.

At resolved $M$, denote the direct vertex shifts by
$\delta a_f(C,M)$ and $\delta v_f(C,M)$:

$$
F_A^f=F_{A,\rm SM}^{f,\mathrm{NNLO+}}(M)+2a_{f,0}(M)\delta a_f(C,M),
$$

$$
s_{\rm eff}^{2,f}=s_{{\rm eff},\rm SM}^{2,f,\mathrm{NNLO+}}(M)
-\frac{\delta v_f/a_{f,0}-v_{f,0}\delta a_f/a_{f,0}^2}{4|Q_f|}.
$$

The displayed FA/SW matching terms are linear **at fixed $M$**: FA has no explicit
$(\delta a_f)^2$, and SW uses the displayed first-order ratio shift rather
than an exact ratio of shifted couplings. In Gmu, $M=M(C)$ is solved at the
requested coefficients; neither SM functions nor direct shifts are expanded
in $M(C)-M(0)$. No separate linear mass derivative is added.

`smeft::directZAxialShift` and `smeft::directZVectorShift` evaluate direct
shifts at the mass supplied by the calculator's automatic `WMassSolution`.

FV uses the same FA/SW values and GRIFFIN's SM imaginary-ratio term:

$$
F_V^f=F_A^f\left[(1-4|Q_f|\operatorname{Re}s_{\rm eff}^{2,f})^2+I_f^2\right],
\qquad
I_f=\frac{a_{f,0}\operatorname{Im}v_{f,1}
-v_{f,0}\operatorname{Im}a_{f,1}}{a_{f,0}^2}.
$$

$I_f$ is evaluated at the resolved mass. Retaining it preserves the
`FV_SMNNLO` construction when all Wilson coefficients vanish.

For each vector/axial channel,

$$
\mathcal A(C)=\mathcal A_{\tt mat\_SMaas}
\bigl(F_A(C),s_{\rm eff}^2(C);M(C)\bigr)+\mathcal A_{4f}(C).
$$

With $s_Z=M_Z^2-iM_Z\Gamma_Z$ and channel/angle labels suppressed,

$$
\mathcal A(C;s)=
\frac{R_{\tt aas}\bigl(F_A(C),s_{\rm eff}^2(C);M(C)\bigr)}{s-s_Z}
+B_{\tt base}\bigl(s;M(C)\bigr)
+\delta B_{\tt aas}\bigl(s;M(C)\bigr)
+\mathcal A_{4f}(C).
$$

$R_{\tt aas}$ is the **complete** `mat_SMaas::coeffR()` residue, not an
addition to `mat_SMNNLO::coeffR()`. $B_{\tt base}$ is
`mat_SMNNLO::resoffZ()` and $\delta B_{\tt aas}$ is
`mat_SMaas::resoffZaas()`. The extension supplies selected
$\alpha\alpha_s$ imaginary-residue and energy-dependent photon/Z
self-energy/mixing contributions. FA/SW already contain mixed EW–QCD terms.

Both full and zero-coefficient results use this assembly at their own
resolved inputs. **The SM calculation repeats the Gmu solve with all
coefficients zero**; it does not reuse the shifted W mass. In the alpha
scheme both use the supplied MW. At zero input $\alpha_s$, the additional
mixed EW–QCD scattering terms vanish.

The massless angular square includes SM–SMEFT interference, contact squares
and vertex/contact products. The finite mass solve and pole/FV construction
also generate nonlinear dependence. These terms are not a complete
higher-power SMEFT expansion.

### Reading the result

| Method | Meaning for observable $O$ |
|---|---|
| `smResult()` | $O(0)$, with zero coefficients and a separate Gmu solve |
| `total()` | $O(C)$, the finite prediction at the requested coefficients |
| `modelCorrection()` | $O(C)-O(0)$, including retained nonlinear dependence |

$O(C)$ means evaluation at $C$, not first order in $C$. The difference is
neither a derivative nor an interference-only quantity; `FinitePrediction`
has no inverse-scale order labels. `smInputValues()` returns the independently
resolved SM inputs, and `smAmplitudes(point)` returns its four amplitudes at
the same requested energy and angle.

## Contact inputs

For $e^+e^-\to\mu^+\mu^-$, the canonical-strength convention gives

$$
g_{LL}=C_{ll}^{1122}+C_{ll}^{1221},\quad
g_{RR}=C_{ee}^{1122},\quad
g_{LR}=C_{le}^{1122},\quad g_{RL}=C_{le}^{2211}.
$$

Cll identifies $1122/2211$ as one value and $1221/2112$ as another.
Cee identifies $1122,2211,1221,2112$ as a single value. Each value is the
canonical total strength of a unique operator. **Setters overwrite; they
do not sum.** Write each canonical value once. `Ceq` is an alias of reordered
`Cqe`. The complete [alias rules and raw-input conversion](SMEFT_IMPLEMENTATION.md#contact-aliases-and-canonical-strengths)
include the procedure for importing older Jackson benchmarks.

## Using the merged API

See [`examples/smeft_example.cc`](../examples/smeft_example.cc) for a complete
program. After constructing and filling its ordinary SM inputs:

```cpp
#include <griffin/SMEFT.h>

using namespace griffin;
// Construct and fill an SMvalGmu smGmu; see the complete example.
smeft::Input coefficients;
coefficients.setLambdaGeV(1000.0);
coefficients.setCphiWB(0.01); // dimensionless; C/Lambda^2 = .01 TeV^-2

prediction::ElectroweakSetup setup(EWInputScheme::GmuAlphaMZ, smGmu);
prediction::SMEFTCalculator calculation(setup, coefficients);
const Fermion initialFermion = Fermion::electron;
const Fermion finalFermion = Fermion::muon;
const double toPicobarns = 0.38937966e9;
prediction::ScatteringPoint point(initialFermion, finalFermion,
                                 90.0*90.0, 0.5, toPicobarns);
auto sigma = calculation.differentialCrossSection(point);
// sigma.smResult(), sigma.total(), sigma.modelCorrection(), all in pb
auto fv = calculation.pseudoObservable(prediction::PseudoObservable::FV,
                                      Fermion::electron);
```

`ScatteringPoint` contains initial/final fermions, $s$ in GeV$^2$,
$\cos\theta$, and an output-unit factor. `0.38937966e9` converts GeV$^{-2}$
to pb; the result is $d\sigma/d\cos\theta$. Reuse a calculator with new
points for an energy/angle scan. It owns snapshots of inputs and coefficients;
construct a new calculator when either changes.

For the alpha scheme, use `SMval` with `EWInputScheme::AlphaMWMZ`.
`ZMassExperiment`/`ZWidthExperiment` are running-width inputs that
`SMval`/`SMvalGmu` convert to `ZMassComplex`/`ZWidthComplex`.
The equations and exact-$s=M_Z^2$ restriction use the **pole mass**,
not the unconverted 91.1876 GeV input. `WMass` likewise returns the W pole
mass; the Gmu solver also supplies its running-width counterpart.

Build and run from the repository root:

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DGRIFFIN_BUILD_EXAMPLES=ON
cmake --build build/release -j 4
ctest --test-dir build/release --output-on-failure
./build/release/examples/smeft_example
```

Tests use bundled fixtures and Python 3's standard library, without network
access. Add `-DGRIFFIN_BUILD_LEGACY_EXAMPLES=ON` for the older SM demos under
`build/release/examples/legacy/`.

### Converting supplied amplitudes

For lower-level work, `<griffin/sm/tools.h>` converts four amplitudes at the
same point into $d\sigma/d\cos\theta$:

```cpp
#include <griffin/sm/tools.h>

// mVV, mAV, mVA and mAA are the four amplitudes at the same energy/angle.
const std::array<griffin::Cplx,4> amplitudes = {{mVV, mAV, mVA, mAA}};
double sigmaPb = griffin::differentialCrossSectionFromAmplitudes(
    amplitudes, griffin::Fermion::electron, griffin::Fermion::muon,
    90.0*90.0, 0.5, 0.38937966e9);
```

It applies the full angular square, color factor and unit conversion; it
does not compute amplitudes or add corrections. Supported final states are
muons, taus and $u,d,s,c,b$, with massless external kinematics. Require finite
inputs, $s>0$, $-1\leq\cos\theta\leq1$ and a positive unit factor.
The full calculator excludes endpoints because its amplitudes are singular
there. Ordinary SMEFT code can use `differentialCrossSection(point)` directly.

## SM accuracy and public class names

`ElectroweakSetup` and `SMEFTCalculator` are actual classes in
`griffin::prediction`. The calculator has **no selectable SM-order option**:
it uses cumulative `FA_SMNNLO`/`SW_SMNNLO`, `dr_SMNNLO` for Gmu, and
`mat_SMaas` for scattering. NNLO+ includes selected three- and four-loop SM
terms; it is not a uniform complete order for every amplitude component.

Current scattering goes beyond the historical Jackson `mat_SMNNLO` setup;
tests retain that archive and check the current assembly separately.
Shared-kernel checks are not independent theory calculations. See the
[testing overview](../README.md#testing) for provenance and results.

### Upstream v1.2 SM amplitude

Direct SM callers can use `mat_SMNNLO2` from `<griffin/sm/xscnnlo2.h>`
(upstream `ca2d523`). It inherits typed `Fermion`/`Current` constructors;
see [testmatel2.cc](../examples/testmatel2.cc) for borrowed FA/SW inputs.
The conversion function above accepts its four amplitudes. Numeric-value
construction owns its fixed form factors; object-based construction borrows
them. Keep the SM input and any borrowed form factors alive while using the
matrix element.

That example compares `mat_SMaas` with `mat_SMNNLO2` for
$e^+e^-\to d\bar d$ at $\cos\theta=0.5$ and $\sqrt{s}=30,50,\ldots,190$ GeV.
It prints the VV amplitude, off-resonance remainder, added vacuum-polarization
term and unpolarized $d\sigma/d\cos\theta$ in nb, using all four V/A channels
and the quark color factor. It demonstrates the active fermionic two-loop
off-resonance additions; it does not assert agreement with reference values.
The automated comparison is [`test_sm_nnlo2.cc`](../tests/test_sm_nnlo2.cc).

`result()` includes the inherited residue/background plus `resoffZ2f()`,
`resoffZfb()` and `resvpol()`. Do not add `mat_SMaas::result()` again.
In this snapshot, `compoffZfb(s)` is disabled in `resoffZfb()`, and
`resaaa` is disabled in `resvpol()`. The added fermion–boson box sector
returns zero with a warning for identical flavors, bottom quarks and
neutrinos. The class therefore does not provide complete NNLO coverage.

Supply strictly positive charged-fermion masses from electron through top:
`resvpol()` takes their logarithms. Grid extrapolation and near-pole box
growth need separate validation; see
[known issues](KNOWN_ISSUES.md#upstream-v12-sm-class).
`SMEFTCalculator` continues to use corrected `mat_SMaas`; this standalone
class introduces no SMEFT loop accuracy or new calculator option.

### Public headers and implementation

Use `<griffin/SMEFT.h>` for ordinary SMEFT work. Individual SM, observable,
scattering and matching headers are listed in the
[header mapping](MIGRATION.md#header-path-mapping).

Use `InputPar`, `Fermion` and `Current` labels throughout; integer-label
overloads and old macros are removed. Pole and running-width mass keys are
explicit and independent of header order. See the
[migration guide](MIGRATION.md#typed-labels) for the mass/width mapping.

## Scope and limitations

The calculator supports WMass, FA, SW and FV; neutral-fermion SW and finite
DeltaR requests are rejected. Its full scattering supports unpolarized
electron annihilation to muons, taus and $u,d,s,c,b$ with massless external
kinematics. Bhabha/neutrino final states, non-electron initial states,
scalar/pseudoscalar currents and endpoint angles are rejected.

Nonphysical mass solutions, nonpositive FA, missing/nonfinite inputs and
nonfinite outputs raise exceptions. These checks do not establish EFT or
interpolation validity. Exact $s=M_Z^2$ remains rejected; pole EWPOs remain
available. See [KNOWN_ISSUES.md](KNOWN_ISSUES.md) for accuracy, process,
near-pole and low-level API limitations.
