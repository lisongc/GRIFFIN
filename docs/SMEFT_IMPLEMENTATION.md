<!--
File: SMEFT_IMPLEMENTATION.md
Purpose: Document Wilson-coefficient conventions, matching and the finite SMEFT calculation.
Authors (merged-version documentation): Jackson Wallace and Lisong Chen
Last updated: 2026-11-16
-->

# SMEFT implementation and input conventions

`prediction::SMEFTCalculator` combines Lisong's typed coefficients with
Jackson's finite mass prescription and full cross-section square. It uses
NNLO+ FA/SW and `mat_SMaas`, with no alternate linearized path or SM-order
selector. The imported `ca2d523` kernels include v1.1.1/v1.2 stabilization;
`mat_SMNNLO2` remains a separate SM class.

This import leaves Wilson normalization, finite Gmu matching and the
FA/SW/FV/full-square prescription unchanged. The
[user guide](USER_GUIDE.md#exact-prescription) contains the complete
prediction equations and usage; this document specifies inputs and matching.

## Operator normalization and scope

$$
\mathcal L_{\rm SMEFT}=\mathcal L_{\rm SM}
+\sum_i\frac{C_i}{\Lambda^2}Q_i.
$$

Every coefficient is real and dimensionless; Lambda is in GeV. Cll/Cee
aliases denote the total strength of one unique operator, not separate
Lagrangian summands. Direct vertices and contacts are tree-level
dimension-six insertions; explicit vertex/pole-input shifts are linear at
fixed resolved MW. Nonlinear dependence from the mass solution, pole
construction and full square is retained. This prescription does not supply
every higher-power SMEFT term or a complete dimension-eight calculation.

## Input object and operator coverage

`smeft::Input` is separate from `inval`. It does not select an electroweak
scheme or calculate MW. Coefficients start at zero, while
`setLambdaGeV()` must be called explicitly before evaluation.

| Category | Stored coefficients |
|---|---|
| Bosonic | `CphiWB`, `CphiD` |
| Lepton currents | `Cphil1`, `Cphil3`, `Cphie` |
| Quark currents | `Cphiq1`, `Cphiq3`, `Cphiu`, `Cphid` |
| Four-fermion | `Cll`, `Clq1`, `Clq3`, `Clu`, `Cld`, `Cqe`, `Ceu`, `Ced`, `Cle`, `Cee` |

`Generation::first`, `second` and `third` are one-based flavor labels.
The API accepts all $3^2$ two-index and $3^4$ four-index label tuples, with
equivalent Cll/Cee labels sharing canonical storage as described below.
Storage does not imply that every coefficient contributes to the supported
flavor-diagonal observables. The container does not impose Hermiticity or a
flavor ansatz. Nonfinite coefficients, invalid generation indices and
nonpositive/nonfinite Lambda are rejected.

`Cqe` is the canonical spelling. `Ceq[p,r,s,t]` is an alias for
`Cqe[s,t,p,r]`, not another coefficient. This pre-existing alias is unchanged.

Complex coefficients, CP-odd directions, running, matching between scales
and RGE evolution are not implemented.

## Contact aliases and canonical strengths

For all generation labels, the supported operator identities are

$$
Q_{ll}^{prst}=Q_{ll}^{stpr},
\qquad
Q_{ee}^{prst}=Q_{ee}^{stpr}=Q_{ee}^{ptsr}=Q_{ee}^{srpt}.
$$

The first identity exchanges the two currents. Cee additionally uses the
four-dimensional same-chirality Fierz identity. The corresponding labels
share one coefficient value in the input API:

| Coefficient | Equivalent labels in the electron–muon example |
|---|---|
| Cll, first independent operator | `1122`, `2211` |
| Cll, second independent operator | `1221`, `2112` |
| Cee, one independent operator | `1122`, `2211`, `1221`, `2112` |

The two Cll rows are independent. Do not extend the Cee Fierz identity
to Cll, which has weak-doublet contractions; see the [Warsaw-basis paper,
Eqs. (4.1)–(4.2) and Section 7](https://arxiv.org/html/1008.4884).
No extra Hermiticity, arbitrary permutations or other coefficient relations
are imposed.

Setters are assignments, not accumulators: the last write through any alias
wins, and every alias getter returns that value. Existing canonical inputs
retain their normalization. To import independently specified old raw
contributions, sum them first and set the canonical value once. For example,
old independent `Cll1221=a` and `Cll2112=b` become

```cpp
const auto one = smeft::Generation::first;
const auto two = smeft::Generation::second;
coefficients.setCll(one, two, two, one, a + b);
```

Two alias assignments retain only the last value. Do not import all 81 raw
tuples including zeros: a later zero can erase a nonzero canonical strength.
Start with zero-initialized `Input` or `clearWilsonCoefficients()` (preserves
Lambda), aggregate raw contributions, then assign each unique operator once.
No public importer performs this aggregation.

### Importing raw Jackson benchmarks

Aggregate Jackson's independently populated raw entries before assigning:

$$
C_{ll,L}^{1122}=C_{ll,J}^{1122}+C_{ll,J}^{2211},\qquad
C_{ll,L}^{1221}=C_{ll,J}^{1221}+C_{ll,J}^{2112},
$$

$$
C_{ee,L}^{1122}=C_{ee,J}^{1122}+C_{ee,J}^{1221}
+C_{ee,J}^{2211}+C_{ee,J}^{2112}.
$$

Cle entries are unchanged and unrelated coefficients are zero. Assign each
canonical sum once; setting an alias later to zero or to one raw piece
overwrites the total. Start with a new zero-initialized `Input`, or call
`clearWilsonCoefficients()` to clear coefficients while preserving Lambda.
Do not import by looping over all 81 raw tuples and assigning unused aliases
zero. No public importer performs this aggregation.

The Gmu contribution remains
$-\tfrac12 C_{ll,L}^{1221}
=-\tfrac12(C_{ll,J}^{1221}+C_{ll,J}^{2112})$;
no benchmark rescaling is needed. Already-canonical inputs are unchanged.

## Input scheme and mass treatment

The scheme belongs to `prediction::ElectroweakSetup`:

- `GmuAlphaMZ` requires `SMvalGmu`/`invalGmu`. `WMassSolution` solves
  Jackson's finite mass equation separately for the SM input values (zero
  coefficients) and the requested coefficients, re-evaluating SM Delta r at each
  iterate. The supplied SM input is not mutated.
- `AlphaMWMZ` requires a directly supplied `SMval`/`inval`. The supplied
  masses and widths remain unchanged.

The Gmu input-matching combination is

$$
\Delta_\mu=C_{\phi l,11}^{(3)}+C_{\phi l,22}^{(3)}
-\frac12C_{ll}^{1221}.
$$

The canonical Cll strength is included once. Since `1221` and `2112`
access the same value, summing both getters would double-count it. With an
imported value $C_{ll}^{1221}=a+b$, this gives the unchanged old raw
contribution $-(a+b)/2$.

The dimension-six input-matching term and finite solver conventions are
given in [USER_GUIDE.md](USER_GUIDE.md#exact-prescription). This is the only
SMEFT W-mass treatment. No linear displacement, SM-response derivative or
additional induced mass correction is added to the vertices.

`SMEFTCalculator` copies both inputs; caller mutations do not change an
existing calculation. Construct a new calculator for a new fit point.
Its `massSolution()` exposes the finite solution without mutating the caller's
SM object. For SM-only solving, use `SMvalGmu`/`invalGmu`.

Use `<griffin/SMEFT.h>` for the public interface; the
[header map](MIGRATION.md#header-path-mapping) describes individual headers.
`ElectroweakSetup` pairs a scheme with ordinary SM inputs; it does not
perform SMEFT matching.

## Direct Z-fermion shifts

`<griffin/smeft/ff0_smeft.h>` exposes

```cpp
double da = smeft::directZAxialShift(fermion, resolvedSM, coefficients);
double dv = smeft::directZVectorShift(fermion, resolvedSM, coefficients);
```

The particle argument is a `Fermion`. These functions return
only the direct tree-level dimension-six current shifts at the supplied
alpha, MW and MZ. They have no scheme selector and never change MW.
The calculator passes its already-resolved mass point to them.

Matching privately selects the generation/category for charged leptons,
their neutrinos, up/charm and down/strange/bottom; unsupported species are
rejected. This lookup is not a public API.

The calculator's private helpers in `src/smeft/EWPO_smeft.cc` use the direct
shifts in

$$
F_A^f=F_{A,\mathrm{SM}}^{f,\mathrm{NNLO+}}(M_W)+2a_{f,0}\delta a_f,
$$

$$
s_{\mathrm{eff}}^{2,f}
=s_{\mathrm{eff,SM}}^{2,f,\mathrm{NNLO+}}(M_W)
-\frac{1}{4|Q_f|}
\left(\frac{\delta v_f}{a_{f,0}}
-\frac{v_{f,0}\delta a_f}{a_{f,0}^2}\right).
$$

All terms use the same resolved mass. NNLO+ means `FA_SMNNLO`/`SW_SMNNLO`
also include selected three- and four-loop SM terms. At fixed MW, FA has no
explicit $(\delta a_f)^2$, and SW is not the exact ratio of shifted couplings.
Evaluation at finite $M_W(C)$ does not linearize its mass dependence.

FV uses the same FA/SW and the existing SM imaginary-ratio term, with the
[exact equation](USER_GUIDE.md#exact-prescription). Public observables and
scattering pole inputs share this evaluation. WMass, FA, SW and FV are
supported; neutral-fermion SW and finite DeltaR requests are rejected.
`EWPO_smeft.cc` contains these private helpers, `SMEFTCalculator.cc` assembles
scattering, and `xsc_smeft.cc` supplies the contacts below.

## Local four-fermion matching

`<griffin/smeft/xsc_smeft.h>` exposes the standalone function
`smeft::fourFermionContactAmplitude(initial, final, initialCurrent, finalCurrent,
coefficients)`, with `Fermion` particle and `Current` current labels. It returns the local reduced
contact amplitude, not a cross section or another SMEFT prediction
prescription. Its chiral coefficients are:

Every Cll/Cee entry below denotes its canonical total strength; aliases must
not be added again to these formulas.

| Final state of generation $i$ | $g_{LL}$ | $g_{LR}$ | $g_{RL}$ | $g_{RR}$ |
|---|---|---|---|---|
| Charged lepton, $i=2,3$ | $C_{ll}^{11ii}+C_{ll}^{1ii1}$ | $C_{le}^{11ii}$ | $C_{le}^{ii11}$ | $C_{ee}^{11ii}$ |
| Up-type quark | $C_{lq}^{(1)11ii}-C_{lq}^{(3)11ii}$ | $C_{lu}^{11ii}$ | $C_{qe}^{ii11}$ | $C_{eu}^{11ii}$ |
| Down-type quark | $C_{lq}^{(1)11ii}+C_{lq}^{(3)11ii}$ | $C_{ld}^{11ii}$ | $C_{qe}^{ii11}$ | $C_{ed}^{11ii}$ |
| $\nu_e$ contact only | $2C_{ll}^{1111}$ | $0$ | $C_{le}^{1111}$ | $0$ |
| $\nu_i$, $i=2,3$, contact only | $C_{ll}^{11ii}$ | $0$ | $C_{le}^{ii11}$ | $0$ |

The projection onto GRIFFIN's vector/axial convention is

$$
\begin{aligned}
4\Lambda^2\mathcal M^{4f}_{VV}&=g_{LL}+g_{LR}+g_{RL}+g_{RR},\\
4\Lambda^2\mathcal M^{4f}_{VA}&=g_{LL}-g_{LR}+g_{RL}-g_{RR},\\
4\Lambda^2\mathcal M^{4f}_{AV}&=g_{LL}+g_{LR}-g_{RL}-g_{RR},\\
4\Lambda^2\mathcal M^{4f}_{AA}&=g_{LL}-g_{LR}-g_{RL}+g_{RR}.
\end{aligned}
$$

The first current labels the initial state. The reduced contact is
independent of $s$, so no separate derivative function is needed. The
calculator adds it once to `mat_SMaas` before squaring, never again as a
Laurent term. Energy dependence remains through normalization and the other
amplitudes. The SM amplitude extension changes neither contact matching nor
the standalone FA/SW/FV and Gmu formulas.

## Cross-section scope

The function `differentialCrossSectionFromAmplitudes()`, declared in
`<griffin/sm/tools.h>` and implemented in `src/sm/tools.cc`, accepts four
amplitudes in the order `VV, AV, VA, AA` and applies the angular factors,
color multiplicity and normalization. It introduces no new amplitude
contributions; `mat_SMNNLO`, `mat_SMaas` and the separate `mat_SMNNLO2`
provide SM amplitudes. The finite calculator still selects `mat_SMaas`.
The calculator uses the function's full square, including SM–SMEFT interference,
contact-squared and vertex/contact products. Nonlinear terms also arise
through the pole construction and finite mass dependence, so the result is
not simply an SM cross section plus a squared dimension-six amplitude.

For every observable, `smResult()` is the independently evaluated
zero-coefficient SM result, `total()` is the requested-coefficient result, and
`modelCorrection()` is their finite difference. In Gmu the SM and
full results have separately solved W masses; in alpha both use the supplied
mass. The finite difference is not a derivative or an interference-only
quantity. See [the prediction equations](USER_GUIDE.md#exact-prescription)
for the complete residue/background decomposition. `smInputValues()` returns
the zero-coefficient electroweak inputs, and `smAmplitudes(point)` returns
the corresponding four amplitudes at the requested energy and angle.

Supported complete processes are unpolarized electron annihilation into
muons, taus and $u,d,s,c,b$. Massless external kinematics do not remove
masses from SM loop functions. Bhabha/neutrino final states, non-electron
initial states and scalar/pseudoscalar currents are rejected. Standalone
neutrino contacts are not complete cross sections; electron-neutrino
production also needs t-channel W exchange.

The conversion function requires finite kinematics, $s>0$, a positive unit factor
and $-1\leq\cos\theta\leq1$. Its default unit factor returns GeV$^{-2}$.
The full calculator requires $-1<\cos\theta<1$ because its amplitudes have
endpoint singularities. Exact $s=M_Z^2$ is rejected. Upstream stabilization
changes the old v1.1 behavior but does not establish precision lineshape
accuracy; the regular background is already included.

The additional SM terms are evaluated at resolved inputs, not new SMEFT
loop amplitudes. `mat_SMNNLO2` has its own disabled sectors, positive-mass
requirements and near-pole box growth; see
[KNOWN_ISSUES.md](KNOWN_ISSUES.md) and its
[SM usage section](USER_GUIDE.md#upstream-v12-sm-amplitude).

## Verification

Tests cover input/alias conventions, direct currents, contacts, finite Gmu
solutions, zero-coefficient SM limits, full squares, scale rescaling, lifetime behavior
and domain errors. Installed-package checks cover the interface used by external users.

The unchanged Jackson archive retains the earlier `mat_SMNNLO` baseline.
Tests reconstruct that baseline and separately check the current
`mat_SMaas` assembly. Shared kernels and resolved inputs make these assembly
checks, not independent SM derivations. See
[testing overview](../README.md#testing) for the independent
upstream references and recorded results.

```bash
./scripts/run_checks.sh
```

Run this command from the project root. The script builds the library and
examples, runs the tests and executes `smeft_example`. Historical comparison documents remain
available as records of earlier prescriptions, not as current API examples.
