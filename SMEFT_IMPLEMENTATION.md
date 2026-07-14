# SMEFT implementation in GRIFFIN

## Scope and conventions

This version adds a tree-level, dimension-six SMEFT layer to the existing
GRIFFIN Standard Model calculation. The convention is

\[
\mathcal L_{\rm SMEFT}=\mathcal L_{\rm SM}
 +\sum_i \frac{C_i}{\Lambda^2}Q_i,
\]

where every `C_i` stored by `smeft::Input` is a real, dimensionless Wilson
coefficient and `Lambda` is supplied in GeV. The implemented SMEFT prediction
is truncated at one dimension-six insertion, namely
\(O(1/\Lambda^2)\). Dimension-six squared terms and dimension-eight terms are
not included.

The current implementation provides:

- tree-level SMEFT shifts of the neutral-current Z-fermion couplings;
- linear shifts of `SW`, `FA`, `FV`, `DeltaR`, and the derived W mass;
- local four-fermion contributions to massless
  \(e^+e^-\rightarrow f\bar f\);
- strict tree-SM/tree-SMEFT interference for differential cross sections;
- the `{Gmu, alpha, MZ}` and `{alpha, MW, MZ}` electroweak input schemes;
- both the original cumulative observable classes and the new order-aware
  prediction-engine interface.

SMEFT NLO corrections are not implemented. The supported cross-section
process domain is listed below and should be read before using the result.

## Input object and operator coverage

`smeft::Input` is deliberately independent of the Standard Model input. It
does not inherit from `inval`, does not compute a W mass, and does not know
which electroweak input scheme will be used. All coefficients start at zero;
`setLambdaGeV()` must be called before any prediction is evaluated.

The stored coefficients are:

| Category | Coefficients |
|---|---|
| Bosonic | `CphiWB`, `CphiD` |
| Two-index lepton currents | `Cphil1`, `Cphil3`, `Cphie` |
| Two-index quark currents | `Cphiq1`, `Cphiq3`, `Cphiu`, `Cphid` |
| Four-fermion | `Cll`, `Clq1`, `Clq3`, `Clu`, `Cld`, `Cqe`, `Ceu`, `Ced`, `Cle`, `Cee` |

`Generation::first`, `second`, and `third` represent the one-based Warsaw
basis labels. Two-index coefficients have all \(3^2\) slots and
four-fermion coefficients all \(3^4\) slots. The container does not impose
Hermiticity or a flavor ansatz.

`Cqe` is the canonical Warsaw/WCxf spelling. The old `Ceq[p,r,s,t]` API is a
current-reordered alias for `Cqe[s,t,p,r]`, not an additional coefficient.
For the right-handed charged-lepton contraction relevant here,
`Cee[1,i,i,1]` and `Cee[1,1,i,i]` are mapped to one canonical storage slot so
that a Fierz-equivalent contribution cannot be counted twice.

The implementation currently accepts real coefficients only. Complex
coefficients, CP-odd directions, running, matching, and RGE evolution are
outside its scope.

## Electroweak input schemes

The scheme is a property of the complete calculation and is stored in
`prediction::CalculationContext`.

### `{Gmu, alpha, MZ}`

Use `SMvalGmu`. The W mass in that object is the purely SM, loop-improved
solution derived from \(G_\mu\), \(\alpha\), and \(M_Z\). Wilson
coefficients never mutate this SM input object.

The muon-decay coefficient is centralized as

\[
\Delta_\mu=C_{\varphi l}^{(3)11}+C_{\varphi l}^{(3)22}
-\frac12\left(C_{ll}^{1221}+C_{ll}^{2112}\right).
\]

At tree level the dimension-six input transformation is

\[
\Delta r^{(6)}=\epsilon_6\left[
2\frac{c_W}{s_W}C_{\varphi WB}
+\frac{c_W^2}{2s_W^2}C_{\varphi D}+\Delta_\mu\right],
\qquad
\epsilon_6=\frac{1}{\sqrt2G_\mu\Lambda^2}.
\]

`derivedMWShiftD6()` converts this raw input shift into a displacement of the
loop-improved SM W-mass solution. It linearizes the existing SM fixed-point
equation,

\[
\delta M_W^2=
\frac{M_Z^2\,\delta x_6}{1-dF_{\rm SM}/dM_W^2},
\]

where the SM response derivative is evaluated numerically at the stored SM
solution. The returned `DerivedMWShiftD6` exposes the raw shift, response
factor, final \(\delta M_W^2/M_Z^2\), and \(\delta M_W\). A shifted W mass
must not be written back into the same input and then passed to SMEFT Z
vertices: that would apply the universal input shift twice.

### `{alpha, MW, MZ}`

Use a direct `SMval`, with W and Z masses supplied independently. The weak
angle is defined by \(c_W=M_W/M_Z\), and

\[
v_\alpha=\frac{2M_Ws_W}{\sqrt{4\pi\alpha}}.
\]

Here \(M_W\) is an input, so its SMEFT prediction shift is zero by definition.
`DeltaR` is instead a prediction of the muon-decay relation. The same
centralized \(\Delta_\mu\) and \(\Delta r^{(6)}\) implementation is used,
with \(\epsilon_6=v_\alpha^2/\Lambda^2\).

`CalculationContext` rejects `SMval` in the Gmu scheme and rejects
`SMvalGmu` in the direct-W scheme. Reference and model providers must also
share the exact same context and SM input object.

## Z-fermion currents and pseudo-observables

`smeft/ff0_smeft.cc` computes signed tree-level shifts
\(\delta a_f\) and \(\delta v_f\). The direct vertex part is evaluated with
`alpha`, `MW`, and `MZ` held fixed. In the Gmu scheme, the analytic derivative
of the ordinary SM tree vertex with respect to
\(x=M_W^2/M_Z^2\) adds the induced derived-W-mass contribution.

The flavor mapper in `SMEFTFlavor.h` is shared by leptons and quarks, so
vertex and contact calculations cannot silently use different generation
assignments.

The pseudo-observables are expanded explicitly to one SMEFT insertion:

\[
\delta F_A^{(6)}=2a_0\,\delta a,
\qquad
\delta F_V^{(6)}=2v_0\,\delta v,
\]

and, for a charged fermion,

\[
\delta s_{\mathrm{eff},f}^2=-\frac{1}{4|Q_f|}
\left(\frac{\delta v}{a_0}
-\frac{v_0\delta a}{a_0^2}\right).
\]

No absolute values and no products of two SMEFT shifts are used. `SW` is
undefined for a neutral fermion because its conventional definition divides
by \(|Q_f|\); both the old SMEFT class and the prediction engine throw a
`std::domain_error` for this request.

The legacy cumulative classes have the following meaning:

| Class | `result()` |
|---|---|
| `SW_SMEFTLO` | SM LO `SW` + linear SMEFT tree shift |
| `FA_SMEFTLO` | SM LO `FA` + linear SMEFT tree shift |
| `FV_SMEFTLO` | SM LO `FV` + linear SMEFT tree shift |
| `dr_SMEFTLO` | SM NLO `DeltaR` + linear SMEFT tree shift |

Each class also provides `resD6Tree()` for the correction alone. For new code,
the provider/engine interface is preferable because it keeps the SM reference
and model correction separate and lets the caller select the SM order.

## Four-fermion amplitudes

The local four-fermion matching first constructs the chiral coefficients
\((C_{LL},C_{LR},C_{RL},C_{RR})\). For a charged-lepton final state
\(\ell_i\ne e\), for example,

\[
\begin{aligned}
C_{LL}&=C_{ll}^{11ii}+C_{ll}^{1ii1},\\
C_{LR}&=C_{le}^{11ii},\\
C_{RL}&=C_{le}^{ii11},\\
C_{RR}&=C_{ee}^{11ii}.
\end{aligned}
\]

They are projected onto GRIFFIN's vector/axial convention as

\[
\begin{aligned}
4\Lambda^2\delta S_{VV}&=C_{LL}+C_{LR}+C_{RL}+C_{RR},\\
4\Lambda^2\delta S_{VA}&=C_{LL}-C_{LR}+C_{RL}-C_{RR},\\
4\Lambda^2\delta S_{AV}&=C_{LL}+C_{LR}-C_{RL}-C_{RR},\\
4\Lambda^2\delta S_{AA}&=C_{LL}-C_{LR}-C_{RL}+C_{RR}.
\end{aligned}
\]

Equivalent mappings for up quarks, down quarks, and neutrinos are centralized
in `smeft/xsc_smeft.cc`. A local contact is independent of `s`, so its exact
regular amplitude equals its Laurent `S` coefficient and its `S'` coefficient
is zero. The numerical amplitude adds the contact background once; it is not
also added separately as an `S` term.

`mat_SMEFTLO::coeffRD6Tree()` linearizes the two Z vertices:

\[
\delta R=(\delta z_i)z_f^{(0)}+z_i^{(0)}(\delta z_f),
\]

and deliberately omits \((\delta z_i)(\delta z_f)\). The complete tree
dimension-six amplitude is the shifted Z pole plus the local contact
background.

## Shared differential cross-section builder

`MasslessFermionCrossSectionBuilder` owns the process validation, channel
ordering, color factor, angular functions, and observable-level square or
interference. Both the prediction engine and the legacy SMEFT helper call this
one implementation.

The channel order is `VV`, `AV`, `VA`, `AA`. Given four SM amplitudes
\(M_0\) and four dimension-six amplitudes \(M_6\), the SMEFT result is
constructed only from

\[
\delta\sigma^{(6)}\propto
2\operatorname{Re}\!\left(M_0^*M_6\right),
\]

including the appropriate even and forward-backward angular terms. The
builder never forms \(|M_6|^2\). Quark final states receive the current fixed
color multiplicity of three.

The supported process is massless annihilation with an electron initial state:

```text
e+ e- -> mu+ mu-, tau+ tau-, u ubar, d dbar, s sbar, c cbar, b bbar
```

The following requests are rejected consistently for SM-only and SM+model
predictions:

- Bhabha scattering, which needs crossed s/t channels and scalar/pseudoscalar
  structures;
- every neutrino final state in the general builder;
- a non-electron initial state;
- nonphysical `s`, `cos(theta)`, or unit-conversion values.

The local neutrino contact coefficients can be inspected separately, but
\(e^+e^-\to\nu_e\bar\nu_e\) is not a complete cross section until the
t-channel W amplitude and its SMEFT corrections are implemented.

## Prediction-engine adapter

`SMEFTProvider` converts these SMEFT primitives into tagged perturbative
series:

- pseudo-observables receive `OrderTag::model(LO, 1)`;
- amplitudes receive the same tag for the tree dimension-six contribution;
- `ProviderOutputKind::AdditiveCorrection` tells the engine that these are
  corrections, not complete model predictions.

The intended combination is `CombinationRule::linearEFT()` with an accuracy
such as `Accuracy::SMNNLOPlus_plus_ModelLO()`. The returned
`PredictionResult` contains the selected SM reference, the linear SMEFT
correction, their total, and both underlying tagged series.

A complete user example is in `examples/testprediction.cc`. The architecture
and extension interface are described in `PREDICTION_ENGINE.md`.

## Verification

The automated tests cover:

- input defaults, scale validation, flavor indices, and coefficient aliases;
- both electroweak input schemes and the derived-W response;
- vertex, pseudo-observable, residue, contact, and amplitude formulas;
- sign reversal, coefficient doubling, and zero-input linearity;
- equality between the legacy and prediction-engine SMEFT LO results;
- the shared cross-section normalization and strict D6 interference;
- consistent rejection of Bhabha, neutrino, and invalid-scheme requests;
- rejection of unavailable SMEFT NLO requests.

Run the complete check and the public example with:

```bash
./scripts/run_prediction_architecture_checks.sh
```

Open limitations and the recommended implementation order are maintained in
`KNOWN_ISSUES.md`.
