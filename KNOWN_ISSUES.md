# Known issues and next development steps

This document records limitations of the current SMEFT/prediction-engine
snapshot. It distinguishes confirmed implementation defects from deliberately
unsupported physics. Unsupported requests are rejected where possible; they
must not be interpreted as zero contributions.

## Priority summary

| Priority | Issue | Present behavior | Recommended next step |
|---|---|---|---|
| High | Legacy owned raw-pointer leaks | Repeated observable construction allocates objects that are never released | Replace ownership with RAII and add sanitizer tests |
| High | Fixed-order SM cross-section expansion | Cumulative amplitudes are squared, retaining terms beyond the named order | Build the cross section from explicit amplitude increments and order-by-order interferences |
| High | `e+e- -> nu_e nubar_e` t-channel | Complete cross section is rejected | Add a dedicated neutral+charged-current process builder |
| Medium | QCD order control | `qcdOrder` is stored but ignored by `Accuracy` | Add an explicit QCD accuracy axis and correctly tag provider terms |
| Medium | Bhabha support in the new engine | Generic builder rejects the process | Add a Bhabha-specific channel basis and adapter |
| Medium | Other neutrino final states | Generic builder rejects all neutrinos | Represent neutral-current neutrino couplings without `SW` and add a dedicated process path |
| Medium | Non-owning lifetime hazards | Contexts/providers retain raw addresses | Make lifetime/immutability explicit or use shared immutable state |
| Medium | Limited SMEFT accuracy and coefficient domain | Tree D6, real coefficients only | Add capabilities, complex coefficients, and higher orders deliberately |
| Medium | W-mass accuracy is outside order tags | Requested loop accuracy does not control how the context W mass was derived | Record derivation accuracy in the context/result provenance |
| Low | Legacy input failures terminate the process | `inval::get/set` call `exit(1)` | Replace with exceptions and typed validation |
| Low | No continuous integration | Tests run only when invoked locally | Add compiler, sanitizer, and platform CI jobs |

## 1. Definite raw-pointer memory leaks

There are two kinds of raw pointer in the code, and they should not be
confused.

### Confirmed owning-pointer leaks

`FV_SMLO` allocates `FA_SMLO` and `SW_SMLO` with `new` but has no destructor.
`FV_SMNLO` and `FV_SMNNLO` overwrite those pointers with newly allocated
higher-order objects, leaking the objects allocated by their base
constructors as well. None of the final allocations is released either.

The value-based `matel` constructor similarly creates four `psobsfix` objects
with `new`. The class has no destructor and therefore never releases them.

This is a real heap-memory leak: after the object is destroyed, the allocated
memory is no longer reachable, so a long fit, scan, service, or repeated test
can grow continuously.

Adding a few `delete` statements is not a safe complete fix because:

- the alternate `matel` constructor borrows pointers to caller-owned
  observables and must not delete them;
- default copy construction would copy owning addresses and could cause
  double deletion after a destructor is added;
- `psobs` has virtual methods but no virtual destructor, so deleting a derived
  object through a `psobs*` would be unsafe;
- the NLO/NNLO constructors currently replace base allocations, so ownership
  must be redesigned rather than patched at the end of the hierarchy.

Recommended fix:

1. Add a virtual destructor to `psobs`.
2. Replace owned sub-observables with values or `std::unique_ptr<psobs>`.
3. Represent borrowed dependencies separately, preferably with references or
   an explicit non-owning wrapper.
4. Delete copying or implement deep-copy semantics for owning classes; support
   moves where useful.
5. Avoid allocating an LO object if a derived constructor will immediately
   replace it.
6. Run the full test suite under AddressSanitizer and LeakSanitizer.

### Non-owning raw pointers

`CalculationContext`, `SMProvider`, `SMEFTProvider`, `PredictionEngine`, and
the new SMEFT observable classes store non-owning addresses. These do not leak
memory, because they do not allocate it. They can, however, dangle if the
referenced input, context, or provider is destroyed first. Mutation after
construction can also change a prediction without changing its identity.

For now, users must keep all referenced objects alive and treat inputs as
immutable while an engine is in use. A future API should use immutable shared
state or otherwise encode this lifetime contract in the type system.

## 2. SM differential cross sections are not strictly expanded by loop order

The engine obtains cumulative amplitude values

\[
M_{LO}=M_0,\qquad
M_{NLO}=M_0+M_1,\qquad
M_{NNLO}=M_0+M_1+M_2,
\]

and currently forms `squaredPrediction()` of each cumulative amplitude. Thus
the value labelled NLO contains \(|M_1|^2\), which is formally NNLO, and the
NNLO square contains products above NNLO. This follows the legacy cumulative
numerical convention but is not the strict fixed-order behavior suggested by
the new order tags.

The SMEFT correction does not have this defect: it is explicitly
\(2\operatorname{Re}(M_0^*M_6)\). The issue concerns the SM reference part of
the generic cross section.

The next implementation should construct observable increments directly:

\[
\begin{aligned}
\sigma_{LO}&=|M_0|^2,\\
\delta\sigma_{NLO}&=2\operatorname{Re}(M_0^*M_1),\\
\delta\sigma_{NNLO}&=|M_1|^2+2\operatorname{Re}(M_0^*M_2).
\end{aligned}
\]

Selected higher-order pieces included in `NNLOPlus` then need explicit tags
and a documented combination policy. Regression tests should compare both
the legacy cumulative mode and the new strict mode before choosing a release
default.

## 3. `e+e- -> nu_e nubar_e` needs a dedicated t-channel implementation

This process is not a simple neutral-current annihilation. In addition to
s-channel Z exchange it contains t-channel W exchange and their interference.
At dimension six, a complete result also requires the consistent shifts of
the charged-current vertices, W propagator/input parameters, and relevant
four-fermion structures.

The code contains local electron-neutrino contact coefficients for diagnostic
use, but `mat_SMEFTLO::resultD6Tree()` and the shared massless cross-section
builder reject the complete process. Returning only the Z pole plus contact
term would be physically incomplete and energy/angle dependent in the wrong
way.

Recommended design:

1. Add an `ElectronNeutrinoProcessBuilder` with explicit s- and t-channel
   amplitudes.
2. Use a chiral/helicity basis that represents the charged-current structure
   directly rather than forcing it through `SW`.
3. Implement the SM W-exchange amplitude first and validate it against a
   trusted tree-level calculation.
4. Add SMEFT W-lepton vertex shifts, input-scheme shifts, propagator effects,
   and local contacts with exactly one D6 insertion.
5. Test the angular dependence, the \(\nu_e\) versus \(\nu_\mu\) difference,
   sign reversal, and Ward/gauge-parameter independence where applicable.

## 4. QCD order is metadata, not a controllable accuracy

`OrderTag` has a `qcdOrder` member, but `Accuracy::accepts()` does not inspect
it and providers currently use the default value zero. Existing mixed
electroweak-QCD and higher-QCD terms are folded into the legacy `NLO` or
`NNLOPlus` cumulative result. Users therefore cannot request, for example,
electroweak NLO with or without a particular QCD order.

The fixed quark color factor of three is also not an inclusive QCD correction;
real radiation, jet definitions, and mass effects are not supplied by the
massless builder.

Recommended fix:

- replace the inert integer with a documented QCD-order enum or a general
  multi-index perturbative order;
- add a QCD selection to `Accuracy`;
- tag each SM contribution with its true \((\alpha,\alpha_s)\) content;
- decide how mixed terms and `NNLOPlus` subsets are selected;
- distinguish a partonic two-body result from an inclusive hadronic
  observable.

Until then, `NNLOPlus` should be described as the content of the existing
GRIFFIN high-order classes, not as independently selectable EW and QCD orders.

## 5. Process coverage gaps

### Bhabha scattering

The new massless annihilation builder correctly rejects
\(e^+e^-\to e^+e^-\). Bhabha scattering needs both s- and t-channel
contributions and the scalar/pseudoscalar structures used by the legacy
GRIFFIN Bhabha implementation. The existing legacy SM example does not make
the generic new-engine result complete.

Add a Bhabha-specific observable builder and an SM adapter before advertising
uniform new-engine coverage. SMEFT Bhabha support additionally requires a
dedicated mapping of crossed four-electron operators.

### Muon and tau neutrinos

The current generic builder rejects every neutrino final state, even though
\(\nu_\mu\) and \(\nu_\tau\) do not have the electron-neutrino W-exchange
diagram. Their neutral-current cross sections still cannot use the usual
`SW` pseudo-observable, because `SW` is undefined for \(Q_f=0\). A future
neutral-current process path should consume signed vector/axial or chiral
couplings directly and support these channels without manufacturing a
neutrino weak-mixing-angle observable.

### Massless approximation and observable definition

The shared builder treats final-state muons, taus, and all supported quarks as
massless and uses a simple color multiplicity. It has no beam polarization,
phase-space cuts, real photon/gluon radiation, ISR convolution, fragmentation,
or jet definition. The API's `unitConversion` is an untyped positive `double`,
so a unit mistake is also possible.

These limitations should remain explicit in public examples and eventually
be represented by process/observable configuration types.

## 6. Order and provider capability gaps

`Accuracy::smPlusModel()` can represent model NLO or `1/Lambda^4`, although
the current SMEFT provider cannot supply either. The engine catches these
requests at runtime. A provider capability query would allow earlier and more
informative validation.

For pseudo-observables, `fullModelDifference()` and `exactModel()` have test
coverage through a toy provider. For differential cross sections, only
`linearEFT()` is implemented. Additive and full-model cross sections need a
clear amplitude-versus-observable combination contract before enabling them.

The reference W mass is another special case. Its loop-improved derivation is
embedded in the `SMvalGmu` input object, while the provider emits the stored W
mass with an LO tag and the engine exempts it from the normal requested-order
check. This is practical but prevents the result from reporting the actual
derivation accuracy. The context should carry W-mass provenance explicitly.

## 7. SMEFT scope limitations

The current model is tree-level dimension six with real coefficients. It does
not provide:

- one-loop SMEFT contributions or SMEFT renormalization;
- Wilson-coefficient running and scale evolution;
- dimension-six squared or dimension-eight terms;
- complex Wilson coefficients and general CP violation;
- automatic Hermiticity/flavor-symmetry enforcement or WCxf import/export;
- a full set of Higgs, charged-current, dipole, or bosonic scattering
  observables;
- uncertainty propagation or theory-covariance output.

These are scope boundaries, not terms that may be silently set to zero. A
future provider should advertise its operator and order capabilities, and the
engine should reject a request outside them.

## 8. Legacy API and error handling

The older core still uses integer macros for parameters, fermions, and current
types, global `using namespace std`, and `exit(1)` for invalid or unset
`inval` entries. The new code adds typed enums and exceptions around its main
boundaries, but it still calls the legacy implementation internally.

Migration should be incremental:

- replace macro indices with typed accessors without changing stored layouts;
- make `inval::get/set` throw exceptions;
- validate finite masses, widths, couplings, and kinematics at construction;
- remove translation-unit-dependent macro behavior;
- preserve compatibility overloads only at the public boundary.

## 9. Testing, reproducibility, and release readiness

The current local tests cover SMEFT input, both EW schemes, matrix elements,
linearity, order selection, domain validation, and legacy/new agreement. The
repository still lacks automated CI, sanitizer jobs, a numerical baseline
artifact with documented tolerances, and independent physics validation for
the new process-level results.

Before calling the branch a stable major release:

1. Fix the ownership leaks and run ASan/LSan/UBSan.
2. Decide and document strict versus legacy SM cross-section order semantics.
3. Add CI for at least two compilers and build configurations.
4. Preserve the two-scheme Jackson comparison inputs and processed tables in
   a reproducible, reviewable test or validation package.
5. Mark Bhabha and neutrino coverage accurately in the release notes.
6. Have the input-scheme and contact conventions reviewed independently.
7. Version the public API and provide a migration example from cumulative
   observable classes to `PredictionEngine`.

The current branch is a strong development release for charged massless
annihilation and linear tree-level SMEFT studies, but these items should be
resolved or explicitly scoped before a production-quality `v2.0` claim.
