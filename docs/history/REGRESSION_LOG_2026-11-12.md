# Archived integration and regression log — through 12 November 2026

This is a historical record, preserved when the current documentation was
cleaned on 13 November 2026. Statements about retained engines, old include
paths and earlier benchmark agreement describe their dated stages, not the
present library. See the [current test report](../REGRESSION_TEST_REPORT.md)
and [documentation index](../README.md) for current instructions.

**Latest update:** 12 November 2026. **Result: PASS.** Original integration:
3 November 2026.

The current calculator uses **`mat_SMaas` for both the full SMEFT prediction
and its zero-coefficient SM reference**. This intentionally extends the SM
scattering content beyond the earlier `mat_SMNNLO` comparison with Jackson;
the current cross sections are therefore not claimed to equal those old
cross-section fixtures. The finite SMEFT prescription is retained, without
extracting its linear part.

All **943 benchmark inputs** pass two complementary checks: **22,475
historical numerical comparisons** against the unchanged Jackson archive
using a test-only reconstruction of the former SM content, and **28,290
comparisons** of the extended public output against a separately assembled
`mat_SMaas` calculation. Release and sanitizer builds each pass **11/11
test suites**. Dated entries below preserve the earlier integration and
cleanup results; their byte-identical scattering outputs refer to the
pre-extension stages, not to the current calculation.

## Mixed electroweak–QCD scattering extension — 12 November 2026

The shared amplitude helper now constructs `mat_SMaas` instead of
`mat_SMNNLO`. It includes both the additional pole-residue terms and
`resoffZaas()` off-resonance terms, evaluated with the appropriate finite
inputs separately for the full prediction and the zero-coefficient reference.
Contact amplitudes are still added exactly once before squaring the complete
amplitude. Jackson's finite Gmu mass solver and the FA/SW/FV constructions
are unchanged. The underlying `xscaas.cc` implementation was not modified.

An exact `alphaS == 0` return was added to `getalphas()` to avoid division
by zero in its iterative formula at the zero-coupling limit. This does not
change positive-coupling results. Tests check that the added scattering
terms vanish in this limit.

Verification:

- Release: **11/11 suites pass**, 13.75 seconds in the final rerun.
- AddressSanitizer/UndefinedBehaviorSanitizer: **11/11 suites pass**,
  81.22 seconds, with no reported address/undefined-behavior errors. Leak
  detection remains disabled for the previously documented legacy
  allocations; this is not a leak-freedom claim.
- A new amplitude test scans **980 points** across both input schemes,
  seven supported final-state flavors, five energies, two angles and seven
  zero/vertex/contact/combined coefficient configurations. It checks all
  four amplitude channels, the residue/background decomposition, full and
  reference predictions, and their finite cross-section difference.
  Additional cases cover the zero-strong-coupling limit.
- All **943 historical cases / 22,475 comparisons** pass with unchanged
  fixtures and tolerances, using the explicitly labelled test-only
  `mat_SMNNLO` assembly. All **28,290 extended-output comparisons** pass
  against the separate `mat_SMaas` assembly. These references share the SM
  kernels and the public finite mass solver: this is an independent
  assembly check, not independent validation of those kernels or a new
  set of Jackson-supplied extended benchmarks.
- Separate invariant checks compare masses, EWPOs and contact amplitudes
  directly with the archived values. All available complete before/after
  mass, width, FA, SW, FV and contact outputs are exactly identical across
  the 943 inputs. Scattering amplitudes and cross sections intentionally
  change.
- The nine original regression-harness self-tests remain, supplemented by
  twelve extended-harness tests. These catch omitted extensions, an
  unchanged reference baseline, incorrect finite differences, corrupted
  historical results and malformed/missing/nonfinite outputs.
- The installed-client test still compiles all **26 public headers**
  independently and runs the example and API tests using only the installed
  package. The exact-real-pole rejection and unsupported-process checks
  remain in place.

The full machine-readable test results are in
[results_2026-11-12_aas.json](../../tests/native_regression/results_2026-11-12_aas.json).
The numerical before/after study is in
[aas_impact_2026-11-12.json](../../tests/native_regression/aas_impact_2026-11-12.json).

### Size of the numerical change

For the existing benchmark inputs, the following are zero-coefficient
$d\sigma/d\cos\theta$ values for $e^+e^-\to\mu^+\mu^-$ in the Gmu scheme,
at $\cos\theta=0.5$. Cross sections are in pb; the last column is
$100(\sigma_{\rm new}-\sigma_{\rm old})/\sigma_{\rm old}$.

| Energy [GeV] | Previous `mat_SMNNLO` | Current `mat_SMaas` | Change [%] |
|---:|---:|---:|---:|
| 50 | 14.10596034 | 14.10210341 | -0.027343 |
| 90 | 448.26992780 | 447.59089166 | -0.151479 |
| 95 | 132.59554209 | 132.64568459 | +0.037816 |
| 120 | 9.77048132 | 9.77193638 | +0.014892 |

As one finite SMEFT example, setting only
$C_{\phi WB}/\Lambda^2=+0.01\,\mathrm{TeV}^{-2}$ at 90 GeV changes
$\Delta\sigma\equiv\sigma(C)-\sigma(0)$ from **-0.5688478552 pb** to
**-0.5638142269 pb**. Its magnitude decreases by **0.8849%**; no
linear projection is involved. This is the current public example's
coefficient strength.

Across all 943 scanned benchmark inputs, the largest relative change in a
**full cross section** is **0.15343%**. The largest relative change in the
**finite SMEFT difference**, measured against the absolute old difference,
is **6.2165%**: a combined vertex/contact case changes from -0.1258462224
to -0.1180230383 pb. That larger percentage uses a small, cancellation-sensitive
difference as its denominator; it is not a 6% change of the full cross
section. These are scan results, not bounds for arbitrary coefficients.

The pre-extension source snapshot and all 943 before/after output records
are preserved in `../diagnostics/aas_extension_20261112_D0YLtp/`, relative
to the merged repository root. The original Lisong and Jackson repositories
were not modified. Nothing was committed or pushed.

## Clean public layout and engine removal — 12 November 2026

The generic `PredictionEngine`, providers, order/accuracy metadata,
perturbative series, combination rules and their old examples/tests have
been removed. There is no optional generic engine left in this version.
`SMEFTCalculator` remains the single finite SMEFT evaluator; its C++ class
name and `griffin::prediction` namespace are unchanged.

Public headers now live under `include/griffin/`, compiled implementations
under `src/{core,sm,numerics,smeft}/`, analytic expression fragments under
`src/sm/expressions/`, and numerical tables under `data/grids/`. Current
documentation lives in `docs/`, with the earlier comparison in
`docs/history/`. The main example is `examples/smeft_example.cc`, and the
build/test helper is `scripts/run_checks.sh`. See [MIGRATION.md](../MIGRATION.md)
for the complete header mapping and removed interfaces.

The new `<griffin/SMEFT.h>` convenience header loads the relevant SM physics
definitions before the legacy input headers redefine mass macros. Missing
include guards were added to retained public headers. No SM or SMEFT physics
formula was changed. Independent source review verified that all **20
retained implementation bodies** differ from the snapshot only in include
paths, and all **24 numerical grids plus 13 expression fragments** are
byte-identical. Unused historical tables were preserved as well.

Verification of the reorganized tree:

- Fresh release build: **9/9 suites pass**, 10.33 seconds.
- Fresh AddressSanitizer/UndefinedBehaviorSanitizer build: **9/9 suites
  pass**, 32.65 seconds. No address/undefined-behavior errors were reported;
  leak detection remains disabled for the previously documented legacy
  allocations, so this is not a leak-freedom claim.
- All **943 archived cases / 22,475 scalar comparisons pass**, with
  unchanged fixtures and tolerances. Numerical summary:
  [results_2026-11-12_layout.json](../../tests/native_regression/results_2026-11-12_layout.json).
- Complete before/after JSONL outputs are **byte-for-byte identical**.
  Their common SHA-256 remains
  `614c6407df0f3fc402e39f6b92f9e6ebf7f6faa76d9f467f259f5b5e59e232bd`.
- All six retained examples run successfully. The renamed finite example
  retains its previous W mass, cross sections and FV values.
- A new isolated installation test compiles all **26 public headers**
  independently, including each twice. It then compiles, links and runs
  copies of the finite example and public-API test using only the installed
  package through `find_package(GRIFFIN)` and `GRIFFIN::griffin`. No source
  or data directory is available on that client's include path. This passes
  in both release and sanitizer builds.
- The replacement public-API tests retain setup validation and independent
  full-square/angular normalization checks. They additionally verify that
  umbrella/advanced-header combinations preserve the pole-mass prefactors
  in both SW and Delta r. Existing finite-physics, contact-alias, Gmu
  matching and numerical regression tests remain in place.
- The clean-API guard rejects retired framework symbols, old top-level
  implementation directories and private includes in public headers.

The prior source tree is recoverable from
`../diagnostics/layout_cleanup_20261112_ObbFpc/before_layout.tar.gz`, relative
to the merged repository root. The same directory holds before/after JSONL
records and the former build trees in `old_builds/`. Those relocated CMake
caches are historical evidence, not reusable current build directories.
The original Lisong and Jackson repositories were not modified. Nothing
was committed, pushed or installed outside local test prefixes.

## Interface and documentation cleanup — 12 November 2026

This intermediate cleanup is historical; the subsequent layout cleanup
removes the entire generic engine rather than retaining its SM-only use.

Removed the undefined `dr_SMNLO::compmw()` declaration from `deltar.h`.
That class evaluates Delta r; the existing SM input solvers remain
`invalGmu`/`SMvalGmu`. The finite SMEFT mass solver is unchanged.

Removed `CombinationRule::linearEFT()`, its enum entry and engine branches,
the interference-only model cross-section path, and the now-unused
`linearInterference()`/`angularInterference()` helpers. There are no
compatibility aliases. All model cross-section accuracy requests to the
generic engine now raise an explicit SM-only error; none silently returns
a reference prediction in place of a requested correction. SM-only requests continue to work even
when a model provider is attached. Generic pseudo-observable additive and
full-model combination rules remain available, as do SM order comparisons.

The `prediction/` directory still contains the active `SMEFTCalculator`
interface. Shared observable names and kinematics now live in
`prediction/Observables.h`. The finite calculator includes this header and
`ElectroweakSetup.h` directly; it no longer pulls in the optional provider,
series, order-selection or combination framework. A preprocessor check and
the updated source guard verify that separation. A symbol inspection of
the rebuilt library found none of the removed interfaces.

Current documentation now distinguishes this active interface from the
optional engine, corrects ownership and local-fixture descriptions, and
records the removed APIs. Historical comparison results are not rewritten.
The input-scheme, SM-preset and contact-storage comments were also clarified.

Verification:

- Rebuilt release: **8/8 suites pass**, 2.74 seconds.
- Rebuilt AddressSanitizer/UndefinedBehaviorSanitizer: **8/8 suites pass**,
  24.10 seconds, with no reported address/undefined-behavior errors. Leak
  detection remains disabled because of the documented legacy allocations;
  this is not a leak-freedom claim.
- All **943 archived cases / 22,475 scalar checks pass**, with unchanged
  fixtures and tolerances. The numerical report is
  [results_2026-11-12.json](../../tests/native_regression/results_2026-11-12.json).
- All **943 complete before/after output records are byte-identical**,
  including masses, contact amplitudes, full amplitudes and cross sections.
  The common output SHA-256 is
  `614c6407df0f3fc402e39f6b92f9e6ebf7f6faa76d9f467f259f5b5e59e232bd`.
- The finite calculator implementation itself is byte-identical to the
  pre-cleanup snapshot; only its header dependencies changed. Both public
  examples retain their previous numerical output.
- Architecture tests cover rejection of every retained model-combination
  rule for cross sections in both input schemes, while preserving SM-only
  results, generic pseudo-observable arithmetic and input validation.

The pre-change API/documents snapshot and before/after JSONL evidence are
saved in `../diagnostics/api_cleanup_20261112_SUPkXY/`. Removed interfaces
can be recovered from `before_api_cleanup.tar.gz` there. No original Lisong
or Jackson source repository was modified, and no changes were committed
or pushed.

## Equivalent contact labels — 9 November 2026

`smeft::Input` now stores one coefficient per operator-identical group:

| Coefficient | Equivalent labels in the electron/muon example |
|---|---|
| Cee | `1122`, `1221`, `2211`, `2112` |
| Cll, first independent contraction | `1122`, `2211` |
| Cll, second independent contraction | `1221`, `2112` |

The same rules work for all three generations: Cll exchanges the two
currents; Cee additionally uses the four-dimensional chiral Fierz identity.
The two Cll rows are **not** identified. No extra Hermiticity relations or
changes to other coefficient tensors are introduced. These identities follow
the operator definitions and Fierz discussion in
[Grzadkowski et al., Sections 4 and 7](https://arxiv.org/html/1008.4884).

An ordinary setter replaces the shared value, including when setting zero.
Repeated equivalent assignments do not add contributions; the last value
wins. For genuinely independent raw entries from Jackson's representation,
sum the entries first and set each unique coefficient once. The existing
fixture import already does this and is unchanged.

Crucially, the Gmu matching combination now reads

$$
\Delta_\mu=C_{\phi l,11}^{(3)}+C_{\phi l,22}^{(3)}
-\frac12 C_{ll,\mathrm{unique}}^{1221}.
$$

The unique coefficient equals Jackson's raw `1221+2112` sum. Summing both
new alias getters would double count, so that old accessor sum was removed.
This is an input-representation change, **not a new W-mass prescription**.
The finite mass equation, its iteration and the full amplitude square remain
unchanged. Likewise, the genuine `2*Cll1111` electron-neutrino local-contact
component factor is retained; it is not a repeated-label sum.

All **943 archived cases / 22,475 scalar checks pass** without changing
fixtures or tolerances. Complete before/after runner output is byte-identical,
including contact amplitudes, masses and finite cross sections. The result is
[results_2026-11-09.json](../../tests/native_regression/results_2026-11-09.json);
before/after outputs are in `../diagnostics/contact_aliases_20261109/`.

Both rebuilt configurations pass **8/8 automated suites**: release in
2.54 seconds, AddressSanitizer/UndefinedBehaviorSanitizer in 23.06 seconds.
No address/undefined-behavior errors were reported. Leak detection remains
disabled for the known legacy allocating SM classes; no leak-freedom claim
is made.

New tests exercise all 81 input and getter labels, finding **45 Cll groups
and 36 Cee groups**, with no unintended identifications. They check repeated
assignments, replacement by a different value, zero, clear, copy, invalid
indices and nonfinite values. Contact checks cover muon and tau final states
and preserve the independent LL contractions and neutrino component factor.
End-to-end tests cover all **16 mixed label combinations**, plus isolated and
repeated aliases, in **both schemes at 50, 90, 95 and 120 GeV**. They compare
contact and full amplitudes, SM references, full cross sections and finite
shifts. Separate checks establish the Gmu factor of one half, no mass shift
from the `1122/2211` group, and explicit import of unequal raw coefficients.
Both public prediction examples still reproduce their previous outputs.

Migration details and examples are in
[SMEFT_IMPLEMENTATION.md](../SMEFT_IMPLEMENTATION.md). In particular, do not
load all redundant raw labels sequentially with the ordinary setters: a later
zero-valued alias would overwrite an earlier nonzero assignment. This change
is local; the original Lisong and Jackson repositories were not modified.

## Direct-shift and W-mass cleanup — 8 November 2026

The actual functions in `smeft/ff0_smeft.h` are now
`directZAxialShift` and `directZVectorShift`, with integer and typed-fermion
overloads. Neither takes an input-scheme argument: each evaluates only the
explicit vertex change at the supplied mass. The former names are removed,
not retained as aliases.

The old linear W-mass helpers and their dependent SMEFT observable classes,
provider, linear scattering predictor and examples have been physically
removed. `smeft/xsc_smeft` now supplies local contact matching only. The
generic SM/provider engine remains, but has no supplied SMEFT provider and
does not implement a SMEFT mass response.

Jackson's finite mass equation and its dimension-six input-matching term
remain unchanged. That matching term drives the finite equation; it is not
a separately appended linear response of an observable. The coefficient
container, contact normalization, FA/SW insertion, FV construction and full
cross-section square are unchanged.

Verification specific to this cleanup:

- Fresh release build: **8/8 suites pass** (1.49 seconds). Fresh
  AddressSanitizer/UndefinedBehaviorSanitizer build: **8/8 suites pass**
  (17.03 seconds), with no reported address/undefined-behavior errors.
  Leak checking remains disabled for the pre-existing allocating SM classes;
  this is not a leak-freedom claim.
- All **943 archived cases / 22,475 scalar checks pass**, with the same
  maximum errors listed below. The machine-readable result is
  [results_2026-11-08.json](../../tests/native_regression/results_2026-11-08.json).
- The original pre-cleanup executable and freshly rebuilt executable were
  each run on the complete fixture set. Their **entire JSONL outputs are
  byte-for-byte identical**, including fields beyond the archived assertions.
- The fixture SHA-256 is unchanged. No expected physics numbers were
  adjusted to make the tests pass.
- Independent tests now check fixed-mass direct vertex formulas for all
  supported flavors, contact matching, and the finite Gmu equation and
  residual. They no longer depend on the removed linear prescription.
- A new source-architecture test rejects the retired implementations and
  symbols. A separate symbol inspection of the rebuilt static library found
  the two new direct-shift functions and none of the retired response APIs.

The complete pre-cleanup source snapshot is recoverable from
`../diagnostics/merged_cleanup_20261108/before_cleanup.tar.gz`. Old generated
release/sanitizer build trees were moved into that same diagnostics directory,
outside the merged source; both current configurations were rebuilt afresh.
The before/after JSONL records and pre-cleanup summary are saved there too.
The original Lisong and Jackson source folders were not modified.

## Public class rename — 6 November 2026

The actual public classes and their header files were renamed to
`ElectroweakSetup` and `SMEFTCalculator`; no compatibility aliases or old
headers are retained. Library callers, examples, tests, CMake source lists,
and current documentation use the new names. Migration instructions are in
[USER_GUIDE.md](../USER_GUIDE.md#sm-accuracy-and-public-class-names).
The numerical prescription, coefficient normalization, and benchmark
fixtures were not changed.

After rebuilding the library and all examples, **7/7 suites pass** in both
release and AddressSanitizer/UndefinedBehaviorSanitizer builds. Both runs
pass all **943 cases and 22,475 numerical comparisons**. Leak detection was
disabled as in the original sanitizer run; no new leak-freedom claim is made.
The updated example reproduces its earlier cross-section, mass, and FV
outputs. The release numerical summary is saved as
[results_2026-11-06.json](../../tests/native_regression/results_2026-11-06.json).

The class rename does not add an SM-order selector. `SMEFTCalculator` still
uses the validated fixed higher-order prescription. The generic
engine's `Accuracy` controls do not apply to it.

## Local version and provenance

- Location: `/Users/lisongchen/Documents/griffin_jakcson/merged_griffin`.
- Local branch: `codex/merged-native-smeft`; integration changes are present
  in its working tree and have not been committed or pushed.
- Lisong base: `d5357a2bbed0aff09375aeb8da32a5d419cfbacc`, from
  `/Users/lisongchen/Documents/GRIFFIN`.
- Jackson finite prescription: `a25076d9cdcbe5171cd73e32d4f1656ee82174e8`, from
  `/Users/lisongchen/Documents/griffin_jakcson/jackson_smeft`.
- The original Jackson source remains clean. Lisong's original source retains
  exactly its pre-existing extra blank line in `smeft/ff0_smeft.cc`; neither
  original repository was edited. The isolated copy uses Lisong's committed
  source plus the new integration.

This is a tested integration of Jackson's calculation into Lisong's typed
input architecture, not a wholesale replacement by Jackson's raw-index input
classes. The SM numerical kernels are unchanged. The input coefficient
container now canonicalizes the equivalent Cll/Cee labels described in the
9 November entry while preserving existing canonical strengths. See
[USER_GUIDE.md](../USER_GUIDE.md) for the exact prescription, contact import
rule, and public API.

## Requested behavior verified

1. **Lisong inputs retained.** Named coefficient setters, generation enums,
   explicit Lambda in GeV and zero defaults are unchanged. Scaling all C by
   four and Lambda by two preserves the complete finite result.
2. **Jackson Gmu prescription.** The finite fixed-point equation re-evaluates
   SM Delta r at each shifted mass; the additive dimension-six mass term is
   outside the square root. No separate linear mass derivative is appended
   to vertices or cross sections. The returned mass and width convention
   matches Jackson's default solver.
3. **FV constructed from FA/SW.** Tests independently verify
   $F_V=F_A[(1-4|Q_f|\operatorname{Re}s_{\rm eff}^2)^2+I_f^2]$ at nonzero
   coefficients and recover the original `FV_SMNNLO` at zero coefficients.
   The SM imaginary-part term $I_f$ is retained. Jackson's old standalone
   tree/direct FV is intentionally not a comparison target.
4. **Full amplitude and full square.** Tests verify the angular square
   independently, including quark color multiplicity. Contact contributions
   occur once. A positive/negative contact-coefficient test explicitly
   detects the retained quadratic contribution, so an accidental linear
   projection cannot pass.
5. **Inputs and ownership.** Input mutation, temporary input/setup lifetimes,
   repeated calls, copy/move and assignment are tested for the snapshot-owning
   API. Original input-container coverage is retained. Following the November 8
   removal, tests of the retired linear predictors were replaced with direct
   vertex/contact checks and finite-equation checks; the numerical fixtures
   are unchanged. Removing those public APIs is an intentional source-interface
   change, documented in [USER_GUIDE.md](../USER_GUIDE.md#migration-from-the-removed-linear-smeft-api).

## Archived comparison coverage

| Set | Cases | What is checked |
|---|---:|---|
| Individual-coefficient scan | 520 | 13 coefficients, both schemes, four energies, five strengths |
| Contact convention controls | 80 | 40 raw-population / matched-strength pairs |
| Updated Tables 6/7 finite results | 120 | Vertex, contact and combined populations |
| Earlier Tables 6/7 finite results | 168 | Including smaller positive and negative strengths |
| Updated Table 4 native points | 40 | Electron FA/SW and available W-mass data |
| Earlier Table 4 finite points | 15 | Zero and positive/negative dimensionless C |
| **Total** | **943** | **22,475 scalar checks** |

The 888 cross-section cases separately check the full result, its zero-C SM
baseline, and their difference. This prevents agreement in a large full
cross section from hiding an incorrect small SMEFT shift. Complex amplitudes
are compared component by component in VV, AV, VA, AA order.

The process is massless unpolarized $e^+e^-\to\mu^+\mu^-$, with
$\cos\theta=0.5$, $\sqrt{s}=50,90,95,120$ GeV and results in pb.
Common inputs and all coefficient values are recorded in
[the fixture README](../../tests/native_regression/README.md) and
[manifest.json](../../tests/native_regression/manifest.json).
In particular, Lambda is 246 GeV and the main individual scan uses
$C=0.060516t$, $t=0,\pm0.001,\pm0.01\,{\rm TeV}^{-2}$.
Contact strengths are explicitly mapped to Lisong's convention; no native
Jackson result is changed to make the comparison agree.

### Historical pre-extension differences — release build

| Quantity | Maximum absolute difference |
|---|---:|
| Full differential cross section | $1.1369\times10^{-13}$ pb |
| Zero-coefficient cross section | $2.8422\times10^{-14}$ pb |
| Finite SMEFT change, checked separately | $1.0081\times10^{-13}$ pb |
| All archived W/Z masses and W widths | 0 |
| All archived FA and SW values | 0 |
| Full complex amplitude components | 0 |
| Contact-amplitude components | $2.4815\times10^{-24}\,{\rm GeV}^{-2}$ |

These are floating-point-level differences, including effects of converting
historical nb outputs to pb. “0” means numerical equality in this tested
build, not a promise of bitwise equality on every compiler/platform. The
largest relative difference in a full cross section is $4.32\times10^{-16}$.

As explicit cancellation-sensitive checks, in Gmu at 50 GeV and
$t=0.01\,{\rm TeV}^{-2}$:

| Single coefficient | Archived finite shift, pb | Merged minus archived, pb |
|---|---:|---:|
| $C_{\phi WB}$ | $4.115163969291985\times10^{-4}$ | $-7.08\times10^{-16}$ |
| $C_{\phi l,11}^{(3)}$ | $3.813860667556401\times10^{-5}$ | $+7.77\times10^{-16}$ |
| $C_{\phi l,22}^{(3)}$ | $3.813860667556401\times10^{-5}$ | $+7.77\times10^{-16}$ |

The machine-readable result is
[results_2026-11-03.json](../../tests/native_regression/results_2026-11-03.json).
Its fixture SHA-256 is
`956234c20636f83187518926da5aa0e8b097297d7a269f012eb520f3bad919f1`.
The manifest records a SHA-256 for each original archive. Fixtures contain
only native finite values; derivative/odd-projection comparison values and
the historical standalone FV are excluded.

### Independent fresh Jackson rebuild

As a further check, Jackson's unchanged current source at `a25076d9` was
rebuilt from scratch outside its source directory, and the original archived
individual-coefficient probe was freshly linked against that library. Its
**520 individual records and 80 contact-mapping records reproduce the original
archive byte-for-byte**. Thus the regression does not rely solely on trusting
an old library binary. No linear projection or source modification was used.

The exact stdout hashes were:

- Individual: `ba47569e137aa1d5b42f161cca804645ec7f76e056f2c7c83cad89ddc3e4476d`.
- Mapping: `78ae701f3001e7f32a712f7f3c13291be7e2908fe513a165ca266110a71aaf10`.

Historical evidence and the exact compiler command have been preserved in
`../diagnostics/layout_cleanup_20261112_ObbFpc/old_builds/build-jackson-reference/`
(relative to the merged repository root). Its `fresh_verification.json`
records that earlier independent rebuild. The archived helper script assumes
the old source layout and is not a current build instruction. The normal
943-case merged-code regression below remains self-contained and does not
need the original repositories or this historical build tree.

## Acceptance criteria

Each scalar must satisfy `abs(actual-expected) <= atol + rtol*abs(expected)`.

| Quantity | Absolute tolerance | Relative tolerance |
|---|---:|---:|
| Full/baseline/shift cross sections | $10^{-9}$ pb | $2\times10^{-11}$ |
| Amplitude real/imaginary components | $2\times10^{-15}$ | $2\times10^{-11}$ |
| Contact components | $2\times10^{-20}$ | $2\times10^{-12}$ |
| Pole/running W mass | $5\times10^{-10}$ GeV | $2\times10^{-12}$ |
| Pole Z mass | $5\times10^{-12}$ GeV | $2\times10^{-12}$ |
| W widths | $5\times10^{-11}$ GeV | $2\times10^{-11}$ |
| FA/SW (also FV in the extended assembly check) | $10^{-12}$ | $2\times10^{-11}$ |

The regression driver also rejects missing results/observables, duplicate or
unexpected identifiers, nonfinite values, and malformed arrays. Nine original
harness self-tests and twelve extended-harness tests check these failure
modes, including an incorrect small SMEFT shift while the full cross section
is correct. The original archived tolerances have not been relaxed.

## Reproduction

Environment used: macOS arm64, AppleClang 16.0.0.16000026, C++11. From the
local merged repository:

```bash
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DGRIFFIN_BUILD_EXAMPLES=ON
cmake --build build/release -j 4
ctest --test-dir build/release --output-on-failure
./build/release/examples/smeft_example
```

Latest November 12 release run: **11/11 suites passed**, 13.75 seconds. The
standard CTest run writes `build/release/native_regression_summary.json`. To retain a fresh
standalone numerical report:

```bash
python3 tests/native_regression/compare_aas.py --runner build/release/tests/native_fixture_runner --summary build/release/regression_summary.json
```

Sanitizer run:

```bash
cmake -S . -B build/sanitize -DBUILD_TESTING=ON -DGRIFFIN_BUILD_EXAMPLES=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined'
cmake --build build/sanitize -j 4
ASAN_OPTIONS=detect_leaks=0 UBSAN_OPTIONS=halt_on_error=1 ctest --test-dir build/sanitize --output-on-failure
```

Latest November 12 sanitizer run: **11/11 suites passed**, 81.22 seconds, without reported
address/undefined-behavior errors. **Leak checking was disabled** because the
suite deliberately exercises legacy allocating FV classes; this is not a
claim that all legacy code is leak-free. The new evaluator uses owned input
values and stack pole inputs rather than those legacy allocating constructors.

## Remaining scope boundaries

- The finite prescription is the only supplied SMEFT prediction path,
  through `SMEFTCalculator`. The old linear SMEFT APIs and the entire generic
  engine/provider/series/order framework have been removed. The retained
  SM classes can still be called directly, but there is no separate
  order-selection engine and no SM-order selector on `SMEFTCalculator`.
- The fixed scattering calculation now includes `mat_SMaas` for both full
  and reference predictions. The earlier `mat_SMNNLO` Jackson comparisons
  remain a separate historical validation stage; the extended predictions
  are checked against a separately assembled low-level calculation.
- No Bhabha/neutrino cross-section implementation was added. Cll/Cee now
  identify their exact operator-equivalent labels; this is not a general
  flavor-symmetry, Hermiticity or external tensor-format importer. Unique
  contact strengths and the explicit Jackson import normalization are retained.
- The inherited SM background has a zero-over-zero expression exactly at
  $s=M_Z^2$. The new API explicitly rejects that point; the exception is
  tested and none of the archived comparison points is affected. Pole EWPO
  calculations remain supported. Fixing that SM limit is a separate task.
- Nonphysical mass solutions, invalid solver options, nonconvergence,
  nonpositive FA, missing/nonfinite inputs, contact overflow and unsupported
  process/current requests are covered by exception tests.
- This validates the merged implementation against the agreed finite
  prescription. It does not establish a complete higher-order SMEFT theory
  calculation or validate arbitrary large-coefficient/extrapolated inputs.
