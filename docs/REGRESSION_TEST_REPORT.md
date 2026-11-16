<!--
File: REGRESSION_TEST_REPORT.md
Purpose: Record current test coverage, numerical evidence, tolerances and rerun commands.
Authors: Jackson Wallace and Lisong Chen
Last updated: 2026-11-16
-->

# Regression test report

The tests cover the merged finite SMEFT prescription: Jackson's Gmu mass
solution, FA/SW-based FV, `mat_SMaas` scattering and the complete amplitude
square. They do not establish general precision Z-lineshape accuracy or
complete higher-order SMEFT accuracy. See [known limits](KNOWN_ISSUES.md).

## Latest verification: smeft_v1.1.1 packaging

The collaborator variant is `smeft_v1.1.1`, with numeric CMake package
version `1.1.1`. On 16 November, its pre-commit Release build with explicit
C++11 passed **13/13 suites** in 20.49 seconds. The installed-client test
requires `find_package(GRIFFIN 1.1.1 EXACT CONFIG REQUIRED)` and passed.
All six examples compiled. The 943 archived cases and the independent
near-pole and scattering-equivalence checks remain enabled.

The intended commit includes the complete public headers, implementation,
grouped tests/examples and all 24 grids. Personal `.vscode/` settings are
ignored, not distributed; generated builds, caches and archived historical
records are not part of the package. Version packaging changes no physics
formula or fixture values. The source is still a collaborator integration
candidate, not a general precision-accuracy certification.

## Test migration and package cleanup verification

On 16 November, a fresh AppleClang 16.0.0 Release build after the cleanup
passed **13/13 suites** in 19.18 seconds, including the installed-package
clients. All six examples compiled. The two SM scattering suites also pass
when selected with `ctest -L sm`.

The migrated equivalence test passes **24,850 comparisons at 700 points**;
the maximum relative difference is $4.244\times10^{-15}$ (tolerance
$5\times10^{-12}$). A separate fresh configuration with
`GRIFFIN_BUILD_EXAMPLES=OFF` still builds and passes this test.

All **943 complete public numerical records** are byte-identical before and
after this cleanup, with SHA-256
`a6f8aa7cc344fad581d679de2e987eeeaf40f9aa23348b9ecc49e5eee8406038`.
All **24 grid files** are retained and byte-identical. No test source or
fixture was discarded. The fresh build and saved outputs for this check are
under `/private/tmp/griffin_sm_tests_iK30TP/`; they are not package dependencies.

The preceding upstream-correction check passed 12/12 suites in both Release
and AddressSanitizer/UndefinedBehaviorSanitizer Debug builds, without reported
sanitizer errors. Sanitizers were not rerun for this file grouping; this
cleanup changes no production calculation.

All **943 archived cases** passed: **22,475** historical Jackson numerical
comparisons and **28,290** current extended-assembly comparisons. Fixture
values, integrity hash and numerical acceptance criteria were not changed.

## What the checks establish

| Check | Coverage and interpretation |
|---|---|
| Historical SMNNLO versus Jackson archive | 943 finite-result cases; preserves the original comparison rather than extracting a linear contribution |
| Public SMaas versus separately assembled SMaas | Same 943 cases; checks amplitude assembly, full/SM results and matching using shared SM kernels |
| Inputs and finite SMEFT prescription | Coefficient validation, Cll/Cee aliases without double counting, Gmu residual and separate zero-coefficient solution, direct shifts, FA/SW-based FV, full square and contact terms |
| Extended scattering | 980 points: both schemes, seven final flavors, five energies, two angles and seven coefficient configurations; four channels and zero-strong-coupling limit |
| Near-pole upstream agreement | 120 frozen down-quark values and 16 muon cross sections from an independently evaluated upstream v1.1.1 build |
| Scattering equivalence | 700 points: seven final flavors, 10 energies (10--190 GeV), five angles and two strong-coupling values; manual angular/chiral formulas versus the conversion helper and zero-Wilson calculator |
| Public/installed interface | Input ownership and lifetime, domain errors, 24 standalone headers, five installed clients, private-member access rejection and absence of retired interfaces |
| Comparison harness | Rejects malformed, missing, duplicate or nonfinite results, wrong small differences, missing mixed terms and corrupted archive values |

The separate SMaas assembly shares the calculation's kernels and resolved
mass inputs. It is **not** a new independent Jackson benchmark or derivation
of the SM kernels. Historical Jackson standalone FV is not a golden target;
the merged FV uses its own documented FA/SW construction.

The equivalence test checks shared-kernel API/algebra consistency, not
physical accuracy. At physical strong coupling, the zero-Wilson calculator
is compared with `mat_SMaas`, not `mat_SMNNLO`; at zero strong coupling
the extra mixed terms vanish and both can be compared.

The historical archive's common process is massless, unpolarized
$e^+e^-\to\mu^+\mu^-$ at $\cos\theta=0.5$ and
$\sqrt{s}=50,90,95,120$ GeV, with cross sections in pb. Inputs, coefficient
strengths and provenance are in the
[fixture guide](../tests/native_regression/README.md) and
[manifest](../tests/native_regression/manifest.json).

## Independent near-pole evidence

Both `src/sm/xscaas.cc` and `data/grids/aasoffz.grid` incorporate
Ayres Freitas's upstream v1.1.1 correction,
`a18f53830283a78155750a23a3594ada8de4c074`. The grid payload is byte-identical
to upstream after its provenance comment; the calculation matches upstream
after accounting for relocated includes and comments.

The fixed expectations in [test_near_pole.cc](../tests/sm/test_near_pole.cc)
come from a separate upstream evaluation on 14 November. They cover the
down-quark scan at 91.140--91.159 GeV and muon offsets of
$\pm10,\pm1,\pm0.5,\pm0.1$ MeV from the exact pole. A fresh standalone
diagnostic reproduced all **364 printed numerical entries** of the two
upstream output files exactly at their 15-digit output precision.
Linking the new test against the saved old library fails on its first mixed
cross section, demonstrating sensitivity to the old defect.

For $e^+e^-\to\mu^+\mu^-$, $\cos\theta=0.5$, with the test's inputs and
pole mass $91.1534806191828$ GeV:

| $\sqrt{s}-M_Z$ [MeV] | Old mixed result [nb] | Corrected result [nb] | Independent upstream [nb] |
|---:|---:|---:|---:|
| -1.0 | 0.298181100695 | 0.954936794294 | 0.954936794294 |
| +0.1 | 3373.996469469 | 0.955102955856 | 0.955102955856 |
| +1.0 | 1.784540713518 | 0.955172829841 | 0.955172829841 |

The +0.1 MeV SMNNLO-only result is 0.955029440547 nb; the corrected mixed
extension changes it by **+0.0076977%**. For the complete AA background,
Im $Q(h)$ at $h=0.03,0.01,0.003\ \mathrm{GeV}^2$ becomes respectively
$-4.9062\times10^{-13}$, $-5.4572\times10^{-14}$ and
$-4.9462\times10^{-15}$, rather than the old $6.738\times10^{-4}$.
The [optional diagnostic](../tests/diagnostics/near_pole_background.cc)
defines and prints $Q(h)$.

Across the 943 pre/post-port public records, masses, widths, FA/SW/FV and
contact amplitudes are identical. Scattering changes, including off-pole
changes from the corrected analytic functions and grids; the largest full
cross-section change in that scan is 0.701522 pb, not a universal bound.
No additional overlap subtraction was introduced. Exactly $s=M_Z^2$
remains explicitly rejected, and arbitrary smaller offsets or beam
convolutions have not been validated by these checks.

## Numerical tolerances

Historical/extended scalar comparisons require
`abs(actual-expected) <= atol + rtol*abs(expected)`;
complex arrays are compared componentwise.

| Quantity | Absolute tolerance | Relative tolerance |
|---|---:|---:|
| Full/SM/difference cross sections | $10^{-9}$ pb | $2\times10^{-11}$ |
| Amplitude components | $2\times10^{-15}$ | $2\times10^{-11}$ |
| Contact components | $2\times10^{-20}$ | $2\times10^{-12}$ |
| Pole/running W mass | $5\times10^{-10}$ GeV | $2\times10^{-12}$ |
| Pole Z mass | $5\times10^{-12}$ GeV | $2\times10^{-12}$ |
| W widths | $5\times10^{-11}$ GeV | $2\times10^{-11}$ |
| FA/SW; FV in the extended assembly | $10^{-12}$ | $2\times10^{-11}$ |

The independent near-pole test uses absolute tolerances of $2\times10^{-8}$
nb for cross sections, $5\times10^{-13}$ for real VV backgrounds and
$5\times10^{-14}$ for imaginary VV backgrounds. Its AA constancy check
uses a floating-point conditioning bound proportional to
$\epsilon M_Z^2/|s-M_Z^2|$, accounting for fused versus separately rounded
denominator subtraction; production arithmetic was not changed.
The scattering-equivalence test uses $5\times10^{-12}$ relative and
$10^{-26}$ absolute tolerance.

## Reproduce

From the package root, using fresh build directories:

```bash
./scripts/run_checks.sh build/release
ctest --test-dir build/release --output-on-failure -L sm
cmake --build build/release --target near_pole_background -j 4
./build/release/tests/near_pole_background
```

The SM tests in `tests/sm/` are part of `BUILD_TESTING`, independent of
`GRIFFIN_BUILD_EXAMPLES`. Generated numerical summaries belong in the build
directory; CTest writes `build/release/native_regression_summary.json`.
A separate full regression run can be saved with:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 tests/native_regression/compare_aas.py --runner build/release/tests/native_fixture_runner --summary build/release/regression_summary.json
```

For sanitizer checks:

```bash
cmake -S . -B build/sanitize -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined'
cmake --build build/sanitize -j 4
ctest --test-dir build/sanitize --output-on-failure
```

The fixture SHA-256 remains
`956234c20636f83187518926da5aa0e8b097297d7a269f012eb520f3bad919f1`.
Python 3's standard library is sufficient. No external comparison checkout
or historical run-record file is required.

## Historical records

Superseded reports and dated JSON run records were moved outside the
package to the local `../cleanup_archive/merged_griffin_v1.1.1_2026-11-16/`
archive (relative to the package root). All active fixtures, manifests,
comparison scripts, tests and grids remain. The upstream-port raw outputs
also remain locally under `/private/tmp/griffin_aas_port_YTMxPs/`.
Neither archive is a build or runtime dependency.
