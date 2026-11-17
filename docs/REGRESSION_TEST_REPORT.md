<!--
File: REGRESSION_TEST_REPORT.md
Purpose: Record current test coverage, numerical evidence, tolerances and rerun commands.
Authors: Jackson Wallace and Lisong Chen
Last updated: 2026-11-17
-->

# Regression test report

The tests cover the merged finite SMEFT prescription: Jackson's Gmu mass
solution, FA/SW-based FV, `mat_SMaas` scattering and the complete amplitude
square. They do not establish general precision Z-lineshape accuracy or
complete higher-order SMEFT accuracy. See [known limits](KNOWN_ISSUES.md).

## Version 1.1.2 packaging: verified 17 November 2026

The current package is designated `smeft_v1.1.2`, with numeric CMake version
`1.1.2`, to distinguish the typed-only interfaces and repaired FV/matrix-element
ownership from the published `smeft_v1.1.1` baseline. The README, current user
and migration instructions, and installed-package checks use the new version.
Upstream v1.1.1 correction credits and the historical records below retain
their original version numbers.

The exported CMake version file uses `ExactVersion`, not `SameMajorVersion`:
removing the old-label interfaces must not advertise backward compatibility.
The isolated installed-client check rejects a non-EXACT request for 1.1.1
before accepting an exact request for 1.1.2 from the same fresh installation.

A fresh AppleClang 16.0.0 Release C++11 build with examples enabled passed
**18/18 tests in 20.29 seconds**, including **9/9 installed clients**, the
24 standalone-header compilation checks and the new version-compatibility
checks. Logs and generated package metadata are under
`/private/tmp/griffin_v112_check_IxeiTX/release/`.

This version-designation step changes packaging/tests/documentation only.
C++ sources, public headers, grids and archived numerical fixtures are
byte-identical to the final reviewed implementation described below. The
earlier sanitizer, native-leak and byte-for-byte numerical checks therefore
cover the same calculation/ownership code; they were not rerun for the label
change. No folder/branch rename, commit or push was performed.

## Final pre-push review: verified 17 November 2026

The candidate was exported into a clean source directory containing all 122
tracked files and the four intended new source files. Fresh, separate builds
used AppleClang 16.0.0 on macOS arm64 with C++11 required. This checks the
uncommitted candidate as a complete package; it is not a checkout of a new
commit, and no commit or push was performed during this review.

| Check | Result |
|---|---|
| Fresh Release, examples enabled | **18/18 passed**, 20.54 seconds |
| Fresh Debug, AddressSanitizer + UndefinedBehaviorSanitizer | **18/18 passed**, 65.54 seconds; no reported sanitizer errors |
| Fresh Release, examples disabled | **17/17 passed**, 21.19 seconds |
| Isolated installed-package clients | **9/9 passed** in each configuration; all 24 public headers compile independently |
| Public calculator, separate SMNNLO and SMaas outputs | **943 records in each mode byte-identical** to the published baseline |
| Six original example programs | All complete outputs **byte-identical** to the published baseline |
| Native macOS leak scans | `testtools`, FV ownership, matrix ownership and the teaching example: **0 leaks, 0 leaked bytes**, all exit statuses zero |
| Package and documentation | All local documentation links resolve; no missing intended files; `git diff --check` clean |

The live GitHub branch still points to
`3277856af174cc887fc9a5d20010dbb5e5d9521b`. All 122 files in the baseline
source copy were checked against that commit, and rerunning its binaries
reproduced the stored baseline outputs. The candidate's 2,829 output records
and six original example outputs match those baseline files byte-for-byte.
The archived harness also passes all 22,475 historical and 28,290 extended
numerical comparisons. Fixtures, tolerances and all 24 grids are unchanged.

The intended package excludes `.vscode/`, build products and temporary
verification artifacts. No tracked source/test/grid file was deleted. The
four new source files must all be added when committing: the typed-input test,
both SM ownership tests and the memory-management teaching example. A commit
of tracked changes alone would omit required files.

The README now visibly credits Jackson Wallace and Lisong Chen and explains
the typed-only interfaces, ownership repair, teaching example and complete
client rebuild requirement. Final cleanup removed two trailing-space warnings,
corrected an obsolete ownership comment and enabled the teaching example in
the documented sanitizer command. User-added comment wording in `classes.h`
and the requested `axiaff`/`swtheff` names are retained.

No new blocker was found for a collaborator review branch within the
documented scope. A pre-existing low-level caveat discovered during review is
now explicit in [KNOWN_ISSUES.md](KNOWN_ISSUES.md#5-ownership-and-legacy-error-handling):
undersized custom-capacity SM input objects are unsafe; use the default
constructors. Existing exact-pole, process/domain, borrowed-lifetime and
precision-validation limits remain. This review is not cross-platform or
general physics-accuracy certification.

Local evidence is under `/private/tmp/griffin_final_v111_smC7SE/`: `source/`,
`release/`, `sanitize/`, `default/`, `native-leaks.txt`,
`final-numerical-summary.json`, the three `final-*.jsonl` files and the
example outputs. These are verification artifacts, not package dependencies.
The native scans needed permission to inspect their processes; macOS restricted
displayed memory contents but completed every scan successfully.

## Memory-management teaching example: verified 17 November 2026

The commented [test_memory_management.cc](../examples/sm/test_memory_management.cc)
prints the lifetimes of instrumented FA/SW helpers managed by the actual
`FV_SMLO` constructors and assignments. It checks copying, moving, replacement,
scope exit, a three-point SM input scan and cleanup after a simulated error.
The final count is **14 created, 14 destroyed, zero still alive**. Result
checks also verify that tracing and ownership operations do not change the
LO calculation. The counter is limited to the instrumented helpers; it is
not a count of all application allocations. No intentional leak is executed.

With `GRIFFIN_BUILD_EXAMPLES=ON`, this example is also registered with CTest
under the `sm`, `ownership` and `example` labels. The suite now contains 18
tests with examples enabled; the existing 17 tests remain available without
examples. No library source, physics formula, grid or fixture was changed
to add the teaching example.

| Check | Result |
|---|---|
| Release C++11 build and full suite | **18/18 passed**, 16.86 seconds |
| Debug AddressSanitizer + UndefinedBehaviorSanitizer build and full suite | **18/18 passed**, 57.28 seconds |
| macOS `leaks --atExit` on the new Release example | **0 leaks, 0 leaked bytes**, exit status zero |

The native leak scan required permission to inspect its launched process;
macOS restricted displayed memory contents but completed the scan. These
checks cover the executed paths, not every possible client program. Builds
and CTest logs are in the existing local `after-clean-build/` and
`sanitize-build/` directories under `/private/tmp/griffin_ownership_v111_gAl9aI/`.
See the [user guide](USER_GUIDE.md#headers-and-ownership) for run commands.

## SM object ownership repair: verified 17 November 2026

The ownership changes from `merged_griffin_v1.2` are ported to this v1.1.1
calculation: FV owns only its final-order FA/SW pair, numerical-value `matel`
owns four fixed form factors, and `psobs` supports destruction through base
pointers. Copies duplicate owned objects while preserving their concrete
types; borrowed SM inputs and external form factors remain borrowed. No
v1.2 scattering formula or grid is imported.

The new suites in [test_fv_ownership.cc](../tests/sm/test_fv_ownership.cc) and
[test_matel_ownership.cc](../tests/sm/test_matel_ownership.cc) directly count
allocations and releases and check copy/move behavior, destruction through
`psobs*`, owning/borrowing transitions and simulated allocation/clone failures.
Each FV order constructs and releases exactly two helpers, rather than
allocating unused lower-order pairs. Numerical-value matrix elements release
their four fixed objects; borrowing matrix elements do not delete external
objects. Copies survive destruction of the source and preserve LO/NLO/NNLO
behavior. Failed assignments leave the destination unchanged.

| Check | Result |
|---|---|
| Fresh Release C++11 suite | **17/17 passed**, including both ownership suites |
| Debug AddressSanitizer + UndefinedBehaviorSanitizer suite | **17/17 passed**, no reported sanitizer errors |
| Isolated installed-package clients | **9/9 passed**, including both ownership suites |
| Before/after calculator and separate SMNNLO/SMaas assemblies | **943 complete records in each mode byte-identical** |
| Before/after six example outputs and direct FV/matrix probe | **Byte-identical**, including the probe's 17-digit values |

macOS `/usr/bin/leaks --atExit` independently measured the Release binaries:

| Executable | Before repair | After repair |
|---|---:|---:|
| SM Z-width example `testtools` | 2 leaks, 64 bytes | **0 leaks, 0 bytes** |
| `test_fv_ownership` | New test | **0 leaks, 0 bytes** |
| `test_matel_ownership` | New test | **0 leaks, 0 bytes** |

The pre-fix example identifies the leaked FA/SW NNLO helpers. This measured
example is not a count of every possible leak in the old library. The native
checker needed permission to inspect the launched process; macOS restricted
displayed memory contents, but completed the leak scans. Its fixed-binary
runs all returned success.

The pre-repair source snapshot, fresh builds, numerical outputs and probe are
under `/private/tmp/griffin_ownership_v111_gAl9aI/`; `after-clean-build/` is the
validated Release build and `sanitize-build/` is the sanitizer build. These
are local verification artifacts, not package dependencies. No fixtures,
tolerances, grids or numerical expression sources changed. All earlier
typed-only input changes remain intact. Rebuild all client objects because
the ownership repair changes class layouts and the base-class virtual table.

The earlier passing AddressSanitizer/UndefinedBehaviorSanitizer runs did
**not** establish leak freedom. Those checks are distinct from explicit
allocation-lifetime checks or a successful leak-detector run; the new
allocation counters and native leak scans test the repaired paths directly. Borrowed-input
lifetime and the separate `inval` destructor restriction still apply.

## Typed-only input migration: verified 16 November 2026

The current source removes the parallel short-macro/integer-label interface
and uses `InputPar`, `Fermion` and `Current` throughout SM and SMEFT code.
Decay-width utilities also use `WidthScheme`. Calculation formulas,
Wilson-coefficient normalization, grids and frozen numerical expectations
are unchanged. Both mass-conversion input classes remain: they perform
different physical tasks, not alternative label interfaces.

A pristine local checkout of the previously published commit
`3277856af174cc887fc9a5d20010dbb5e5d9521b` supplied the before-change outputs.
Both versions were built with AppleClang 16.0.0, Release and explicit C++11.

| Check after removing the old labels | Result |
|---|---|
| Release C++11 build and regression suite | **15/15 passed**, all six examples built |
| Debug AddressSanitizer + UndefinedBehaviorSanitizer suite | **15/15 passed**, no reported sanitizer errors |
| Before/after public calculator output | All **943 complete records byte-identical** |
| Before/after separate SMNNLO and SMaas assemblies | All **943 records in each mode byte-identical** |
| Before/after six example programs | All six complete outputs byte-identical |
| Installed-package checks | Seven clients passed, including typed-only checks in both header orders; 24 standalone headers compile |

The three before/after modes contain **84,870 numeric scalar outputs** in
total. Byte-identical means that the full saved outputs match, not merely
that their differences fall below a numerical tolerance. The output SHA-256
hashes, identical for the two builds, are:

- Public calculator: `a6f8aa7cc344fad581d679de2e987eeeaf40f9aa23348b9ecc49e5eee8406038`.
- SMNNLO assembly: `adec2a30c11d1d66777e33c1d207b6b3b05fd05fd6c53c98ea87160d307f50a8`.
- SMaas assembly: `66bee8e93dbacec63020d08d673e795f0ed0a4107d9b9beec3120032c4509a4c`.

The added [typed-input test](../tests/test_typed_inputs.cc) verifies at compile
time that integer-label constructors, integer-label setters/getters and short-label
macros are absent, including inherited SM interfaces and low-level SMEFT
helpers. It is built in both SM-header orders and checks the experimental
and complex-pole mass slots and inline SM mass factors at runtime. Both
variants also run against an isolated installed package with C++11.

All 24 grids, the 943-case archive, its manifest and numerical tolerances
are untouched. No existing test program was removed. A normalized source
comparison also checked that the 13 SM expression fragments, core formulas
and mass-conversion formulas retain their arithmetic and input slots.
The local before/after artifacts are under
`/private/tmp/griffin_typed_baseline_08hYiv/`; the sanitizer build is under
`/private/tmp/griffin_typed_sanitize_build/`. These directories are not
package dependencies. See [migration instructions](MIGRATION.md#use-typed-labels-only)
for updating client calls; this intentionally breaks the retired interface.

## Previous verification: smeft_v1.1.1 packaging

The collaborator variant at that verification was `smeft_v1.1.1`, with numeric
CMake package version `1.1.1`. On 16 November, its pre-commit Release build with explicit
C++11 passed **13/13 suites** in 20.49 seconds. The installed-client test
then required `find_package(GRIFFIN 1.1.1 EXACT CONFIG REQUIRED)` and passed.
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
| SM object ownership | FV and matrix-element allocation balance, base-pointer destruction, independent copies, safe moves, borrowed-object preservation and allocation-failure cleanup |
| Public/installed interface | Input ownership and lifetime, domain errors, 24 standalone headers, nine installed clients, typed-only inputs in both header orders, private-member access rejection and absence of retired interfaces |
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
after accounting for relocated includes, comments and typed-label substitutions.

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
ctest --test-dir build/release --output-on-failure -L ownership
cmake --build build/release --target near_pole_background -j 4
./build/release/tests/near_pole_background
```

The SM tests in `tests/sm/` are part of `BUILD_TESTING`, independent of
`GRIFFIN_BUILD_EXAMPLES`. Generated numerical summaries belong in the build
directory; CTest writes `build/release/native_regression_summary.json`.
A macOS native leak check can be run separately on a Release executable:

```bash
/usr/bin/leaks --atExit -- ./build/release/tests/test_fv_ownership
/usr/bin/leaks --atExit -- ./build/release/tests/test_matel_ownership
```

Permission to inspect the launched process may be needed. Allocation-count
checks remain part of ordinary CTest and do not require this platform tool.
A separate full regression run can be saved with:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 tests/native_regression/compare_aas.py --runner build/release/tests/native_fixture_runner --summary build/release/regression_summary.json
```

For sanitizer checks:

```bash
cmake -S . -B build/sanitize -DBUILD_TESTING=ON -DGRIFFIN_BUILD_EXAMPLES=ON -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined'
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
