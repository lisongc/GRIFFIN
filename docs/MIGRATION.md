<!--
File: MIGRATION.md
Purpose: Help users adapt existing client programs and coefficient inputs to the current package.
Authors: Jackson Wallace and Lisong Chen
Last updated: 2026-11-17
-->

# Updating an existing program

This is a practical upgrade checklist, **not a change archive**. Use it when
adapting an existing GRIFFIN program or importing coefficient inputs.
For a new program, start with the [user guide](USER_GUIDE.md).

This package is identified as `smeft_v1.1.2`; its numeric CMake version is
`1.1.2`, and its library target remains `GRIFFIN::griffin`.
It removes the old integer-label interfaces and repairs FV/matrix-element
ownership. Old-label clients require source changes and a complete rebuild;
this is not a drop-in-compatible patch for `smeft_v1.1.1` programs. Copy/move
constructors and the distinct value-based/object-based constructors remain.

## 1. Rebuild and use the public include paths

Choose a fresh build directory after moving or renaming the package:
CMake caches retain absolute paths. Use a fresh installation prefix too,
because installing over an existing prefix does not remove stale headers.

### SM object ownership

The FV/`matel` ownership repair changes class layouts, and `psobs` now has a
virtual destructor. Rebuild the library and **all** client object files;
do not link old binaries against the repaired headers/library.

Typed constructor calls are unchanged. FV classes now own their FA/SW
pair, and copies duplicate the helpers while preserving the correction order.
Custom FV subclasses must treat protected `fa`/`sw` as owning
`std::unique_ptr` members, not unmanaged pointers. Custom FA/SW subclasses
used as FV helpers must override the protected `clone()` function to copy
their actual derived type. Standard LO/NLO/NNLO classes already do this.

Numerical-value `matel` construction owns four fixed form factors;
object-based construction borrows them. Copies duplicate owned values or
retain external borrowed references, and moves transfer ownership. Reassign
moved-from objects before using them. Keep SM inputs and borrowed form
factors alive; `matel::setinput()` does not mutate external form factors.
The calculation formulas are unchanged from the published `smeft_v1.1.1` baseline.

### Build and installed headers

The package requires C++11 and CMake 3.20 or newer. Prefer the CMake library
target `GRIFFIN::griffin` rather than maintaining a list of source files.
For an installed package:

```cmake
find_package(GRIFFIN 1.1.2 EXACT CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE GRIFFIN::griffin)
```

Supply the installation prefix through `CMAKE_PREFIX_PATH`. When building
directly from source, add the package with `add_subdirectory()` and link
the same target.

The exported version file uses exact-version compatibility. Even a versioned
request for 1.1.1 without `EXACT` will not accept this 1.1.2 package. This avoids
claiming backward compatibility after removal of the old interfaces. An
unversioned `find_package(GRIFFIN CONFIG REQUIRED)` is still allowed, but does
not verify which interface version the client was written for.

Ordinary SMEFT programs can use just:

```cpp
#include <griffin/SMEFT.h>
```

For individual headers, use the following current locations:

| Needed interface | Include |
|---|---|
| SM inputs and amplitudes | `griffin/sm/SMval.h`, `griffin/sm/SMvalG.h`, `griffin/sm/xscnnlo.h`, `griffin/sm/xscaas.h` |
| Electroweak setup and scheme | `griffin/sm/ElectroweakSetup.h` |
| Wilson coefficients | `griffin/smeft/SMEFTInput.h` |
| SMEFT calculator | `griffin/smeft/SMEFTCalculator.h` |
| Observable selections | `griffin/smeft/EWPO_smeft.h` |
| Scattering request and contact amplitudes | `griffin/smeft/xsc_smeft.h` |
| Amplitude-to-cross-section function | `griffin/sm/tools.h` |

Programs using the original flat SM layout must qualify their includes:
for example, `"xscaas.h"` becomes `<griffin/sm/xscaas.h>`.
Shared headers are under `griffin/core/`; numerical headers are under
`griffin/numerics/`. Client include paths should point to `include/`,
not `src/` or `data/`. Grids are compiled into the library; no runtime
grid directory is required.

## 2. Check the input and result conventions

### Use typed labels only

The old short macros and integer-label overloads are removed, not retained
as compatibility aliases. Update both SM-only and SMEFT clients, then rebuild
the library and the clients together.

| Older call or label | Current spelling |
|---|---|
| `sm.set(al, value)` | `sm.set(InputPar::alpha, value)` |
| `sm.set(MZ, value)` after including old `SMval.h` | `sm.set(InputPar::ZMassExperiment, value)` |
| `sm.get(MZc)` | `sm.get(InputPar::ZMassComplex)` |
| `FA_SMNNLO(ELE, sm)` | `FA_SMNNLO(Fermion::electron, sm)` |
| `amplitude.setform(VEC, AXV)` | `amplitude.setform(Current::vector, Current::axial)` |
| `amplitude.setftype(ELE, MUO)` | `amplitude.setftype(Fermion::electron, Fermion::muon)` |
| `COMPPOLESCHEME` in `partzwidth` / `zwidth` | `WidthScheme::ComplexPole` |
| `RUNWIDTHSCHEME` in `partzwidth` / `zwidth` | `WidthScheme::RunningWidth` |

Previously `MW`, `MZ`, `GamW` and `GamZ` could mean different storage slots
depending on whether `SMval.h` had already been included. Do not translate
those four names blindly. Use `WMassExperiment`, `ZMassExperiment`,
`WWidthExperiment` and `ZWidthExperiment` for supplied running-width values;
use the corresponding `*Complex` names for complex-pole values. All are
members of `InputPar`. The physics kernels always use the complex-pole slots,
and `SMval`/`SMvalGmu` perform the same conversions as before.

Matrix-element constructors accept `Fermion` for the initial/final species
and `Current` for their currents. Form-factor constructors, `setftype`,
vertex functions, box functions and contact functions likewise accept the
corresponding typed labels only. Pass these enum values directly; converting
them with `index(...)` is not a substitute for migrating a call.
Wilson coefficient generations continue to use `smeft::Generation`.
Physical Cll/Cee aliases and the reordered Ceq/Cqe convention are unchanged.
The width-output labels are declared in `<griffin/sm/tools.h>`; omitted
width-scheme arguments still default to the complex-pole convention.

Use `SMval` with `EWInputScheme::AlphaMWMZ`, or `SMvalGmu` with
`EWInputScheme::GmuAlphaMZ`. Set Wilson coefficients through
`smeft::Input` and explicitly supply a positive Lambda in GeV.
Coefficients are dimensionless: an imported dimensional coefficient must
be converted to the convention $C_i/\Lambda^2$.

`SMEFTCalculator` solves Jackson's finite Gmu equation at the supplied
coefficients and evaluates the prediction at that mass. It solves the SM
case separately with all coefficients zero. In the alpha scheme the W
mass is supplied directly. The calculator uses fixed SM content, constructs
FV from FA/SW, and squares the complete scattering amplitude.

Read `total()` for the full result, `smResult()` for the zero-coefficient
SM result, and `modelCorrection()` for their finite difference.
For scattering, `ScatteringPoint` takes $s$ in GeV$^2$, not $\sqrt{s}$.
The default cross-section unit is GeV$^{-2}$; use `0.38937966e9` for pb.
The full input example is [smeft_example.cc](../examples/smeft_example.cc).

### Import equivalent contact labels once

Cll identifies `prst/stpr`; Cee identifies `prst/stpr/ptsr/srpt`.
Equivalent labels access **one value**. Setters assign rather than accumulate,
so importing independently specified raw contributions requires summing
their strengths before one assignment per unique operator.

For electron–muon inputs, Cll `1122/2211` and `1221/2112` are two separate
groups; Cee `1122/2211/1221/2112` is one group. Gmu matching uses half the
canonical Cll1221 strength once. Do not loop over every alias, including
zero-valued entries, after setting a nonzero coefficient: a later assignment
can overwrite it. See the [coefficient conventions](SMEFT_IMPLEMENTATION.md#contact-aliases-and-canonical-strengths)
and [import example](USER_GUIDE.md#contact-inputs-and-importing-jackson-benchmarks).

## 3. Compare like-for-like calculations and run the tests

The calculator uses `mat_SMaas`, including selected mixed electroweak–QCD
residue/background contributions and the upstream v1.1.1 correction.
Its zero-coefficient scattering result therefore need not equal a program
using only `mat_SMNNLO`. The [calculation specification](SMEFT_IMPLEMENTATION.md)
defines the retained content; [KNOWN_ISSUES.md](KNOWN_ISSUES.md) lists
unsupported processes and the exact-$s=M_Z^2$ restriction.

From the package root:

```bash
./scripts/run_checks.sh build/current
ctest --test-dir build/current --output-on-failure -L sm
```

SM demonstration sources are in `examples/sm/`; executable names are
unchanged. The SMEFT example remains `examples/smeft_example.cc`.
The checks in `tests/sm/` include the independent upstream near-pole
comparison and `test_scattering_equivalence.cc`, which compares manual
SM formulas with the cross-section helper and zero-coefficient calculator.
They run with `BUILD_TESTING=ON`, even when examples are disabled.
See the [test report](REGRESSION_TEST_REPORT.md) for coverage and tolerances.
