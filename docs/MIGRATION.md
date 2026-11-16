<!--
File: MIGRATION.md
Purpose: Help users adapt existing client programs and coefficient inputs to the current package.
Authors: Jackson Wallace and Lisong Chen
Last updated: 2026-11-16
-->

# Updating an existing program

This is a practical upgrade checklist, **not a change archive**. Use it when
adapting an existing GRIFFIN program or importing coefficient inputs.
For a new program, start with the [user guide](USER_GUIDE.md).

This package is identified as `smeft_v1.1.1`; its numeric CMake version is
`1.1.1`, and its library target remains `GRIFFIN::griffin`.

## 1. Rebuild and use the public include paths

Choose a fresh build directory after moving or renaming the package:
CMake caches retain absolute paths. Use a fresh installation prefix too,
because installing over an existing prefix does not remove stale headers.

The package requires C++11 and CMake 3.20 or newer. Prefer the CMake library
target `GRIFFIN::griffin` rather than maintaining a list of source files.
For an installed package:

```cmake
find_package(GRIFFIN 1.1.1 CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE GRIFFIN::griffin)
```

Supply the installation prefix through `CMAKE_PREFIX_PATH`. When building
directly from source, add the package with `add_subdirectory()` and link
the same target.

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
