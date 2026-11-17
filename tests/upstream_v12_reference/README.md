# Upstream v1.2 amplitude reference

`points.inc` was generated on 2026-11-16 by `generate.cc`, compiled and linked
against the **unmodified upstream** repository at commit
[`ca2d5234d7712063b192e079bd8e3a4346fe9a7e`](https://github.com/lisongc/GRIFFIN/commit/ca2d5234d7712063b192e079bd8e3a4346fe9a7e).
It was not produced by the merged library or by a second implementation that
calls the merged library.

The generation environment was Apple Clang 16.0.0 on arm64 macOS. The upstream
library was built in Release mode (`-O3 -DNDEBUG -std=c++11`), and the generator
used `-O2 -std=c++11`. Floating-point values are emitted in scientific notation
with `std::numeric_limits<double>::max_digits10` digits after the decimal point.

There are 36 amplitude points: muon, up-quark and down-quark final states;
`(sqrt(s), cos(theta))` equal to `(50, -0.45)`, `(90, 0.3)`, and `(160, 0.6)`;
and the four VV, AV, VA, AA current combinations. All charged-fermion masses
are positive. `generate.cc` records the complete input point, including the
running-width masses that `SMval` converts to complex-pole masses.

Each row contains the integer final-state/current labels, energy, cosine,
and five real/imaginary pairs in this order:

1. `mat_SMNNLO2::result()`;
2. a separately constructed `mat_SMaas::result()`;
3. `mat_SMNNLO2::resoffZ2f()`;
4. `mat_SMNNLO2::resoffZfb()`;
5. `mat_SMNNLO2::resvpol()`.

These fixtures verify the port of the active upstream calculation. The named
commit excludes `compoffZfb()` from `resoffZfb()` and excludes `resaaa` from
`resvpol()`. The reference therefore does not establish complete off-resonance
NNLO accuracy. It does not cover the unresolved near-pole behavior, identical
initial/final flavors, bottom/neutrino final states, zero charged-fermion
masses, angular endpoints, or extrapolation beyond the supplied grids. No
failure or NaN in those domains is asserted as a required regression result.

To reproduce using a clean checkout of the exact commit, run from the merged
repository root (substitute the two absolute paths):

```sh
upstream_source=/absolute/path/to/unmodified/GRIFFIN
reference_build=/absolute/path/to/temporary/reference-build
git -C "$upstream_source" rev-parse HEAD
git -C "$upstream_source" status --short
cmake -S "$upstream_source" -B "$reference_build" -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build "$reference_build" -j 4
c++ -std=c++11 -O2 -I "$upstream_source" tests/upstream_v12_reference/generate.cc "$reference_build/libgriffin.a" -o "$reference_build/generate-v12-reference"
"$reference_build/generate-v12-reference" > "$reference_build/points.inc"
diff -u tests/upstream_v12_reference/points.inc "$reference_build/points.inc"
```

The policy-minimum option permits the upstream CMake file's old minimum
version on CMake 4.x. It does not modify the source calculation. Different
compiler/platform combinations can yield small floating-point differences;
`test_sm_nnlo2` uses a relative tolerance of `2e-10` and absolute tolerance
of `2e-17` for the frozen values, plus tighter same-build decomposition checks.
Regenerate reference values only against an explicitly reviewed upstream
commit, not merely to make a failing port test pass.

SHA-256 checksums of the original generation inputs/output:

```text
5aa1e0f311a9bc5e897ea71ce92ff45166478979988bca41e112801349e7cf53  generate.cc
906d5e1e70f903ce45289b58139faaa898d96a0ba42a74a14a016c65ef61db5a  points.inc
b1173dc0ba1398f7ac6e3a98b5508a861dda4fdf39117eecf605a315cece030f  upstream libgriffin.a
```

The archive checksum identifies the local build used for generation, not a
portable build requirement. The archive itself is not stored in this repository.
