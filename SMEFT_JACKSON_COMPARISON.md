# SMEFT implementation comparison: current GRIFFIN and Jackson

## How to read this comparison

This document summarizes a comparison performed locally on 2026-07-14 between
the current GRIFFIN SMEFT implementation and Jackson's SMEFT implementation.
The numerical results were generated through the two implementations' public
interfaces and were not regenerated while preparing this summary.

Two comparisons answer different questions:

- **Native public comparison:** what each implementation returns through its
  own public classes and input objects.
- **Common-reference diagnostic:** whether the primitive formulas agree after
  both are deliberately forced to use the same complex-pole masses and a
  strict linear reconstruction.

The native results can disagree even when the primitive algebra agrees. That
is exactly what happens here.

For new GRIFFIN, every reported SMEFT number is structurally linear in a
dimension-six coefficient. For Jackson, the quoted public coefficient is the
odd projection

```text
[O(C=+0.001) - O(C=-0.001)] / 0.002.
```

This is a derivative only when the tested public result is smooth and linear
at zero. Because the tested Jackson implementation uses `fabs` and nonlinear
products for some observables, several values below are only odd projections.
A zero Jackson entry for `FA` therefore does not mean that the implementation
returns no change.

## Benchmark definition

The comparison uses both electroweak input schemes:

```text
Gmu scheme:    {Gmu, alpha, MZ}
Alpha scheme:  {alpha, MW, MZ}
```

The cross-section benchmark is

```text
process       e+ e- -> mu+ mu-
cos(theta)    0.5
Lambda        246 GeV
conversion    0.38937966e6 nb GeV^2
```

The vertex scenario sets `CphiD`, `CphiWB`, `Cphil1_11`, `Cphil3_11`,
`Cphie_11`, `Cphil1_22`, `Cphil3_22`, and `Cphie_22` to the common benchmark
value `0.060516 * scale`. The contact scenario sets
`Cll_1122=Cll_1221=c`, `Cee_1122=c/3`, `Cle_1122=c/2`, and
`Cle_2211=c/4`, with `c=0.060516*scale`. The combined scenario applies both.

## Executive conclusion

The new implementation and Jackson contain the same signed direct tree-level
Z-current formula. When both are evaluated with identical complex-pole masses,
their direct \(\delta a\) and \(\delta v\) coefficients agree to numerical
precision. Their raw `{Gmu,alpha,MZ}` dimension-six input transformation also
agrees.

The public implementations are nevertheless not equivalent. The main causes
are:

1. Jackson's compiled direct SMEFT current reads different W/Z mass slots from
   the ordinary SM tree current.
2. Jackson constructs `FA` with `fabs(delta a)`, losing the Wilson-coefficient
   sign and making the result nonanalytic at the SM point.
3. Jackson reconstructs `FV` from fully shifted factors without expanding to
   one dimension-six insertion.
4. Jackson's public cross section is not restricted to tree-SM/tree-D6
   interference.
5. Jackson's four-fermion `VV` and `AA` projection signs differ from the
   GRIFFIN amplitude convention, and redundant flavor aliases can be counted
   independently.
6. The two zero-Wilson SM line shapes already differ near the Z pole.

Changing to `{alpha,MW,MZ}` removes only the Wilson-dependent derived-W-mass
feedback. It does not remove Jackson's mass-slot, `FA`, `FV`, contact, or
cross-section construction issues.

## Native Standard Model reference values

The electron EW pseudo-observables are already close in the two branches:

| Scheme | Observable | New GRIFFIN | Jackson | Relative difference |
|---|---|---:|---:|---:|
| Gmu | `SW` | 0.231524072386 | 0.231516529073 | -0.00326% |
| Gmu | `FA` | 0.0344618847451 | 0.0344628085127 | +0.00268% |
| Gmu | `FV` | 0.000193358366928 | 0.000193513229979 | +0.0801% |
| Alpha | `SW` | 0.231172010026 | 0.231164288194 | -0.00334% |
| Alpha | `FA` | 0.0344989583309 | 0.0344999146215 | +0.00277% |
| Alpha | `FV` | 0.000200833992945 | 0.000200995709670 | +0.0805% |

The native SM differential cross sections show a larger, sign-changing
line-shape difference near the Z pole:

| Scheme | `sqrt(s)` [GeV] | New [nb] | Jackson [nb] | Relative difference |
|---|---:|---:|---:|---:|
| Gmu | 50 | 0.0141059603 | 0.0141095332 | +0.0253% |
| Gmu | 90 | 0.448269928 | 0.434383009 | -3.0979% |
| Gmu | 95 | 0.132595542 | 0.134327929 | +1.3065% |
| Gmu | 120 | 0.00977048132 | 0.00977983654 | +0.0957% |
| Alpha | 50 | 0.0141014944 | 0.0141050660 | +0.0253% |
| Alpha | 90 | 0.449797783 | 0.435867767 | -3.0970% |
| Alpha | 95 | 0.132943175 | 0.134681005 | +1.3072% |
| Alpha | 120 | 0.00978359736 | 0.00979298638 | +0.0960% |

The almost identical pattern in both schemes demonstrates that this SM
line-shape difference is not caused by the SMEFT Gmu input transformation.

## Native linear SMEFT response of electron EWPOs

Every entry is the response per unit dimensionless Wilson coefficient at
`Lambda=246 GeV`. `New` is the signed linear term. `Jackson` is the public odd
projection defined above.

### `{alpha,MW,MZ}`

| Coefficient | `SW` new | `SW` Jackson | `FA` new | `FA` Jackson | `FV` new | `FV` Jackson |
|---|---:|---:|---:|---:|---:|---:|
| `CphiD` | -0.403155130 | -0.403426418 | +0.0426398672 | 0 | +0.0120225230 | +0.00839945286 |
| `CphiWB` | -0.431951839 | -0.432242499 | +0.128148392 | 0 | +0.0138436178 | +0.00902170788 |
| `Cphil1_11` | -0.231402725 | -0.231597791 | +0.0686509107 | 0 | +0.00741622234 | +0.00482556881 |
| `Cphil3_11` | -0.231402725 | -0.231597791 | +0.0686509107 | 0 | +0.00741622234 | +0.00482556881 |
| `Cphie_11` | -0.287453768 | -0.287696083 | -0.0686509107 | 0 | +0.00741622234 | +0.00599443216 |

The `SW` difference is only about 0.07--0.08%, matching the isolated
mass-slot effect. Jackson's `FA` odd projection vanishes because his public
term is proportional to `fabs(delta a)` and is therefore even under
\(C\to-C\). The corresponding Jackson even residuals are nonzero—for example
`8.5311e-5` for `CphiD` and `2.56469e-4` for `CphiWB` at the finite benchmark.
The 19--35% `FV` differences remain because Jackson multiplies fully shifted
`FA` and `SW` rather than constructing \(2v_0\delta v\).

### `{Gmu,alpha,MZ}`

| Coefficient | `SW` new | `SW` Jackson | `FA` new | `FA` Jackson | `FV` new | `FV` Jackson |
|---|---:|---:|---:|---:|---:|---:|
| `CphiD` | +0.185762120 | +0.164778262 | -0.0195756509 | -0.0599699410 | -0.00545481253 | -0.00369092648 |
| `CphiWB` | +0.791077997 | +0.747781174 | -0.000949166380 | -0.124592832 | -0.0222936416 | -0.0159846060 |
| `Cphil1_11` | -0.232051624 | -0.232245576 | +0.0686509107 | 0 | +0.00731602195 | +0.00474347496 |
| `Cphil3_11` | +0.0856793373 | +0.0743225982 | +0.0351683560 | -0.0322294506 | -0.00201399014 | -0.00169423388 |
| `Cphie_11` | -0.287409832 | -0.287650046 | -0.0686509107 | 0 | +0.00731602195 | +0.00587507765 |

The direct-only `SW` directions remain close. Directions that also move the
derived W mass differ by roughly 5--13% because Jackson feeds a
Wilson-dependent W mass through cumulative SM objects, whereas new GRIFFIN
keeps the SM reference fixed and adds the analytically linear response.

## Native SMEFT differential cross section

The following table reports the combined vertex+contact benchmark. It is the
coefficient of the benchmark `scale`, in nb.

| Scheme | `sqrt(s)` [GeV] | New linear term [nb] | Jackson odd projection [nb] | Relative difference |
|---|---:|---:|---:|---:|
| Alpha | 50 | -0.000877465 | -0.001085393 | +23.70% |
| Alpha | 90 | +0.579089155 | +0.142783721 | -75.34% |
| Alpha | 95 | +0.137067771 | +0.026007802 | -81.03% |
| Alpha | 120 | +0.006772256 | -0.000759449 | -111.21% |
| Gmu | 50 | +0.000923988 | +0.000608701 | -34.12% |
| Gmu | 90 | -0.091797693 | -0.420740509 | +358.34% |
| Gmu | 95 | -0.012830177 | -0.107989404 | +741.68% |
| Gmu | 120 | +0.001267153 | -0.005749350 | -553.72% |

At 90 GeV, separating the two sources makes the disagreement easier to
locate:

| Scheme | Scenario | New [nb] | Jackson [nb] |
|---|---|---:|---:|
| Alpha | Vertex only | +0.584873172 | +0.138368805 |
| Alpha | Contact only | -0.005784017 | +0.004412627 |
| Gmu | Vertex only | -0.132288942 | -0.464073801 |
| Gmu | Contact only | +0.040491249 | +0.043297368 |

The Alpha scheme removes the derived-W feedback but leaves large vertex and
contact differences. The public cross sections therefore must not be treated
as two numerical implementations of the same truncated observable.

## What the common-reference diagnostic establishes

The common diagnostic converts the same running-width inputs once and uses:

| Parameter | Common complex-pole value [GeV] |
|---|---:|
| `MW` | 80.3499709226283 |
| `MZ` | 91.1534806191828 |
| `GammaW` | 2.08429885879891 |
| `GammaZ` | 2.49426637877282 |

Representative signed direct-current coefficients then agree exactly:

| Coefficient | Quantity | New GRIFFIN | Jackson formula with common masses |
|---|---|---:|---:|
| `CphiD` | `delta a` | -0.117223941478 | -0.117223941478 |
| `CphiD` | `delta v` | -0.305956479149 | -0.305956479149 |
| `CphiWB` | `delta a` | -0.352300806953 | -0.352300806953 |
| `Cphil1_11` | `delta a` | -0.188732537671 | -0.188732537671 |

Calling Jackson's literal compiled current instead gives shifts of about
0.0366--0.0843% because it reads different mass slots. In the Gmu scheme the
raw dimension-six `DeltaR`/W-mass input transformation also agrees
coefficient by coefficient. Only the SM fixed-point response differs slightly:

```text
new GRIFFIN response factor   1.09111501347680
Jackson response factor       1.09105595754966
relative difference          -0.00541244%
```

For example, `CphiD` gives `delta MW=-30.4737823135 GeV per C` in new GRIFFIN
and `-30.4721329388 GeV per C` in Jackson at `Lambda=246 GeV`. This small
response difference cannot explain the large native `FA`, `FV`, and
cross-section disagreements.

When both branches are forced to the common masses and `SW`, `FA`, `FV`, and
the cross section are reconstructed externally with exactly one D6 insertion,
the Alpha-scheme results agree to numerical precision. This is the strongest
evidence that the primary disagreement is in Jackson's public observable
construction and conventions, not in the underlying signed direct-current
matching.

## Pinpointed issues

- **Mass slotting:** the direct SMEFT Z-current and the SM tree current use
  different W/Z mass slots in the tested Jackson implementation. This causes
  the small residual difference that remains after derived-W feedback is
  removed.
- **Nonlinear SMEFT construction:** `FA` uses `fabs(delta a)`, while `FV` and
  the public cross section are formed from fully shifted quantities instead
  of being expanded to exactly one dimension-six insertion. This loses sign
  information in `FA` and introduces terms beyond `1/Lambda^2`.
- **Contact matching conventions:** the `VV` and `AA` projections have a
  different sign pattern, and symmetry-related `Cll` and `Cee` flavor entries
  can be counted as separate inputs.
- **Matrix-element integration:** the direct SMEFT Z-current shift depends on
  separately supplied form-factor objects rather than being completely owned
  by the SMEFT matrix-element path.

## Bottom line

The signed direct-current formulas agree after a common-mass translation, but
the two public implementations are not yet numerically interchangeable. The
current GRIFFIN implementation enforces an arbitrary-scale, strictly linear
dimension-six construction. The tested Jackson public `FA`, `FV`, contact,
and cross-section definitions follow different mass, sign, and truncation
conventions. Those differences should be resolved or documented jointly
before either public result is used as a reference for validating the other.
