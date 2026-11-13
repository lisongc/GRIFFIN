# Documentation

These guides describe the local merged GRIFFIN implementation. Start with
`<griffin/SMEFT.h>` and the [working example](../examples/smeft_example.cc).

**Review-branch warning:** the unchanged `mat_SMaas` scattering background
has unresolved near-real-Z-pole behaviour. Do not treat this branch as a
validated precision Z-lineshape calculation. The [full-background warning](KNOWN_ISSUES.md#1-unresolved-near-pole-scattering-background)
includes a small reproducible diagnostic; passing ordinary regression tests
does not certify that region.

## What calculation is implemented?

**GRIFFIN's NNLO+ SM form factors and mixed electroweak–QCD-extended
scattering amplitude, with dimension-six corrections inserted using a
finite, unlinearized prescription.**

- **NNLO+** means NNLO form factors plus the selected higher-order SM terms
  already supplied by GRIFFIN.
- **Extended scattering amplitude** means `mat_SMaas`: its selected
  fermionic mixed electroweak–QCD terms affect both the Z-pole residue and
  the off-resonance background. It is not complete off-resonance NNLO.
- **Finite, unlinearized** means evaluating at the supplied coefficient
  values, solving the Gmu mass equation without a mass-derivative
  approximation, and squaring the resulting amplitude without discarding
  the nonlinear products. It does not mean that every possible higher-power
  SMEFT term is included. The explicit dimension-six pole-input shifts
  remain linear at fixed mass; see the [precise formulas](USER_GUIDE.md).

`reference()` is the result with all Wilson coefficients zero, `total()`
is the full result at the requested coefficients, and `modelCorrection()` is
`total() - reference()`. The shift is a finite difference, not a derivative
or a projection onto the term linear in the coefficients.

## Where to read next

| Your question | Document |
|---|---|
| How do I supply inputs, run a scan and read results? | [User guide](USER_GUIDE.md) |
| What are the coefficient conventions and exact formulas? | [SMEFT implementation](SMEFT_IMPLEMENTATION.md) |
| How do the classes and directories fit together? | [Architecture](ARCHITECTURE.md) |
| Which processes and kinematics are unsupported? | [Known limitations](KNOWN_ISSUES.md) |
| How do I update an older client program? | [Migration guide](MIGRATION.md) |
| What was tested, and how can I rerun it? | [Current regression report](REGRESSION_TEST_REPORT.md) |

The files in [history/](history/README.md) preserve older comparisons and
dated development checks. Their old class names, linear prescriptions and
benchmark numbers are historical evidence, not instructions for the current
calculator. In particular, the old Jackson cross-section benchmarks used
`mat_SMNNLO`; they are tested separately from today's extended `mat_SMaas`
calculation.
