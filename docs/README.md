<!--
File: README.md
Purpose: Route readers to the current guides and historical evidence.
Authors (merged-version documentation): Jackson Wallace and Lisong Chen
Last updated: 2026-11-16
-->

# Documentation

Start with the [build instructions](../README.md), `<griffin/SMEFT.h>` and
the [SMEFT example](../examples/smeft_example.cc).

`SMEFTCalculator` uses the finite SMEFT prescription with `mat_SMaas`.
The imported `mat_SMNNLO2` class is available separately for direct SM use.
Read the [current pole restrictions](KNOWN_ISSUES.md#1-corrected-mixed-background-and-remaining-pole-restriction)
before using either calculation for a lineshape study.

| Your question | Document |
|---|---|
| How do I supply inputs, run a scan and read results? | [User guide](USER_GUIDE.md) |
| What are the coefficient conventions and exact formulas? | [SMEFT implementation](SMEFT_IMPLEMENTATION.md) |
| How do the classes, lifetimes and directories fit together? | [Architecture](ARCHITECTURE.md) |
| Which processes and kinematics are unsupported? | [Known limitations](KNOWN_ISSUES.md) |
| How do external users update an existing program? | [Migration guide](MIGRATION.md) |
| What is tested, and how can I run it? | [Testing](../README.md#testing) |

[Historical records](history/README.md) preserve earlier comparisons and
development checks. Their old APIs, linear prescriptions and benchmark
numbers describe earlier implementations; use the guides above for the
current calculator. The historical Jackson `mat_SMNNLO` baseline and the
current `mat_SMaas`/upstream-reference checks remain distinct.
