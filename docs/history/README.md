# Historical records

These documents describe earlier versions. They are retained for provenance;
do not use their class names, build commands or accuracy statements as
instructions for the current library.

- [July implementation comparison](SMEFT_JACKSON_COMPARISON.md): older
  branches, including the removed linear SMEFT implementation and
  derivative/odd-projection diagnostics. Its reported disagreements do
  not describe the merged calculator.
- [Integration and regression log through 12 November 2026](REGRESSION_LOG_2026-11-12.md):
  dated merge, API, contact-alias, layout and amplitude-extension checks.
  Intermediate statements about a retained prediction engine apply only
  to those earlier stages; the engine was subsequently removed.
- [Pre-import regression chronology, 13–16 November 2026](REGRESSION_LOG_2026-11-16_PRE_IMPORT.md):
  source/interface cleanup checks, before/after hashes, fresh-checkout audit,
  historical amplitude-impact study and the last pre-import sanitizer run.
  Old v1.1 near-pole observations and header counts belong to those snapshots;
  the upstream v1.1.1/v1.2 import supersedes their current-behavior claims.
- [Regression audit, 16 November 2026](REGRESSION_REPORT_2026-11-16.md):
  upstream import, independent-reference comparisons, numerical tolerances,
  diagnostics and the subsequent example-layout verification.

Dated JSON reports and diagnostic outputs are stored in `validation/`.
Their contents retain the commands and paths used at the time of each run.

For current behavior use the [documentation index](../README.md) and
[user guide](../USER_GUIDE.md). Current test commands are in the
[repository README](../../README.md#testing).
