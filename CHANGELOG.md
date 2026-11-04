# Changelog

## [1.0.0] - 2026-11-02

- First public release as `gap-arrangement-dynamics`.
- Include Python simulations, synthetic initial data, frozen reference results,
  independent numerical checks and complete reproduction instructions.
- Provide software citation information without an article DOI or publication claim.
- Align package, citation and release versions at 1.0.0.
- Configure cross-platform tests and full clean-environment reproduction in GitHub Actions.
- Remove upload instructions and retain a distribution-only repository policy.
- Read release documentation explicitly as UTF-8 on all supported platforms.
- Archive the repository after successful release verification to prevent comments.

## [0.5.1]

### Fixed

- Validate every CSV field and value against the numerical JSON results.
- Regenerate figures and check PNG pixels and PDF content, detecting swapped,
  edited or incomplete outputs.
- Add regression tests for corrupted CSVs, changed images and embedded metadata.
- Document the response coefficient and the rational observation-time bound.

### Updated

- Identify Vangelis D. Karalis as the code creator and software citation author.
- Acknowledge Georgia Karali (Department of Mathematics, National and Kapodistrian
  University of Athens) as a contributor to the study.
- Remove embedded authoring metadata from PDF/PNG figures, runtime inventories
  from public JSON, and redundant environment/validation snapshots.
- Deliver a source-only ZIP without duplicated Git bundles or release manifests.
- Declare Pillow explicitly for export validation; use modern MIT license metadata.

## [0.5.0]

- Initial standalone particle and cell simulations.
- Synthetic input generation, independent numerical audits and frozen reference data.
- Python 3.12 environment specification, clean reproduction script and CI workflows.
