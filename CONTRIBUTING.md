# Repository policy

This repository distributes the simulation code and data for reproducibility.
External contributions, support requests and discussions are not accepted.
Issues, pull requests and discussions are disabled. No external collaborators
are granted write access.

The repository is archived after release verification to make code and comments
read-only. Only the owner can unarchive it for a future update.

The MIT license permits reuse of the software. For local modifications, run
`python -m unittest discover -s tests -v` and the full
`python -m gap_arrangement reproduce` pipeline. Keep the frozen reference data
unchanged when checking results.
