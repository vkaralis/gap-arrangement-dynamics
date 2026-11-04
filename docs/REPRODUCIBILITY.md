# Reproducibility and provenance

## Scope and source

The v1.0.0 pipeline recomputes manuscript Figures 1–2 and Tables 1–2 from
analytic initial conditions and compares the results with the supplied
reference data. Figures 3–5, Tables 3–5 and the Table 5 refinement record
(`simulations/wasserstein_results.json`) belong to the separate manuscript
supplement. See the [manuscript figure and table guide](../README.md#manuscript-figures-and-tables).

The complete pipeline runs the audit before generating Figure 1, because its
relaxation panels use audit results. Missing audit files cause an explicit error.
Reference plots can be regenerated with the `figures` command without running
the ODE solvers; this is documented as rendering existing data, not recomputation.

## Environment

The supported interpreter series is **Python 3.12 only** (`>=3.12,<3.13`).
The tested interpreter is Python 3.12.14. Packaging, the source entry point,
`.python-version`, documentation and GitHub Actions agree on this support range.
The direct numerical dependencies are **NumPy 2.3.5**, **SciPy 1.17.1** and
**Matplotlib 3.10.9**, with **Pillow 12.3.0** for figure export and image checks.
`requirements-lock.txt` pins all runtime dependencies;
`requirements-build.txt` pins the build tools. `requirements.txt` lists only the
four direct dependencies, consistently with `pyproject.toml`.

`requirements-lock.txt` is the explicit runtime install specification. It is
used with `requirements-build.txt` for a fully pinned installation. Platform-specific
wheels can differ, so cross-platform bitwise equality is not assumed. Actual versions
are recorded in locally generated run manifests, not embedded in public figures.

NumPy FFTs and the solvers use CPU arrays. The code defaults OpenBLAS and OpenMP
to one thread when their environment variables have not already been set.
Figures use Matplotlib's noninteractive `Agg` backend and DejaVu fonts.
The only random values are transport test cases, with fixed seeds `20261102`
and `7713`; the scientific trajectories are deterministic.

## Validation layers

1. Fast tests compare independent fields, circular transport against linear
   programming, conservative gap and position evolution, affine exact solutions,
   zero-diffusion closure, quadrature and the rational certificate. CLI tests
   cover execution outside the checkout, missing inputs, protected references
   and disabled-assertion mode.
2. The complete audit checks all main particle sizes and frequencies, cell
   discretizations, time/spatial refinements, conservation, positivity and
   theoretical identities. Acceptance thresholds are enforced in
   `simulation.validate_results` and `validation.validate_results`.
3. `verify` compares the scientific outputs and every saved curve with the
   frozen synthetic reference results. JSON comparisons use relative tolerance `2e-5`
   and absolute tolerance `2e-10`; curves use `2e-5` and `2e-11`. Near-zero
   solver discrepancies are checked against acceptance thresholds rather than
   required to match their original roundoff values. The displayed LaTeX tables
   must match the reference text. Every CSV header, row and numeric cell is
   checked against its validated JSON source, using relative tolerance `1e-12`
   and zero absolute tolerance. Missing rows, non-finite values and extra columns
   fail verification. Both figures are regenerated in a temporary directory from
   the validated saved data: PNG pixels and metadata-free PDF bytes must match
   these freshly rendered outputs. A swapped figure or changed pixel fails.
   This compares within the current rendering environment, rather than requiring
   identical font rendering across operating systems. To verify results made
   with a different rendering environment, regenerate their figures first.

These regression tolerances are acceptance criteria, not bounds on the error
relative to an exact continuum solution. If a check fails, inspect the equations,
settings, environment, and refinement results before changing any tolerance.

The numerical outputs record the actual values of all reported audit checks.
The certificate uses exact rational arithmetic for the specific analytical bound
encoded there; it does not turn the floating-point simulations into interval
certificates.

## Clean-environment procedure

Run `python scripts/clean_reproduce.py --work-dir clean-run` with Python 3.12
and a destination that does not already exist. The script:

1. Copies the source and supplied data to a fresh work directory, excluding
   existing environments, generated outputs, caches and Git metadata.
2. Creates a venv with `include-system-site-packages = false`, removes inherited
   Python path overrides and disables the user package directory.
3. Downloads all pinned runtime/build dependencies from PyPI using
   `--isolated --no-cache-dir`, then installs a newly built package wheel.
4. Runs `pip check` and asserts that every distribution and numerical module
   comes from the new venv, with no external site-packages directories.
5. Runs all tests and the complete pipeline from outside the copied source using
   `python -I`. The supplied initial arrays are checked before integration.
6. Records tests, input agreement, curve comparison, source hashes, timestamps
   and the actual environment in `clean_reproduction.json`.

The logs remain in the selected work directory. A failed step stops the run and
records failure. An existing venv is never accepted as a new clean test. The
base Python interpreter and operating system are shared with the host; this is
Python-package isolation, not a container or a claim of testing another OS.

## Portable execution and run records

`python -m gap_arrangement reproduce` is a full recomputation. It does not copy
reference results into generated outputs. It writes a manifest with actual UTC
timestamps, package versions, source hashes and output hashes. A failed or
interrupted command never leaves a previous successful manifest as the current
run record: the status starts as `running`, then becomes `passed` or `failed`.
A hard process kill can leave `running`. The manifest's status applies to its
named command; only `reproduce` runs the entire pipeline. Rerunning `verify`
checks existing files without reintegrating the equations.

The source-checkout runner and an installed package use the same code. The
installed wheel contains the Python implementation; the repository also supplies
the synthetic inputs, frozen reference data, tests, documentation and workflows.
Retain `inputs/` and `reference/` when using an installed wheel for reproduction;
pass `--inputs-dir` and `--reference-dir` from another working directory.

## Continuous integration

The `Numerical checks` GitHub Actions workflow runs the fast checks on Linux,
macOS and Windows with Python 3.12 after pushes or manual dispatch. The
`Full reproduction` workflow runs the clean-environment script for version tags
or when manually triggered in the Actions tab, and uploads the outputs and
verification record. See the
[GitHub Actions runs](https://github.com/vkaralis/gap-arrangement-dynamics/actions)
for the status of a specific commit.

Workflow syntax follows the official
[checkout](https://github.com/actions/checkout),
[setup-python](https://github.com/actions/setup-python), and
[upload-artifact](https://github.com/actions/upload-artifact/tree/v4)
documentation.

## Public files and local run records

Public PDF figures omit Creator, Producer, creation/modification dates and other
optional document information. PNG figures contain pixel data without text,
EXIF, software, color-profile or resolution metadata. Reference JSON files retain
scientific parameters/results and omit runtime timings and machine inventories.
The redundant environment inventory and historical machine-specific validation
record have been removed. The published source contains only the repository files;
Git history, virtual environments, caches and duplicated release bundles are
excluded. Numerical data, tests, licensing and installation files are retained.

A new run still generates its own environment, timestamps and hashes locally for
reproducibility, under ignored `outputs/` or `clean-run/` directories. Those files
are not part of the upload. The package version, dependency requirements, MIT
notice and intentional software citation are retained because they identify and
make the software usable.

## Software credit and release

The code creator and software citation author is **Vangelis D. Karalis**.
**Georgia Karali**, Department of Mathematics, National and Kapodistrian
University of Athens, is acknowledged as a contributor to the study. This
acknowledgment does not assign her code authorship or a more specific role.

`CITATION.cff`, `pyproject.toml`, the package and `CHANGELOG.md` describe version
**1.0.0**, with tag **v1.0.0** and license MIT. The repository is
[vkaralis/gap-arrangement-dynamics](https://github.com/vkaralis/gap-arrangement-dynamics).
No software DOI has been assigned, and the associated manuscript has not yet
been submitted. The citation identifies the software only.
GitHub describes software citation support in
[About CITATION files](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files).
