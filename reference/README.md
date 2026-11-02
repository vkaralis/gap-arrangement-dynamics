# Frozen numerical reference

These synthetic numerical results are frozen regression inputs. A complete run computes new results in `outputs/`;
it does not substitute these files for integration.

| File | Interpretation |
| --- | --- |
| `data/simulation_results.json` | Parameters, particle-error rows, response rows and refinement diagnostics |
| `data/simulation_curves.npz` | Main time-dependent observables and circular distances |
| `data/audit_results.json` | Independent solvers, conservation/identity checks, quadrature, refinements and relaxation summaries |
| `data/audit_curves.npz` | Relaxation errors and normalized-cell variance |
| `data/sanity_results.json` | Supplementary affine experiment (`beta=0`, `T=0.01`) |
| `figures/` | Two reference figures, PDF and PNG |
| `tables/` | Two rounded reference LaTeX tables and full-precision CSV exports |

NPZ files are NumPy compressed archives. Load them with
`numpy.load(path, allow_pickle=False)`. Every curve is a one-dimensional array
aligned with the archive's `times` array. These are dimensionless model variables
on a circle of length one; no physical-unit calibration is implied.

`simulation_curves.npz` keys:

- `times`: the union of the 101 uniform output times and additional small times.
- `cell_w1`: circular distance between nonlinear `k=1` and `k=2` cell maps.
- `cell_observable`, `affine_observable`: `J2-J1` for the two diffusion laws.
- `w1_n8`, `w1_n16`, `w1_n32`: nonlinear particle circular distances, `N=n^3`.
- `zero_w1_n8`, `zero_w1_n16`, `zero_w1_n32`: corresponding zero-diffusion controls.

`audit_curves.npz` keys:

- `times`: the shared output times.
- `relaxation_error_L`: `max_a|X-Xbar| + max_a|q-qbar|`, for
  `L = 1,2,4,8,16,32` and cell diffusivity `lambda=L`.
- `variance_L`: `max_a integral (G/q-1)^2 dtheta` for the same diffusivities.

The JSON keys `max_*`, `*_error`, `*_difference`, and `*_residual` refer to the
specific sampled comparisons described in [the methods](../docs/NUMERICAL_METHODS.md).
Runtime timings and machine inventories are omitted from the public reference
JSON. They may be recorded locally by a new run and are not scientific regression
targets. CSV files preserve computed precision; the LaTeX
tables round values for presentation. The coefficient CSV contains the theoretical
zero-time values separately from the positive-time response data.
