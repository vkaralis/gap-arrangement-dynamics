"""Deterministic synthetic initial conditions and experiment specifications."""
from hashlib import sha256
from pathlib import Path
import json

import numpy as np

from . import simulation as sim


def specification():
    """Machine-readable description of this fixed experiment, not a sweep API."""
    return {
        "schema_version": 1,
        "data_kind": "synthetic; analytic initial conditions; no observational data",
        "units": "dimensionless; circle length one",
        "model": {"alpha": sim.ALPHA, "delta": sim.DELTA, "eta": sim.ETA,
                  "beta": sim.BETA, "m0": 1., "gamma": 1., "kappa": 1.,
                  "K": "1 + eta*cos(2*pi*x)", "W": "-kappa*cos(2*pi*x)",
                  "Phi": "exp(-g)", "A": "g + beta*g**2"},
        "initial": {"rho0": "1 + alpha*cos(2*pi*x)",
                    "quantile": "X0 + alpha*sin(2*pi*X0)/(2*pi) = a",
                    "q0": "1/rho0(X0)",
                    "G0": "q0(a)*(1 + delta*cos(2*pi*k*theta))",
                    "particle_gaps": "G0(i/N, i/p) divided by its discrete mean",
                    "particle_lifts": "x0=0; successive cumulative gaps divided by N"},
        "main": {"n": [8, 16, 32], "N": [512, 4096, 32768], "p": [8, 16, 32],
                 "k": [1, 2], "lambda": 1., "control_lambda": 0., "T": sim.T,
                 "nu_N": "lambda*(p/N)**2", "uniform_times": sim.GRID.tolist(),
                 "response_times": sim.SMALL.tolist(), "saved_times": sim.TIMES.tolist(),
                 "cell_grid": [512, 64], "coarse_cell_grid": [256, 32],
                 "cell_dt": 5e-6, "coarse_cell_dt": 2e-5,
                 "particle_dt": 1e-5, "particle_refined_dt": 5e-6,
                 "affine_labels": 1024, "w1_labels": 8192, "w1_refined_labels": 16384},
        "audit": {"integrator": "DOP853", "rtol": 2e-12, "atol": 2e-14,
                  "particle_step_cap": "1/(4*lambda*p**2) if lambda>0",
                  "cell_step_cap": "1/(lambda*(pi*m)**2) if lambda>0",
                  "initial_step_cap": 1e-6, "cell_grid": [64, 64],
                  "label_refinement": [128, 64], "theta_refinement": [64, 128],
                  "finite_difference_m": [16, 32, 64, 128], "finite_difference_a": 32,
                  "early_grid": [128, 32],
                  "early_times": [0., 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3],
                  "relaxation_lambda": [1., 2., 4., 8., 16., 32.],
                  "relaxation_grid": [64, 32], "relaxation_refined_grid": [64, 64],
                  "transport_seed": 7713, "transport_cases": 60},
        "supplementary_affine": {"beta": 0., "T": .01, "n": [4, 8, 16],
                                  "reference_labels": 4096, "refined_labels": 8192,
                                  "dt": 2e-5, "refined_dt": 1e-5},
        "production_transport_seed": 20261102,
    }


def initial_arrays():
    """Generate the numerical preparations used by the production and audit runs."""
    from . import validation as audit
    arrays = {"times": sim.TIMES, "uniform_times": sim.GRID,
              "response_times": sim.SMALL}
    for n in [4, 8, 16, 32]:
        count = n**3
        _, q = sim.initial(count)
        for k in [1, 2]:
            g = q * (1 + sim.DELTA * np.cos(2*np.pi*k*np.arange(count)/n))
            g /= g.mean()
            arrays[f"particle_n{n}_k{k}_g0"] = g
            arrays[f"particle_n{n}_k{k}_x0"] = np.r_[0., np.cumsum(g[:-1])/count]
    for na, m in [(512, 64), (256, 32)]:
        x, q = sim.initial(na)
        for k in [1, 2]:
            g = q[:, None] * (1 + sim.DELTA*np.cos(2*np.pi*k*np.arange(m)/m))
            arrays[f"cell_a{na}_m{m}_k{k}"] = np.c_[x, g]
    cases = [(64, 64, 1), (64, 64, 2), (128, 64, 2), (64, 128, 2),
             (64, 65, 1), (64, 65, 2), (128, 32, 1), (128, 32, 2), (64, 32, 1)]
    cases += [(32, m, 2) for m in [16, 32, 64, 128]]
    for na, m, k in cases:
        arrays[f"audit_cell_a{na}_m{m}_k{k}"] = audit.init(na, k, m)
    for na in [64, 1024, 4096, 8192]:
        x, q = sim.initial(na)
        arrays[f"mean_a{na}_x0_q0"] = np.c_[x, q]
    return arrays


def write_inputs(destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "experiment.json").write_text(
        json.dumps(specification(), indent=2) + "\n", encoding="utf-8")
    arrays = initial_arrays()
    np.savez_compressed(destination / "initial_conditions.npz", **arrays)
    return {"arrays": len(arrays), "kind": "deterministic synthetic initial conditions"}


def validate_inputs(directory):
    directory = Path(directory)
    actual_specification = json.loads((directory / "experiment.json").read_text(encoding="utf-8"))
    if actual_specification != specification():
        raise AssertionError("Synthetic experiment settings differ from the implemented experiment.")
    expected = initial_arrays()
    with np.load(directory / "initial_conditions.npz", allow_pickle=False) as actual:
        assert set(actual.files) == set(expected), "Different synthetic initial-condition keys"
        maximum = 0.
        for key, values in expected.items():
            assert actual[key].shape == values.shape, f"Synthetic array shape: {key}"
            np.testing.assert_allclose(actual[key], values, rtol=2e-13, atol=2e-14, err_msg=key)
            maximum = max(maximum, float(np.max(abs(actual[key] - values))))
    return {"status": "passed", "arrays": len(expected),
            "max_absolute_difference": maximum,
            "sha256": {name: sha256((directory / name).read_bytes()).hexdigest()
                       for name in ["experiment.json", "initial_conditions.npz"]}}
