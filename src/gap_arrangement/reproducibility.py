"""Numerical regression checks and portable run provenance."""
from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
import json
import platform
import sys

import numpy as np

from . import __version__


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def compare_numeric(actual, expected, label, rtol=2e-5, atol=2e-10):
    """Compare scientific values, permitting small platform roundoff differences."""
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise AssertionError(f"{label}: different keys")
        for key in expected:
            compare_numeric(actual[key], expected[key], f"{label}.{key}", rtol, atol)
    elif isinstance(expected, list):
        if len(actual) != len(expected):
            raise AssertionError(f"{label}: different lengths")
        for index, (a, e) in enumerate(zip(actual, expected)):
            compare_numeric(a, e, f"{label}[{index}]", rtol, atol)
    elif isinstance(expected, (int, float)):
        if not np.isfinite(actual) or not np.isclose(actual, expected, rtol=rtol, atol=atol):
            raise AssertionError(f"{label}: got {actual}, expected {expected}")
    elif actual != expected:
        raise AssertionError(f"{label}: values differ")


def scientific_values(simulation, audit, affine):
    # Error estimates close to machine precision are checked against acceptance
    # thresholds, not against their exact values on the original machine.
    particle_keys = ["p", "N", "scale", "max_mesh_E_plus_D", "scaled_error",
                     "W1_initial", "W1_final", "W1_final_error",
                     "zero_W1_initial", "zero_W1_final", "min_gap", "max_gap"]
    affine_keys = ["n", "N", "scale", "max_mesh_E_plus_D", "scaled_error"]
    return {
        "parameters": simulation["parameters"],
        "C2": simulation["C2"],
        "C2_over_2": simulation["C2_over_2"],
        "affine_C2_over_2": simulation["affine_C2_over_2"],
        "final_cell_W1": simulation["final_cell_W1"],
        "particles": [{k: row[k] for k in particle_keys} for row in simulation["particle_rows"]],
        "response": simulation["observable_rows"],
        "coefficients": {beta: {k: value[k] for k in
                          ["C2", "C2_over_2", "W1_leading_coefficient"]}
                         for beta, value in audit["coefficients"].items()},
        "early_time": audit["early_time"],
        "finite_difference": audit["finite_difference"],
        "relaxation": audit["relaxation"],
        "affine_parameters": affine["parameters"],
        "affine_C2": affine["C2"],
        "affine_particles": [{k: row[k] for k in affine_keys} for row in affine["particle_rows"]],
        "affine_response": [{k: row[k] for k in ["t", "ratio"]} for row in affine["observable_rows"]],
    }


def verify(output_dir, reference_dir):
    from .simulation import validate_results as validate_simulation
    from .validation import validate_results as validate_audit
    from .certificate import certificate

    output_dir, reference_dir = Path(output_dir), Path(reference_dir)
    output, reference = output_dir / "data", reference_dir / "data"
    names = ["simulation_results.json", "audit_results.json", "sanity_results.json"]
    actual = [read_json(output / n) for n in names]
    expected = [read_json(reference / n) for n in names]
    validate_simulation(actual[0])
    validate_audit(actual[1])
    affine = actual[2]
    assert affine["reference_spatial_refinement_max"] < 1e-9
    assert all(r["time_refinement_difference"] < 1e-7 for r in affine["particle_rows"])
    assert abs(affine["observable_rows"][-1]["ratio"] / affine["C2_over_2"] - 1) < .01
    compare_numeric(scientific_values(*actual), scientific_values(*expected), "results")
    max_curve_difference = {}
    for name in ["simulation_curves.npz", "audit_curves.npz"]:
        with np.load(output / name, allow_pickle=False) as a, np.load(reference / name, allow_pickle=False) as b:
            assert set(a.files) == set(b.files), f"{name}: different curve names"
            maximum = 0.0
            for key in a.files:
                assert a[key].shape == b[key].shape, f"{name}/{key}: different shape"
                assert np.isfinite(a[key]).all(), f"{name}/{key}: non-finite values"
                np.testing.assert_allclose(a[key], b[key], rtol=2e-5, atol=2e-11,
                                           err_msg=f"{name}/{key}")
                maximum = max(maximum, float(np.max(abs(a[key] - b[key]))))
            max_curve_difference[name] = maximum
    assert read_json(output / "certificate.json") == certificate()
    from .exports import verify_exports
    export_checks = verify_exports(output_dir, actual[0], output)
    # The displayed tables are the rounded scientific reference outputs.
    for name in ["table1.tex", "table2.tex"]:
        assert (output_dir / "tables" / name).read_text() == (reference_dir / "tables" / name).read_text(), name
    result = {"status": "passed", "json_rtol": 2e-5, "json_atol": 2e-10,
              "curve_rtol": 2e-5, "curve_atol": 2e-11,
              "max_curve_absolute_difference": max_curve_difference,
              "export_checks": export_checks,
              "note": "Finite numerical checks, not interval certification or a proof."}
    (output / "verification.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def write_manifest(output_dir, command, started, reference_dir, inputs_dir):
    root = Path(output_dir)
    sources = Path(__file__).resolve().parent
    environment = {"python": sys.version, "platform": platform.platform(),
                   "machine": platform.machine()}
    for package in ["numpy", "scipy", "matplotlib", "contourpy", "cycler",
                    "fonttools", "kiwisolver", "packaging", "pillow",
                    "pyparsing", "python-dateutil", "six"]:
        environment[package] = version(package)
    outputs = {p.relative_to(root).as_posix(): sha256(p.read_bytes()).hexdigest()
               for folder in ["data", "figures", "tables", "inputs"]
               for p in sorted((root / folder).glob("*")) if p.is_file()}
    reference = {p.relative_to(reference_dir).as_posix(): sha256(p.read_bytes()).hexdigest()
                 for p in sorted((Path(reference_dir) / "data").glob("*")) if p.is_file()}
    manifest = {"software": "gap-arrangement-dynamics", "version": __version__,
                "command": command, "status": "passed", "started_utc": started,
                "completed_utc": datetime.now(timezone.utc).isoformat(),
                "environment": environment,
                "source_sha256": {p.name: sha256(p.read_bytes()).hexdigest()
                                  for p in sorted(sources.glob("*.py"))},
                "output_sha256": outputs, "reference_sha256": reference,
                "input_sha256": {p.name: sha256(p.read_bytes()).hexdigest()
                                 for p in sorted(Path(inputs_dir).glob("*"))
                                 if p.is_file() and p.suffix in [".json", ".npz"]},
                "scope": "Status applies to the named command; full pipeline is 'reproduce'."}
    (root / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
