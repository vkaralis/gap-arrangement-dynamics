import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from gap_arrangement import simulation as sim
from gap_arrangement.synthetic import initial_arrays, write_inputs, validate_inputs


class SyntheticInputTests(unittest.TestCase):
    def test_saved_particle_preparation_matches_solver(self):
        arrays = initial_arrays()
        positions = sim.particles(4, 1)
        np.testing.assert_array_equal(positions[0], arrays["particle_n4_k1_x0"])
        gaps = 64 * np.diff(np.r_[positions[0], positions[0, 0] + 1])
        np.testing.assert_allclose(gaps, arrays["particle_n4_k1_g0"], atol=2e-14, rtol=2e-14)

    def test_inputs_regenerate_and_detect_modified_parameters(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            write_inputs(path)
            self.assertEqual(validate_inputs(path)["max_absolute_difference"], 0.)
            config = json.loads((path / "experiment.json").read_text())
            config["model"]["beta"] = 1.
            (path / "experiment.json").write_text(json.dumps(config))
            with self.assertRaises(AssertionError):
                validate_inputs(path)

    def test_inputs_detect_modified_array(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            write_inputs(path)
            with np.load(path / "initial_conditions.npz", allow_pickle=False) as archive:
                arrays = {key: archive[key] for key in archive.files}
            arrays["particle_n8_k1_g0"][0] += .01
            np.savez_compressed(path / "initial_conditions.npz", **arrays)
            with self.assertRaises(AssertionError):
                validate_inputs(path)
