"""Fast checks against independent formulations; full audits use `reproduce`."""
from fractions import Fraction
import unittest

import numpy as np

from gap_arrangement import simulation as sim
from gap_arrangement import validation as audit
from gap_arrangement.certificate import certificate
from gap_arrangement.reproducibility import compare_numeric


class NumericalTests(unittest.TestCase):
    def test_moment_fields_against_direct_pairwise_sum(self):
        rng = np.random.default_rng(19)
        x, h = rng.random(23), rng.random(23)
        np.testing.assert_allclose(sim.fields(x, h), audit.dense_fields(x, h), rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(audit.complex_fields(x, h), audit.dense_fields(x, h), rtol=1e-13, atol=1e-13)

    def test_circle_transport_against_linear_programming(self):
        result = audit.transport_checks()
        self.assertEqual(result["cases"], 60)
        for key in ["max_LP_error", "max_rotation_error", "max_symmetry_error"]:
            self.assertLess(result[key], 1e-11)
        self.assertAlmostEqual(sim.circle_w1(np.array([.99]), np.array([.01])), .02, places=14)

    def test_position_ode_against_conservative_gap_ode(self):
        x, g = audit.independent_particles(4, 1, 1.)
        reference = sim.particles(4, 1)
        reference_gaps = 64 * np.diff(np.c_[reference, reference[:, :1] + 1], axis=1)
        np.testing.assert_allclose(x, reference, rtol=0, atol=1e-9)
        np.testing.assert_allclose(g, reference_gaps, rtol=0, atol=1e-7)
        np.testing.assert_allclose(g.mean(1), 1, rtol=0, atol=1e-12)
        self.assertGreater(g.min(), 0)

    def test_cell_solver_against_affine_heat_formula(self):
        y = audit.independent_cell(32, 32, 2, beta=0)
        exact = sim.affine(32, 2)
        expected = exact[:, :, 1, None] * (
            1 + .25 * np.exp(-16 * np.pi**2 * sim.TIMES[:, None, None]) *
            np.cos(4 * np.pi * np.arange(32)[None, None, :] / 32))
        np.testing.assert_allclose(y[:, :, 0], exact[:, :, 0], rtol=0, atol=1e-9)
        np.testing.assert_allclose(y[:, :, 1:], expected, rtol=0, atol=1e-8)

    def test_zero_diffusion_equal_laws(self):
        times = np.array([0., .001, .002])
        y1 = audit.independent_cell(32, 33, 1, lam=0, times=times)
        y2 = audit.independent_cell(32, 33, 2, lam=0, times=times)
        np.testing.assert_allclose(y1[:, :, 0], y2[:, :, 0], rtol=0, atol=1e-11)
        np.testing.assert_allclose(y1[:, :, 1:].mean(2), y2[:, :, 1:].mean(2), rtol=0, atol=1e-11)

    def test_second_derivative_coefficient_independent_quadrature(self):
        for beta in [0., .5]:
            coefficient, _ = audit.coefficient_quad(beta)
            self.assertAlmostEqual(coefficient, sim.coefficient(beta), places=10)

    def test_exact_rational_time_bound(self):
        result = certificate()
        self.assertLess(Fraction(result["L3_upper_rational"]), 10**9)
        self.assertEqual(Fraction(result["certified_time_rational"]), Fraction(729, 8192 * 10**9))

    def test_regression_check_rejects_changed_or_nonfinite_results(self):
        for value in [.8, float("nan"), float("inf")]:
            with self.assertRaises(AssertionError):
                compare_numeric({"coefficient": value}, {"coefficient": .741984411175}, "test")


if __name__ == "__main__":
    unittest.main()
