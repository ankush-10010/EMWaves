"""Unit and validation tests for Rectangular Waveguide TM mode physics."""

import unittest
import math
import numpy as np

from core.waveguide import RectangularWaveguide, eps_0, mu_0
C_0 = 299792458.0
from core.fields import TMFieldCalculator
from core.cutoff import compute_cutoff_table, cutoff_frequency, get_dominant_tm_mode


class TestWaveguidePhysics(unittest.TestCase):
    """Test waveguide mathematical and physical calculations."""

    def setUp(self):
        # WR-90 waveguide dimensions: a = 22.86 mm, b = 10.16 mm
        self.a = 0.02286
        self.b = 0.01016
        self.wg = RectangularWaveguide(self.a, self.b, eps_r=1.0, mu_r=1.0)

    def test_dominant_tm_mode(self):
        """TM mode requires m>=1 and n>=1. The lowest TM mode is TM11."""
        dom_mode, fc = get_dominant_tm_mode(self.a, self.b)
        self.assertEqual(dom_mode, (1, 1))
        # fc for TM11 in WR-90 is ~16.145 GHz
        self.assertAlmostEqual(fc / 1e9, 16.145, places=2)

    def test_invalid_modes(self):
        """TM00, TM01, TM10 must raise ValueError."""
        with self.assertRaises(ValueError):
            self.wg.cutoff_frequency(0, 0)
        with self.assertRaises(ValueError):
            self.wg.cutoff_frequency(1, 0)
        with self.assertRaises(ValueError):
            self.wg.cutoff_frequency(0, 1)
        with self.assertRaises(ValueError):
            self.wg.cutoff_frequency(-1, 2)

    def test_invalid_dimensions(self):
        """Negative or zero dimensions/permittivities must raise ValueError."""
        with self.assertRaises(ValueError):
            RectangularWaveguide(-0.01, 0.02)
        with self.assertRaises(ValueError):
            RectangularWaveguide(0.02, 0.01, eps_r=0)
        with self.assertRaises(ValueError):
            RectangularWaveguide(0.02, 0.01, mu_r=-1)

    def test_cutoff_frequency_formula(self):
        """Verify fc = (c / (2*pi*sqrt(mu*eps))) * sqrt((m*pi/a)^2 + (n*pi/b)^2)."""
        m, n = 2, 1
        expected_kc = math.sqrt((2 * math.pi / self.a)**2 + (1 * math.pi / self.b)**2)
        expected_fc = expected_kc / (2 * math.pi * math.sqrt(mu_0 * eps_0))
        computed_fc = self.wg.cutoff_frequency(m, n)
        self.assertAlmostEqual(computed_fc, expected_fc, places=4)

    def test_dielectric_loading(self):
        """Dielectric loading with eps_r = 4 should halve the cutoff frequency."""
        wg_dielectric = RectangularWaveguide(self.a, self.b, eps_r=4.0, mu_r=1.0)
        fc_air = self.wg.cutoff_frequency(1, 1)
        fc_dielectric = wg_dielectric.cutoff_frequency(1, 1)
        self.assertAlmostEqual(fc_dielectric, fc_air / 2.0, places=4)

    def test_propagation_and_attenuation(self):
        """Operating above cutoff -> propagating (real beta); below -> evanescent (imaginary beta)."""
        fc = self.wg.cutoff_frequency(1, 1)
        # Below cutoff
        f_below = fc * 0.8
        self.assertFalse(self.wg.is_propagating(1, 1, f_below))
        beta_below = self.wg.propagation_constant(1, 1, f_below)
        self.assertIsInstance(beta_below, complex)
        self.assertEqual(beta_below.real, 0.0)
        self.assertGreater(beta_below.imag, 0.0)

        # Above cutoff
        f_above = fc * 1.2
        self.assertTrue(self.wg.is_propagating(1, 1, f_above))
        beta_above = self.wg.propagation_constant(1, 1, f_above)
        self.assertIsInstance(beta_above, float)
        self.assertGreater(beta_above, 0.0)

    def test_phase_and_group_velocity(self):
        """vp * vg must equal c_medium^2 for propagating modes."""
        fc = self.wg.cutoff_frequency(1, 1)
        f = fc * 1.5
        vp = self.wg.phase_velocity(1, 1, f)
        vg = self.wg.group_velocity(1, 1, f)
        c_med = self.wg.c_medium
        self.assertGreater(vp, c_med)  # vp > c
        self.assertLess(vg, c_med)     # vg < c
        self.assertAlmostEqual(vp * vg, c_med**2, places=3)

    def test_boundary_conditions(self):
        """Ez must vanish identically along PEC walls (x=0, x=a, y=0, y=b)."""
        calc = TMFieldCalculator(self.a, self.b, 1, 1, 1.0, 1.0, 20e9)
        xy = calc.compute_xy_fields(nx=100, ny=100, z=0, t=0)
        Ez = xy['Ez']

        # Check boundaries: bottom (y=0), top (y=b), left (x=0), right (x=a)
        np.testing.assert_allclose(Ez[0, :], 0.0, atol=1e-12)   # y = 0
        np.testing.assert_allclose(Ez[-1, :], 0.0, atol=1e-12)  # y = b
        np.testing.assert_allclose(Ez[:, 0], 0.0, atol=1e-12)   # x = 0
        np.testing.assert_allclose(Ez[:, -1], 0.0, atol=1e-12)  # x = a

    def test_mode_table_sorting(self):
        """Mode table must be sorted monotonically by cutoff frequency."""
        table = compute_cutoff_table(self.a, self.b, 1.0, 1.0, m_max=4, n_max=4)
        frequencies = [item['fc'] for item in table]
        self.assertEqual(frequencies, sorted(frequencies))
        self.assertEqual((table[0]['m'], table[0]['n']), (1, 1))


if __name__ == '__main__':
    unittest.main()
