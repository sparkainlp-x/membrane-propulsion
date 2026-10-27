"""Tests for the receiver sensor. Run: python test_sensor.py"""

import math
import unittest

from membrane_propulsion import MembranePropulsion
from sensor import ReceiverSensor


def locked_model(seed=42, drive_k=1, steps=600):
    m = MembranePropulsion(seed=seed, drive_k=drive_k)
    m.run(steps)
    return m


class TestSensor(unittest.TestCase):
    def test_sensor_is_read_only(self):
        m = locked_model()
        before_theta = list(m.theta)
        before_probe = m.probe
        before_t = m.t
        s = ReceiverSensor(m, position=10.0)
        s.sample()
        self.assertEqual(m.theta, before_theta)
        self.assertEqual(m.probe, before_probe)
        self.assertEqual(m.t, before_t)

    def test_phase_at_integer_position_equals_oscillator(self):
        m = locked_model()
        s = ReceiverSensor(m, position=7.0)
        r = s.sample()
        self.assertAlmostEqual(r["phase"], m.theta[7], places=9)

    def test_signal_pair_is_unit(self):
        m = locked_model()
        s = ReceiverSensor(m, position=10.5)
        r = s.sample()
        self.assertAlmostEqual(r["signal_cos"] ** 2 + r["signal_sin"] ** 2,
                               1.0, places=9)

    def test_forward_wave_k_est_and_omega(self):
        m = locked_model(drive_k=1)
        s = ReceiverSensor(m, position=10.0)
        s.sample()
        m.step()
        r = s.sample()
        # Locked k=+1 wave: gradient ~ 2*pi/N, omega ~ -drive_omega (-1.0)
        self.assertAlmostEqual(r["k_est"], 1.0, places=1)
        self.assertAlmostEqual(r["omega"], -1.0, places=1)
        self.assertGreater(r["coherence"], 0.9)

    def test_reverse_wave_k_est(self):
        m = locked_model(drive_k=-1)
        s = ReceiverSensor(m, position=10.0)
        s.sample()
        m.step()
        r = s.sample()
        self.assertAlmostEqual(r["k_est"], -1.0, places=1)
        self.assertAlmostEqual(r["omega"], -1.0, places=1)

    def test_uniform_drive_k_est_zero(self):
        m = locked_model(drive_k=0)
        s = ReceiverSensor(m, position=10.0)
        s.sample()
        m.step()
        r = s.sample()
        self.assertAlmostEqual(r["k_est"], 0.0, places=1)
        self.assertAlmostEqual(r["omega"], -1.0, places=1)

    def test_noise_lowers_coherence(self):
        clean = locked_model(drive_k=1)
        s1 = ReceiverSensor(clean, position=10.0)
        s1.sample()
        clean.step()
        r_clean = s1.sample()
        noisy = MembranePropulsion(seed=42, drive_k=1, noise=2.0)
        noisy.run(600)
        s2 = ReceiverSensor(noisy, position=10.0)
        s2.sample()
        noisy.step()
        r_noisy = s2.sample()
        self.assertLess(r_noisy["coherence"], r_clean["coherence"])

    def test_sensor_does_not_know_drive(self):
        # The sensor class takes no drive parameters; it infers
        # k and omega from the vibration alone.
        import inspect
        params = inspect.signature(ReceiverSensor.__init__).parameters
        self.assertNotIn("drive_k", params)
        self.assertNotIn("drive_omega", params)


if __name__ == "__main__":
    unittest.main(verbosity=2)
