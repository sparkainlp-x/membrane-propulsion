"""Tests for membrane-propulsion. Run: python tests.py

Displacement is steady-state displacement: 200 warmup steps are
discarded first, so transient synchronisation drift is not counted
as propulsion. See DESIGN.md for the twisted-state control caveat.
"""

import math
import unittest

from membrane_propulsion import MembranePropulsion, run_case, wrap


class TestPropulsion(unittest.TestCase):
    def test_wrap(self):
        self.assertAlmostEqual(wrap(math.pi + 0.1), -math.pi + 0.1, places=9)
        self.assertAlmostEqual(wrap(-math.pi - 0.1), math.pi - 0.1, places=9)
        self.assertAlmostEqual(wrap(0.5), 0.5, places=9)

    def test_uniform_ring_probe_stationary(self):
        # All phases equal: gradient 0, probe must not move.
        m = MembranePropulsion(seed=1, drive_amp=0.0, coupling=0.0)
        m.theta = [0.3] * m.n
        disp = m.run(100)
        self.assertAlmostEqual(disp, 0.0, places=9)

    def test_forward_drive_drifts_forward(self):
        disp, wave_order, _ = run_case(seed=42, steps=400, drive_k=1)
        self.assertAlmostEqual(disp, 20.0, places=1)
        self.assertGreater(wave_order, 0.9)

    def test_reverse_drive_drifts_backward(self):
        disp, wave_order, _ = run_case(seed=42, steps=400, drive_k=-1)
        self.assertAlmostEqual(disp, -20.0, places=1)
        self.assertGreater(wave_order, 0.9)

    def test_forward_drive_sign_controlled_across_seeds(self):
        # Direction is controlled by the drive, across seeds.
        for seed in (0, 1, 2, 3, 42):
            disp, _, _ = run_case(seed=seed, steps=400, drive_k=1)
            self.assertGreater(disp, 10.0)
            disp, _, _ = run_case(seed=seed, steps=400, drive_k=-1)
            self.assertLess(disp, -10.0)

    def test_uniform_drive_no_drift(self):
        disp, _, _ = run_case(seed=42, steps=400, drive_k=0)
        self.assertLess(abs(disp), 1.0)

    def test_no_drive_sign_not_controlled(self):
        # Control, honestly stated: without a drive the ring can fall
        # into a twisted state (seed 2: -28.1 sites) and drift on its
        # own. What the drive controls is the SIGN, consistently
        # across seeds (see the test above); no-drive signs vary.
        disps = [run_case(seed=s, steps=400, drive_amp=0.0)[0]
                 for s in (42, 0, 1, 2, 3)]
        self.assertTrue(any(d > 0 for d in disps))
        self.assertTrue(any(d < 0 for d in disps))

    def test_determinism(self):
        d1, _, _ = run_case(seed=7, steps=200, drive_k=1)
        d2, _, _ = run_case(seed=7, steps=200, drive_k=1)
        self.assertEqual(d1, d2)

    def test_strong_noise_collapses_drift_seed42(self):
        clean, _, _ = run_case(seed=42, steps=400, drive_k=1, noise=0.0)
        noisy, _, _ = run_case(seed=42, steps=400, drive_k=1, noise=2.0)
        self.assertLess(abs(noisy), abs(clean) * 0.7)

    def test_probe_wraps_in_range(self):
        m = MembranePropulsion(seed=3, drive_k=1)
        m.run(400)
        self.assertGreaterEqual(m.probe, 0.0)
        self.assertLess(m.probe, m.n)


if __name__ == "__main__":
    unittest.main(verbosity=2)
