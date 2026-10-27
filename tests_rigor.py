"""Tests for the rigor experiments (rigor.py, render_results.py), the Ring helper
and the packaging metadata. Standard library only. Run: python tests_rigor.py
"""

import json
import math
import os
import re
import unittest

import render_results
import rigor
from membrane_propulsion import MembranePropulsion
from sensor import ReceiverSensor

HERE = os.path.dirname(os.path.abspath(__file__))


def read(name):
    with open(os.path.join(HERE, name), encoding="utf-8") as fh:
        return fh.read()


class TestStatistics(unittest.TestCase):
    def test_wilson_known_values(self):
        lo, hi = rigor.wilson(0, 10)
        self.assertEqual(lo, 0.0)
        self.assertAlmostEqual(hi, 0.2775, places=3)
        lo, hi = rigor.wilson(50, 100)       # textbook: 0.404 to 0.596
        self.assertAlmostEqual(lo, 0.4038, places=3)
        self.assertAlmostEqual(hi, 0.5962, places=3)
        self.assertTrue(math.isnan(rigor.wilson(0, 0)[0]))

    def test_auc_known_values(self):
        self.assertEqual(rigor.auc([3, 4], [1, 2]), 1.0)
        self.assertEqual(rigor.auc([1, 2], [3, 4]), 0.0)
        self.assertEqual(rigor.auc([1, 2], [1, 2]), 0.5)

    def test_quantile_mean_sd(self):
        self.assertEqual(rigor.mean([1, 2, 3]), 2.0)
        self.assertAlmostEqual(rigor.sd([1, 2, 3]), 1.0)
        self.assertEqual(rigor.quantile([1, 2, 3, 4, 5], 0.5), 3)

    def test_linfit_and_pearson(self):
        a, b, r2 = rigor.linfit([0, 1, 2, 3], [1, 3, 5, 7])
        self.assertAlmostEqual(b, 2.0)
        self.assertAlmostEqual(a, 1.0)
        self.assertAlmostEqual(r2, 1.0)
        self.assertAlmostEqual(rigor.pearson([1, 2, 3], [3, 2, 1]), -1.0)

    def test_bootstrap_ci_deterministic_and_brackets_mean(self):
        xs = [float(i % 7) for i in range(50)]
        a = rigor.boot_ci(xs, 300, "t")
        b = rigor.boot_ci(xs, 300, "t")
        self.assertEqual(a, b)
        self.assertLessEqual(a[0], rigor.mean(xs))
        self.assertGreaterEqual(a[1], rigor.mean(xs))

    def test_stream_is_reproducible_and_tagged(self):
        self.assertEqual(rigor.stream("a", 1).random(), rigor.stream("a", 1).random())
        self.assertNotEqual(rigor.stream("a", 1).random(), rigor.stream("b", 1).random())

    def test_clean_rounds_and_removes_negative_zero(self):
        c = rigor.clean({"x": -0.0, "y": 1.23456789012345, "z": [float("nan")]})
        self.assertEqual(json.dumps(c["x"]), "0.0")
        self.assertEqual(c["y"], 1.2345678901)
        self.assertIsNone(c["z"][0])


class TestRing(unittest.TestCase):
    def test_ring_is_bitwise_identical_to_the_module(self):
        for kw in (dict(), dict(drive_k=-1), dict(drive_k=0), dict(drive_amp=0.0),
                   dict(noise=1.0), dict(n=32, coupling=3.0)):
            a = MembranePropulsion(seed=5, **kw)
            b = rigor.Ring(seed=5, **kw)
            for _ in range(150):
                a.step()
                b.step()
            self.assertEqual(a.theta, b.theta, kw)
            self.assertEqual(a.probe, b.probe, kw)
            self.assertEqual(a.probe_unwrapped, b.probe_unwrapped, kw)

    def test_winding_of_a_twisted_state(self):
        r = rigor.Ring(seed=1, init=("twist", 3), drive_amp=0.0)
        self.assertEqual(rigor.winding(r.theta), 3)

    def test_probe_law_zero_and_const(self):
        z = rigor.Ring(seed=1, probe_mode="zero")
        z.run(100)
        self.assertEqual(z.probe_unwrapped, 0.0)
        c = rigor.Ring(seed=1, probe_mode="const", const_v=2.0)
        c.run(100)
        self.assertAlmostEqual(c.probe_unwrapped, 2.0 * 100 * c.dt, places=9)

    def test_unknown_init_rejected(self):
        with self.assertRaises(ValueError):
            rigor.Ring(init="nope")


class TestExperiments(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a = rigor.run_all(rigor.QUICK)

    def test_run_is_deterministic(self):
        b = rigor.run_all(rigor.QUICK)
        self.assertEqual(json.dumps(self.a, sort_keys=True), json.dumps(b, sort_keys=True))

    def test_sign_control_by_construction(self):
        ms = self.a["multiseed"]
        self.assertEqual(ms["fwd"]["positive"]["rate"], 1.0)
        self.assertEqual(ms["rev"]["negative"]["rate"], 1.0)
        self.assertAlmostEqual(ms["fwd"]["stats"]["mean"], 20.0, places=2)
        self.assertAlmostEqual(ms["rev"]["stats"]["mean"], -20.0, places=2)
        self.assertLess(ms["uniform"]["max_abs"], 0.01)

    def test_nodrive_control_is_not_zero(self):
        # the honest control: without drive the ring still drifts in some seeds
        nd = self.a["multiseed"]["nodrive"]
        self.assertGreater(nd["abs_stats"]["max"], 5.0)

    def test_k_fit_is_linear_by_construction(self):
        fit = self.a["k_mu"]["k_fit"]
        self.assertAlmostEqual(fit["slope"], 20.0, places=2)
        self.assertGreater(fit["r2"], 0.99999)

    def test_reproduction_includes_sensor_table(self):
        s = self.a["reproduction"]["sensor_table"]
        self.assertAlmostEqual(s["forward k=+1"]["phase"], 3.077, places=3)
        self.assertAlmostEqual(s["reverse k=-1"]["phase"], 1.114, places=3)
        self.assertAlmostEqual(s["uniform k=0"]["phase"], 2.096, places=3)
        for v in s.values():
            self.assertAlmostEqual(v["omega"], -1.0, places=6)

    def test_readme_seed_label_is_seed_2_not_3(self):
        nd = self.a["reproduction"]["nodrive_by_seed"]
        self.assertAlmostEqual(nd["2"], -28.1, places=1)
        self.assertGreater(nd["3"], -5.0)

    def test_sensor_recovers_k_and_is_read_only(self):
        rec = self.a["sensor"]["recovery"]
        for row in rec["rows"]:
            self.assertLess(row["max_abs_k_err"], 1e-6)
            self.assertLess(row["max_abs_omega_err"], 1e-6)
        ro = self.a["sensor"]["readonly"]
        self.assertEqual(ro["snapshot_violations"], 0)
        self.assertEqual(ro["trajectory_violations"], 0)

    def test_sensor_coherence_is_geometric(self):
        for row in self.a["sensor"]["recovery"]["k_rows"]:
            self.assertAlmostEqual(row["coherence_mean"], row["geometric_coherence"], places=6)

    def test_sensor_wrap_and_seam(self):
        pos = self.a["sensor"]["position"]
        self.assertLess(pos["max_abs_k_err"], 1e-6)
        self.assertLess(pos["interp_max_err_vs_short_arc_midpoint"], 1e-9)
        self.assertGreater(pos["straddling_the_pi_wrap"], 0)

    def test_sensor_omega_moved_between_samples_is_contaminated(self):
        pos = self.a["sensor"]["position"]
        self.assertGreater(abs(pos["omega_moved_between_samples_mean"] + 1.0), 1.0)

    def test_renderer_uses_only_stored_numbers(self):
        text = render_results.render(self.a)
        self.assertIn("not propulsion", text.lower())
        self.assertIn("by construction", text.lower())
        self.assertIn("seed 2", text)


class TestStoredResults(unittest.TestCase):
    def test_results_json_matches_results_md(self):
        with open(os.path.join(HERE, "results", "results.json"), encoding="utf-8") as fh:
            res = json.load(fh)
        self.assertEqual(read("RESULTS.md"), render_results.render(res))

    def test_results_json_has_every_section(self):
        with open(os.path.join(HERE, "results", "results.json"), encoding="utf-8") as fh:
            res = json.load(fh)
        for key in ("config", "reproduction", "multiseed", "anatomy", "k_mu", "locking", "noise",
                    "warmup", "dt", "N", "nulls", "sensor"):
            self.assertIn(key, res)
        self.assertEqual(res["multiseed"]["seeds"], 200)


class TestPackaging(unittest.TestCase):
    def test_versions_agree(self):
        v = "0.2.0"
        self.assertRegex(read("pyproject.toml"), r'version = "%s"' % re.escape(v))
        self.assertRegex(read("CITATION.cff"), r'(?m)^version: %s$' % re.escape(v))
        self.assertEqual(json.loads(read(".zenodo.json"))["version"], v)
        self.assertIn("## [0.2.0]", read("CHANGELOG.md"))
        self.assertIn("## [0.1.0]", read("CHANGELOG.md"))

    def test_license_author_orcid_and_verified_dois(self):
        z = json.loads(read(".zenodo.json"))
        self.assertEqual(z["license"], "AGPL-3.0-only")
        self.assertEqual(z["creators"][0]["orcid"], "0009-0000-9778-5374")
        self.assertEqual(z["creators"][0]["affiliation"], "Spark AI NLP, Fredericton, New Brunswick, Canada")
        self.assertNotIn("doi", z)
        cff = read("CITATION.cff")
        self.assertIn("AGPL-3.0-only", cff)
        self.assertIn("0009-0000-9778-5374", cff)
        # DOIs are only the two read from the public Zenodo record of the v0.2.0 release
        self.assertRegex(cff, r'(?m)^doi: "10\.5281/zenodo\.23687826"$')
        self.assertIn("10.5281/zenodo.23687825", cff)
        self.assertNotIn("10.5281", read(".zenodo.json"))
        self.assertIn("GNU AFFERO GENERAL PUBLIC LICENSE", read("LICENSE").upper())

    def test_ci_covers_all_python_versions(self):
        ci = read(os.path.join(".github", "workflows", "tests.yml"))
        for v in ("3.10", "3.11", "3.12", "3.13"):
            self.assertIn('"%s"' % v, ci)

    def test_honesty_statement_present(self):
        for name in ("README.md", "RESULTS.md", ".zenodo.json", "CITATION.cff"):
            self.assertRegex(read(name).lower(), r"not (a )?(propulsion|thrust)|synthetic", name)

    def test_core_module_unchanged_marker(self):
        # the core file must still define the documented public API
        import membrane_propulsion as mp
        for name in ("MembranePropulsion", "run_case", "wrap", "TWO_PI"):
            self.assertTrue(hasattr(mp, name))
        s = ReceiverSensor(MembranePropulsion(seed=1))
        self.assertTrue(hasattr(s, "sample"))


if __name__ == "__main__":
    unittest.main()
