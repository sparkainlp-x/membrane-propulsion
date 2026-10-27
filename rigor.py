"""Rigor experiments for membrane-propulsion v0.1.0 (stdlib only, seeded).

SYNTHETIC toy. Nothing here is physics, propulsion, thrust, force, energy or a
property of any real system. "Displacement" is the drift of a probe whose
velocity is, BY ASSUMPTION, proportional to a local phase gradient.
Every number in RESULTS.md is produced by ``run_all`` below and stored in
``results/results.json``.

Run:  python rigor.py --write        # full run, writes results/results.json and RESULTS.md
      python rigor.py --render-only  # re-render RESULTS.md from the JSON
      python rigor.py --quick        # small run, prints a summary only

membrane_propulsion.py is imported unchanged. The ``Ring`` subclass below only
adds (a) a history-free step that is bit-for-bit identical to the original in
its default mode (tested) and (b) the control variants (other probe laws,
other initial conditions).
"""

import argparse
import json
import math
import os
import random
import sys
import zlib

from membrane_propulsion import MembranePropulsion, TWO_PI, run_case, wrap
from sensor import ReceiverSensor

HERE = os.path.dirname(os.path.abspath(__file__))
STEP_TIME = 0.05          # default dt
T_WARM = 200 * STEP_TIME  # default warm-up, in time units (10.0)
T_MEAS = 400 * STEP_TIME  # default measurement window (20.0)

FULL = {
    "seeds": 200, "boot": 5000, "perm": 2000,
    "anatomy_seeds": 100, "anatomy_long_seeds": 40, "anatomy_warm_time": [10, 50, 200, 1000],
    "k_values": [-3, -2, -1, 0, 1, 2, 3], "k_seeds": 30, "mu_values": [0.0, 0.5, 1.0, 2.0],
    "A_values": [0.0, 0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0, 3.0], "lock_seeds": 40,
    "tongue_A": [0.0, 0.5, 1.0, 1.5, 2.0, 3.0], "tongue_omega": [0.0, 0.5, 1.0, 1.5, 2.0, 3.0],
    "tongue_seeds": 10, "K_values": [0.0, 0.5, 1.0, 2.0, 4.0, 8.0], "K_seeds": 30,
    "noise_values": [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0], "noise_seeds": 100,
    "warm_steps": [0, 20, 50, 100, 200, 400, 800, 2000], "warm_seeds": 60,
    "dt_values": [0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625], "dt_seeds": 20,
    "N_values": [16, 32, 64, 128, 256], "N_seeds": 40,
    "init_noise_seeds": 40,
    "sens_seeds": 200, "sens_k_values": [-5, -3, -2, -1, 0, 1, 2, 3, 5], "sens_k_seeds": 40,
    "sens_noise": [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0], "sens_noise_seeds": 100,
    "sens_pos_seeds": 20, "sens_intervals": [1, 10, 30, 50, 60, 70, 100],
    "sens_thresholds": [0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999],
}
QUICK = {
    "seeds": 8, "boot": 200, "perm": 200,
    "anatomy_seeds": 6, "anatomy_long_seeds": 3, "anatomy_warm_time": [10, 50],
    "k_values": [-1, 0, 1], "k_seeds": 3, "mu_values": [0.0, 1.0],
    "A_values": [0.0, 1.5], "lock_seeds": 3,
    "tongue_A": [0.0, 1.5], "tongue_omega": [0.0, 1.0], "tongue_seeds": 2,
    "K_values": [0.0, 2.0], "K_seeds": 3,
    "noise_values": [0.0, 2.0], "noise_seeds": 4,
    "warm_steps": [0, 200], "warm_seeds": 4,
    "dt_values": [0.1, 0.05], "dt_seeds": 2,
    "N_values": [16, 64], "N_seeds": 3,
    "init_noise_seeds": 3,
    "sens_seeds": 6, "sens_k_values": [-1, 0, 1], "sens_k_seeds": 3,
    "sens_noise": [0.0, 2.0], "sens_noise_seeds": 6,
    "sens_pos_seeds": 2, "sens_intervals": [1, 70],
    "sens_thresholds": [0.9, 0.99],
}


# ---------------------------------------------------------------- statistics
def stream(tag, seed=0):
    return random.Random(zlib.crc32(tag.encode()) * 1_000_003 + int(seed))


def wilson(k, n, z=1.959964):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, centre - half), min(1.0, centre + half))


def mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def sd(xs):
    xs = list(xs)
    if len(xs) < 2:
        return 0.0
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def quantile(xs, q):
    s = sorted(xs)
    if not s:
        return float("nan")
    pos = q * (len(s) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (pos - lo)


def boot_ci(values, boot, tag="boot", stat=mean):
    """Percentile bootstrap 95% interval of a statistic of `values`."""
    values = list(values)
    if len(values) < 2:
        v = stat(values) if values else float("nan")
        return (v, v)
    rng = stream(tag, len(values))
    n = len(values)
    reps = sorted(stat([values[rng.randrange(n)] for _ in range(n)]) for _ in range(boot))
    return (reps[int(0.025 * boot)], reps[min(boot - 1, int(0.975 * boot))])


def pearson(a, b):
    ma, mb = mean(a), mean(b)
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    if sa == 0 or sb == 0:
        return float("nan")
    return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / (sa * sb)


def linfit(x, y):
    """Least-squares y = a + b x; returns (a, b, r2)."""
    mx, my = mean(x), mean(y)
    sxx = sum((u - mx) ** 2 for u in x)
    if sxx == 0:
        return (my, 0.0, float("nan"))
    b = sum((u - mx) * (v - my) for u, v in zip(x, y)) / sxx
    a = my - b * mx
    ss_tot = sum((v - my) ** 2 for v in y)
    ss_res = sum((v - (a + b * u)) ** 2 for u, v in zip(x, y))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return (a, b, r2)


def clean(x):
    """JSON-safe copy: inf/nan -> None; floats rounded to 10 decimals so stored
    results do not depend on float summation details (sum() became compensated
    in Python 3.12)."""
    if isinstance(x, float):
        return round(x, 10) + 0.0 if math.isfinite(x) else None
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    return x


def rate(k, n, boot=None, tag="rate"):
    lo, hi = wilson(k, n)
    return {"k": k, "n": n, "rate": k / n if n else float("nan"), "wilson": [lo, hi]}


# ---------------------------------------------------------------- the ring
class Ring(MembranePropulsion):
    """MembranePropulsion with a history-free step and control variants.

    probe_mode: 'gradient' (the original law), 'zero' (v = 0), 'const'
    (v = const_v whatever the ring does), 'mean_gradient' (v from the ring's
    mean wrapped gradient instead of the gradient under the probe),
    'random_site' (gradient at a uniformly random site each step).
    init: 'random' (the original), 'zero', 'lock' (theta = drive phase),
    ('twist', w) (winding number w, theta_i = 2 pi w i / N).
    In the default modes ``step`` performs exactly the arithmetic of the
    original (tests compare bit for bit).
    """

    def __init__(self, *args, probe_mode="gradient", init="random", const_v=1.0, **kw):
        super().__init__(*args, **kw)
        self.probe_mode = probe_mode
        self.const_v = const_v
        n = self.n
        if init == "zero":
            self.theta = [0.0] * n
        elif init == "lock":
            self.theta = [wrap(self._drive_phase(i)) for i in range(n)]
        elif isinstance(init, tuple) and init[0] == "twist":
            self.theta = [wrap(TWO_PI * init[1] * i / n) for i in range(n)]
        elif init != "random":
            raise ValueError("unknown init %r" % (init,))
        self._site_rng = random.Random(zlib.crc32(b"site") + 12345)

    def step(self):
        n = self.n
        dt = self.dt
        th = self.theta
        t = self.t
        sin = math.sin
        pi = math.pi
        k_ = self.drive_k
        om = self.drive_omega
        A = self.drive_amp
        K = self.coupling
        w0 = self.omega0
        noise = self.noise
        sq = math.sqrt(dt)
        gauss = self.rng.gauss
        new = []
        for i in range(n):
            ct = (sin(th[(i - 1) % n] - th[i]) + sin(th[(i + 1) % n] - th[i])) / 2.0
            dp = k_ * TWO_PI * i / n - om * t
            dr = A * sin(dp - th[i])
            nz = 0.0
            if noise:
                nz = noise * sq * gauss(0.0, 1.0)
            x = th[i] + dt * (w0 + K * ct + dr) + nz
            new.append(x if -pi < x <= pi else wrap(x))
        self.theta = new
        self.t += dt
        mode = self.probe_mode
        if mode == "gradient":
            idx = int(self.probe) % n
            nxt = (idx + 1) % n
            g = wrap(new[nxt] - new[idx])
            v = self.mu * g * n / TWO_PI
        elif mode == "zero":
            v = 0.0
        elif mode == "const":
            v = self.const_v
        elif mode == "mean_gradient":
            v = self.mu * self.mean_gradient() * n / TWO_PI
        elif mode == "random_site":
            idx = self._site_rng.randrange(n)
            g = wrap(new[(idx + 1) % n] - new[idx])
            v = self.mu * g * n / TWO_PI
        else:
            raise ValueError(mode)
        self.probe = (self.probe + v * dt) % n
        self.probe_unwrapped += v * dt
        return v


def winding(theta):
    """Integer winding number of a phase ring: sum of wrapped neighbour differences / 2 pi."""
    n = len(theta)
    s = 0.0
    for i in range(n):
        s += wrap(theta[(i + 1) % n] - theta[i])
    return int(round(s / TWO_PI))


def mean_wave_order(m, steps, every=10):
    """Run `steps` steps; return (displacement, mean wave-order over the window, final wave-order)."""
    start = m.probe_unwrapped
    acc, cnt = 0.0, 0
    for s in range(steps):
        m.step()
        if (s + 1) % every == 0:
            acc += m.wave_order()
            cnt += 1
    return m.probe_unwrapped - start, (acc / cnt if cnt else m.wave_order()), m.wave_order()


def trial(seed, steps=400, warmup=200, **kw):
    """One run on the fast ring. kw: MembranePropulsion arguments + probe_mode, init, const_v."""
    m = Ring(seed=seed, **kw)
    w0 = winding(m.theta)
    if warmup:
        m.run(warmup)
    w1 = winding(m.theta)
    disp, wo_mean, wo_end = mean_wave_order(m, steps)
    return {"disp": disp, "wo": wo_end, "wo_mean": wo_mean, "w_init": w0, "w_warm": w1,
            "w_end": winding(m.theta)}


def disp_only(seed, steps=400, warmup=200, **kw):
    m = Ring(seed=seed, **kw)
    if warmup:
        m.run(warmup)
    return m.run(steps)


def summ(xs, boot, tag):
    xs = list(xs)
    return {"n": len(xs), "mean": mean(xs), "sd": sd(xs), "min": min(xs), "max": max(xs),
            "q05": quantile(xs, 0.05), "q25": quantile(xs, 0.25), "median": quantile(xs, 0.5),
            "q75": quantile(xs, 0.75), "q95": quantile(xs, 0.95),
            "mean_ci": list(boot_ci(xs, boot, tag))}


# ---------------------------------------------------------------- experiments
def e_reproduction():
    """The README / DESIGN numbers, recomputed (module code, not the fast ring)."""
    rows = {}
    for label, kw in [("forward k=+1", dict(drive_k=1)), ("reverse k=-1", dict(drive_k=-1)),
                      ("uniform k=0", dict(drive_k=0)), ("no drive", dict(drive_amp=0.0)),
                      ("forward + noise 2.0", dict(drive_k=1, noise=2.0))]:
        d, wo, _ = run_case(seed=42, steps=400, **kw)
        rows[label] = {"disp": d, "wave_order": wo}
    order = (42, 0, 1, 2, 3)
    nodrive = {str(s): run_case(seed=s, steps=400, drive_amp=0.0)[0] for s in order}
    noisy = {str(s): run_case(seed=s, steps=400, drive_k=1, noise=2.0)[0] for s in order}
    fwd = {str(s): run_case(seed=s, steps=400, drive_k=1)[0] for s in (0, 1, 2, 3, 42)}
    rev = {str(s): run_case(seed=s, steps=400, drive_k=-1)[0] for s in (0, 1, 2, 3, 42)}
    sensor_rows = {}
    for label, k in (("forward k=+1", 1), ("reverse k=-1", -1), ("uniform k=0", 0)):
        m = MembranePropulsion(seed=42, drive_k=k)
        m.run(600)
        sn = ReceiverSensor(m, position=10.0)
        sn.sample()
        m.step()
        r = sn.sample()
        sensor_rows[label] = {x: r[x] for x in ("phase", "omega", "k_est", "coherence")}
    return {"sensor_table": sensor_rows, "table": rows, "nodrive_by_seed": nodrive, "noisy_by_seed": noisy,
            "fwd_by_seed": fwd, "rev_by_seed": rev,
            "uniform_no_warmup": run_case(seed=42, steps=400, warmup=0, drive_k=0)[0],
            "uniform_with_warmup": run_case(seed=42, steps=400, warmup=200, drive_k=0)[0],
            "T_measure": 400 * STEP_TIME, "T_warmup": 200 * STEP_TIME}


def e_multiseed(cfg):
    """Sign control, spread and the no-drive control over many seeds."""
    S, boot = cfg["seeds"], cfg["boot"]
    cond = {"fwd": dict(drive_k=1), "rev": dict(drive_k=-1), "uniform": dict(drive_k=0),
            "nodrive": dict(drive_amp=0.0), "noise2": dict(drive_k=1, noise=2.0)}
    raw = {c: [trial(s, **kw) for s in range(S)] for c, kw in cond.items()}
    T = 400 * STEP_TIME
    out = {"seeds": S, "T": T, "raw": {c: [{"disp": r["disp"], "wo": r["wo"], "wo_mean": r["wo_mean"],
                                              "w_end": r["w_end"], "w_init": r["w_init"]} for r in v]
                                        for c, v in raw.items()}}
    f = [r["disp"] for r in raw["fwd"]]
    r_ = [r["disp"] for r in raw["rev"]]
    u = [r["disp"] for r in raw["uniform"]]
    nd = [r["disp"] for r in raw["nodrive"]]
    n2 = [r["disp"] for r in raw["noise2"]]
    out["fwd"] = {"stats": summ(f, boot, "ms-fwd"), "positive": rate(sum(d > 0 for d in f), S),
                  "gt10": rate(sum(d > 10 for d in f), S),
                  "max_abs_dev_from_mu_k_T": max(abs(d - 1.0 * 1 * T) for d in f),
                  "wave_order_min": min(r["wo"] for r in raw["fwd"]),
                  "locked_wo_gt_0.9": rate(sum(r["wo_mean"] > 0.9 for r in raw["fwd"]), S)}
    out["rev"] = {"stats": summ(r_, boot, "ms-rev"), "negative": rate(sum(d < 0 for d in r_), S),
                  "lt_minus10": rate(sum(d < -10 for d in r_), S),
                  "max_abs_dev_from_mu_k_T": max(abs(d + T) for d in r_),
                  "locked_wo_gt_0.9": rate(sum(r["wo_mean"] > 0.9 for r in raw["rev"]), S)}
    out["both_signs_right"] = rate(sum(1 for a, b in zip(f, r_) if a > 10 and b < -10), S)
    out["uniform"] = {"stats": summ(u, boot, "ms-uni"), "abs_lt_1": rate(sum(abs(d) < 1 for d in u), S),
                      "max_abs": max(abs(d) for d in u)}
    abs_nd = [abs(d) for d in nd]
    out["nodrive"] = {
        "stats": summ(nd, boot, "ms-nd"), "abs_stats": summ(abs_nd, boot, "ms-nd-abs"),
        "abs_ge_10": rate(sum(a >= 10 for a in abs_nd), S),
        "abs_ge_20": rate(sum(a >= 20 for a in abs_nd), S),
        "positive": rate(sum(d > 0 for d in nd), S),
        "nonzero_end_winding": rate(sum(r["w_end"] != 0 for r in raw["nodrive"]), S),
        "mean_abs_over_driven": mean(abs_nd) / T,
        "wave_order_mean": mean(r["wo"] for r in raw["nodrive"]),
        "max_abs": max(abs_nd)}
    out["noise2"] = {"stats": summ(n2, boot, "ms-n2"), "positive": rate(sum(d > 0 for d in n2), S),
                     "gt10": rate(sum(d > 10 for d in n2), S),
                     "wave_order_mean": mean(r["wo"] for r in raw["noise2"])}
    return out


def e_anatomy(cfg):
    """What the no-drive 'twisted state' is: winding numbers, relaxation, local sampling."""
    S = cfg["anatomy_seeds"]
    w_init = [winding(Ring(seed=s, drive_amp=0.0).theta) for s in range(S)]
    dist = {}
    for w in w_init:
        dist[str(w)] = dist.get(str(w), 0) + 1
    out = {"seeds": S, "init_winding_counts": dist, "init_winding_zero": rate(sum(w == 0 for w in w_init), S),
           "relaxation_time_estimate": 64 ** 2 / (math.pi ** 2 * 2.0), "by_warmup_time": []}
    for T in cfg["anatomy_warm_time"]:
        seeds = S if T <= 200 else cfg["anatomy_long_seeds"]
        warm = int(round(T / STEP_TIME))
        rows = [trial(s, warmup=warm, drive_amp=0.0) for s in range(seeds)]
        ds = [r["disp"] for r in rows]
        ew = [r["w_end"] for r in rows]
        # drift expected if the probe sat in a uniformly twisted state of winding w: mu * w * T
        pred = [1.0 * w * 400 * STEP_TIME for w in ew]
        out["by_warmup_time"].append({
            "warmup_time": T, "seeds": seeds,
            "nonzero_winding": rate(sum(w != 0 for w in ew), seeds),
            "abs_ge_10": rate(sum(abs(d) >= 10 for d in ds), seeds),
            "mean_abs_disp": mean(abs(d) for d in ds), "max_abs_disp": max(abs(d) for d in ds),
            "mean_abs_w_end": mean(abs(w) for w in ew),
            "corr_disp_vs_w_T": pearson(ds, pred) if seeds > 2 else float("nan"),
            "mean_abs_gap_to_uniform_prediction": mean(abs(d - p) for d, p in zip(ds, pred)),
            "winding_counts": {str(w): ew.count(w) for w in sorted(set(ew))}})
    return out


def e_k_mu(cfg):
    """Displacement against k and mu: linear BY CONSTRUCTION; this verifies it."""
    T = 400 * STEP_TIME
    rows = []
    for k in cfg["k_values"]:
        ds = [disp_only(s, drive_k=k) for s in range(cfg["k_seeds"])]
        rows.append({"k": k, "mean": mean(ds), "sd": sd(ds), "min": min(ds), "max": max(ds),
                     "predicted_mu_k_T": 1.0 * k * T})
    a, b, r2 = linfit([r["k"] for r in rows], [r["mean"] for r in rows])
    mus = []
    for mu in cfg["mu_values"]:
        ds = [disp_only(s, drive_k=1, mu=mu) for s in range(cfg["k_seeds"])]
        mus.append({"mu": mu, "mean": mean(ds), "sd": sd(ds), "predicted_mu_k_T": mu * 1 * T})
    a2, b2, r22 = linfit([r["mu"] for r in mus], [r["mean"] for r in mus])
    return {"T": T, "k_rows": rows, "k_fit": {"intercept": a, "slope": b, "r2": r2, "slope_expected": T},
            "mu_rows": mus, "mu_fit": {"intercept": a2, "slope": b2, "r2": r22, "slope_expected": T}}


def e_locking(cfg):
    """Locking against drive strength, drive frequency and coupling."""
    S = cfg["lock_seeds"]
    rows = []
    for A in cfg["A_values"]:
        rs = [trial(s, drive_k=1, drive_amp=A) for s in range(S)]
        ds = [r["disp"] for r in rs]
        rows.append({"A": A, "locked": rate(sum(r["wo_mean"] > 0.9 for r in rs), S),
                     "mean_wo": mean(r["wo_mean"] for r in rs), "mean_disp": mean(ds),
                     "mean_abs_disp": mean(abs(d) for d in ds),
                     "positive": rate(sum(d > 0 for d in ds), S)})
    tongue = []
    TS = cfg["tongue_seeds"]
    for A in cfg["tongue_A"]:
        row = []
        for om in cfg["tongue_omega"]:
            rs = [trial(s, drive_k=1, drive_amp=A, drive_omega=om) for s in range(TS)]
            row.append({"omega": om, "locked_frac": sum(r["wo_mean"] > 0.9 for r in rs) / TS,
                        "mean_wo": mean(r["wo_mean"] for r in rs)})
        tongue.append({"A": A, "cells": row})
    ks = []
    for A in (0.5, 1.5):
        for K in cfg["K_values"]:
            rs = [trial(s, drive_k=1, drive_amp=A, coupling=K) for s in range(cfg["K_seeds"])]
            ks.append({"A": A, "K": K, "locked": rate(sum(r["wo_mean"] > 0.9 for r in rs), cfg["K_seeds"]),
                       "mean_wo": mean(r["wo_mean"] for r in rs)})
    return {"seeds": S, "A_rows": rows, "tongue": tongue, "tongue_seeds": TS, "K_rows": ks,
            "predicted_threshold_A_eq_abs_omega_minus_omega0": 1.0,
            "predicted_lag_rad_at_A1p5": math.asin(1.0 / 1.5)}


def e_noise(cfg):
    S = cfg["noise_seeds"]
    rows = []
    for sg in cfg["noise_values"]:
        f = [trial(s, drive_k=1, noise=sg) for s in range(S)]
        r_ = [trial(s, drive_k=-1, noise=sg) for s in range(S)]
        fd = [x["disp"] for x in f]
        rd = [x["disp"] for x in r_]
        rows.append({"sigma": sg,
                     "fwd_positive": rate(sum(d > 0 for d in fd), S),
                     "rev_negative": rate(sum(d < 0 for d in rd), S),
                     "fwd_gt10": rate(sum(d > 10 for d in fd), S),
                     "both_gt10": rate(sum(1 for a, b in zip(fd, rd) if a > 10 and b < -10), S),
                     "fwd_mean": mean(fd), "fwd_sd": sd(fd), "rev_mean": mean(rd),
                     "mean_wo": mean(x["wo_mean"] for x in f),
                     "fwd_ci": list(boot_ci(fd, 2000, "noise-%s" % sg))})
    return {"seeds": S, "rows": rows}


def e_warmup(cfg):
    S = cfg["warm_seeds"]
    rows = []
    for w in cfg["warm_steps"]:
        row = {"warmup_steps": w, "warmup_time": w * STEP_TIME}
        for name, kw in (("fwd", dict(drive_k=1)), ("uniform", dict(drive_k=0)),
                         ("nodrive", dict(drive_amp=0.0))):
            ds = [disp_only(s, warmup=w, **kw) for s in range(S)]
            row[name] = {"mean": mean(ds), "mean_abs": mean(abs(d) for d in ds),
                         "max_abs": max(abs(d) for d in ds)}
        rows.append(row)
    return {"seeds": S, "rows": rows}


def e_dt(cfg):
    """Euler step convergence at fixed physical durations (warm-up 10, measure 20)."""
    S = cfg["dt_seeds"]
    rows = []
    for dt in cfg["dt_values"]:
        wu, st = int(round(T_WARM / dt)), int(round(T_MEAS / dt))
        row = {"dt": dt, "steps": st}
        for name, kw in (("fwd", dict(drive_k=1)), ("rev", dict(drive_k=-1)), ("uniform", dict(drive_k=0)),
                         ("nodrive", dict(drive_amp=0.0))):
            rs = [trial(s, steps=st, warmup=wu, dt=dt, **kw) for s in range(S)]
            ds = [r["disp"] for r in rs]
            row[name] = {"mean": mean(ds), "sd": sd(ds), "min": min(ds), "max": max(ds),
                         "mean_wo": mean(r["wo_mean"] for r in rs),
                         "finite": all(math.isfinite(d) for d in ds)}
        rows.append(row)
    # same seeds, nodrive: how much does the individual trajectory move with dt (vs dt=0.00625 or smallest)?
    ref_dt = min(cfg["dt_values"])
    ref = [disp_only(s, steps=int(round(T_MEAS / ref_dt)), warmup=int(round(T_WARM / ref_dt)), dt=ref_dt,
                     drive_amp=0.0) for s in range(S)]
    for row in rows:
        dt = row["dt"]
        cur = [disp_only(s, steps=int(round(T_MEAS / dt)), warmup=int(round(T_WARM / dt)), dt=dt,
                         drive_amp=0.0) for s in range(S)]
        row["nodrive_mean_abs_diff_to_finest"] = mean(abs(a - b) for a, b in zip(cur, ref))
        row["nodrive_sign_agree_with_finest"] = sum((a > 0) == (b > 0) for a, b in zip(cur, ref)) / S
    return {"seeds": S, "rows": rows, "finest_dt": ref_dt, "T_warm": T_WARM, "T_meas": T_MEAS}


def e_N(cfg):
    S = cfg["N_seeds"]
    rows = []
    for n in cfg["N_values"]:
        f = [trial(s, n=n, drive_k=1) for s in range(S)]
        nd = [trial(s, n=n, drive_amp=0.0) for s in range(S)]
        fd = [x["disp"] for x in f]
        nds = [x["disp"] for x in nd]
        rows.append({"N": n, "fwd_mean": mean(fd), "fwd_min": min(fd), "fwd_max": max(fd),
                     "fwd_locked": rate(sum(x["wo_mean"] > 0.9 for x in f), S),
                     "nodrive_abs_ge_10": rate(sum(abs(d) >= 10 for d in nds), S),
                     "nodrive_mean_abs": mean(abs(d) for d in nds),
                     "nodrive_nonzero_winding": rate(sum(x["w_end"] != 0 for x in nd), S),
                     "nodrive_init_nonzero_winding": rate(sum(x["w_init"] != 0 for x in nd), S)})
    return {"seeds": S, "rows": rows}


def e_nulls(cfg, multiseed):
    S = cfg["seeds"]
    rng = stream("null-shuffle", 0)
    f = [r["disp"] for r in multiseed["raw"]["fwd"]]
    r_ = [r["disp"] for r in multiseed["raw"]["rev"]]
    nd = [r["disp"] for r in multiseed["raw"]["nodrive"]]
    real = [(1, d) for d in f] + [(-1, d) for d in r_]
    agree_real = sum(1 for lab, d in real if d * lab > 0) / len(real)
    perm = []
    labels = [lab for lab, _ in real]
    for _ in range(cfg["perm"]):
        sh = labels[:]
        rng.shuffle(sh)
        perm.append(sum(1 for lab, (_, d) in zip(sh, real) if d * lab > 0) / len(real))
    perm.sort()
    # no-drive displacements scored against random +/-1 labels (no drive -> no label to follow)
    nd_labels_rate = []
    for _ in range(cfg["perm"]):
        labs = [rng.choice((1, -1)) for _ in nd]
        nd_labels_rate.append(sum(1 for lab, d in zip(labs, nd) if d * lab > 0) / len(nd))
    nd_labels_rate.sort()
    out = {"sign_agreement_real": agree_real, "sign_agreement_shuffled_mean": mean(perm),
           "sign_agreement_shuffled_q025": perm[int(0.025 * len(perm))],
           "sign_agreement_shuffled_q975": perm[min(len(perm) - 1, int(0.975 * len(perm)))],
           "nodrive_vs_random_labels_mean": mean(nd_labels_rate),
           "nodrive_vs_random_labels_q025": nd_labels_rate[int(0.025 * len(nd_labels_rate))],
           "nodrive_vs_random_labels_q975": nd_labels_rate[min(len(nd_labels_rate) - 1, int(0.975 * len(nd_labels_rate)))]}
    # probe-law variants
    PS = min(S, 40)
    laws = []
    for mode in ("gradient", "zero", "const", "mean_gradient", "random_site"):
        row = {"probe_mode": mode}
        for name, kw in (("fwd", dict(drive_k=1)), ("rev", dict(drive_k=-1)), ("nodrive", dict(drive_amp=0.0))):
            ds = [disp_only(s, probe_mode=mode, **kw) for s in range(PS)]
            row[name] = {"mean": mean(ds), "min": min(ds), "max": max(ds)}
        row["sign_follows_drive"] = rate(
            sum(1 for s in range(PS) if disp_only(s, probe_mode=mode, drive_k=1) > 1
                and disp_only(s, probe_mode=mode, drive_k=-1) < -1), PS)
        laws.append(row)
    out["probe_laws"] = laws
    # initial-condition variants (noise 0: deterministic given the init) and with noise 1.0
    inits = [("random", "random"), ("zero", "zero"), ("lock", "lock"), ("twist -1", ("twist", -1)),
             ("twist +2", ("twist", 2)), ("twist -5", ("twist", -5))]
    ir = []
    NS = cfg["init_noise_seeds"]
    for name, init in inits:
        d0 = disp_only(0, drive_k=1, init=init)
        wo0 = trial(0, drive_k=1, init=init)["wo_mean"]
        dn = [disp_only(s, drive_k=1, init=init, noise=1.0) for s in range(NS)]
        ir.append({"init": name, "noise0_disp": d0, "noise0_wave_order": wo0,
                   "noise1_positive": rate(sum(d > 0 for d in dn), NS), "noise1_mean": mean(dn)})
    out["initial_conditions"] = ir
    return out



# ---------------------------------------------------------------- sensor
def auc(pos, neg):
    """P(pos > neg) + 0.5 P(tie) (Mann-Whitney), 0.5 = no separation."""
    if not pos or not neg:
        return float("nan")
    gt = eq = 0
    for p in pos:
        for q in neg:
            if p > q:
                gt += 1
            elif p == q:
                eq += 1
    return (gt + 0.5 * eq) / (len(pos) * len(neg))


def sensor_read(model, position=10.0, window=2):
    """README procedure: sample, one step, sample; returns the second reading."""
    s = ReceiverSensor(model, position=position, window=window)
    s.sample()
    model.step()
    return s.sample()


def state_snapshot(m):
    return (list(m.theta), m.probe, m.probe_unwrapped, m.t, m.rng.getstate(), list(m.history))


def e_sensor_recovery(cfg):
    """k_est and omega against the TRUE k and Omega, locked rings, many seeds."""
    S, boot = cfg["sens_seeds"], cfg["boot"]
    rows = []
    for k in (-1, 0, 1):
        errk, erro, coh = [], [], []
        for s in range(S):
            m = MembranePropulsion(seed=s, drive_k=k)
            m.run(600)
            r = sensor_read(m)
            errk.append(r["k_est"] - k)
            erro.append(r["omega"] - (-1.0))
            coh.append(r["coherence"])
        rows.append({"k": k, "k_err": summ(errk, boot, "sk%d" % k), "omega_err": summ(erro, boot, "so%d" % k),
                     "max_abs_k_err": max(abs(x) for x in errk), "max_abs_omega_err": max(abs(x) for x in erro),
                     "k_err_lt_0p1": rate(sum(abs(x) < 0.1 for x in errk), S),
                     "coherence": summ(coh, boot, "sc%d" % k)})
    # omega is -Omega for other Omega values too (frequency locking)
    om_rows = []
    for om in (0.25, 0.5, 1.0, 1.4):
        errs = []
        for s in range(min(S, 40)):
            m = MembranePropulsion(seed=s, drive_k=1, drive_omega=om)
            m.run(600)
            errs.append(sensor_read(m)["omega"] + om)
        om_rows.append({"Omega": om, "mean_err": mean(errs), "max_abs_err": max(abs(e) for e in errs)})
    # larger |k|
    kr = []
    for k in cfg["sens_k_values"]:
        errk, coh = [], []
        for s in range(cfg["sens_k_seeds"]):
            m = MembranePropulsion(seed=s, drive_k=k)
            m.run(600)
            r = sensor_read(m)
            errk.append(r["k_est"] - k)
            coh.append(r["coherence"])
        kr.append({"k": k, "max_abs_k_err": max(abs(e) for e in errk), "mean_k_err": mean(errk),
                   "coherence_mean": mean(coh), "coherence_min": min(coh),
                   "geometric_coherence": math.sin(5 * math.pi * k / 64) / (5 * math.sin(math.pi * k / 64)) if k else 1.0})
    return {"seeds": S, "rows": rows, "omega_rows": om_rows, "k_rows": kr}


def e_sensor_readonly(cfg):
    """Sampling must not change the model: state snapshot and trajectory comparisons."""
    S = cfg["sens_seeds"]
    snap_violations = traj_violations = 0
    cases = 0
    for s in range(S):
        for kw in (dict(drive_k=1), dict(drive_amp=0.0), dict(drive_k=-1, noise=2.0)):
            m = MembranePropulsion(seed=s, **kw)
            m.run(100)
            before = state_snapshot(m)
            sen = ReceiverSensor(m, position=13.37)
            sen.sample()
            sen.follow_probe()
            sen.sample()
            after = state_snapshot(m)
            snap_violations += before != after
            # trajectory with a sensor sampling every step vs none
            a = MembranePropulsion(seed=s, **kw)
            b = MembranePropulsion(seed=s, **kw)
            sb = ReceiverSensor(b, position=3.0)
            for i in range(120):
                a.step()
                b.step()
                sb.sample()
                if i % 7 == 0:
                    sb.follow_probe()
            traj_violations += (a.theta != b.theta or a.probe != b.probe or a.rng.getstate() != b.rng.getstate())
            cases += 1
    return {"cases": cases, "snapshot_violations": snap_violations, "trajectory_violations": traj_violations}


def e_sensor_states(cfg):
    """Sensor on states that are not a clean locked wave: no drive (twisted/relaxing), unlocked drive."""
    S = cfg["sens_seeds"]
    boot = cfg["boot"]
    nd = []
    for s in range(S):
        m = MembranePropulsion(seed=s, drive_amp=0.0)
        m.run(600)
        w = winding(m.theta)
        ks = []
        sen = ReceiverSensor(m, position=0.0)
        for p in range(m.n):
            sen.position = float(p)
            ks.append(sen.sample()["gradient"] * m.n / TWO_PI)
        r = sensor_read(m)
        nd.append({"w": w, "k_est_10": r["k_est"], "k_mean_all_sites": mean(ks), "k_sd_sites": sd(ks),
                   "omega": r["omega"], "coherence": r["coherence"]})
    ndw = [r["w"] for r in nd]
    out = {"nodrive": {
        "seeds": S, "mean_all_sites_equals_winding_max_gap": max(abs(r["k_mean_all_sites"] - r["w"]) for r in nd),
        "k_est_single_site_abs_err_vs_winding": summ([abs(r["k_est_10"] - r["w"]) for r in nd], boot, "snd-err"),
        "k_est_single_site_within_0p5_of_winding": rate(sum(abs(r["k_est_10"] - r["w"]) < 0.5 for r in nd), S),
        "k_est_single_site_rounds_to_winding": rate(sum(round(r["k_est_10"]) == r["w"] for r in nd), S),
        "nonzero_winding": rate(sum(w != 0 for w in ndw), S),
        "omega": summ([r["omega"] for r in nd], boot, "snd-om"),
        "coherence": summ([r["coherence"] for r in nd], boot, "snd-coh"),
        "coherence_gt_0p99": rate(sum(r["coherence"] > 0.99 for r in nd), S)}}
    # no-drive vs driven: can single readings tell them apart?
    drv = []
    for s in range(S):
        m = MembranePropulsion(seed=s, drive_k=1)
        m.run(600)
        drv.append(sensor_read(m))
    out["driven_vs_nodrive"] = {
        "auc_coherence_driven_higher": auc([r["coherence"] for r in drv], [r["coherence"] for r in nd]),
        "auc_abs_omega_driven_higher": auc([abs(r["omega"]) for r in drv], [abs(r["omega"]) for r in nd]),
        "nodrive_abs_omega_gt_0p5": rate(sum(abs(r["omega"]) > 0.5 for r in nd), S)}
    # unlocked drive (A < Omega): what does the sensor read?
    ul = []
    for A in (0.25, 0.5, 0.9, 1.1):
        rows = []
        for s in range(min(S, 40)):
            m = MembranePropulsion(seed=s, drive_k=1, drive_amp=A)
            m.run(600)
            rows.append(sensor_read(m))
        ul.append({"A": A, "omega_mean": mean(r["omega"] for r in rows), "omega_sd": sd(r["omega"] for r in rows),
                   "k_est_mean": mean(r["k_est"] for r in rows), "k_est_sd": sd(r["k_est"] for r in rows),
                   "coherence_mean": mean(r["coherence"] for r in rows)})
    out["partial_drive"] = ul
    return out


def e_sensor_noise(cfg):
    """Noise sweep: reading quality, and can coherence tell locked from unlocked?"""
    S = cfg["sens_noise_seeds"]
    rows = []
    pool = {}
    nd_ref = []
    for s in range(S):
        m = MembranePropulsion(seed=s, drive_amp=0.0)
        m.run(600)
        nd_ref.append(sensor_read(m)["coherence"])
    for sg in cfg["sens_noise"]:
        rd, wos, om_avg = [], [], []
        for s in range(S):
            m = MembranePropulsion(seed=s, drive_k=1, noise=sg)
            m.run(600)
            sen = ReceiverSensor(m, position=10.0)
            sen.sample()
            for _ in range(20):
                m.step()
                r = sen.sample()
            rd.append(r)
            wos.append(m.wave_order())
            om_avg.append(mean(h[2] for h in sen.history[1:]))
        coh = [r["coherence"] for r in rd]
        pool[sg] = coh
        locked = [c for c, w in zip(coh, wos) if w > 0.9]
        unl = [c for c, w in zip(coh, wos) if w <= 0.9]
        rows.append({"sigma": sg, "coherence_mean": mean(coh), "coherence_sd": sd(coh),
                     "k_err_sd": sd(r["k_est"] - 1 for r in rd), "k_err_mean": mean(r["k_est"] - 1 for r in rd),
                     "k_est_rounds_to_1": rate(sum(round(r["k_est"]) == 1 for r in rd), S),
                     "omega_single_mean": mean(r["omega"] for r in rd), "omega_single_sd": sd(r["omega"] for r in rd),
                     "omega_avg20_mean": mean(om_avg), "omega_avg20_sd": sd(om_avg),
                     "locked_n": len(locked), "unlocked_n": len(unl),
                     "auc_coherence_locked_vs_unlocked": auc(locked, unl) if locked and unl else None,
                     "auc_coherence_vs_nodrive": auc(coh, nd_ref)})
    # threshold sensitivity: classify "driven (sigma-level) vs no-drive" by coherence > t
    thr = []
    for t in cfg["sens_thresholds"]:
        row = {"threshold": t, "nodrive_false_positive": mean(c > t for c in nd_ref)}
        for sg in (0.0, 1.0, 2.0):
            if sg in pool:
                row["tpr_sigma_%s" % sg] = mean(c > t for c in pool[sg])
        thr.append(row)
    return {"seeds": S, "rows": rows, "thresholds": thr}


def e_sensor_position(cfg):
    """Position sweep incl. the seam (63->0) and the +-pi wrap; sampling-interval aliasing."""
    S = cfg["sens_pos_seeds"]
    positions = [i * 0.5 for i in range(128)]
    kerr, om, coh, moved = [], [], [], []
    wrap_straddle = interp_err = naive_err = 0
    wrap_total = 0
    per_site_phase_err = []
    seam_k = []
    for s in range(S):
        m = MembranePropulsion(seed=s, drive_k=1)
        m.run(600)
        sen = ReceiverSensor(m, position=0.0)
        sen.sample()
        m.step()
        for p in positions:
            sen.position = p
            r = sen.sample()
            kerr.append(r["k_est"] - 1.0)
            coh.append(r["coherence"])
            if p % 1 == 0.5:
                i = int(p) % m.n
                j = (i + 1) % m.n
                a, b = m.theta[i], m.theta[j]
                wrap_total += 1
                straddles = abs(a - b) > math.pi
                wrap_straddle += straddles
                # reference: midpoint along the SHORT arc
                ref = wrap(a + 0.5 * wrap(b - a))
                interp_err = max(interp_err, abs(wrap(r["phase"] - ref)))
                if straddles:
                    naive_err = max(naive_err, abs(wrap(0.5 * (a + b) - ref)))
            if 63.0 <= p < 64.0:
                seam_k.append(r["k_est"])
        # omega at fixed positions: fresh sensor, two samples one step apart
        for p in positions[::4]:
            sq = ReceiverSensor(m, position=p)
            sq.sample()
            m.step()
            om.append(sq.sample()["omega"])
        # weakness probe: moving the sensor BETWEEN two samples contaminates omega
        sm = ReceiverSensor(m, position=0.0)
        sm.sample()
        m.step()
        sm.position = 8.0
        moved.append(sm.sample()["omega"])
    # sampling interval aliasing
    alias = []
    for every in cfg["sens_intervals"]:
        vals = []
        for s in range(min(S, 10)):
            m = MembranePropulsion(seed=s, drive_k=1)
            m.run(600)
            sen = ReceiverSensor(m, position=10.0)
            sen.sample()
            for _ in range(every):
                m.step()
            vals.append(sen.sample()["omega"])
        alias.append({"every_steps": every, "interval_time": every * STEP_TIME,
                      "omega_est_mean": mean(vals), "omega_est_sd": sd(vals),
                      "omega_times_interval": -1.0 * every * STEP_TIME})
    return {"seeds": S, "positions": len(positions), "k_err": summ(kerr, cfg["boot"], "pos-k"),
            "max_abs_k_err": max(abs(x) for x in kerr),
            "coherence_min": min(coh), "coherence_max": max(coh),
            "half_site_positions": wrap_total, "straddling_the_pi_wrap": wrap_straddle,
            "interp_max_err_vs_short_arc_midpoint": interp_err,
            "naive_linear_average_max_err_when_straddling": naive_err,
            "seam_k_est_min": min(seam_k), "seam_k_est_max": max(seam_k),
            "omega_fixed_positions_n": len(om), "omega_all_mean": mean(om),
            "omega_all_min": min(om), "omega_all_max": max(om),
            "omega_moved_between_samples_mean": mean(moved),
            "omega_moved_between_samples_min": min(moved),
            "omega_moved_between_samples_max": max(moved), "aliasing": alias}


def e_sensor(cfg):
    return {"recovery": e_sensor_recovery(cfg), "readonly": e_sensor_readonly(cfg),
            "states": e_sensor_states(cfg), "noise": e_sensor_noise(cfg), "position": e_sensor_position(cfg)}


def run_all(cfg):
    ms = e_multiseed(cfg)
    res = {"config": cfg, "reproduction": e_reproduction(), "multiseed": ms,
           "anatomy": e_anatomy(cfg), "k_mu": e_k_mu(cfg), "locking": e_locking(cfg),
           "noise": e_noise(cfg), "warmup": e_warmup(cfg), "dt": e_dt(cfg), "N": e_N(cfg),
           "nulls": e_nulls(cfg, ms), "sensor": e_sensor(cfg)}
    return clean(res)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--write", action="store_true", help="full run; write results/results.json and RESULTS.md")
    ap.add_argument("--render-only", action="store_true")
    ap.add_argument("--quick", action="store_true", help="small run; print a summary")
    ap.add_argument("--out-dir", default=HERE)
    args = ap.parse_args(argv)
    import render_results
    res_dir = os.path.join(args.out_dir, "results")
    json_path = os.path.join(res_dir, "results.json")
    if args.render_only:
        with open(json_path, encoding="utf-8") as fh:
            res = json.load(fh)
        render_results.write(res, os.path.join(args.out_dir, "RESULTS.md"))
        return 0
    cfg = QUICK if args.quick and not args.write else FULL
    res = run_all(cfg)
    if args.write:
        os.makedirs(res_dir, exist_ok=True)
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump(res, fh, indent=1, sort_keys=True, ensure_ascii=False)
            fh.write("\n")
        render_results.write(res, os.path.join(args.out_dir, "RESULTS.md"))
        print("wrote", json_path)
    else:
        ms = res["multiseed"]
        print("quick run: seeds=%d fwd mean %.2f rev mean %.2f nodrive mean|d| %.2f" % (
            ms["seeds"], ms["fwd"]["stats"]["mean"], ms["rev"]["stats"]["mean"],
            ms["nodrive"]["abs_stats"]["mean"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
