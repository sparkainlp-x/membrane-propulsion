# membrane-propulsion v0.2.0

**SYNTHETIC toy: a travelling wave on a phase-oscillator membrane drifts
a probe particle, and a passive receiver sensor inside the membrane
perceives the vibration mathematically. Classical simulation only.
Not a propulsion device, not a force, not hardware, not a measurement
of anything real.**

## What it is

A ring of 64 Kuramoto phase oscillators (the membrane) is driven by a
travelling external wave, mode number `k`. A probe particle rides the
membrane: its velocity is proportional to the local phase gradient
between the two oscillators it sits between.

That probe law is a modelling assumption, stated in the code and in
DESIGN.md — not a physical force law. "Propulsion" here means net drift
of the probe in the toy, in sites per time unit.

## Measured (seed 42, 200 warmup + 400 steps, N=64)

| Case | Displacement (sites) | Wave-order |
|---|---|---|
| Forward drive, k=+1 | +20.0 | 1.000 |
| Reverse drive, k=-1 | -20.0 | 1.000 |
| Uniform drive, k=0 | -0.0 | 1.000 |
| No drive | +1.8 | 0.418 |
| Forward + strong noise (sigma=2.0) | +1.4 | 0.384 |

Direction is controlled by the drive sign, consistently across seeds
0, 1, 2, 3, 42 (all +k > +10, all -k < -10). Uniform drive gives ~0.

**Honest control that does NOT go away:** without a drive, the ring
can fall into a twisted state and drift on its own — seed 2 gives
-28.1 sites steady-state with no drive at all. So the drive controls
the *sign*, it does not create drift from nothing. That twisted-state
result is kept in the tests, not hidden.

## Receiver sensor (new in v0.2.0)

`sensor.py` — a passive sensor at a fixed position inside the ring.
It reads, it does not write: sampling changes no phase, no probe, no
drive (test-pinned, like the Yang Observer Shield).

Mathematics, for a sensor at position p between oscillators i and i+1,
fraction f = p - i:

- Interpolated phase: `phase = arg((1-f)*exp(i*theta_i) + f*exp(i*theta_{i+1}))`
- Vibration signal: `s(t) = cos(phase)`, quadrature `sin(phase)`
- Local frequency: `omega = wrap(phase(t) - phase(t-dt)) / dt`
- Gradient: `g = wrap(theta_{i+1} - theta_i)`, wave number `k_est = g*N/(2*pi)`
- Local coherence (window +/-2 sites): `r = |mean exp(i*theta_j)|`, spread `1 - r`

Measured (seed 42, locked, position 10):

| Case | phase | omega | k_est | coherence |
|---|---|---|---|---|
| Forward, k=+1 | +3.077 | -1.000 | +1.00 | 0.990 |
| Reverse, k=-1 | +1.114 | -1.000 | -1.00 | 0.990 |
| Uniform, k=0 | +2.096 | -1.000 | +0.00 | 1.000 |

Note the sign prediction that the tests pin: a locked oscillator
follows drive phase `k*2*pi*i/N - Omega*t`, so at fixed i its phase
velocity is **-Omega** (-1.000 here), in all three cases. The sensor
does not know k or Omega in advance — it infers them from the
vibration alone (test-pinned: no drive parameters in its constructor).
Strong noise lowers coherence below the locked value (test-pinned).

## Run

```
python demo.py          # drift table + sensor readings above
python tests.py         # 10/10 tests (drift)
python test_sensor.py   # 8/8 tests (sensor)
```

## Files

- `membrane_propulsion.py` — Kuramoto ring + probe (stdlib only)
- `sensor.py` — passive receiver sensor (stdlib only)
- `demo.py`, `tests.py`, `test_sensor.py`, `DESIGN.md`, `README.md`, `LICENSE` (AGPL-3.0-only)

## Limits

- No newtons, no thrust, no fluid, no air, no space drive.
- No holographic memory is used for the drift: the memory and the
  propulsion layers are deliberately kept separate (microtubule-stars
  showed a drive degrades recall, so claiming both at once would be
  dishonest).
- The global Kuramoto order parameter is ~0 for a locked travelling
  wave by design (phases spread around the circle); locking is measured
  with `wave_order()`, in the drive frame, instead.
- Noise results are seed-dependent: strong noise collapses drift on
  seed 42, but across seeds the noisy drift varies in sign — the test
  pins only seed 42, and the README says so.

## Multi-seed controls (added in the packaged 0.2.0)

`rigor.py` runs 200-seed controls (standard library only); the results are in
[`RESULTS.md`](RESULTS.md), generated from `results/results.json`. Run
`python3 rigor.py --write` to reproduce (about 10 minutes) or `python3 rigor.py --quick` for a smoke run.
Short version, all measured, none of it physics:

- With the drive locked, the sign and size of the displacement (`mu*k*T`) are **by construction**.
- The no-drive ring drifts on its own in many seeds (about 70% end in a twisted state; 31.5% drift 10 sites
  or more in either direction), so single-seed no-drive numbers are not representative.
- The sign control under noise sigma = 2 is close to a coin flip across seeds.
- One documentation correction: the -28.1 no-drive drift is **seed 2**, not seed 3 (see the DESIGN addendum).
- The sensor's `k_est` and `omega` equal the drive's values on a locked ring by construction, and its
  coherence cannot tell a driven ring from an undriven one.

This is a synthetic toy; it is not propulsion, thrust, force, or a space drive.
`tests_rigor.py` tests the controls and the packaging metadata; `CHANGELOG.md` lists the changes.
