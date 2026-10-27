# DESIGN.md — membrane-propulsion

## Lineage

Reuses, in spirit, the sparkainlp-x membrane frame: a Kuramoto phase
ring (microtubule-stars), travelling enveloped waves (membrane-heart /
membrane-brain), and a seeded, reproducible, stdlib-only toy with
tests. It deliberately does NOT reuse the holographic memory or the
refusal gate for the drift itself.

## Design decisions

1. **Probe law is assumed, not derived.** `v = mu * g * N / (2*pi)`,
   with `g` the wrapped neighbour phase difference. A locked k-wave
   then gives `v ~ mu*k` by construction. The experiment is therefore
   not "does a wave push?" — it is "does the ring lock to the drive,
   and does the sign follow k, while the controls stay near zero?"
   Locking (wave-order 1.000) and sign control across seeds are the
   measured claims; the probe law itself is an input, and is labelled
   as such in every deliverable.

2. **Steady-state measurement only.** 200 warmup steps are discarded.
   Counting transient synchronisation drift as propulsion would be a
   measurement error: the uniform-drive case goes from +7.6 sites
   (no warmup) to -0.0 (with warmup). All reported displacements are
   steady-state.

3. **Neighbour coupling, not mean-field.** Neighbour coupling lets the
   ring hold a travelling wave and, crucially, lets it fall into a
   twisted state. Mean-field coupling would have hidden that. The
   twisted state is the most interesting control result in the toy:
   seed 2 drifts -28.1 sites with no drive at all.

4. **Wave-order, not global order, for locking.** A k=1 travelling
   wave spreads phases uniformly around the circle, so the classic
   Kuramoto `r` is ~0 even when perfectly locked. The first test
   version used `r` and failed for exactly this reason — a test bug,
   fixed by measuring locking in the drive frame.

## Falsified / rejected designs

- **"Drive creates drift from nothing."** Falsified by the no-drive
  twisted state (seed 2, -28.1). Revised claim: drive controls sign.
- **"Global order measures locking."** Falsified as above; replaced
  by wave-order.
- **"Count all displacement from t=0."** Falsified by the uniform-drive
  transient (+7.6 -> -0.0 with warmup). Replaced by steady-state.
- **"Noise always collapses drift."** Only pinned on seed 42
  (+20.0 -> +1.4). Across seeds, noisy drift varies in sign
  (+1.4, -17.6, +10.2, +11.2, -12.2). Not claimed generally.
- **Holographic memory as the propulsion mechanism.** Rejected: it
  would conflate recall (which a drive degrades, per
  microtubule-stars) with drift, and no mechanism for it is in the
  code. Memory is out of scope for v0.1.0.

## v0.2.0 — receiver sensor

The sensor is read-only by construction: `sample()` touches no model
state, and a test pins that (theta, probe, t unchanged). That is not
a nicety — a sensor that steered the membrane would make every
reading circular (it would measure its own push).

Circular interpolation (not linear phase averaging) is used for the
interpolated phase, so the +/-pi wrap does not corrupt the reading
when the sensor sits across the wrap.

Frequency is a measured derivative from the sensor's own history,
not a parameter read back: the sensor constructor takes no drive
parameters at all. The -Omega sign is therefore a prediction the
tests check, not a value copied from the model.

## Open questions

- Does a 2D membrane (not a ring) show the same sign control, or do
  defects / vortices break it? Untested.
- Does a probe that feeds back onto the oscillators (two-way
  coupling) still drift, or does it kill the wave? Untested.
- Is there a probe law derived from a real wave equation (acoustic
  radiation force, etc.) that gives the same drift without assuming
  it? That would be a different, physical model — not this toy.

## Addendum (packaging, 0.2.0): what the multi-seed controls show

Added after running `rigor.py` (see `RESULTS.md`). The text above is left unchanged; where the data
contradict it, this addendum is the correction.

1. **Seed label (corrected in this PR).** The v0.2.0 upload said "seed 3 drifts -28.1 sites with no drive at all"
   in this file, the README, `demo.py` and a comment in `tests.py`; that was a mislabel, now fixed to seed 2.
   Measured no-drive
   displacements: seed 0 -2.2, seed 1 +0.7, seed 2 **-28.1**, seed 3 -1.0, seed 42 +1.8. The -28.1 is seed 2.
   Only the label text changed; no code or behaviour changed (the test itself does not depend on the label).
2. **"Revised claim: drive controls sign."** Over 200 seeds the claim holds for the locked drive (sign right
   in 200/200, 100.0%, Wilson 95% 98.1-100.0%) but it is true **by construction**: the probe law is assumed
   (`v = mu * local phase gradient`) and a locked wave has gradient `2*pi*k/N` everywhere, so
   displacement = `mu*k*T`. It is not evidence of any mechanism. A probe law of zero gives zero motion, and
   a constant law gives drive-independent motion (null controls in `RESULTS.md`).
3. **The no-drive control is not small.** 69.5% of no-drive rings end in a twisted state; 31.5% drift 10 sites
   or more, 14.0% drift 20 sites or more (as far as the driven ring), and the maximum is 98.1 sites. The sign
   is a coin flip (47.5% positive). The README's "+1.8" is an unrepresentative single seed.
4. **The default warmup is short.** The twisted-state relaxation time (~N^2/(pi^2 K), about 208 time units)
   is far longer than the default 10 time units; no-drive drift persists even after a warmup of 2000 steps.
5. **Noise.** The claim "strong noise collapses drift" is only a seed-42 statement. At sigma = 2 across
   100 seeds the forward sign survives in 57% (Wilson 95% 47.2-66.3%).
6. **dt.** The driven result is converged for dt <= 0.2; dt = 0.4 is wrong. Individual no-drive displacements
   are not converged at dt = 0.05.
7. **Sensor.** `k_est`, `omega` and `coherence` on a locked ring are by construction (adjacent-phase
   difference; coherence is a geometric function of |k|). Coherence does not separate a driven from an
   undriven ring (AUC 0.585) and ranks noisy driven rings below undriven ones; `omega` does separate them
   (AUC 1.0) but is meaningless if the sensor is moved between samples and aliases above a sampling
   interval of pi/Omega. These are usage limits, not bugs.
8. **Still not claimed:** any physical interpretation. The open questions above remain open.
