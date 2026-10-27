# Changelog

## [0.2.0] - local packaging 2026-10-26 (not published, not released)

Base: the supplied v0.2.0 upload. `membrane_propulsion.py` is byte-identical to the
supplied v0.1.0 and v0.2.0 files (no behaviour change).

### Added (supplied v0.2.0)
- `sensor.py` and `test_sensor.py` (8 tests): a passive, read-only receiver sensor
  (interpolated phase, omega, k_est, coherence). Updated README, DESIGN and `demo.py`.

### Added (packaging and controls)
- `rigor.py`: stdlib-only, seeded multi-seed experiments with Wilson and bootstrap intervals:
  sign control (200 seeds), the no-drive control (twisted states and their drift), k and mu sweeps,
  locking threshold and Arnold-tongue grid, noise sweep, warmup, dt and N sweeps, null controls
  (shuffled labels, probe-law ablation, initial conditions), and sensor checks
  (k_est recovery, omega, read-only guarantee, twisted and no-drive states, noise and
  coherence ROC, position sweep, the +-pi wrap and the seam, omega aliasing).
- `render_results.py` and `RESULTS.md` (generated from `results/results.json`).
- `tests_rigor.py`: statistics helpers, Ring-versus-module bitwise equality, determinism,
  stored-results consistency, packaging metadata.
- `pyproject.toml`, CI workflow (Python 3.10 to 3.13), `CITATION.cff`, `.zenodo.json`
  (no DOI; none has been issued), `.gitignore`.
- README section pointing to RESULTS.md; DESIGN.md addendum (original text untouched).

### Documentation findings (no code change)
- The v0.2.0 README, DESIGN, `demo.py` and a `tests.py` comment attribute the no-drive drift of -28.1 to
  "seed 3"; it is seed 2. The label text is corrected (a separate commit); no code or behaviour changed.
  See the DESIGN addendum and RESULTS.md section 1.

## [0.1.0] - supplied upload

Initial supplied version: `membrane_propulsion.py` (Kuramoto ring, travelling drive, probe),
`demo.py`, `tests.py` (10 tests), README, DESIGN, LICENSE (AGPL-3.0-only).
