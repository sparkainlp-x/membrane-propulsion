"""membrane-propulsion: a synthetic travelling-wave drift toy.

A ring of Kuramoto phase oscillators (the membrane) is driven by a
travelling external wave. A probe particle rides the membrane: its
velocity is proportional to the local phase gradient between the two
oscillators it sits between.

This is a modelling assumption, stated plainly, not a physical force
law. "Propulsion" here means net drift of the probe in the toy, in
sites per time unit. There is no force in newtons, no thrust, no
hardware, no fluid, and no claim about the real world.

Model
-----
Oscillator i (i = 0..N-1, ring neighbours) obeys, Euler step dt::

    dtheta_i = omega0
               + K * (sin(theta_{i-1} - theta_i)
                      + sin(theta_{i+1} - theta_i)) / 2
               + A * sin(drive_phase_i - theta_i)
               + sigma * sqrt(dt) * noise

with drive_phase_i = k * 2*pi*i / N - Omega * t, where k is an integer
mode number (the travelling-wave direction and winding), A the drive
amplitude, Omega the drive frequency.

Probe, continuous position p in [0, N), unwrapped position pu::

    g = wrap(theta_{idx+1} - theta_{idx})   # local phase gradient
    v = mu * g * N / (2*pi)                 # so a locked k-wave gives v ~ mu*k
    p = (p + v * dt) % N
    pu += v * dt

If the oscillators lock to the drive, g ~ 2*pi*k / N and v ~ mu*k:
k = +1 drifts forward, k = -1 drifts backward, k = 0 (uniform drive)
has g ~ 0 and no drift. With no drive, coupling synchronises the ring
to a common phase, so g -> 0 and net drift stays small. Strong noise
breaks locking and collapses the drift. Those three controls are what
the tests pin; without them the drift would be a tautology.

Standard library only. AGPL-3.0-only.
"""

import math
import random

TWO_PI = 2.0 * math.pi


def wrap(angle):
    """Wrap an angle to (-pi, pi]."""
    while angle > math.pi:
        angle -= TWO_PI
    while angle <= -math.pi:
        angle += TWO_PI
    return angle


class MembranePropulsion:
    """Kuramoto ring + phase-gradient probe."""

    def __init__(self, n=64, coupling=2.0, omega0=0.0, drive_amp=1.5,
                 drive_k=1, drive_omega=1.0, mu=1.0, noise=0.0,
                 dt=0.05, seed=42, probe0=0.0):
        self.n = n
        self.coupling = coupling
        self.omega0 = omega0
        self.drive_amp = drive_amp
        self.drive_k = drive_k
        self.drive_omega = drive_omega
        self.mu = mu
        self.noise = noise
        self.dt = dt
        self.t = 0.0
        self.rng = random.Random(seed)
        # Random initial phases, seeded (deterministic per seed).
        self.theta = [self.rng.uniform(-math.pi, math.pi) for _ in range(n)]
        self.probe = probe0 % n
        self.probe_unwrapped = probe0
        self.history = []  # (t, probe_unwrapped, mean_phase_gradient)

    def _drive_phase(self, i):
        return self.drive_k * TWO_PI * i / self.n - self.drive_omega * self.t

    def mean_gradient(self):
        """Mean wrapped neighbour phase difference, in rad per site."""
        total = 0.0
        for i in range(self.n):
            total += wrap(self.theta[(i + 1) % self.n] - self.theta[i])
        return total / self.n

    def order_parameter(self):
        """Kuramoto order r = |mean exp(i theta)|, 0 (incoherent) to 1.

        Note: a locked travelling wave with winding k != 0 spreads
        phases around the whole circle, so this r is ~0 by design even
        when the wave is perfectly locked. Use wave_order() for that.
        """
        re = sum(math.cos(th) for th in self.theta) / self.n
        im = sum(math.sin(th) for th in self.theta) / self.n
        return math.hypot(re, im)

    def wave_order(self):
        """Locking to the drive frame: |mean exp(i (theta_i - drive_i))|.

        ~1 when oscillators ride the travelling drive, ~0 when not.
        """
        re = 0.0
        im = 0.0
        for i in range(self.n):
            d = self.theta[i] - self._drive_phase(i)
            re += math.cos(d)
            im += math.sin(d)
        return math.hypot(re / self.n, im / self.n)

    def step(self):
        n = self.n
        dt = self.dt
        new_theta = []
        for i in range(n):
            coupling_term = (
                math.sin(self.theta[(i - 1) % n] - self.theta[i])
                + math.sin(self.theta[(i + 1) % n] - self.theta[i])
            ) / 2.0
            drive_term = self.drive_amp * math.sin(
                self._drive_phase(i) - self.theta[i])
            noise_term = 0.0
            if self.noise:
                noise_term = self.noise * math.sqrt(dt) * self.rng.gauss(0.0, 1.0)
            new_theta.append(
                self.theta[i]
                + dt * (self.omega0 + self.coupling * coupling_term + drive_term)
                + noise_term)
        self.theta = [wrap(th) for th in new_theta]
        self.t += dt
        # Probe rides the local gradient.
        idx = int(self.probe) % n
        nxt = (idx + 1) % n
        g = wrap(self.theta[nxt] - self.theta[idx])
        v = self.mu * g * n / TWO_PI
        self.probe = (self.probe + v * dt) % n
        self.probe_unwrapped += v * dt
        self.history.append((self.t, self.probe_unwrapped, self.mean_gradient()))
        return v

    def run(self, steps):
        """Run `steps` Euler steps; return unwrapped probe displacement."""
        start = self.probe_unwrapped
        for _ in range(steps):
            self.step()
        return self.probe_unwrapped - start


def run_case(seed=42, steps=400, warmup=200, **kwargs):
    """Build, warm up, then measure steady-state drift.

    The warmup matters: displacement during the initial transient
    (ring synchronising from random phases) is not propulsion, so it
    is discarded and only steady-state drift is reported. Returns
    (displacement, final wave-order, model).
    """
    model = MembranePropulsion(seed=seed, **kwargs)
    if warmup:
        model.run(warmup)
    start = model.probe_unwrapped
    model.run(steps)
    return model.probe_unwrapped - start, model.wave_order(), model
