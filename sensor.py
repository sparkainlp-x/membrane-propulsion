"""Receiver sensor inside the membrane.

A passive, read-only sensor that sits at a fixed position inside the
Kuramoto ring and perceives the local vibration mathematically. It
never writes to the membrane: sampling does not change any phase,
the probe, or the drive. That is deliberate — like the Yang Observer
Shield in membrane-arena-3x3, the sensor observes, it does not steer.

Mathematics (all quantities are toy metrics, not physical units)
-----------------------------------------------------------------
Let the sensor sit at continuous position p in [0, N), between
oscillator i = floor(p) and i+1, with fraction f = p - i.

1. Interpolated phase (circular interpolation, so the wrap at
   +/-pi does not break the reading)::

       z = (1 - f) * exp(i * theta_i) + f * exp(i * theta_{i+1})
       phase = arg(z)

2. Vibration signal (the analog waveform the sensor "hears")::

       s(t) = cos(phase(t))

   plus its quadrature sin(phase(t)), so phase and amplitude are
   both recoverable from the pair (cos, sin).

3. Local frequency: phase velocity between two successive samples::

       omega = wrap(phase(t) - phase(t - dt)) / dt

   For a ring locked to a travelling drive with drive frequency
   Omega (drive phase = k*2*pi*i/N - Omega*t), a locked oscillator
   obeys theta_i ~ drive_phase_i + const, so omega ~ -Omega.
   That sign is a testable prediction, not a convention: the drive
   phase decreases in time at fixed i.

4. Local phase gradient (same definition as the probe)::

       g = wrap(theta_{i+1} - theta_i)

   and the estimated wave number::

       k_est = g * N / (2*pi)

   A locked k=+1 wave gives k_est ~ +1, k=-1 gives ~ -1, a uniform
   drive (k=0) gives ~ 0.

5. Local coherence over a window of +/-2 sites (5 oscillators)::

       r_local = |mean_j exp(i * theta_j)|,  j in window

   r_local = 1 means the window is perfectly in step; lower means
   the vibration is spread / incoherent there.

6. Vibration spread (a toy distortion metric, not energy)::

       spread = 1 - r_local

Sample returns all of these in one dict, plus the raw interpolated
phase. The sensor keeps its own history (t, phase, omega) so
frequency is a measured derivative, not a parameter read back from
the model: the sensor does not know drive_k or drive_omega, it
infers them from the vibration alone.

Standard library only. AGPL-3.0-only.
"""

import cmath
import math

TWO_PI = 2.0 * math.pi


def wrap(angle):
    while angle > math.pi:
        angle -= TWO_PI
    while angle <= -math.pi:
        angle += TWO_PI
    return angle


class ReceiverSensor:
    """Passive vibration sensor at a fixed position inside the ring."""

    def __init__(self, model, position=0.0, window=2):
        self.model = model
        self.position = position % model.n
        self.window = window
        self.history = []  # (t, phase, omega)
        self._prev_phase = None
        self._prev_t = None

    def _interpolated_phase(self):
        n = self.model.n
        i = int(self.position) % n
        j = (i + 1) % n
        f = self.position - int(self.position)
        z = (1.0 - f) * cmath.exp(1j * self.model.theta[i]) \
            + f * cmath.exp(1j * self.model.theta[j])
        return cmath.phase(z)

    def _gradient(self):
        n = self.model.n
        i = int(self.position) % n
        j = (i + 1) % n
        return wrap(self.model.theta[j] - self.model.theta[i])

    def _local_coherence(self):
        n = self.model.n
        i = int(self.position) % n
        total = 0j
        count = 0
        for d in range(-self.window, self.window + 1):
            total += cmath.exp(1j * self.model.theta[(i + d) % n])
            count += 1
        return abs(total) / count

    def sample(self):
        """Read the vibration now. Read-only: changes no model state."""
        phase = self._interpolated_phase()
        gradient = self._gradient()
        coherence = self._local_coherence()
        t = self.model.t
        if self._prev_phase is None:
            omega = 0.0
        else:
            dt = t - self._prev_t
            omega = wrap(phase - self._prev_phase) / dt if dt else 0.0
        reading = {
            "t": t,
            "position": self.position,
            "phase": phase,
            "signal_cos": math.cos(phase),
            "signal_sin": math.sin(phase),
            "omega": omega,
            "gradient": gradient,
            "k_est": gradient * self.model.n / TWO_PI,
            "coherence": coherence,
            "spread": 1.0 - coherence,
        }
        self.history.append((t, phase, omega))
        self._prev_phase = phase
        self._prev_t = t
        return reading

    def follow_probe(self):
        """Move the sensor to the probe's position (still read-only)."""
        self.position = self.model.probe
