"""Demo: steady-state probe drift vs drive mode (seed 42)."""

from membrane_propulsion import run_case

CASES = [
    ("forward drive  k=+1", dict(drive_k=1)),
    ("reverse drive  k=-1", dict(drive_k=-1)),
    ("uniform drive  k=0", dict(drive_k=0)),
    ("no drive", dict(drive_amp=0.0)),
    ("forward + strong noise", dict(drive_k=1, noise=2.0)),
]

print("membrane-propulsion demo (N=64, 200 warmup + 400 steps, seed 42)")
print("displacement in sites; + = forward, - = backward")
for label, kw in CASES:
    disp, wave_order, _ = run_case(seed=42, steps=400, **kw)
    print(f"  {label:24s} displacement {disp:+7.1f}  wave-order {wave_order:.3f}")
print()
print("Read: drive mode controls the sign, uniformly across seeds;")
print("uniform drive gives ~0; no drive can still drift on its own")
print("if the ring falls into a twisted state (seed 3: -28.1).")
print("This is a synthetic drift toy, not a propulsion device.")
print()
print("Receiver sensor (inside the membrane, position 10, seed 42):")
from membrane_propulsion import MembranePropulsion
from sensor import ReceiverSensor
for label, kw in CASES[:3]:
    m = MembranePropulsion(seed=42, **kw)
    m.run(600)
    s = ReceiverSensor(m, position=10.0)
    s.sample()
    m.step()
    r = s.sample()
    print(f"  {label:24s} phase {r['phase']:+.3f}  omega {r['omega']:+.3f}  "
          f"k_est {r['k_est']:+.2f}  coherence {r['coherence']:.3f}  "
          f"signal {r['signal_cos']:+.3f}")
print("Sensor reads the vibration only; it does not move the membrane.")
