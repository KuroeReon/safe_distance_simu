import numpy as np
from velocity_track import velocity_track
from config import GRAVITY, PETBOTTLE_VOLUME

water_volumes_m3 = np.arange(0.000010, PETBOTTLE_VOLUME, 0.000010)
petbottle_masses = np.round(np.arange(0.3, 0.71, 0.1), 2)

print(f"{'Mass (kg)':>12}  {'Safe Distance (m)':>18}")
print("-" * 33)

for mass in petbottle_masses:
    best_sd = 0.0
    for vol in water_volumes_m3:
        bv = velocity_track(vol, mass)
        if bv is None:
            continue
        sd = bv ** 2 / (-GRAVITY)
        if sd > best_sd:
            best_sd = sd
    print(f"{mass:>12.2f}  {best_sd:>18.2f}")
