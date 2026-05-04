import numpy as np
from velocity_track import velocity_track
from config import GRAVITY, PETBOTTLE_VOLUME

water_volumes = np.arange(0.000010, PETBOTTLE_VOLUME, 0.000010)
petbottle_masses = np.round(np.arange(0.3, 0.71, 0.01), 2)

print(f"{'Mass (kg)':>12}  {'Safe Distance (m)':>18}")
print("-" * 33)

for dry_mass in petbottle_masses:
    best_safe_distance = 0.0
    for water_volume in water_volumes:
        burnout_velocity = velocity_track(water_volume, dry_mass)
        if burnout_velocity is None:
            continue
        safe_distance = burnout_velocity ** 2 / (-GRAVITY)
        if safe_distance > best_safe_distance:
            best_safe_distance = safe_distance
    print(f"{dry_mass:>12.2f}  {best_safe_distance:>18.2f}")
