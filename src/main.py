import numpy as np
from velocity_track import velocity_track
from config import GRAVITY, PETBOTTLE_VOLUME

petbottle_mass = float(input("ペットボトルの質量 [kg]: "))
max_safe_distance = 0

for water_volume in np.arange(0, PETBOTTLE_VOLUME / 2, 0.00001):
    burnout_velocity = velocity_track(water_volume, petbottle_mass)
    safe_distance = (burnout_velocity) ** 2 / (-GRAVITY)
    if safe_distance > max_safe_distance:
        max_safe_distance = safe_distance

print(f"保安距離: {max_safe_distance:.3f} m")