import numpy as np
import matplotlib.pyplot as plt
from thrust_mass_track import thrust_mass_track_main
from config import PETBOTTLE_MASS, PETBOTTLE_VOLUME, DT

# シミュレーション条件
water_volume = 0.000500 # ボトル容量の35%を水として使用
petbottle_mass = PETBOTTLE_MASS

total_mass, thrusts = thrust_mass_track_main(water_volume, petbottle_mass)

time = np.arange(len(thrusts)) * DT

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(time, thrusts, color="tab:blue", linewidth=2)
ax.set_xlabel("Time (s)")
ax.set_ylabel("Thrust (N)")
ax.set_title(
    f"Thrust vs Time\n"
    f"(water volume: {water_volume * 1000000:.0f} mL, "
    f"bottle mass: {petbottle_mass * 1000:.0f} g)"
)
ax.grid(True, linestyle="--", alpha=0.6)
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
plt.tight_layout()
plt.show()
