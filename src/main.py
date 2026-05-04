import numpy as np
from thrust_mass_track import thrust_mass_track_main
from config import GRAVITY, PETBOTTLE_VOLUME, DT

water_volumes = np.arange(0.000010, PETBOTTLE_VOLUME, 0.000010)
petbottle_masses = np.round(np.arange(0.3, 0.71, 0.01), 2)

best_safe_distances = np.zeros(len(petbottle_masses))
best_burnout_velocities = np.zeros(len(petbottle_masses))

# 推力はpetbottle_massに依存しないため、water_volumeごとに1回だけ計算する。
# total_mass(t) = propellant_mass(t) + petbottle_mass なので、
# petbottle_mass=0で計算したtotal_massがそのままpropellant_massesになる。
for water_volume in water_volumes:
    propellant_masses, thrusts = thrust_mass_track_main(water_volume, 0.0)

    n = len(thrusts)
    thrusts_arr = np.array(thrusts)           # shape (n,)
    prop_mass_arr = np.array(propellant_masses[:n])  # shape (n,)

    # dv/dt = thrust(t) / mass(t) は速度に依存しないため、RK4はシンプソン則に帰着する。
    # k2 = k3 = f(t + DT/2) を隣接ステップの線形補間で評価。
    # shape: (n-1,)
    thrust_a, thrust_b = thrusts_arr[:-1], thrusts_arr[1:]
    prop_a,   prop_b   = prop_mass_arr[:-1], prop_mass_arr[1:]
    thrust_mid = (thrust_a + thrust_b) / 2
    prop_mid   = (prop_a   + prop_b)   / 2

    # shape: (n-1, n_masses)
    mass_a   = prop_a[:, np.newaxis]   + petbottle_masses[np.newaxis, :]
    mass_b   = prop_b[:, np.newaxis]   + petbottle_masses[np.newaxis, :]
    mass_mid = prop_mid[:, np.newaxis] + petbottle_masses[np.newaxis, :]

    k1 = thrust_a[:, np.newaxis]   / mass_a
    k2 = thrust_mid[:, np.newaxis] / mass_mid
    k4 = thrust_b[:, np.newaxis]   / mass_b

    dv = (k1 + 4 * k2 + k4) / 6 * DT  # shape: (n-1, n_masses)
    burnout_velocities = np.sum(dv, axis=0)  # shape (n_masses,)
    safe_distances = burnout_velocities ** 2 / (-GRAVITY)

    improved = safe_distances > best_safe_distances
    best_safe_distances = np.where(improved, safe_distances, best_safe_distances)
    best_burnout_velocities = np.where(improved, burnout_velocities, best_burnout_velocities)

print(f"{'Mass (kg)':>12}  {'Safe Distance (m)':>18}  {'Burnout Velocity (m/s)':>18}")
print("-" * 56)
for dry_mass, best_safe_dist, best_vel in zip(petbottle_masses, best_safe_distances, best_burnout_velocities):
    print(f"{dry_mass:>12.2f}  {best_safe_dist:>18.2f}  {best_vel:>18.4f}")
