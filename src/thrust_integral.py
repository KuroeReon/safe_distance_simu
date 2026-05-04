from dataclasses import dataclass
from thrust_mass_track import thrust_mass_track_main
import numpy as np
from config import DT, MAX_T


class State:
    def __init__(self, impulse):
        self.impulse = impulse

    def __add__(self, other):
        return State(self.impulse + other.impulse)

    def __mul__(self, other):
        return State(self.impulse * other)

    def __truediv__(self, other):
        return State(self.impulse / other)

    def __rmul__(self, other):
        return self * other

    def __str__(self):
        return f"State(impulse={self.impulse})"


def get_thrust_at_time(t, dt, thrust_list):
    times = np.arange(len(thrust_list)) * dt
    return np.interp(t, times, thrust_list, left=0.0, right=0.0)


def derivative(state, thrust_list, t, dt):
    thrust = get_thrust_at_time(t, dt, thrust_list)
    return State(thrust)


def runge_kutta(state, thrust_list, t, dt):
    k1 = derivative(state, thrust_list, t, dt)
    k2 = derivative(state + k1 * dt / 2, thrust_list, t + dt / 2, dt)
    k3 = derivative(state + k2 * dt / 2, thrust_list, t + dt / 2, dt)
    k4 = derivative(state + k3 * dt, thrust_list, t + dt, dt)
    return state + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6


def run_simulation(initial_state, thrust_list, dt, max_t):
    state = initial_state
    impulses = [state.impulse]
    times = [0.0]
    t = 0.0
    while t < max_t:
        state = runge_kutta(state, thrust_list, t, dt)
        t += dt
        impulses.append(state.impulse)
        times.append(t)
    return times, impulses


def thrust_integral_main(water_volume, petbottle_mass):
    total_mass, thrusts = thrust_mass_track_main(water_volume, petbottle_mass)

    initial_state = State(impulse=0.0)
    times, impulses = run_simulation(initial_state, thrusts, DT, MAX_T)

    total_impulse = impulses[-1]
    return total_mass, total_impulse