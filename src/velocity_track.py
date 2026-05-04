from dataclasses import dataclass
from thrust_mass_track import thrust_mass_track_main
import matplotlib.pyplot as plt
import numpy as np
from config import GRAVITY, CROSS_SECTIONAL_AREA, DRAG_COEFFICIENT, AIR_DENSITY, DT, MAX_T

class State:
    def __init__(self, velocity):
        self.velocity = velocity

    def __add__(self, other):
        return State(self.velocity + other.velocity)

    def __mul__(self, other):
        return State(self.velocity * other)

    def __truediv__(self, other):
        return State(self.velocity / other)

    def __rmul__(self, other):
        return self * other

    def __str__(self):
        return (f"velocity={self.velocity})")

@dataclass
class Config:
    gravity: float
    cross_sectional_area: float
    drag_coefficient: float
    air_density: float

def get_value_at_time(t, dt, data):
    times = np.arange(len(data)) * dt
    return np.interp(t, times, data)

def derivative(state, config, mass_list, thrust_list, t, dt):
    mass = get_value_at_time(t, dt, mass_list)
    thrust = get_value_at_time(t, dt, thrust_list)
    force = (thrust + config.gravity * mass - config.air_density * config.drag_coefficient * config.cross_sectional_area * state.velocity
             * abs(state.velocity) * 0.5)
    d_velocity = force / mass
    return State(d_velocity)

def runge_kutta(derivative, state, config, t, dt, mass, thrust):
    k1 = derivative(state, config, mass, thrust, t, dt)
    k2 = derivative(state + k1 * dt / 2, config, mass, thrust, t + dt / 2, dt)
    k3 = derivative(state + k2 * dt / 2, config, mass, thrust, t + dt / 2, dt)
    k4 = derivative(state + k3 * dt, config, mass, thrust, t + dt, dt)
    return state + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6

def run_simulation(config, initial_state, t, dt, max_t, states, thrust, mass):
    state = initial_state
    states.append(state)
    while t < max_t:
        state = runge_kutta(derivative, states[-1], config, t, dt, mass, thrust)
        if state.velocity < states[-1].velocity:
            return state.velocity
        states.append(state)
        t += dt
    return states[-1].velocity

def velocity_track(water_volume, petbottle_mass):
    total_mass, thrusts = thrust_mass_track_main(water_volume, petbottle_mass)

    config = Config(
        gravity=GRAVITY,
        cross_sectional_area=CROSS_SECTIONAL_AREA,
        drag_coefficient=DRAG_COEFFICIENT,
        air_density=AIR_DENSITY,
    )
    initial_state = State(velocity=0)
    states = []

    burnout_velocity = run_simulation(config, initial_state, 0, DT, MAX_T, states, thrusts, total_mass)

    return burnout_velocity