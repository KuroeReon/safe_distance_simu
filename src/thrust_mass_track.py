from dataclasses import dataclass
import matplotlib.pyplot as plt
from config import (
    PETBOTTLE_MASS, PETBOTTLE_VOLUME, WATER_VOLUME, INITIAL_BOTTLE_PRESSURE,
    NOZZLE_AREA, WATER_DENSITY, ATMOS_PRESSURE, SPECIFIC_HEAT_RATIO,
    GAS_CONSTANT, ABSOLUTE_TEMPERATURE, AIR_MOLAR_MASS, DT, MAX_T,
)

class State:
    def __init__(self, water_mass, bottle_pressure, air_volume, total_mass):
        self.water_mass = water_mass
        self.bottle_pressure = bottle_pressure
        self.air_volume = air_volume
        self.total_mass = total_mass

    def __add__(self, other):
        return State(self.water_mass + other.water_mass, self.bottle_pressure + other.bottle_pressure,
                     self.air_volume + other.air_volume, self.total_mass + other.total_mass)

    def __mul__(self, other):
        return State(self.water_mass * other, self.bottle_pressure * other, self.air_volume * other, self.total_mass * other)

    def __truediv__(self, other):
        return State(self.water_mass / other, self.bottle_pressure / other, self.air_volume / other, self.total_mass / other)

    def __rmul__(self, other):
        return self * other

    def __str__(self):
        return (f"State(water_mass={self.water_mass}, bottle_pressure={self.bottle_pressure}, "
                f"air_volume={self.air_volume}, total_mass={self.total_mass})")

def derivative_water(state, config):
    exit_velocity = (2 * (state.bottle_pressure - config.atmos_pressure) / config.water_density) ** 0.5
    d_water_mass = -config.water_density * config.nozzle_area * exit_velocity
    d_bottle_pressure = -config.specific_heat_ratio * state.bottle_pressure * config.nozzle_area * exit_velocity / state.air_volume
    d_air_volume = config.nozzle_area * exit_velocity
    d_total_mass = -config.water_density * config.nozzle_area * exit_velocity
    return State(d_water_mass, d_bottle_pressure, d_air_volume, d_total_mass)

def basic_math_air(state, config):
    judge_pressure = ((2 / (config.specific_heat_ratio + 1)) ** (config.specific_heat_ratio / (config.specific_heat_ratio - 1))) * state.bottle_pressure
    if judge_pressure >= config.atmos_pressure:
        exit_pressure = judge_pressure
    else:
        exit_pressure = config.atmos_pressure
    if exit_pressure > state.bottle_pressure:
        exit_pressure = state.bottle_pressure
    exit_velocity = ((2 * config.specific_heat_ratio / (config.specific_heat_ratio - 1) * state.bottle_pressure /
                     ((state.bottle_pressure / config.const) ** (1 / config.specific_heat_ratio)) *
                     (1 - ((exit_pressure / state.bottle_pressure) ** (1 - 1 / config.specific_heat_ratio)))) ** 0.5)
    d_air_mass = ((exit_pressure / config.const) ** (1 / config.specific_heat_ratio)) * config.nozzle_area * exit_velocity
    return exit_pressure, exit_velocity, d_air_mass

def derivative_air(state, config):
    _, _, d_air_mass = basic_math_air(state, config)
    d_bottle_pressure = (-config.specific_heat_ratio * state.bottle_pressure / state.air_volume * d_air_mass *
                         ((state.bottle_pressure / config.const) ** (-1 / config.specific_heat_ratio)))
    d_water_mass = 0
    d_air_volume = 0
    d_total_mass = -d_air_mass
    return State(d_water_mass, d_bottle_pressure, d_air_volume, d_total_mass)

@dataclass
class Config:
    water_density: float
    nozzle_area: float
    atmos_pressure: float
    specific_heat_ratio: float
    const: float
    gas_constant: float
    abusolute_temperture: float
    air_molar_mass: float

def runge_kutta(derivative, state, config, dt):
    k1 = derivative(state, config)
    k2 = derivative(state + k1 * dt / 2, config)
    k3 = derivative(state + k2 * dt / 2, config)
    k4 = derivative(state + k3 * dt, config)
    return state + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6

def run_simulation_water(config, initial_state, t, dt, max_t, states, thrusts, exit_velocities, mass):
    state = initial_state
    states.append(state)
    while t < max_t and state.water_mass > 0:
        t += dt
        thrust = 2 * config.nozzle_area * (state.bottle_pressure - config.atmos_pressure)
        exit_velocity = (2 * (state.bottle_pressure - config.atmos_pressure) / config.water_density) ** 0.5
        state = runge_kutta(derivative_water, state, config, dt)
        mass.append(state.total_mass)
        thrusts.append(thrust)
        exit_velocities.append(exit_velocity)
        states.append(state)
    return t

def run_simulation_air(config, initial_state, t, dt, max_t, states, thrusts, exit_velocities, mass):
    state = initial_state
    while t < max_t and state.bottle_pressure > config.atmos_pressure:
        t += dt
        state = runge_kutta(derivative_air, state, config, dt)
        exit_pressure, exit_velocity, d_air_mass = basic_math_air(state, config)
        thrust = d_air_mass * exit_velocity + config.nozzle_area * (exit_pressure - config.atmos_pressure)
        thrusts.append(thrust)
        states.append(state)
        exit_velocities.append(exit_velocity)
        mass.append(state.total_mass)
    return states, thrusts

def thrust_mass_track_main(water_volume, petbottle_mass=None):
    if petbottle_mass is None:
        petbottle_mass = PETBOTTLE_MASS

    initial_air_volume = PETBOTTLE_VOLUME - water_volume
    initial_water_mass = water_volume * WATER_DENSITY
    initial_air_mass = (AIR_MOLAR_MASS * (INITIAL_BOTTLE_PRESSURE - ATMOS_PRESSURE) * initial_air_volume
                        / GAS_CONSTANT / ABSOLUTE_TEMPERATURE)

    config = Config(
        water_density=WATER_DENSITY,
        nozzle_area=NOZZLE_AREA,
        atmos_pressure=ATMOS_PRESSURE,
        specific_heat_ratio=SPECIFIC_HEAT_RATIO,
        const=0,
        gas_constant=GAS_CONSTANT,
        abusolute_temperture=ABSOLUTE_TEMPERATURE,
        air_molar_mass=AIR_MOLAR_MASS,
    )

    thrusts = []
    exit_velocities = []
    states = []
    total_mass = [initial_air_mass + initial_water_mass + petbottle_mass]

    initial_state_1 = State(
        water_mass=initial_water_mass,
        bottle_pressure=INITIAL_BOTTLE_PRESSURE,
        air_volume=initial_air_volume,
        total_mass=initial_air_mass + initial_water_mass + petbottle_mass,
    )
    initial_t_2 = run_simulation_water(config, initial_state_1, 0, DT, MAX_T, states, thrusts, exit_velocities, total_mass)

    initial_state_2 = states[-1]
    config.const = (GAS_CONSTANT * ABSOLUTE_TEMPERATURE) ** SPECIFIC_HEAT_RATIO
    run_simulation_air(config, initial_state_2, initial_t_2, DT, MAX_T, states, thrusts, exit_velocities, total_mass)

    return total_mass, thrusts