import math

# ---- 内圧、水の量など ----
PETBOTTLE_MASS = 0.470          # kg → 入力がもしなかったら
PETBOTTLE_VOLUME = 0.002000     # m^3 
INITIAL_BOTTLE_PRESSURE = 500000  # Pa
NOZZLE_DIAMETER = 0.0089        # m

# ---- 機体 ----
CROSS_SECTIONAL_DIAMETER = 0.09  # m

# ---- 物理定数 ----
GRAVITY = -9.807                # m/s^2
WATER_DENSITY = 998             # kg/m^3
ATMOS_PRESSURE = 101325.0       # Pa
SPECIFIC_HEAT_RATIO = 1.4
GAS_CONSTANT = 8.3144           # J/(mol·K)
ABSOLUTE_TEMPERATURE = 303      # K
AIR_MOLAR_MASS = 0.02897        # kg/mol
AIR_DENSITY = 1.25              # kg/m^3
DRAG_COEFFICIENT = 0.5

# ---- Derived ----
NOZZLE_AREA = math.pi * (NOZZLE_DIAMETER / 2) ** 2
CROSS_SECTIONAL_AREA = math.pi * (CROSS_SECTIONAL_DIAMETER / 2) ** 2

# ---- Simulation ----
DT = 0.001   # s
MAX_T = 5.0  # s
