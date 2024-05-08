import numpy as np

import rv_from_r0v0_ta

mu = 398600

R0 = np.array([8182.4, -6865.9, 0.0])
V0 = np.array([0.47572, 8.8116, 0.0])
dt = 120.0

R, V = rv_from_r0v0_ta.rv_from_r0v0_ta(R0, V0, dt, mu)

print(f"{R=}")
print(f"{V=}")
