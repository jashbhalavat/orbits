import rk45

import numpy as np
import matplotlib.pyplot as plt
import math

mu = 398600
minutes = 60


def rates(t, f):
    x = f[0][0]
    dx = f[1][0]
    d2x = -mu / pow(x, 2)
    dfdt = np.array([[dx], [d2x]])
    return dfdt


x0 = 6500
v0 = 7.8
y0 = [x0, v0]
t0 = 0.0
tf = 70 * minutes

[t, f] = rk45.rk45(rates, [t0, tf], y0)

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(t, f[:, 0], label="rk45")
plt.title("Position")
plt.subplot(2, 1, 1)
plt.plot(t, f[:, 1], label="rk45")
plt.title("Velocity")
plt.savefig("example_1_20.jpg")
