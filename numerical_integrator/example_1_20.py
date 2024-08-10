import rkf45

import numpy as np
import matplotlib.pyplot as plt

mu = 398600
minutes = 60

x0 = 6500
v0 = 7.8
y0 = [x0, v0]
t0 = 0
tf = 70 * minutes

def rates(t, f):
    x = f[0][0]
    dx = f[1][0]
    d2x = -mu/pow(x, 2)
    dfdt = np.array([[dx], [d2x]])
    return dfdt

[t, f] = rkf45.rkf45(rates, [t0, tf], y0)
t = [ti/minutes for ti in t]

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(t, f[:,0])
plt.xlabel("Time [min]")
plt.ylabel("Position [km]")

plt.subplot(2, 1, 2)
plt.plot(t, f[:,1])
plt.xlabel("Time [min]")
plt.ylabel("Velocity [km/s]")

plt.suptitle("Position and Velocity vs Time, Example 1.20")

plt.show()