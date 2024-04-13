import rk

import numpy as np
import matplotlib.pyplot as plt
import math

# Inputs
m = 1.0
z = 0.03
wn = 1.0
F0 = 1.0
w = 0.4 * wn


def rates(t, f):
    x = f[0][0]
    dx = f[1][0]
    d2x = F0 / m * np.sin(w * t[0]) - 2 * z * wn * dx - pow(wn, 2) * x

    dfdt = np.array([[dx], [d2x]])
    return dfdt


x0 = 0.0
x_dot0 = 0.0
f0 = [x0, x_dot0]

t0 = 0.0
tf = 110.0
tspan = [t0, tf]

# RK1
h = 0.01
[t1, f1] = rk.rk1_4(rates, tspan, f0, h, 4)


# Exact soln
wd = wn * np.sqrt(1 - pow(z, 2))
den = pow((pow(wn, 2) - pow(w, 2)), 2) + pow(2 * w * wn * z, 2)
C1 = (pow(wn, 2) - pow(w, 2)) / den * F0 / m
C2 = -2 * w * wn * z / den * F0 / m
A = (
    x0 * wn / wd
    + x_dot0 / wd
    + (pow(w, 2) + (2 * pow(z, 2) - 1) * pow(wn, 2)) / den * w / wd * F0 / m
)
B = x0 + 2 * w * wn * z / den * F0 / m

t = np.linspace(t0, tf, 5000)
x = np.zeros_like(t)
for i in range(len(t)):
    x[i] = (
        (A * np.sin(wd * t[i]) + B * np.cos(wd * t[i])) * math.exp(-wn * z * t[i])
        + C1 * np.sin(w * t[i])
        + C2 * np.cos(w * t[i])
    )


plt.figure()
plt.plot(t1, f1, label="rk")
plt.plot(t, x, label="true")
plt.legend()
plt.savefig("example_1_18.jpg")
