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
    f = f.flatten()
    x = f[0]
    dx = np.array([f[1]])
    d2x = F0 / m * np.sin(w * t) - 2 * z * wn * dx - pow(wn, 2) * x

    # dfdt = np.array([[dx], [d2x]])
    dfdt = np.concatenate((dx, d2x), axis=0)
    return dfdt


x0 = np.array([0.0])
x_dot0 = np.array([0.0])
f0 = np.concatenate((x0, x_dot0), axis=None)

t0 = 0.0
tf = 110.0
tspan = [t0, tf]

# RK1
h = 0.01
[t1, f1] = rk.rk1_4(rates, tspan, f0, h, 3)


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
plt.plot(t1/np.max(t1), f1[0,:]/np.max(f1[0,:]), label="rk")
plt.plot(t/np.max(t), x/np.max(x), label="true")
plt.legend()
plt.savefig("example_1_18.jpg")
