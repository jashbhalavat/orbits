import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, '../numerical_integrator')
import rk

# Universal gravitational constant [km^3/(kg*s^2)]
G = 6.67259e-20


def rates(t, y):
    y = y.flatten()
    # print(f"rates {y.shape=}")
    R1 = np.array([y[0], y[1], y[2]])
    R2 = np.array([y[3], y[4], y[5]])

    V1 = np.array([y[6], y[7], y[8]])
    V2 = np.array([y[9], y[10], y[11]])

    r = np.linalg.norm(R1-R2)

    A1 = G * m2 *(R2 - R1) / pow(r, 3)
    A2 = G * m1 *(R1 - R2) / pow(r, 3)

    dydt = np.concatenate((V1, V2, V1, A2), axis=0)
    # print(f"rates {dydt.shape=}")
    return dydt


# Input data
m1 = 1e26
m2 = 1e26
t0 = 0.0
tf = 480.0
tspan = [t0, tf]
R1_0 = np.array([0.0, 0.0, 0.0])
R2_0 = np.array([3000.0, 0.0, 0.0])
V1_0 = np.array([10.0, 20.0, 30.0])
V2_0 = np.array([0.0, 40.0, 0.0])

y0 = np.concatenate((R1_0, R2_0, V1_0, V2_0), axis=None)

h = 0.01

[t, y] = rk.rk1_4(rates, tspan, y0, h, 4)

X1 = y[0, :]
X1 = X1[np.logical_not(np.isnan(X1))]
Y1 = y[1, :]
Z1 = y[2, :]
X2 = y[3, :]
Y2 = y[4, :]
Z2 = y[5, :]

XG = np.empty([len(X1)])
YG = np.empty([len(Y1)])
ZG = np.empty([len(Z1)])

for i in range(len(X1)):
    XG[i] = (m1*X1[i] + m2*X2[i])/(m1 + m2)
    YG[i] = (m1*Y1[i] + m2*Y2[i])/(m1 + m2)
    ZG[i] = (m1*Z1[i] + m2*Z2[i])/(m1 + m2)

print(X1)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(X1, Y1[0:len(X1)], Z1[0:len(X1)], 'red')
ax.plot3D(X2, Y2[0:len(X2)], Z2[0:len(X2)], 'green')
ax.plot3D(XG, YG[0:len(XG)], ZG[0:len(XG)], 'blue')
plt.show()