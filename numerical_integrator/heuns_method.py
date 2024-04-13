import numpy as np


def heuns_method(ode_function, tspan, y0, h):
    """Algorithm 1.2 in Orbital Mechanics for Engineering Students (3rd edition)

    This function uses Heun's Metohd of first-order differential equations dy/dt = f(t,y)

    Inputs:
        ode_function - function where derivates f are computed
        tspan - the vector [t0 tf] giving the time interval for the solution
        y0 - column vector of initial values of the vector y
        h - time step
        rk - =1 for RK1, =2 for RK2, =3 for RK3, =4 for RK4

    """
    tol = 1e-6
    eps = 1e-12
    itermax = 100

    t0 = tspan[0]
    tf = tspan[1]

    t = t0
    y = np.empty([len(y0), 1])
    for i in range(len(y0)):
        y[i] = y0[i]
    tout = [t]
    yout = [y0[0]]  # Doesn't work if y is 3d

    while t < tf:
        h = min(h, tf - 1)
        t1 = [t]
        y1 = y
        f1 = ode_function(t1, y1)
        y2 = y1 + f1 * h
        t2 = [t1[0] + h]
        err = tol + 1
        iteration = 0
        while err > tol and iteration <= itermax:
            y2p = y2
            f2 = ode_function(t2, y2p)
            favg = (f1 + f2) / 2
            y2 = y1 + favg * h
            for i in range(len(y2)):
                temp = abs(y2[i] - y2p[i]) / (y2[i] + eps)
            err = max(temp)
            iteration += 1

        if iteration > itermax:
            print(f"Maximum number of iterations {itermax}")
            print(f" exceeded at at time = {t}")
            print(f" in function heun.")

        t = t + h
        y = y2
        tout.append(t)
        yout.append(y[0][0])

    return [tout, yout]
