import numpy as np


def rk1_4(ode_function, tspan, y0, h, rk):
    """Algorithm 1.1 in Orbital Mechanics for Engineering Students (3rd edition)

    This function uses a selected Runge-Kutta procedure to integrate a system of first-order differential equations dy/dt = f(t,y)

    Inputs:
        ode_function - function where derivates f are computed
        tspan - the vector [t0 tf] giving the time interval for the solution
        y0 - column vector of initial values of the vector y
        h - time step
        rk - =1 for RK1, =2 for RK2, =3 for RK3, =4 for RK4

    """
    if rk == 1:
        n_stages = 1
        a = np.array([[0]])
        b = 0
        c = np.array([[1]])
    elif rk == 2:
        n_stages = 2
        a = np.array([[0], [1]])
        b = np.array([[0], [1]])
        c = np.array([[1 / 2], [1 / 2]])
    elif rk == 3:
        n_stages = 3
        a = np.array([[0], [1 / 2], [1]])
        b = np.array([[0, 0], [1 / 2, 0], [-1, 2]])
        c = np.array([[1 / 6], [2 / 3], [1 / 6]])
    elif rk == 4:
        n_stages = 4
        a = np.array([[0], [1 / 2], [1 / 2], [1]])
        b = np.array([[0, 0, 0], [1 / 2, 0, 0], [0, 1 / 2, 0], [0, 0, 1]])
        c = np.array([[1 / 6], [1 / 3], [1 / 3], [1 / 6]])
    else:
        print(f"rk must be 1, 2, 3, or 4, but it was {rk}.")

    t0 = tspan[0]
    tf = tspan[1]

    t = t0
    y = np.empty([len(y0), 1])
    for i in range(len(y0)):
        y[i] = y0[i]
    tout = [t]
    yout = [y0[0]]  # Doesn't work if y is 3d

    f = np.empty([len(y0), rk])

    while t < tf:
        ti = t
        yi = y

        for i in range(n_stages):
            t_inner = ti + a[i] * h
            y_inner = yi

            for j in range(i - 1):
                y_inner = y_inner + h * b[i, j] * f[:, j]

            f[:, i] = ode_function(t_inner, y_inner).reshape(
                2,
            )

        h = min(h, tf - t)
        t = t + h
        y = yi + h * np.matmul(f, c)
        tout.append(t)
        yout.append(y[0][0])

    return [tout, yout]
