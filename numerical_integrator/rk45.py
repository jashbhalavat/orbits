import numpy as np


def rk45(ode_function, tspan, y0, tolerance=1e-4):
    """Algorithm 1.3 in Orbital Mechanics for Engineering Students (3rd edition)

    This function uses Runge-Kutta with variable step size

    Inputs:
        ode_function - function where derivates f are computed
        tspan - the vector [t0 tf] giving the time interval for the solution
        y0 - column vector of initial values of the vector y
        tolerance - desired tolerance
        rk - =1 for RK1, =2 for RK2, =3 for RK3, =4 for RK4

    """
    a = np.array([0, 1 / 4, 3 / 8, 12 / 13, 1.0, 1 / 2])
    b = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [1 / 4, 0.0, 0.0, 0.0, 0.0],
            [3 / 32, 9 / 32, 0.0, 0.0, 0.0],
            [1932 / 2197, -7200 / 2197, 7296 / 2197, 0.0, 0.0],
            [439 / 216, -8, 3680 / 513, -845 / 4104, 0.0],
            [-8 / 27, 2, -3544 / 2565, 1859 / 4104, -11 / 40],
        ]
    )
    c4 = np.array([25 / 216, 0.0, 1408 / 2565, 2197 / 4104, -1 / 5, 0])
    c5 = np.array([16 / 315, 0.0, 6656 / 12825, 28561 / 56430, -9 / 50, 2 / 55])

    t0 = tspan[0]
    tf = tspan[1]

    t = t0
    y = np.empty([len(y0), 1])
    for i in range(len(y0)):
        y[i] = y0[i]
    tout = [t]
    yout = [y0]

    h = (tf - t0) / 100.0

    f = np.empty([len(y0), 6])

    while t < tf:
        hmin = 160 * t
        ti = t
        yi = y
        # Evaluate the time derivative at six points within the current interval
        for i in range(6):
            t_inner = ti + a[i] * h
            y_inner = yi
            for j in range(i - 1):
                temp = h * b[i, j] * f[:, j]
                y_inner[0] += temp[0]
                y_inner[1] += temp[1]
            temp = ode_function(t_inner, y_inner)
            f[0, i] = temp[0]
            f[1, i] = temp[1]

        te = h * f * [c4 - c5]
        te_max = np.max(abs(te))

        ymax = np.max(abs(y))
        te_allowed = tolerance * max(ymax, 1.0)

        delta = pow((te_allowed / (te_max)), 1 / 5)

        h = min(h, tf - t)
        t = t + h
        mat_product = np.matmul(f, c5)
        y[0] = yi[0] + h * mat_product[0]
        y[1] = yi[1] + h * mat_product[1]
        tout.append(t)
        temp_y = [y[0][0], y[1][0]]
        yout.append(temp_y)

        h = min(delta * h, 4 * h)
        if h < hmin:
            print(
                f"Warning: Step size fell below its minimum allowable value {hmin} at time {t}"
            )

    return [tout, yout]