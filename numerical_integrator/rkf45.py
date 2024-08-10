import numpy as np


def rkf45(ode_function, tspan, y0, tolerance=1e-8):
    """Algorithm 1.3 in Orbital Mechanics for Engineering Students (3rd edition)

    This function uses the Runge-Kutta-Fehlberg 4(5) algorithm to integrate a system of first-order differential equations dy/dt = f(t,y).

    Inputs:
        ode_function - function where derivates f are computed
        tspan - the vector [t0 tf] giving the time interval for the solution
        y0 - column vector of initial values of the vector y
        tolerance - truncation error cannot exceed this tolerance

    """
    a = np.array([[0], [1/4], [3/8], [12/13], [1], [1/2]])
    b = np.array([
        [0, 0, 0, 0, 0],
        [1/4, 0, 0, 0, 0],
        [3/32, 9/32, 0, 0, 0],
        [1932/2197, -7200/2197, 7296/2197, 0, 0],
        [439/216, -8, 3680/513, -845/4014, 0],
        [-8/27, 2, -3544/2656, 1859/4104, -11/40]
    ])
    c4 = np.array([[25/216], [0], [1408/2565], [2197/4104], [-1/5], [0]])
    c5 = np.array([[16/135], [0], [6656/12825], [28561/56430], [-9/50], [2/55]])

    eps = 2.2204e-16

    t0 = tspan[0]
    tf = tspan[1]

    t = t0
    y = np.empty([len(y0), 1])
    for i in range(len(y0)):
        y[i] = y0[i]
    tout = [t]
    yout = np.array(y0)  # Doesn't work if y is 3d
    yout = yout.reshape(1, len(y0))
    h = (tf - t0)/100

    f = np.empty([len(y0), 6])

    while t < tf:
        hmin = 16 * np.spacing(t)
        ti = t
        yi = y

        for i in range(6):
            t_inner = ti + a[i] * h
            y_inner = yi

            for j in range(i):
                y_inner_temp = np.zeros_like(y_inner)
                for k in range(len(y_inner)):
                    y_inner_temp[k] = y_inner[k] + h * b[i, j] * f[k, j]
                y_inner = y_inner_temp
            
            # print(ode_function(t_inner, y_inner).reshape(2))
            # print(f)
            f[:, i] = ode_function(t_inner, y_inner).reshape(
                2,
            )

        te = h *  np.matmul(f, (c4 - c5))
        te_max = max(abs(te))

        ymax = max(abs(y))
        te_allowed = tolerance * max(ymax, 1.0)

        delta = pow((te_allowed / (te_max + eps)), 1/5)

        if te_max <= te_allowed:
            h = min(h, tf - t)
            t = t + h
            y = yi + h * np.matmul(f, c5)
            tout.append(t[0])
            yout = np.append(yout, np.transpose(y), axis=0)
        
        h = min(delta*h, 4*h)
        if h < hmin:
            print(f"\n\nWarning: Step size fell below its minimum allowable value {hmin} at time {t}.\n\n")

    return [tout, yout]
