import numpy as np


def auc_trapezoid(t, c):
    t, c = np.asarray(t, dtype=float), np.asarray(c, dtype=float)
    return float(np.sum(np.diff(t) * (c[1:] + c[:-1]) / 2))


def terminal_slope(t, c, n_points=3):
    # log-linear regression on the last n points
    slope, _ = np.polyfit(t[-n_points:], np.log(c[-n_points:]), 1)
    return -slope


def nca(t, c, dose, n_points=3):
    t, c = np.asarray(t, dtype=float), np.asarray(c, dtype=float)
    i = int(np.argmax(c))
    ke = terminal_slope(t, c, n_points)
    auc_last = auc_trapezoid(t, c)
    auc_inf = auc_last + c[-1] / ke
    cl = dose / auc_inf
    return {"Cmax": c[i], "Tmax": t[i], "AUC_0_last": auc_last, "AUC_0_inf": auc_inf,
            "ke": ke, "t_half": np.log(2) / ke, "CL_F": cl, "V_F": cl / ke}
