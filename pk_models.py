"""1-cmt oral absorption models. Units: h, mg, mg/L, L (V is V/F, F unknown)."""

import numpy as np


def one_compartment_oral(t, dose, ka, ke, V):
    t = np.asarray(t, dtype=float)
    c = dose * ka / (V * (ka - ke)) * (np.exp(-ke * t) - np.exp(-ka * t))
    return np.where(t < 0, 0.0, c)  # nothing before the dose


def one_compartment_oral_lag(t, dose, ka, ke, V, tlag):
    # same curve shifted by tlag, so C = 0 until absorption starts
    return one_compartment_oral(np.asarray(t, dtype=float) - tlag, dose, ka, ke, V)


def multiple_dose(t, dose, tau, n_doses, ka, ke, V):
    # superposition, ok because the model is linear
    t = np.asarray(t, dtype=float)
    c = np.zeros_like(t)
    for i in range(n_doses):
        c = c + one_compartment_oral(t - i * tau, dose, ka, ke, V)
    return c
