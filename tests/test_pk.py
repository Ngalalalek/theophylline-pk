import numpy as np
import pytest

from nca import auc_trapezoid
from pk_models import multiple_dose, one_compartment_oral, one_compartment_oral_lag

DOSE = 320.0
KA, KE, V = 1.5, 0.09, 32.0


def test_lag():
    t = np.array([0.0, 0.3, 1.0, 5.0])
    c = one_compartment_oral_lag(t, DOSE, KA, KE, V, tlag=0.5)
    assert c[0] == 0 and c[1] == 0
    assert c[3] == pytest.approx(one_compartment_oral(4.5, DOSE, KA, KE, V))


def test_steady_state_trough():
    # 叠加 40 次以后应该等于稳态谷浓度的公式
    tau, n = 12.0, 40
    trough = multiple_dose(n * tau, DOSE, tau, n, KA, KE, V)
    expected = DOSE * KA / (V * (KA - KE)) * (
        np.exp(-KE * tau) / (1 - np.exp(-KE * tau)) - np.exp(-KA * tau) / (1 - np.exp(-KA * tau)))
    assert trough == pytest.approx(expected, rel=1e-6)


def test_auc():
    # 直线 y = x，0-4 的面积是 8
    assert auc_trapezoid([0, 2, 4], [0, 2, 4]) == pytest.approx(8.0)
