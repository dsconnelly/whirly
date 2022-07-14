from functools import partial

import numpy as np

from whirly.fourier import FourierField
from whirly.solvers import NavierStokesSolver

def test_taylor_green():
    m, p = 101, 1
    tau, T = 0.01, 0.5
    nu = 1e-2

    def zeta_func(x, y, t):
        decay = 2 * np.pi * np.exp(-4 * (np.pi ** 2) * nu * t)
        oscillation = np.cos(2 * np.pi * x) + np.cos(2 * np.pi * y)

        return decay * oscillation

    zeta_i = FourierField.from_func(partial(zeta_func, t=0), m, p)
    zeta_f = FourierField.from_func(partial(zeta_func, t=T), m, p)

    solver = NavierStokesSolver(tau, m, p, Re=(1 / nu))
    solution = solver.solve(zeta_i, T)[-1]

    assert np.allclose(zeta_f.real, solution.real)
