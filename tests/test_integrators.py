import numpy as np

from whirly.integrators import IFRK4Integrator

def test_fourth_order():
    nonlinear = lambda y: y ** 2
    y_func = lambda t: np.exp(t + 1) / (1 - np.exp(t + 1))

    taus, errors = [0.01, 0.001], []
    for tau in taus:
        integrator = IFRK4Integrator(tau, 1, nonlinear)
        y = y_func(0)

        n_steps = int(1 / tau)
        for _ in range(n_steps):
            y = integrator.step(y)

        errors.append(abs(y - y_func(1)))

    p = np.log(errors[0] / errors[1]) / np.log(taus[0] / taus[1])

    assert round(p) == 4
