import numpy as np

class IFRK4Integrator:
    """
    Implementation of fourth-order Runge-Kutta integrating factor time stepping.

    We follow equations (18) and (19) of Yang et al. (2021), implementing an
    integrating factor scheme based on the classic RK4 scheme.

    """

    def __init__(self, tau, L, nonlinear):
        """
        Initializes an integrator given parameters of the ODE system.

        Parameters
        ----------
        tau : float
            The time step. Corresponds to h in Yang et al. (2021).
        L : float or np.ndarray
            The scalar or array of scalars representing the linear part of the
            ODE system.
        nonlinear : callable
            Function returning the nonlinear part of the ODE system. Corresponds
            to N in Yang et al. (2021).

        """

        self.tau, self.nonlinear = tau, nonlinear
        self.E = np.exp(self.tau * L)
        self.H = np.sqrt(self.E)

    def step(self, y):
        """
        Steps an array of function values forward one time step.

        Parameters
        ----------
        y : np.ndarray
            The array of function values at the current time step.

        Returns
        -------
        y_next : np.ndarray
            The array of function values at the next time step.

        """

        N = {1 : self.nonlinear(y)}
        N[2] = self.nonlinear(self.H * (y + 0.5 * self.tau * N[1]))
        N[3] = self.nonlinear(self.H * y + 0.5 * self.tau * N[2])
        N[4] = self.nonlinear(self.E * y + self.tau * self.H * N[3])

        return self.E * y + self.tau * (
            (1 / 6) * self.E * N[1] +
            (1 / 3) * self.H * N[2] +
            (1 / 3) * self.H * N[3] +
            (1 / 6) * N[4]
        )
