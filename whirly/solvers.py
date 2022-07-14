import numpy as np

from whirly.integrators import IFRK4Integrator
from whirly.utils import make_wavenumbers

class PseudospectralSolver:
    """
    Abstract class for pseudospectral solutions to doubly-periodic PDEs.

    Subclasses should
        - define an __init__ method that accepts a time step tau, grid number m,
            and domain size p, and pass all to PseudospectralSolver.__init__
        - define an array L with the coefficient for each Fourier mode in the
            linear part of the PDE
        - define a nonlinear method that accepts a FourierField representing
            the vorticity and returns a FourierField representing the nonlinear
            part of the PDE

    Instances of PseudospectralSolver should not be created directly.

    """

    def __init__(self, tau, m, p):
        """
        Initialize a PseudospectralSolver.

        Note that PseudospectralSolver is meant to be subclassed and instances
        should not be created directly.

        Parameters
        ----------
        tau : float
            The time step.
        m : int
            The number of grid nodes.
        p : float
            The length of a side of the domain.

        """

        self.tau = tau
        self.m = m
        self.p = p

    def solve(self, zeta_initial, T, output_tau=None):
        """
        Solve a nonlinear PDE on a doubly-periodic square domain.

        This method can be used by classes that subclass PseudospectralSolver
        and relies on subclasses implementing nonlinear and defining L.

        Parameters
        ----------
        zeta_initial : whirly.fourier.FourierField
            A FourierField representing the initial solution field.
        T : float
            The final time to solve until.
        output_tau : float
            The time frequency which the solution field should be output. Must
            be smaller than the effective time step, which is calculated by
            adjusting the tau attribute of the instance such that it divides
            evenly into T. If None, the solution will be output every time step.

        Returns
        -------
        solution : [whirly.fourier.FourierField]
            A list of FourierField instances representing the solution field at
            time frequency governed by output_tau

        """

        n_steps = round(T / self.tau)
        tau = T / n_steps

        if output_tau is None:
            output_tau = tau

        skip = round(output_tau / tau)
        zeta = zeta_initial
        outputs = [zeta]

        integrator = IFRK4Integrator(tau, self.L, self.nonlinear)
        for i in range(1, n_steps + 1):
            zeta = integrator.step(zeta)
            if i % skip == 0:
                outputs.append(zeta)

        return outputs

class NavierStokesSolver(PseudospectralSolver):
    """PseudospectralSolver for two-dimensional Navier-Stokes."""

    def __init__(self, tau, m, p, Re):
        """
        Initializes a NavierStokesSolver.

        Parameters
        ----------
        tau : float
            The time step.
        m : int
            The number of grid nodes.
        p : float
            The length of a side of the domain.
        Re : float
            The Reynolds number.

        """

        super().__init__(tau, m, p)

        self.k, self.ell = make_wavenumbers(self.m, self.p)
        K_sq = (self.k ** 2) + (self.ell ** 2)
        self.L = -(1 / Re) * K_sq

        mask = K_sq != 0
        self.inverse_laplacian = np.zeros_like(K_sq)
        self.inverse_laplacian[mask] = -1 / K_sq[mask]

    def nonlinear(self, zeta):
        """
        Calculates the advection of vorticity.

        First, the streamfunction is calculated by solving Laplace's equation,
        and then u and v are found through spectral differentiation and used to
        compute the (negative) advection term.

        Parameters
        ----------
        zeta : whirly.fourier.FourierField
            The current vorticity field.

        Returns
        -------
        advection : whirly.fourier.FourierField
            The nonlinear advection term -u dot grad(zeta).

        """

        psi = self.inverse_laplacian * zeta
        u = -1j * self.ell * psi
        v = 1j * self.k * psi

        zeta_x = 1j * self.k * zeta
        zeta_y = 1j * self.ell * zeta

        return -(u * zeta_x + v * zeta_y)
