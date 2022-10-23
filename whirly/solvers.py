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

    def solve(self, q_initial, T, output_tau=None):
        """
        Solve a nonlinear PDE on a doubly-periodic square domain.

        This method can be used by classes that subclass PseudospectralSolver
        and relies on subclasses implementing nonlinear and defining L.

        Parameters
        ----------
        q_initial : whirly.fourier.FourierField
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
        q = q_initial
        outputs = [q]

        integrator = IFRK4Integrator(tau, self.L, self.nonlinear)
        for i in range(1, n_steps + 1):
            q = integrator.step(q)
            if i % skip == 0:
                outputs.append(q)

        return outputs

class VorticitySolver(PseudospectralSolver):
    """
    Abstract PseudospectralSolver subclass for vorticity equations.

    VorticitySolver subclasses assume that the solution field is a vorticity q
    that can be inverted to find a streamfunction, and consequently to give the
    advecting velocity. Subclasses should implement a make_operator method
    returning an array such that in Fourier space operator * psi = q holds.

    """

    def __init__(self, tau, m, p, Re):
        """
        Initializes a NavierStokesSolver.

        This method calls the make_operator method defined by subclasses and
        uses it to set the inverse_operator attribute of the object, which is
        then used to invert the vorticity in the nonlinear method.

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
        self.K_sq = (self.k ** 2) + (self.ell ** 2)
        self.L = -(1 / Re) * self.K_sq

        operator = self.make_operator()
        mask = operator != 0

        self.inverse_operator = np.zeros_like(operator)
        self.inverse_operator[mask] = 1 / operator[mask]

    def nonlinear(self, q):
        """
        Calculates the advection of vorticity.

        First, the streamfunction is calculated using the inverse_operator
        attribute, and then u and v are found through spectral differentiation
        and used to compute the (negative) advection term.

        Parameters
        ----------
        q : whirly.fourier.FourierField
            The current vorticity field.

        Returns
        -------
        advection : whirly.fourier.FourierField
            The nonlinear advection term -u dot grad(q).

        """

        psi = self.inverse_operator * q
        u = -1j * self.ell * psi
        v = 1j * self.k * psi

        q_x = 1j * self.k * q
        q_y = 1j * self.ell * q

        return -(u * q_x + v * q_y)

class NavierStokesSolver(VorticitySolver):
    """PseudospectralSolver for two-dimensional Navier-Stokes."""

    def make_operator(self):
        """The vorticity is simply the Laplacian of the streamfunction."""
        return -self.K_sq
