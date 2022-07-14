from itertools import product

import numpy as np

import whirly.utils as utils

class FourierField:
    """
    Array of Fourier coefficients with automatic spectral transforms.

    FourierField objects represent functions given by their Fourier coefficients
    on a given doubly-periodic square domain. They define arithmetic operations
    in real and spectral space, antialiasing when necessary, and can also be
    resampled to different grid resolutions.

    """

    # Prevent numpy from broadcasting instead of using arithmetic defined here.
    __array_ufunc__ = None

    def __init__(self, hats, p):
        """
        Initializes a FourierField given Fourier coefficients and a domain size.

        Note that in general FourierField objects should be created using one of
        the class constructors defined below, which in turn call this method.
        This method also checks that the grid has an odd number of points, to
        prevent differentiation errors from the unmatched mode.

        Parameters
        ----------
        hats : np.ndarray
            Two-dimensional array of Fourier coefficients, in the order used by
            numpy FFT routines.
        p : float
            The length of a side of the domain.

        Raises
        ------
        ValueError
            If the grid does not have an odd number of points.

        """

        self.hats = hats
        self.p = p

        if self.m % 2 == 0:
            raise ValueError('grid size must be odd')

    def __add__(self, other):
        """
        Adds one FourierField to another or to a real function.

        Parameters
        ----------
        other : int or float or complex or np.ndarray or FourierField
            Value to add to. If a numeric type or array, other represents a
            function by its values in real space, and addition is performed in
            real space after an IFFT. If a FourierField, other represents a
            function by its Fourier coefficients, and linearity of the FFT is
            leveraged to perform the addition in spectral space.

        Returns
        -------
        field : FourierField
            Spectral representation of the sum of self and other.

        """

        if isinstance(other, (int, float, complex, np.ndarray)):
            return type(self).from_data(self.real + other, self.p)

        return type(self)(self.hats + other.hats, self.p)

    def __sub__(self, other):
        """Subtraction is just negation followed by addition."""
        return self + -other

    def __mul__(self, other):
        """
        Multiplies one FourierField by another or by numerical values.

        Parameters
        ----------
        other : int or float or complex or np.ndarray or FourierField
            Value to multiply by. If a numeric type or array, other represents a
            scaling factor or array of factors for the Fourier coefficients of
            self, and multiplication is performed in spectral space. If other is
            a FourierField, multiplication is performed in real space after an
            IFFT, with both fields temporarily resampled at higher resolution to
            avoid aliasing effects.

        Returns
        -------
        field : FourierField
            Spectral representation of the product of self and other.

        """

        if isinstance(other, (int, float, complex, np.ndarray)):
            return type(self)(other * self.hats, self.p)

        m_aa = 3 * (self.m + 1) // 2
        if m_aa % 2 == 0:
            m_aa = m_aa + 1

        data = self.resample(m_aa).real * other.resample(m_aa).real

        return type(self).from_real(data, self.p).resample(self.m)

    def __rmul__(self, other):
        """Makes FourierField multiplication is always commutative."""
        return self * other

    def __neg__(self):
        """Allows FourierField objects to be negated directly."""
        return -1 * self

    @property
    def m(self):
        """The number of nodes in the spectral representation."""
        return self.hats.shape[0]

    @property
    def real(self):
        """The real-space representation of the function, computed via IFFT."""
        return utils.ifft(self.hats)

    def resample(self, m):
        """
        Resample a FourierField to a new spatial resolution.

        FourierField objects can be resampled to finer resolutions, in which
        case the array of coefficients is padded with zeros, or to coarser
        resolutions, in which case the coefficients are truncated.

        Parameters
        ----------
        m : int
            The new spatial resolution.

        Returns
        -------
        field : FourierField
            The spectral representation of the function resampled at the new
            spatial resolution.

        """

        if m == self.m:
            return self

        c = (min(m, self.m) + 1) // 2
        slices = [slice(c), slice(-(c - 1), None)]

        hats = np.zeros((m, m), dtype=self.hats.dtype)
        for idx, jdx in product(slices, slices):
            hats[idx, jdx] = self.hats[idx, jdx]

        hats = ((m / self.m) ** 2) * hats

        return type(self)(hats, self.p)

    @classmethod
    def from_func(cls, func, m, p=1):
        """
        Construct a FourierField from a real-valued function.

        Parameters
        ----------
        func : callable
            Real-valued function of two spatial coordinates. Should have the
            signature func(x, y), where x and y are numpy arrays.
        m : int
            The desired grid resolution.
        p : float, optional
            The size of the square domain. Default is 1, in which case the
            domain is the doubly-periodic unit square.

        Returns
        -------
        field : FourierField
            Spectral representation of the function applied to [0, p]^2.

        """

        x, y = utils.make_grid(m, p)
        hats = utils.fft(func(x, y))

        return cls(hats, p)

    @classmethod
    def from_real(cls, real, p=1):
        """
        Construct a FourierField from gridded real data on a square domain.

        Parameters
        ----------
        real : np.ndarray
            Real-space function values on a square grid.
        p : float, optional
            The size of the square domain. Default is 1, in which case the
            domain is the doubly-periodic unit square.

        Returns
        -------
        field : FourierField
            Spectral representation of the function.

        """

        return cls(utils.fft(real), p)
