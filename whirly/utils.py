import numpy as np

# Provide short aliases for Fourier operations, taking real parts as necessary.
fft = np.fft.fft2
ifft = lambda a: np.fft.ifft2(a).real

def make_grid(m, p):
    """Makes a grid with m nodes on the square domain [0, p]^2."""
    x, y = np.mgrid[0:p:((m + 1) * 1j), 0:p:((m + 1) * 1j)]

    return x[:-1, :-1], y[:-1, :-1]

def make_wavenumbers(m, p):
    """Makes arrays of Fourier wavenumbers for m nodes on [0, p]^2."""
    freqs = 2 * np.pi * np.fft.fftfreq(m, d=(p / m))
    ell, k = np.meshgrid(freqs, freqs)

    return k, ell
