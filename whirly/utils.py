import numpy as np

# Provide short aliases for Fourier operations, taking real parts as necessary.
fft = np.fft.fft2
ifft = lambda a: np.fft.ifft2(a).real

def make_grid(m, L):
    """Makes a grid with m nodes on the square domain [L, L]."""

    x, y = np.mgrid[0:L:((m + 1) * 1j), 0:L:((m + 1) * 1j)]

    return x[:-1, :-1], y[:-1, :-1]
