import matplotlib.pyplot as plt
import numpy as np

from whirly.utils import make_grid

def show(field, ax=None, cmax=None):
    """
    Plot a FourierField object on its domain.

    Parameters
    ----------
    field : whirly.fourier.FourierField
        The function to be plotted.
    ax : matplotlib.pyplot.axis, optional
        The axis on which to make the plot. If None, an axis will be created.
    cmax : float, optional
        The maximum value to use for the symmetric colormap. If None, cmax will
        be computed as the infinity norm of field.

    Returns
    -------
    ax : matplotlib.pyplot.axis
        The axis containing the plot. If ax is passed as an argument, the
        returned axis is the same object.
    mesh : matplotlib.collections.QuadMesh
        The result of the call to pcolormesh.

    """

    if ax is None:
        fig, ax = plt.subplots()
        fig.set_size_inches(6, 6)

    ax.set_xlim(0, field.p)
    ax.set_ylim(0, field.p)

    ticks = np.linspace(0, field.p, 6)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    data = field.real
    if cmax is None:
        cmax = abs(data).max()

    x, y = make_grid(field.m, field.p)
    mesh = ax.pcolormesh(
        x, y, data,
        vmin=(-cmax),
        vmax=cmax,
        cmap='RdBu_r',
        shading='nearest'
    )

    return ax, mesh
