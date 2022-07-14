import matplotlib.pyplot as plt
import numpy as np

from matplotlib.animation import FuncAnimation

from whirly.fourier import FourierField
from whirly.utils import make_grid

def make_animation(outputs, output_tau, m=None):
    """
    Animate the solution to a two-dimensional PDE.

    Parameters
    ----------
    outputs : [whirly.fourier.FourierField]
        The set of snapshots of the solution, most likely the output of a call
        to the solve method of a PseudospectralSolver subclass.
    output_tau : float
        The time step between snapshots in outputs.
    m : int, optional
        If m is not None, snapshots will be resampled to m grid points.

    Returns
    -------
    animation : matplotlib.animation.FuncAnimation
        The animated solution.

    """

    p = outputs[0].p
    if m is None:
        m = outputs[0].m

    fields = [field.resample(m).real for field in outputs]
    vmax = max(abs(field).max() for field in fields)

    fig, ax = plt.subplots()
    fig.set_size_inches(6, 6)

    zero = FourierField.from_real(np.zeros((m, m)), p)
    _, mesh = plot(zero, ax, vmax)

    def update(i):
        mesh.set_array(fields[i].ravel())
        ax.set_title(f't = {(i * output_tau):.3f}')

        return mesh

    return FuncAnimation(fig, update, np.arange(len(fields)), interval=20)

def plot(field, ax=None, vmax=None):
    """
    Plot a FourierField object on its domain.

    Parameters
    ----------
    field : whirly.fourier.FourierField
        The function to be plotted.
    ax : matplotlib.pyplot.axis, optional
        The axis on which to make the plot. If None, an axis will be created.
    vmax : float, optional
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
    if vmax is None:
        vmax = abs(data).max()

    x, y = make_grid(field.m, field.p)
    mesh = ax.pcolormesh(
        x, y, data,
        vmin=(-vmax),
        vmax=vmax,
        cmap='RdBu_r',
        shading='nearest'
    )

    return ax, mesh
