#!/usr/bin/env python
"""
Shared preview-plot rendering for flash2openpmd's CLI and GUI, so both
draw the density array exactly the same way.
"""

import matplotlib
from matplotlib.colors import LogNorm

try:
    import scienceplots  # noqa: F401
    matplotlib.style.use(["science", "ieee", "no-latex"])
    # The "ieee" style sets figure.dpi=600 for print-resolution export,
    # which blows up interactive/embedded (Tk canvas, plt.show()) figures
    # to thousands of pixels. Keep the style's colors/ticks/fonts, but use
    # a screen-appropriate dpi for on-screen display.
    matplotlib.rcParams["figure.dpi"] = 100
except ImportError:
    pass


_AXIS_LABELS_3D = ("z", "y", "x")
_AXIS_LABELS_2D = ("y", "x")


def axis_labels_for(ndim):
    """Axis labels, in array-dimension order, for a density array of this rank."""

    if ndim == 3:
        return _AXIS_LABELS_3D
    if ndim == 2:
        return _AXIS_LABELS_2D
    raise ValueError(f"Unsupported density.ndim={ndim}; expected 2 or 3")


def default_slice_indices(density):
    """
    The slice that reproduces the tool's original preview behavior: no
    slicing (full 2D plot) for a 2D array, a central y-slice for a 3D one.
    """

    if density.ndim == 3:
        labels = axis_labels_for(3)
        return {"y": density.shape[labels.index("y")] // 2}
    return {}


def make_preview_figure(density, log_scale=True, cmap="jet", title="", fig=None, slice_indices=None):
    """
    Render a preview of `density` into a matplotlib Figure, slicing along
    zero or more axes to reduce it to something plottable.

    Parameters
    ----------
    density : ndarray
        2D (axes "y","x") or 3D (axes "z","y","x") array.
    log_scale : bool
        Use log color/y scaling (falls back to linear if no positive values).
    cmap : string
        Matplotlib colormap name, used for 2D plots.
    title : string
        Figure title.
    fig : matplotlib.figure.Figure or None
        Figure to draw into (cleared first). A new Figure is created if None.
    slice_indices : dict or None
        Maps an axis label to the integer index to slice it at; axes not
        present are kept as plotted dimensions. `None` or `{}` uses
        `default_slice_indices(density)`. For a 3D array: one key leaves a
        2D plot, two keys leave a 1D line plot. For a 2D array: no keys
        leaves a 2D plot, one key leaves a 1D line plot. Slicing every axis
        (leaving nothing to plot) raises ValueError.

    Returns
    -------
    matplotlib.figure.Figure
    """

    if fig is None:
        from matplotlib.figure import Figure
        fig = Figure(figsize=(8, 6))
    else:
        fig.clf()

    labels = axis_labels_for(density.ndim)

    if not slice_indices:
        slice_indices = default_slice_indices(density)

    unknown = set(slice_indices) - set(labels)
    if unknown:
        raise ValueError(f"Unknown slice axis {unknown} for a {density.ndim}D array; expected any of {labels}")

    remaining_ndim = density.ndim - len(slice_indices)
    if remaining_ndim < 1:
        raise ValueError("Cannot slice every axis; leave at least one axis unsliced to plot.")
    if remaining_ndim > 2:
        raise ValueError(f"{remaining_ndim} axes would remain unsliced; slice enough axes to leave 1 or 2.")

    data = density
    remaining_labels = list(labels)
    clamped_indices = {}
    for label, idx in slice_indices.items():
        dim = remaining_labels.index(label)
        idx = max(0, min(int(idx), data.shape[dim] - 1))
        clamped_indices[label] = idx
        data = data.take(idx, axis=dim)
        remaining_labels.pop(dim)

    slice_desc = ", ".join(f"{label}={clamped_indices[label]}" for label in labels if label in clamped_indices)
    full_title = f"{title} ({slice_desc})".strip() if slice_desc else title

    ax = fig.add_subplot(111)

    if data.ndim == 1:
        ax.plot(data)
        if log_scale and bool((data > 0).any()):
            ax.set_yscale("log")
        ax.set_xlabel(remaining_labels[0] if remaining_labels else "index")
        ax.set_ylabel("value")
    else:
        norm = None
        if log_scale:
            positive = data[data > 0]
            if positive.size:
                norm = LogNorm(vmin=positive.min(), vmax=data.max())
        img = ax.imshow(data.T, origin="lower", norm=norm, cmap=cmap, aspect="auto")
        fig.colorbar(img, ax=ax, orientation="vertical")
        ax.set_xlabel(remaining_labels[1])
        ax.set_ylabel(remaining_labels[0])

    ax.set_title(full_title)

    return fig
