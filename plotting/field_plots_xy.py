"""2D field plots in the x-y transverse cross-section for TM modes."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm

from .utils import (
    apply_style, draw_waveguide_boundary_xy,
    get_field_colormap, add_mode_label,
)


def plot_E_field_xy(field_data, a, b, m, n, fig=None, ax=None):
    """Plot electric field in the x-y plane (transverse cross-section).

    Shows E_z as filled contour with (E_x, E_y) transverse arrows overlaid.

    Parameters
    ----------
    field_data : dict
        Output from TMFieldCalculator.compute_xy_fields().
    a, b : float
        Waveguide dimensions (m).
    m, n : int
        Mode indices.
    fig, ax : matplotlib objects, optional
        Existing figure/axes. Created if None.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes.
    """
    apply_style()
    if fig is None or ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    X = field_data['x']
    Y = field_data['y']
    Ez = field_data['Ez']
    Ex = field_data['Ex']
    Ey = field_data['Ey']

    # Ez filled contour (diverging colormap)
    vmax = np.max(np.abs(Ez)) if np.max(np.abs(Ez)) > 0 else 1.0
    try:
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    except Exception:
        norm = Normalize(vmin=-vmax, vmax=vmax)

    cf = ax.contourf(X, Y, Ez, levels=50, cmap=get_field_colormap('diverging'),
                     norm=norm)
    cbar = fig.colorbar(cf, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label(r'$E_z$ (V/m)', fontsize=11)

    # Transverse E-field arrows (subsample for clarity)
    skip = max(1, X.shape[0] // 15)
    E_mag = np.sqrt(Ex**2 + Ey**2)
    E_mag_max = np.max(E_mag) if np.max(E_mag) > 0 else 1.0

    ax.quiver(
        X[::skip, ::skip], Y[::skip, ::skip],
        Ex[::skip, ::skip], Ey[::skip, ::skip],
        color='black', alpha=0.7, scale=E_mag_max * 20,
        width=0.004, headwidth=4, headlength=5,
        zorder=5,
    )

    draw_waveguide_boundary_xy(ax, a, b)
    ax.set_xlabel('x (m)', fontsize=12)
    ax.set_ylabel('y (m)', fontsize=12)
    ax.set_title(f'$\\mathbf{{E}}$-field, $TM_{{{m}{n}}}$ mode (x-y plane)',
                 fontsize=14, fontweight='bold')
    ax.set_aspect('equal')
    fig.tight_layout()
    return fig, ax


def plot_H_field_xy(field_data, a, b, m, n, fig=None, ax=None):
    """Plot magnetic field in the x-y plane using streamlines.

    Parameters
    ----------
    field_data : dict
        Output from TMFieldCalculator.compute_xy_fields().
    a, b : float
        Waveguide dimensions (m).
    m, n : int
        Mode indices.
    fig, ax : matplotlib objects, optional

    Returns
    -------
    fig, ax
    """
    apply_style()
    if fig is None or ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    X = field_data['x']
    Y = field_data['y']
    Hx = field_data['Hx']
    Hy = field_data['Hy']

    H_mag = np.sqrt(Hx**2 + Hy**2)
    H_mag_max = np.max(H_mag) if np.max(H_mag) > 0 else 1.0

    # Use streamplot for magnetic field lines
    x1d = field_data.get('x1d', X[0, :])
    y1d = field_data.get('y1d', Y[:, 0])

    try:
        strm = ax.streamplot(
            x1d, y1d, Hx, Hy,
            color=H_mag,
            cmap=get_field_colormap('magnitude'),
            linewidth=1.5,
            density=1.5,
            arrowsize=1.5,
            zorder=3,
        )
        cbar = fig.colorbar(strm.lines, ax=ax, shrink=0.85, pad=0.02)
        cbar.set_label(r'$|\mathbf{H}_t|$ (A/m)', fontsize=11)
    except Exception:
        # Fallback to quiver if streamplot fails
        skip = max(1, X.shape[0] // 15)
        ax.quiver(
            X[::skip, ::skip], Y[::skip, ::skip],
            Hx[::skip, ::skip], Hy[::skip, ::skip],
            H_mag[::skip, ::skip],
            cmap=get_field_colormap('magnitude'),
            alpha=0.8, zorder=3,
        )

    draw_waveguide_boundary_xy(ax, a, b)
    ax.set_xlabel('x (m)', fontsize=12)
    ax.set_ylabel('y (m)', fontsize=12)
    ax.set_title(f'$\\mathbf{{H}}$-field, $TM_{{{m}{n}}}$ mode (x-y plane)',
                 fontsize=14, fontweight='bold')
    ax.set_aspect('equal')
    fig.tight_layout()
    return fig, ax


def plot_Ez_contour_xy(field_data, a, b, m, n, fig=None, ax=None):
    """Plot E_z as a standalone contour plot.

    Parameters
    ----------
    field_data : dict
    a, b : float
    m, n : int
    fig, ax : optional

    Returns
    -------
    fig, ax
    """
    apply_style()
    if fig is None or ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    X = field_data['x']
    Y = field_data['y']
    Ez = field_data['Ez']

    vmax = np.max(np.abs(Ez)) if np.max(np.abs(Ez)) > 0 else 1.0
    try:
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    except Exception:
        norm = Normalize(vmin=-vmax, vmax=vmax)

    cf = ax.contourf(X, Y, Ez, levels=50, cmap=get_field_colormap('diverging'), norm=norm)
    cs = ax.contour(X, Y, Ez, levels=15, colors='k', linewidths=0.5, alpha=0.5)
    cbar = fig.colorbar(cf, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label(r'$E_z$ (V/m)', fontsize=11)

    draw_waveguide_boundary_xy(ax, a, b)
    ax.set_xlabel('x (m)', fontsize=12)
    ax.set_ylabel('y (m)', fontsize=12)
    ax.set_title(f'$E_z$, $TM_{{{m}{n}}}$ mode (x-y plane)',
                 fontsize=14, fontweight='bold')
    ax.set_aspect('equal')
    fig.tight_layout()
    return fig, ax
