"""2D field plots in the y-z longitudinal cross-section for TM modes."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm

from .utils import (
    apply_style, draw_waveguide_walls_yz,
    get_field_colormap,
)


def plot_E_field_yz(field_data, b, m, n, fig=None, ax=None):
    """Plot electric field in the y-z plane (longitudinal section).

    Shows (E_y, E_z) as quiver arrows colored by magnitude.

    Parameters
    ----------
    field_data : dict
        Output from TMFieldCalculator.compute_yz_fields().
    b : float
        Waveguide height (m).
    m, n : int
        Mode indices.
    fig, ax : optional

    Returns
    -------
    fig, ax
    """
    apply_style()
    if fig is None or ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    Y = field_data['y']
    Z = field_data['z']
    Ey = field_data['Ey']
    Ez = field_data['Ez']

    E_mag = np.sqrt(Ey**2 + Ez**2)
    E_mag_max = np.max(E_mag) if np.max(E_mag) > 0 else 1.0

    # Background contour of |E|
    cf = ax.contourf(Z, Y, E_mag, levels=30, cmap=get_field_colormap('magnitude'),
                     alpha=0.4)

    # Quiver arrows
    skip_y = max(1, Y.shape[0] // 20)
    skip_z = max(1, Y.shape[1] // 25)

    # In quiver: x-axis is z, y-axis is y, so U=Ez, V=Ey
    q = ax.quiver(
        Z[::skip_y, ::skip_z], Y[::skip_y, ::skip_z],
        Ez[::skip_y, ::skip_z], Ey[::skip_y, ::skip_z],
        E_mag[::skip_y, ::skip_z],
        cmap=get_field_colormap('magnitude'),
        alpha=0.85, scale=E_mag_max * 25,
        width=0.003, headwidth=4, headlength=5,
        zorder=5,
    )
    cbar = fig.colorbar(q, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label(r'$|\mathbf{E}|$ (V/m)', fontsize=11)

    z_min, z_max = Z.min(), Z.max()
    draw_waveguide_walls_yz(ax, b, z_min, z_max)

    ax.set_xlabel('z (m)', fontsize=12)
    ax.set_ylabel('y (m)', fontsize=12)
    ax.set_title(f'$\\mathbf{{E}}$-field, $TM_{{{m}{n}}}$ mode (y-z plane)',
                 fontsize=14, fontweight='bold')
    ax.set_ylim(-0.05 * b, 1.05 * b)
    fig.tight_layout()
    return fig, ax


def plot_H_field_yz(field_data, b, m, n, fig=None, ax=None):
    """Plot magnetic field in the y-z plane.

    For TM modes, H_x is perpendicular to the y-z plane (out-of-plane).
    Shows H_x as filled contour and H_y as arrows if significant.

    Parameters
    ----------
    field_data : dict
        Output from TMFieldCalculator.compute_yz_fields().
    b : float
        Waveguide height (m).
    m, n : int
        Mode indices.
    fig, ax : optional

    Returns
    -------
    fig, ax
    """
    apply_style()
    if fig is None or ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    Y = field_data['y']
    Z = field_data['z']
    Hx = field_data['Hx']
    Hy = field_data.get('Hy', np.zeros_like(Hx))

    # H_x contour (out-of-plane component)
    vmax = np.max(np.abs(Hx)) if np.max(np.abs(Hx)) > 0 else 1.0
    try:
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    except Exception:
        norm = Normalize(vmin=-vmax, vmax=vmax)

    cf = ax.contourf(Z, Y, Hx, levels=50, cmap=get_field_colormap('diverging'),
                     norm=norm)
    cbar = fig.colorbar(cf, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label(r'$H_x$ (A/m) [out-of-plane]', fontsize=11)

    # Overlay Hy as arrows if non-negligible
    Hy_max = np.max(np.abs(Hy))
    if Hy_max > 1e-15:
        skip_y = max(1, Y.shape[0] // 15)
        skip_z = max(1, Y.shape[1] // 20)
        # Arrows in y-z plane: U=direction along z (Hy component along z is 0 for Hz=0)
        # Actually Hy is along y-hat, so in y-z plane plot: U=0, V=Hy
        ax.quiver(
            Z[::skip_y, ::skip_z], Y[::skip_y, ::skip_z],
            np.zeros_like(Hy[::skip_y, ::skip_z]), Hy[::skip_y, ::skip_z],
            color='black', alpha=0.5, scale=Hy_max * 20,
            width=0.003, zorder=5,
        )

    z_min, z_max = Z.min(), Z.max()
    draw_waveguide_walls_yz(ax, b, z_min, z_max)

    ax.set_xlabel('z (m)', fontsize=12)
    ax.set_ylabel('y (m)', fontsize=12)
    ax.set_title(f'$\\mathbf{{H}}$-field, $TM_{{{m}{n}}}$ mode (y-z plane)',
                 fontsize=14, fontweight='bold')
    ax.set_ylim(-0.05 * b, 1.05 * b)
    fig.tight_layout()
    return fig, ax
