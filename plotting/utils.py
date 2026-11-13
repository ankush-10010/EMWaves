"""Shared plotting utilities and styling for waveguide field visualization."""

import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


# Publication-quality style configuration
STYLE_CONFIG = {
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.figsize': (8, 6),
    'figure.dpi': 100,
    'lines.linewidth': 1.5,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linewidth': 0.5,
}


def apply_style():
    """Apply publication-quality style to matplotlib."""
    plt.rcParams.update(STYLE_CONFIG)


def draw_waveguide_boundary_xy(ax, a, b, **kwargs):
    """Draw the rectangular waveguide cross-section boundary.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.
    a : float
        Waveguide width (m).
    b : float
        Waveguide height (m).
    """
    defaults = dict(linewidth=2.5, edgecolor='black', facecolor='none', zorder=10)
    defaults.update(kwargs)
    rect = patches.Rectangle((0, 0), a, b, **defaults)
    ax.add_patch(rect)
    ax.set_xlim(-0.05 * a, 1.05 * a)
    ax.set_ylim(-0.05 * b, 1.05 * b)


def draw_waveguide_walls_yz(ax, b, z_min, z_max, **kwargs):
    """Draw waveguide walls in the y-z plane.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    b : float
        Waveguide height (m).
    z_min, z_max : float
        Range of z-axis.
    """
    defaults = dict(color='black', linewidth=2.5, zorder=10)
    defaults.update(kwargs)
    ax.plot([z_min, z_max], [0, 0], **defaults)
    ax.plot([z_min, z_max], [b, b], **defaults)


def draw_waveguide_walls_xz(ax, a, z_min, z_max, **kwargs):
    """Draw waveguide walls in the x-z plane."""
    defaults = dict(color='black', linewidth=2.5, zorder=10)
    defaults.update(kwargs)
    ax.plot([z_min, z_max], [0, 0], **defaults)
    ax.plot([z_min, z_max], [a, a], **defaults)


def add_mode_label(ax, m, n, position='title'):
    """Add TM_mn mode label to plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    m, n : int
        Mode indices.
    position : str
        'title' to set as title suffix, 'annotation' to add text box.
    """
    label = f'$TM_{{{m}{n}}}$'
    if position == 'title':
        current_title = ax.get_title()
        if current_title:
            ax.set_title(f'{current_title}')
        else:
            ax.set_title(label)
    else:
        ax.annotate(
            label,
            xy=(0.95, 0.95),
            xycoords='axes fraction',
            fontsize=14,
            fontweight='bold',
            ha='right',
            va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8),
        )


def save_plot(fig, filename, directory=None):
    """Save figure as PNG and PDF.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    filename : str
        Base filename (without extension).
    directory : str or None
        Output directory. Defaults to report/figures/.
    """
    if directory is None:
        directory = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'report', 'figures'
        )
    os.makedirs(directory, exist_ok=True)

    png_path = os.path.join(directory, f'{filename}.png')
    pdf_path = os.path.join(directory, f'{filename}.pdf')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(pdf_path, bbox_inches='tight')
    return png_path, pdf_path


def get_field_colormap(field_type='diverging'):
    """Return a colormap name for the given field visualization type.

    Parameters
    ----------
    field_type : str
        'diverging' for bipolar fields (Ez), 'magnitude' for |E|/|H|,
        'phase' for phase plots.

    Returns
    -------
    str
        Matplotlib colormap name.
    """
    cmaps = {
        'diverging': 'RdBu_r',
        'magnitude': 'viridis',
        'sequential': 'plasma',
        'phase': 'hsv',
    }
    return cmaps.get(field_type, 'viridis')
