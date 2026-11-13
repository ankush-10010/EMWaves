"""Cutoff frequency chart and dispersion diagram for TM modes."""

import math
import numpy as np
import matplotlib.pyplot as plt

from .utils import apply_style, get_field_colormap

# Physical constants
C_0 = 299792458.0
MU_0 = 4.0 * math.pi * 1e-7
EPS_0 = 1.0 / (MU_0 * C_0**2)


def plot_cutoff_chart(mode_table, operating_freq=None, fig=None, ax=None):
    """Plot horizontal bar chart of cutoff frequencies for TM modes.

    Parameters
    ----------
    mode_table : list of dict
        List with keys 'm', 'n', 'fc_ghz' (or 'fc' in Hz).
    operating_freq : float or None
        Operating frequency in Hz. If given, draws a vertical line and
        colors bars green (propagating) or red (evanescent).
    fig, ax : optional

    Returns
    -------
    fig, ax
    """
    apply_style()
    if fig is None or ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(9, 6))

    # Limit to first 15 modes
    modes = mode_table[:15]
    if not modes:
        ax.text(0.5, 0.5, 'No TM modes to display', ha='center', va='center',
                transform=ax.transAxes, fontsize=14)
        return fig, ax

    labels = [f'$TM_{{{d["m"]}{d["n"]}}}$' for d in modes]
    fc_ghz = [d.get('fc_ghz', d.get('fc', 0) / 1e9) for d in modes]

    # Color by propagation status
    if operating_freq is not None:
        f_op_ghz = operating_freq / 1e9
        colors = ['#2ecc71' if fc < f_op_ghz else '#e74c3c' for fc in fc_ghz]
    else:
        colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(fc_ghz)))

    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, fc_ghz, color=colors, edgecolor='black',
                   linewidth=0.5, height=0.6, zorder=3)

    # Add frequency labels on bars
    for i, (bar, fc) in enumerate(zip(bars, fc_ghz)):
        ax.text(bar.get_width() + max(fc_ghz) * 0.02, bar.get_y() + bar.get_height() / 2,
                f'{fc:.2f} GHz', va='center', fontsize=9)

    # Operating frequency line
    if operating_freq is not None:
        f_op_ghz = operating_freq / 1e9
        ax.axvline(x=f_op_ghz, color='blue', linestyle='--', linewidth=2,
                   label=f'Operating: {f_op_ghz:.2f} GHz', zorder=5)
        ax.legend(fontsize=10, loc='lower right')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Cutoff Frequency (GHz)', fontsize=12)
    ax.set_title('Cutoff Frequencies for TM Modes', fontsize=14, fontweight='bold')
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.3)
    fig.tight_layout()
    return fig, ax


def plot_dispersion_diagram(a, b, eps_r, mu_r, m_max=3, n_max=3,
                            freq_max_ghz=30.0, fig=None, ax=None):
    """Plot dispersion diagram (β vs frequency) for multiple TM modes.

    Parameters
    ----------
    a, b : float
        Waveguide dimensions (m).
    eps_r, mu_r : float
        Relative material parameters.
    m_max, n_max : int
        Maximum mode indices to plot.
    freq_max_ghz : float
        Maximum frequency for the plot (GHz).
    fig, ax : optional

    Returns
    -------
    fig, ax
    """
    apply_style()
    if fig is None or ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(9, 6))

    epsilon = EPS_0 * eps_r
    mu = MU_0 * mu_r
    v_medium = 1.0 / math.sqrt(mu * epsilon)

    freq_hz = np.linspace(0.01e9, freq_max_ghz * 1e9, 1000)
    freq_ghz = freq_hz / 1e9
    omega = 2.0 * np.pi * freq_hz

    # Light line in medium
    k_medium = omega / v_medium
    ax.plot(freq_ghz, k_medium, 'k--', linewidth=1.5, alpha=0.6,
            label=f'Light line ($v = c/\\sqrt{{\\epsilon_r \\mu_r}}$)')

    colors = plt.cm.tab10(np.linspace(0, 1, m_max * n_max))
    idx = 0
    for m_val in range(1, m_max + 1):
        for n_val in range(1, n_max + 1):
            kc = math.sqrt((m_val * math.pi / a)**2 + (n_val * math.pi / b)**2)
            fc = kc / (2.0 * math.pi * math.sqrt(mu * epsilon))

            # Only plot above cutoff
            mask = freq_hz > fc
            if not np.any(mask):
                idx += 1
                continue

            k = omega[mask] * math.sqrt(mu * epsilon)
            beta = np.sqrt(k**2 - kc**2)

            ax.plot(freq_ghz[mask], beta, color=colors[idx], linewidth=2,
                    label=f'$TM_{{{m_val}{n_val}}}$ ($f_c$={fc/1e9:.1f} GHz)')
            idx += 1

    ax.set_xlabel('Frequency (GHz)', fontsize=12)
    ax.set_ylabel(r'$\beta$ (rad/m)', fontsize=12)
    ax.set_title('Dispersion Diagram — TM Modes', fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left', ncol=2)
    ax.set_xlim(0, freq_max_ghz)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig, ax
