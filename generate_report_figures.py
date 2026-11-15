#!/usr/bin/env python3
"""Script to generate all publication-quality figures for the LaTeX report."""

import os
import sys

# Ensure project root is in sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from core.waveguide import RectangularWaveguide
from core.fields import TMFieldCalculator
from core.cutoff import compute_cutoff_table
from plotting.field_plots_xy import plot_E_field_xy, plot_H_field_xy, plot_Ez_contour_xy
from plotting.field_plots_yz import plot_E_field_yz, plot_H_field_yz
from plotting.cutoff_chart import plot_cutoff_chart, plot_dispersion_diagram
from plotting.utils import apply_style, save_plot

# Standard WR-90 dimensions
A = 0.02286  # 22.86 mm
B = 0.01016  # 10.16 mm
EPS_R = 1.0
MU_R = 1.0

FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'report', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

def generate_mode_plots(m, n, freq_ghz=None):
    """Generate xy and yz plots for a given mode TM_mn."""
    wg = RectangularWaveguide(A, B, EPS_R, MU_R)
    fc = wg.cutoff_frequency(m, n)
    fc_ghz = fc / 1e9

    if freq_ghz is None:
        freq_ghz = fc_ghz * 1.25  # 25% above cutoff for good propagation

    freq = freq_ghz * 1e9
    T = 1.0 / freq
    calc = TMFieldCalculator(A, B, m, n, EPS_R, MU_R, freq)

    # Compute fields with phase quadrature for combined view
    xy_ez = calc.compute_xy_fields(nx=60, ny=60, z=0, t=0)
    xy_trans = calc.compute_xy_fields(nx=60, ny=60, z=0, t=T/4)
    xy_data = {
        'Ex': xy_trans['Ex'], 'Ey': xy_trans['Ey'],
        'Ez': xy_ez['Ez'],
        'Hx': xy_trans['Hx'], 'Hy': xy_trans['Hy'],
        'x': xy_ez['x'], 'y': xy_ez['y'],
        'x1d': xy_ez['x1d'], 'y1d': xy_ez['y1d'],
    }

    # 1. E-field xy
    fig, _ = plot_E_field_xy(xy_data, A, B, m, n)
    save_plot(fig, f'E_field_xy_TM{m}{n}', FIG_DIR)
    plt.close(fig)

    # 2. H-field xy
    fig, _ = plot_H_field_xy(xy_data, A, B, m, n)
    save_plot(fig, f'H_field_xy_TM{m}{n}', FIG_DIR)
    plt.close(fig)

    # 3. Ez contour xy
    fig, _ = plot_Ez_contour_xy(xy_data, A, B, m, n)
    save_plot(fig, f'Ez_contour_xy_TM{m}{n}', FIG_DIR)
    plt.close(fig)

    # 4. Y-Z fields
    yz_data = calc.compute_yz_fields(ny=40, nz=120, x=A/2, t=0, n_wavelengths=3.0)
    fig, _ = plot_E_field_yz(yz_data, B, m, n)
    save_plot(fig, f'E_field_yz_TM{m}{n}', FIG_DIR)
    plt.close(fig)

    fig, _ = plot_H_field_yz(yz_data, B, m, n)
    save_plot(fig, f'H_field_yz_TM{m}{n}', FIG_DIR)
    plt.close(fig)

    print(f"Generated plots for TM_{m}{n} (fc = {fc_ghz:.2f} GHz, f_op = {freq_ghz:.2f} GHz)")

def main():
    apply_style()
    print("Generating report figures...")

    # Generate plots for TM11 (dominant TM mode)
    generate_mode_plots(1, 1, freq_ghz=20.0)

    # Generate plots for TM21
    generate_mode_plots(2, 1, freq_ghz=25.0)

    # Generate plots for TM12
    generate_mode_plots(1, 2, freq_ghz=35.0)

    # Cutoff chart
    table = compute_cutoff_table(A, B, EPS_R, MU_R, m_max=5, n_max=5)
    fig, _ = plot_cutoff_chart(table, operating_freq=20.0e9)
    save_plot(fig, 'cutoff_chart', FIG_DIR)
    plt.close(fig)
    print("Generated cutoff chart")

    # Dispersion diagram
    fig, _ = plot_dispersion_diagram(A, B, EPS_R, MU_R, m_max=3, n_max=3, freq_max_ghz=40.0)
    save_plot(fig, 'dispersion_diagram', FIG_DIR)
    plt.close(fig)
    print("Generated dispersion diagram")

    print(f"All figures saved to {FIG_DIR}")

if __name__ == '__main__':
    main()
