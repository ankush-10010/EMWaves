#!/usr/bin/env python3
"""Command-Line Interface (CLI) for Rectangular Waveguide TM Mode Analyzer.

Allows fast headless computation, parameter sweeps, and plot generation
without opening the graphical interface.

Usage:
    python cli.py --a 0.02286 --b 0.01016 --m 1 --n 1 --freq 20.0 --table --plot
"""

import argparse
import sys
import os

# Add project root to sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from core.waveguide import RectangularWaveguide
from core.fields import TMFieldCalculator
from core.cutoff import compute_cutoff_table, format_cutoff_table
from plotting.field_plots_xy import plot_E_field_xy, plot_H_field_xy
from plotting.field_plots_yz import plot_E_field_yz, plot_H_field_yz
from plotting.cutoff_chart import plot_cutoff_chart, plot_dispersion_diagram
from plotting.utils import apply_style, save_plot

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rectangular Waveguide TM Mode Analyzer CLI",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--a", type=float, default=0.02286, help="Waveguide width in meters")
    parser.add_argument("--b", type=float, default=0.01016, help="Waveguide height in meters")
    parser.add_argument("--m", type=int, default=1, help="Mode index along x (m >= 1)")
    parser.add_argument("--n", type=int, default=1, help="Mode index along y (n >= 1)")
    parser.add_argument("--eps_r", type=float, default=1.0, help="Relative permittivity")
    parser.add_argument("--mu_r", type=float, default=1.0, help="Relative permeability")
    parser.add_argument("--freq", type=float, default=20.0, help="Operating frequency in GHz")
    parser.add_argument("--m-max", type=int, default=5, help="Max mode index m for mode table")
    parser.add_argument("--n-max", type=int, default=5, help="Max mode index n for mode table")
    parser.add_argument("--table", action="store_true", help="Print mode cutoff frequency table")
    parser.add_argument("--plot", action="store_true", help="Generate and save publication plots")
    parser.add_argument("--out-dir", type=str, default="./output_plots", help="Output directory for plots")
    return parser.parse_args()


def main():
    args = parse_args()
    freq_hz = args.freq * 1e9

    print("=" * 65)
    print("      RECTANGULAR WAVEGUIDE TM MODE ANALYZER")
    print("=" * 65)

    wg = RectangularWaveguide(args.a, args.b, eps_r=args.eps_r, mu_r=args.mu_r)
    calc = TMFieldCalculator(args.a, args.b, args.m, args.n, args.eps_r, args.mu_r, freq_hz)

    fc = wg.cutoff_frequency(args.m, args.n)
    fc_ghz = fc / 1e9
    kc = wg.cutoff_wavenumber(args.m, args.n)

    print(f"Waveguide Dimensions : a = {args.a*1e3:.2f} mm, b = {args.b*1e3:.2f} mm")
    print(f"Dielectric Medium    : eps_r = {args.eps_r:.2f}, mu_r = {args.mu_r:.2f}")
    print(f"Operating Frequency  : {args.freq:.4f} GHz")
    print("-" * 65)
    print(f"Mode Analyzed        : TM_{args.m}{args.n}")
    print(f"Cutoff Wavenumber kc : {kc:.4f} rad/m")
    print(f"Cutoff Frequency fc  : {fc_ghz:.4f} GHz ({fc:.4e} Hz)")

    if calc.is_propagating:
        beta = wg.propagation_constant(args.m, args.n, freq_hz)
        lg = wg.guide_wavelength(args.m, args.n, freq_hz)
        vp = wg.phase_velocity(args.m, args.n, freq_hz)
        vg = wg.group_velocity(args.m, args.n, freq_hz)
        ztm = wg.wave_impedance(args.m, args.n, freq_hz)

        print(f"Propagation Status   : PROPAGATING (f > fc)")
        print(f"Phase Constant beta  : {beta:.4f} rad/m")
        print(f"Guide Wavelength lg  : {lg*1e3:.4f} mm")
        print(f"Phase Velocity vp    : {vp:.4e} m/s ({vp/3e8:.3f} c)")
        print(f"Group Velocity vg    : {vg:.4e} m/s ({vg/3e8:.3f} c)")
        print(f"Wave Impedance Z_TM  : {ztm:.4f} Ohms")
    else:
        print(f"Propagation Status   : EVANESCENT (f < fc)")
        print(f"Attenuation Const a  : {calc.alpha:.4f} Np/m ({calc.alpha*8.686:.2f} dB/m)")

    print("=" * 65)

    if args.table or args.plot:
        table = compute_cutoff_table(args.a, args.b, args.eps_r, args.mu_r,
                                     m_max=args.m_max, n_max=args.n_max)
        if args.table:
            print("\nCutoff Frequency Table:")
            print(format_cutoff_table(table[:15]))
            print()

    if args.plot:
        apply_style()
        os.makedirs(args.out_dir, exist_ok=True)
        print(f"Generating field plots in '{args.out_dir}'...")

        T = 1.0 / freq_hz
        xy_ez = calc.compute_xy_fields(nx=50, ny=50, z=0, t=0)
        xy_trans = calc.compute_xy_fields(nx=50, ny=50, z=0, t=T/4)
        xy_data = {
            'Ex': xy_trans['Ex'], 'Ey': xy_trans['Ey'],
            'Ez': xy_ez['Ez'],
            'Hx': xy_trans['Hx'], 'Hy': xy_trans['Hy'],
            'x': xy_ez['x'], 'y': xy_ez['y'],
            'x1d': xy_ez['x1d'], 'y1d': xy_ez['y1d'],
        }
        yz_data = calc.compute_yz_fields(ny=40, nz=100, x=args.a/2, t=0)

        # Plot 1: E-field xy
        fig, _ = plot_E_field_xy(xy_data, args.a, args.b, args.m, args.n)
        save_plot(fig, f"E_field_xy_TM{args.m}{args.n}", args.out_dir)
        plt.close(fig)

        # Plot 2: H-field xy
        fig, _ = plot_H_field_xy(xy_data, args.a, args.b, args.m, args.n)
        save_plot(fig, f"H_field_xy_TM{args.m}{args.n}", args.out_dir)
        plt.close(fig)

        # Plot 3: E-field yz
        fig, _ = plot_E_field_yz(yz_data, args.b, args.m, args.n)
        save_plot(fig, f"E_field_yz_TM{args.m}{args.n}", args.out_dir)
        plt.close(fig)

        # Plot 4: H-field yz
        fig, _ = plot_H_field_yz(yz_data, args.b, args.m, args.n)
        save_plot(fig, f"H_field_yz_TM{args.m}{args.n}", args.out_dir)
        plt.close(fig)

        # Plot 5: Cutoff chart
        fig, _ = plot_cutoff_chart(table, operating_freq=freq_hz)
        save_plot(fig, "cutoff_chart", args.out_dir)
        plt.close(fig)

        # Plot 6: Dispersion
        fig, _ = plot_dispersion_diagram(args.a, args.b, args.eps_r, args.mu_r)
        save_plot(fig, "dispersion_diagram", args.out_dir)
        plt.close(fig)

        print(f"Saved 6 publication figures to '{args.out_dir}'.")


if __name__ == '__main__':
    main()
