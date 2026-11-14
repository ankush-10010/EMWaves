"""Main Tkinter GUI application for Rectangular Waveguide TM Mode Analyzer."""

import sys
import os
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from core.waveguide import RectangularWaveguide
from core.fields import TMFieldCalculator
from core.cutoff import compute_cutoff_table
from plotting.field_plots_xy import plot_E_field_xy, plot_H_field_xy
from plotting.field_plots_yz import plot_E_field_yz, plot_H_field_yz
from plotting.cutoff_chart import plot_cutoff_chart, plot_dispersion_diagram
from plotting.utils import apply_style, save_plot


class WaveguideApp(tk.Tk):
    """Main application window for the Rectangular Waveguide TM Mode Analyzer."""

    def __init__(self):
        super().__init__()

        self.title("Rectangular Waveguide TM Mode Analyzer")
        self.geometry("1300x800")
        self.minsize(1100, 700)

        # Apply matplotlib style
        apply_style()

        # Configure ttk style
        self.style = ttk.Style(self)
        try:
            self.style.theme_use('clam')
        except Exception:
            pass

        self.style.configure('Title.TLabel', font=('Helvetica', 16, 'bold'))
        self.style.configure('Header.TLabel', font=('Helvetica', 11, 'bold'))
        self.style.configure('Result.TLabel', font=('Consolas', 10))
        self.style.configure('Status.TLabel', font=('Helvetica', 11, 'bold'))
        self.style.configure('Analyze.TButton', font=('Helvetica', 11, 'bold'))

        self._create_menu()
        self._create_layout()

        # Store figures for export
        self.figures = {}

    def _create_menu(self):
        """Create menu bar."""
        menubar = tk.Menu(self)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Export All Plots", command=self._export_plots)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.quit)
        menubar.add_cascade(label="File", menu=file_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=self._show_about)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.config(menu=menubar)

    def _create_layout(self):
        """Create the main layout: left panel (inputs) | right panel (plots + results)."""
        # Main paned window
        main_pane = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        main_pane.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # ===== LEFT PANEL =====
        left_frame = ttk.Frame(main_pane, width=320)
        main_pane.add(left_frame, weight=0)

        # Title
        ttk.Label(left_frame, text="⚡ TM Mode Analyzer",
                  style='Title.TLabel').pack(pady=(10, 15))

        # --- Waveguide Parameters ---
        params_frame = ttk.LabelFrame(left_frame, text=" Waveguide Parameters ",
                                       padding=10)
        params_frame.pack(fill=tk.X, padx=10, pady=5)

        self.var_a = self._add_entry(params_frame, "Width a (m):", "0.02286", 0)
        self.var_b = self._add_entry(params_frame, "Height b (m):", "0.01016", 1)

        # --- Mode Selection ---
        mode_frame = ttk.LabelFrame(left_frame, text=" Mode Selection ", padding=10)
        mode_frame.pack(fill=tk.X, padx=10, pady=5)

        ttk.Label(mode_frame, text="Mode index m:").grid(row=0, column=0, sticky='w', pady=3)
        self.spin_m = ttk.Spinbox(mode_frame, from_=1, to=20, width=8)
        self.spin_m.set(1)
        self.spin_m.grid(row=0, column=1, padx=5, pady=3)

        ttk.Label(mode_frame, text="Mode index n:").grid(row=1, column=0, sticky='w', pady=3)
        self.spin_n = ttk.Spinbox(mode_frame, from_=1, to=20, width=8)
        self.spin_n.set(1)
        self.spin_n.grid(row=1, column=1, padx=5, pady=3)

        # --- Medium Properties ---
        medium_frame = ttk.LabelFrame(left_frame, text=" Medium Properties ", padding=10)
        medium_frame.pack(fill=tk.X, padx=10, pady=5)

        self.var_eps_r = self._add_entry(medium_frame, "Rel. Permittivity εr:", "1.0", 0)
        self.var_mu_r = self._add_entry(medium_frame, "Rel. Permeability μr:", "1.0", 1)

        # --- Operating Frequency ---
        freq_frame = ttk.LabelFrame(left_frame, text=" Operating Frequency ", padding=10)
        freq_frame.pack(fill=tk.X, padx=10, pady=5)

        self.var_freq = self._add_entry(freq_frame, "Frequency (GHz):", "15.0", 0)

        # --- Buttons ---
        btn_frame = ttk.Frame(left_frame)
        btn_frame.pack(fill=tk.X, padx=10, pady=15)

        analyze_btn = ttk.Button(btn_frame, text="🔬  Analyze",
                                  command=self.analyze, style='Analyze.TButton')
        analyze_btn.pack(fill=tk.X, pady=3, ipady=5)

        reset_btn = ttk.Button(btn_frame, text="↺  Reset", command=self._reset)
        reset_btn.pack(fill=tk.X, pady=3)

        export_btn = ttk.Button(btn_frame, text="💾  Export Plots",
                                 command=self._export_plots)
        export_btn.pack(fill=tk.X, pady=3)

        # --- Results ---
        results_frame = ttk.LabelFrame(left_frame, text=" Results ", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        self.results_text = tk.Text(results_frame, height=12, font=('Consolas', 9),
                                     wrap=tk.WORD, bg='#f8f9fa', relief=tk.FLAT,
                                     borderwidth=1)
        results_scroll = ttk.Scrollbar(results_frame, orient=tk.VERTICAL,
                                        command=self.results_text.yview)
        self.results_text.configure(yscrollcommand=results_scroll.set)
        self.results_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        results_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        self.results_text.insert(tk.END, "Click 'Analyze' to compute results.\n")
        self.results_text.config(state=tk.DISABLED)

        # ===== RIGHT PANEL =====
        right_frame = ttk.Frame(main_pane)
        main_pane.add(right_frame, weight=1)

        # Tabbed notebook for plots
        self.notebook = ttk.Notebook(right_frame)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # Create tabs
        self.tab_names = [
            "E-field (x-y)", "H-field (x-y)",
            "E-field (y-z)", "H-field (y-z)",
            "Cutoff Chart", "Dispersion"
        ]
        self.tab_frames = {}
        self.tab_canvases = {}
        self.tab_toolbars = {}
        self.tab_figures = {}

        for name in self.tab_names:
            frame = ttk.Frame(self.notebook)
            self.notebook.add(frame, text=name)
            self.tab_frames[name] = frame

            fig = Figure(figsize=(8, 6), dpi=100)
            canvas = FigureCanvasTkAgg(fig, master=frame)
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

            toolbar_frame = ttk.Frame(frame)
            toolbar_frame.pack(fill=tk.X)
            toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
            toolbar.update()

            self.tab_figures[name] = fig
            self.tab_canvases[name] = canvas
            self.tab_toolbars[name] = toolbar

    def _add_entry(self, parent, label, default, row):
        """Add a labeled entry field."""
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky='w', pady=3)
        var = tk.StringVar(value=default)
        entry = ttk.Entry(parent, textvariable=var, width=15)
        entry.grid(row=row, column=1, padx=5, pady=3)
        return var

    def _get_params(self):
        """Read and validate all input parameters."""
        try:
            a = float(self.var_a.get())
            b = float(self.var_b.get())
            m = int(self.spin_m.get())
            n = int(self.spin_n.get())
            eps_r = float(self.var_eps_r.get())
            mu_r = float(self.var_mu_r.get())
            freq_ghz = float(self.var_freq.get())
            freq = freq_ghz * 1e9
        except ValueError as e:
            raise ValueError(f"Invalid input: {e}")

        if a <= 0 or b <= 0:
            raise ValueError("Waveguide dimensions a, b must be positive.")
        if m < 1 or n < 1:
            raise ValueError("TM mode indices m, n must be ≥ 1.")
        if eps_r <= 0 or mu_r <= 0:
            raise ValueError("εr, μr must be positive.")
        if freq <= 0:
            raise ValueError("Frequency must be positive.")

        return a, b, m, n, eps_r, mu_r, freq

    def analyze(self):
        """Run the analysis: compute fields, update plots, display results."""
        try:
            a, b, m, n, eps_r, mu_r, freq = self._get_params()
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            return

        try:
            # Create waveguide and field calculator
            wg = RectangularWaveguide(a, b, eps_r, mu_r)
            calc = TMFieldCalculator(a, b, m, n, eps_r, mu_r, freq)

            # Compute fields
            # For x-y plane: use t such that transverse fields are maximized
            # At z=0, t=0: Ez ~ cos(0)=1 (max), Ex,Ey,Hx,Hy ~ sin(0)=0
            # At z=0, t=T/4: Ez ~ cos(pi/2)=0, Ex,Ey,Hx,Hy ~ sin(pi/2)=max
            # We compute both and use appropriate one, OR use amplitude patterns
            T = 1.0 / freq if freq > 0 else 1.0
            # For E-field xy: show Ez contour (t=0) with transverse arrows (t=T/4)
            xy_data_ez = calc.compute_xy_fields(nx=50, ny=50, z=0, t=0)
            xy_data_trans = calc.compute_xy_fields(nx=50, ny=50, z=0, t=T/4)
            # Merge: use Ez from t=0 snapshot, transverse from t=T/4
            xy_data = {
                'Ex': xy_data_trans['Ex'], 'Ey': xy_data_trans['Ey'],
                'Ez': xy_data_ez['Ez'],
                'Hx': xy_data_trans['Hx'], 'Hy': xy_data_trans['Hy'],
                'x': xy_data_ez['x'], 'y': xy_data_ez['y'],
                'x1d': xy_data_ez['x1d'], 'y1d': xy_data_ez['y1d'],
            }
            yz_data = calc.compute_yz_fields(ny=40, nz=100, x=a/2, t=0)

            # Compute mode table
            m_max = max(m + 2, 5)
            n_max = max(n + 2, 5)
            mode_table = compute_cutoff_table(a, b, eps_r, mu_r, m_max, n_max)

            # Update all plots
            self._update_plot("E-field (x-y)", plot_E_field_xy, xy_data, a, b, m, n)
            self._update_plot("H-field (x-y)", plot_H_field_xy, xy_data, a, b, m, n)
            self._update_plot("E-field (y-z)", plot_E_field_yz, yz_data, b, m, n)
            self._update_plot("H-field (y-z)", plot_H_field_yz, yz_data, b, m, n)
            self._update_cutoff_chart(mode_table, freq)
            self._update_dispersion(a, b, eps_r, mu_r)

            # Update results
            self._update_results(wg, calc, m, n, freq, mode_table)

        except Exception as e:
            messagebox.showerror("Computation Error", f"Error during analysis:\n{e}")
            import traceback
            traceback.print_exc()

    def _update_plot(self, tab_name, plot_func, *args):
        """Clear and redraw a plot tab."""
        fig = self.tab_figures[tab_name]
        fig.clear()
        ax = fig.add_subplot(111)

        if tab_name in ("E-field (y-z)", "H-field (y-z)"):
            plot_func(*args, fig=fig, ax=ax)
        else:
            plot_func(*args, fig=fig, ax=ax)

        self.tab_canvases[tab_name].draw()
        self.figures[tab_name] = fig

    def _update_cutoff_chart(self, mode_table, freq):
        """Update cutoff frequency chart."""
        fig = self.tab_figures["Cutoff Chart"]
        fig.clear()
        ax = fig.add_subplot(111)
        plot_cutoff_chart(mode_table, operating_freq=freq, fig=fig, ax=ax)
        self.tab_canvases["Cutoff Chart"].draw()
        self.figures["Cutoff Chart"] = fig

    def _update_dispersion(self, a, b, eps_r, mu_r):
        """Update dispersion diagram."""
        fig = self.tab_figures["Dispersion"]
        fig.clear()
        ax = fig.add_subplot(111)
        plot_dispersion_diagram(a, b, eps_r, mu_r, fig=fig, ax=ax)
        self.tab_canvases["Dispersion"].draw()
        self.figures["Dispersion"] = fig

    def _update_results(self, wg, calc, m, n, freq, mode_table):
        """Update the results text panel."""
        self.results_text.config(state=tk.NORMAL)
        self.results_text.delete(1.0, tk.END)

        fc = wg.cutoff_frequency(m, n)
        fc_ghz = fc / 1e9
        freq_ghz = freq / 1e9

        lines = []
        lines.append(f"═══ TM{m}{n} Mode Analysis ═══\n")
        lines.append(f"Waveguide: {wg.a*1e3:.3f} × {wg.b*1e3:.3f} mm")
        lines.append(f"Medium: εr = {wg.eps_r}, μr = {wg.mu_r}")
        lines.append(f"f_op = {freq_ghz:.4f} GHz\n")

        lines.append(f"── Cutoff ──")
        lines.append(f"f_c = {fc_ghz:.4f} GHz")
        lines.append(f"k_c = {wg.cutoff_wavenumber(m, n):.4f} rad/m\n")

        if calc.is_propagating:
            lines.append(f"✅ Status: PROPAGATING")
            beta = wg.propagation_constant(m, n, freq)
            lines.append(f"β = {beta:.4f} rad/m")
            lg = wg.guide_wavelength(m, n, freq)
            if lg:
                lines.append(f"λ_g = {lg*1e3:.4f} mm")
            vp = wg.phase_velocity(m, n, freq)
            vg = wg.group_velocity(m, n, freq)
            if isinstance(vp, (int, float)):
                lines.append(f"v_p = {vp:.4e} m/s")
                lines.append(f"v_g = {vg:.4e} m/s")
                lines.append(f"v_p/c = {vp/3e8:.4f}")
            z_tm = wg.wave_impedance(m, n, freq)
            if isinstance(z_tm, (int, float)):
                lines.append(f"Z_TM = {z_tm:.4f} Ω")
        else:
            lines.append(f"❌ Status: EVANESCENT")
            lines.append(f"α = {calc.alpha:.4f} Np/m")

        lines.append(f"\n── Mode Table (sorted by f_c) ──")
        lines.append(f"{'Mode':<8} {'f_c (GHz)':>12} {'Status':<12}")
        lines.append(f"{'─'*8} {'─'*12} {'─'*12}")
        for mode in mode_table[:15]:
            mode_label = f"TM{mode['m']}{mode['n']}"
            status = "Propagating" if freq > mode['fc'] else "Evanescent"
            marker = "✅" if freq > mode['fc'] else "❌"
            lines.append(f"{mode_label:<8} {mode['fc_ghz']:>12.4f} {marker} {status}")

        self.results_text.insert(tk.END, '\n'.join(lines))
        self.results_text.config(state=tk.DISABLED)

    def _reset(self):
        """Reset all inputs to defaults."""
        self.var_a.set("0.02286")
        self.var_b.set("0.01016")
        self.spin_m.set(1)
        self.spin_n.set(1)
        self.var_eps_r.set("1.0")
        self.var_mu_r.set("1.0")
        self.var_freq.set("15.0")

        self.results_text.config(state=tk.NORMAL)
        self.results_text.delete(1.0, tk.END)
        self.results_text.insert(tk.END, "Click 'Analyze' to compute results.\n")
        self.results_text.config(state=tk.DISABLED)

        for name in self.tab_names:
            fig = self.tab_figures[name]
            fig.clear()
            self.tab_canvases[name].draw()

    def _export_plots(self):
        """Export all current plots to the report/figures directory."""
        if not self.figures:
            messagebox.showinfo("Export", "No plots to export. Run analysis first.")
            return

        out_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'report', 'figures'
        )
        os.makedirs(out_dir, exist_ok=True)

        try:
            m = int(self.spin_m.get())
            n = int(self.spin_n.get())
        except ValueError:
            m, n = 1, 1

        name_map = {
            "E-field (x-y)": f"E_field_xy_TM{m}{n}",
            "H-field (x-y)": f"H_field_xy_TM{m}{n}",
            "E-field (y-z)": f"E_field_yz_TM{m}{n}",
            "H-field (y-z)": f"H_field_yz_TM{m}{n}",
            "Cutoff Chart": "cutoff_chart",
            "Dispersion": "dispersion_diagram",
        }

        saved = []
        for tab_name, fig in self.figures.items():
            fname = name_map.get(tab_name, tab_name.replace(' ', '_'))
            png_path, pdf_path = save_plot(fig, fname, out_dir)
            saved.append(fname)

        messagebox.showinfo(
            "Export Complete",
            f"Saved {len(saved)} plots to:\n{out_dir}\n\nFiles: {', '.join(saved)}"
        )

    def _show_about(self):
        """Show about dialog."""
        messagebox.showinfo(
            "About",
            "Rectangular Waveguide TM Mode Analyzer\n\n"
            "Project 13 — Electromagnetic Waves\n\n"
            "Analyzes TM_mn modes in a rectangular waveguide\n"
            "with perfectly conducting walls.\n\n"
            "Features:\n"
            "• Cutoff frequency computation\n"
            "• 2D E-field & H-field visualization\n"
            "• Dispersion diagram\n"
            "• Publication-quality plots\n\n"
            "Built with Python, NumPy, Matplotlib, Tkinter"
        )


def main():
    """Launch the application."""
    app = WaveguideApp()
    app.mainloop()


if __name__ == '__main__':
    main()
