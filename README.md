# Rectangular Waveguide TM Mode Analyzer

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

> **Project 13: Electromagnetic Waves & Radiating Systems**  
> An advanced, interactive computational suite and publication-grade reporting framework for analyzing **Transverse Magnetic ($\text{TM}_{mn}$)** modes in perfectly conducting rectangular waveguides.

---

## 🌟 Key Features

- **Physics Engine (`core/`):**
  - Full Maxwell-Helmholtz modal solution for $\text{TM}_{mn}$ modes ($H_z = 0, m \ge 1, n \ge 1$).
  - Exact computation of cutoff wavenumber $k_{c,mn}$, cutoff frequency $f_{c,mn}$, propagation constant $\beta$ (real) / attenuation constant $\alpha$ (imaginary), guide wavelength $\lambda_g$, phase velocity $v_p$, group velocity $v_g$, and wave impedance $Z_{\text{TM}}$.
  - Comprehensive mode hierarchy sorting and filtering (propagating vs. evanescent).
  - Handles arbitrary waveguide dimensions ($a, b$) and lossless dielectric media ($\epsilon_r, \mu_r$).

- **Visualization Suite (`plotting/`):**
  - **Transverse ($x$-$y$) Plane:** $E_z$ filled contour overlaid with $(E_x, E_y)$ vector quiver arrows, alongside magnetic $(H_x, H_y)$ streamline loops circulating around axial displacement currents.
  - **Longitudinal ($y$-$z$) Plane:** Vector flow of $(E_y, E_z)$ along multiple guide wavelengths and out-of-plane magnetic field $H_x$ contour.
  - **Cutoff Spectrum Bar Chart:** Ranked cutoff frequencies with color-coded operational states.
  - **Dispersion Diagram ($\beta$ vs. $f$):** Modal dispersion curves compared against the light line.

- **Graphical User Interface (`gui/`):**
  - Clean Tkinter desktop application with modern ttk widgets.
  - Interactive parameter input (waveguide dimensions, mode indices, permittivity, permeability, operating frequency).
  - Multi-tab visualization notebook with integrated Matplotlib navigation toolbars (pan, zoom, home).
  - Monospace telemetry console showing real-time physical properties and mode ranking tables.
  - One-click batch export of publication-quality figures (300 DPI PNG and lossless vector PDF).

- **Command-Line Interface (`cli.py`):**
  - Fast headless calculation, parametric sweeps, mode table generation, and automated plot generation without GUI.

- **Comprehensive LaTeX Report (`report/`):**
  - Complete theoretical derivation from first-principles Maxwell's equations.
  - Modular structure (`introduction.tex`, `theory.tex`, `implementation.tex`, `results.tex`, `conclusion.tex`).
  - Standard IEEE citations and bibliography (`references.bib`).
  - Embedded high-resolution diagrams generated directly from the codebase.

---

## 📂 Repository Structure

```text
EmWavesProject/
├── main.py                     # Primary GUI application launcher
├── cli.py                      # Headless Command-Line Interface
├── generate_report_figures.py  # Batch generator for all report figures
├── requirements.txt            # Python dependencies (NumPy, Matplotlib)
├── core/                       # Physics and electrodynamics engine
│   ├── __init__.py
│   ├── waveguide.py            # RectangularWaveguide class (geometry & dispersion)
│   ├── fields.py               # TMFieldCalculator class (vectorized 2D field tensors)
│   └── cutoff.py               # Cutoff frequency calculations & mode tables
├── plotting/                   # Visualization subroutines
│   ├── __init__.py
│   ├── field_plots_xy.py       # Transverse cross-section field visualizers
│   ├── field_plots_yz.py       # Longitudinal cross-section field visualizers
│   ├── cutoff_chart.py         # Modal cutoff bar charts & dispersion curves
│   └── utils.py                # Publication styling, colormaps & boundary drawing
├── gui/                        # Desktop GUI package
│   ├── __init__.py
│   └── app.py                  # WaveguideApp Tkinter application
├── tests/                      # Validation test suite
│   └── test_waveguide.py       # 9 comprehensive unit & physics validation tests
```

---

## 🚀 Quick Start

### 1. Installation
Clone the repository and install the standard dependencies:
```bash
cd EmWavesProject
pip install -r requirements.txt
```

### 2. Launch the Graphical Interface (GUI)
Run the application launcher:
```bash
python main.py
```
1. Enter waveguide width $a$ and height $b$ (defaults to standard WR-90: $a=22.86\text{ mm}, b=10.16\text{ mm}$).
2. Select mode indices $m \ge 1$ and $n \ge 1$ (e.g., $m=1, n=1$ for the dominant $\text{TM}_{11}$ mode).
3. Specify relative permittivity $\epsilon_r$, relative permeability $\mu_r$, and excitation frequency $f$ in GHz.
4. Click **🔬 Analyze** to calculate modal dispersion metrics and update all plot tabs.
5. Click **💾 Export Plots** to save all plots to `report/figures/`.

### 3. Using the Command-Line Interface (CLI)
Compute parameters and export plots headlessly:
```bash
# Analyze TM11 mode in WR-90 at 20 GHz, print cutoff table, and save plots
python cli.py --a 0.02286 --b 0.01016 --m 1 --n 1 --freq 20.0 --table --plot --out-dir ./output_plots

# Analyze higher order mode TM21 with dielectric loading (Teflon: eps_r = 2.1)
python cli.py --a 0.02286 --b 0.01016 --m 2 --n 1 --eps_r 2.1 --freq 18.0 --table
```

### 4. Running Validation Unit Tests
Execute the test suite to verify analytical compliance and boundary conditions:
```bash
python -m unittest tests/test_waveguide.py
```

---

## 📐 Mathematical Summary for $\text{TM}_{mn}$ Modes

For Transverse Magnetic ($\text{TM}$) modes, $H_z = 0$. The longitudinal electric field is given by:
$$E_z(x, y, z) = E_0 \sin\left(\frac{m\pi x}{a}\right) \sin\left(\frac{n\pi y}{b}\right) e^{-j\beta z}, \quad (m \ge 1, n \ge 1)$$

Transverse components:
$$E_x = -\frac{j\beta}{k_c^2} \left(\frac{m\pi}{a}\right) E_0 \cos\left(\frac{m\pi x}{a}\right) \sin\left(\frac{n\pi y}{b}\right) e^{-j\beta z}$$
$$E_y = -\frac{j\beta}{k_c^2} \left(\frac{n\pi}{b}\right) E_0 \sin\left(\frac{m\pi x}{a}\right) \cos\left(\frac{n\pi y}{b}\right) e^{-j\beta z}$$
$$H_x = \frac{j\omega\epsilon}{k_c^2} \left(\frac{n\pi}{b}\right) E_0 \sin\left(\frac{m\pi x}{a}\right) \cos\left(\frac{n\pi y}{b}\right) e^{-j\beta z}$$
$$H_y = -\frac{j\omega\epsilon}{k_c^2} \left(\frac{m\pi}{a}\right) E_0 \cos\left(\frac{m\pi x}{a}\right) \sin\left(\frac{n\pi y}{b}\right) e^{-j\beta z}$$

Cutoff wavenumber and cutoff frequency:
$$k_{c,mn} = \sqrt{\left(\frac{m\pi}{a}\right)^2 + \left(\frac{n\pi}{b}\right)^2}$$
$$f_{c,mn} = \frac{c}{2\sqrt{\epsilon_r\mu_r}} \sqrt{\left(\frac{m}{a}\right)^2 + \left(\frac{n}{b}\right)^2}$$

Propagation constant and dispersion:
$$\beta_{mn} = \sqrt{k^2 - k_{c,mn}^2} = \frac{2\pi f}{c}\sqrt{\epsilon_r\mu_r} \sqrt{1 - \left(\frac{f_{c,mn}}{f}\right)^2}$$
$$v_p \cdot v_g = \frac{c^2}{\epsilon_r\mu_r}, \quad Z_{\text{TM}} = \eta \sqrt{1 - \left(\frac{f_{c,mn}}{f}\right)^2}$$

---

## 👥 Authors
Developed for **Project 13: Electromagnetic Waves Course**.
