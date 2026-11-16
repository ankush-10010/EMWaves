# Project 13: Rectangular Waveguide TM Mode Analyzer — Implementation Plan

> **Deadline:** November 17, 2026 (midnight) · **Demo:** November 18, 2026 (Wednesday afternoon)

---

## 1. Project Overview

Build a **Python GUI application** that analyzes TM<sub>mn</sub> modes in a perfectly conducting rectangular waveguide, computes cutoff frequencies, and renders publication-quality 2D field plots. Deliver alongside a **LaTeX report** with full theoretical derivation, embedded diagrams, and results.

---

## 2. Physics & Theory Summary

### 2.1 Governing Equations (TM<sub>mn</sub> Modes)

For TM modes, **H<sub>z</sub> = 0** everywhere. The z-component of the electric field satisfies:

$$E_z(x, y, z, t) = E_0 \sin\!\left(\frac{m\pi x}{a}\right) \sin\!\left(\frac{n\pi y}{b}\right) e^{j(\omega t - \beta z)}$$

where $m \geq 1$, $n \geq 1$.

### 2.2 Cutoff Frequency

$$f_{c,mn} = \frac{1}{2\pi\sqrt{\mu\epsilon}} \sqrt{\left(\frac{m\pi}{a}\right)^2 + \left(\frac{n\pi}{b}\right)^2} = \frac{c}{2\sqrt{\epsilon_r \mu_r}} \sqrt{\left(\frac{m}{a}\right)^2 + \left(\frac{n}{b}\right)^2}$$

### 2.3 Transverse Field Components

From $E_z$, all other components are derived using:

$$E_x = -\frac{j\beta}{k_c^2}\frac{\partial E_z}{\partial x}, \quad E_y = -\frac{j\beta}{k_c^2}\frac{\partial E_z}{\partial y}$$

$$H_x = \frac{j\omega\epsilon}{k_c^2}\frac{\partial E_z}{\partial y}, \quad H_y = -\frac{j\omega\epsilon}{k_c^2}\frac{\partial E_z}{\partial x}$$

where $k_c^2 = \left(\frac{m\pi}{a}\right)^2 + \left(\frac{n\pi}{b}\right)^2$ and $\beta = \sqrt{k^2 - k_c^2}$ (propagation constant).

---

## 3. Deliverables

| # | Deliverable | Format |
|---|-------------|--------|
| 1 | Python GUI Application | `.py` files (Tkinter-based) |
| 2 | Cutoff frequency table for multiple (m,n) | Displayed in GUI + Report |
| 3 | E-field plots in x-y plane | Matplotlib (quiver + contour) |
| 4 | E-field plots in y-z plane | Matplotlib (quiver + contour) |
| 5 | H-field plots in x-y plane | Matplotlib (quiver + streamplot) |
| 6 | H-field plots in y-z plane | Matplotlib (quiver + streamplot) |
| 7 | LaTeX Report (PDF) | Full theory + results + diagrams |

---

## 4. Project Structure

```
EmWavesProject/
├── main.py                     # Entry point — launches GUI
├── gui/
│   ├── __init__.py
│   ├── app.py                  # Main Tkinter application window
│   ├── input_panel.py          # Parameter input frame (a, b, m, n, εr, μr, freq)
│   ├── results_panel.py        # Cutoff frequency table display
│   └── plot_panel.py           # Tabbed notebook for embedded Matplotlib plots
├── core/
│   ├── __init__.py
│   ├── waveguide.py            # RectangularWaveguide class — all physics
│   ├── fields.py               # Field computation on grids (Ex, Ey, Ez, Hx, Hy)
│   └── cutoff.py               # Cutoff frequency calculation & mode table
├── plotting/
│   ├── __init__.py
│   ├── field_plots_xy.py       # E & H field plots in x-y cross-section
│   ├── field_plots_yz.py       # E & H field plots in y-z longitudinal section
│   ├── cutoff_chart.py         # Bar chart / table of cutoff frequencies
│   └── utils.py                # Color maps, styling, annotations
├── report/
│   ├── main.tex                # Master LaTeX document
│   ├── sections/
│   │   ├── introduction.tex
│   │   ├── theory.tex          # Full derivation from Maxwell's equations
│   │   ├── implementation.tex  # Code architecture & algorithm description
│   │   ├── results.tex         # Plots, tables, discussion
│   │   └── conclusion.tex
│   ├── figures/                # Auto-generated plots saved here
│   ├── references.bib
│   └── Makefile                # pdflatex build
├── requirements.txt            # numpy, matplotlib, tkinter (stdlib)
└── README.md
```

---

## 5. Implementation Phases

### Phase 1: Core Physics Engine (`core/`)
**Priority: CRITICAL · Estimated: 2–3 hours**

1. **`waveguide.py` — `RectangularWaveguide` class**
   - Constructor takes `a, b, eps_r, mu_r`
   - Methods: `cutoff_frequency(m, n)`, `propagation_constant(m, n, freq)`, `wave_impedance(m, n, freq)`
   - Validation: `m >= 1`, `n >= 1` for TM modes

2. **`fields.py` — Field grid computation**
   - `compute_Ez(x, y, z, t, m, n, ...)` → 2D array
   - `compute_Ex(...)`, `compute_Ey(...)` from Ez derivatives
   - `compute_Hx(...)`, `compute_Hy(...)` from Ez derivatives
   - Use `numpy` meshgrid for vectorized computation
   - For x-y plane: fix `z = 0`, compute all transverse components
   - For y-z plane: fix `x = x0` (e.g., `x = a/2`), compute relevant components

3. **`cutoff.py` — Mode table**
   - Generate cutoff frequencies for modes `(m, n)` where `m = 1..M`, `n = 1..N`
   - Sort by ascending frequency
   - Return as structured data (list of dicts or DataFrame)

### Phase 2: Plotting Module (`plotting/`)
**Priority: CRITICAL · Estimated: 3–4 hours**

4. **`field_plots_xy.py` — Cross-section plots (z = const)**
   - **E-field in x-y:** Quiver plot of `(Ex, Ey)` + filled contour of `|E_z|` as background
   - **H-field in x-y:** Streamplot of `(Hx, Hy)` with colorbar for magnitude
   - Waveguide boundary drawn as rectangle
   - Proper axis labels, title with mode `TM_{mn}`, units

5. **`field_plots_yz.py` — Longitudinal plots (x = const)**
   - **E-field in y-z:** Quiver plot of `(Ey, Ez)` showing propagation pattern
   - **H-field in y-z:** Quiver/streamplot of `(Hy, 0)` — note Hz=0 for TM
   - Show multiple wavelengths along z-axis
   - Waveguide walls at `y = 0` and `y = b`

6. **`cutoff_chart.py` — Frequency visualization**
   - Horizontal bar chart of first ~10 modes sorted by cutoff frequency
   - Operating frequency line overlay (if above cutoff → propagating)

7. **Styling (`utils.py`)**
   - Consistent `rcParams` for publication quality
   - Custom colormaps (diverging for fields, sequential for magnitudes)
   - LaTeX rendering in labels: `plt.rc('text', usetex=False)` with mathtext

### Phase 3: GUI Application (`gui/`)
**Priority: HIGH · Estimated: 3–4 hours**

8. **`app.py` — Main window**
   - Title: "Rectangular Waveguide TM Mode Analyzer"
   - Layout: Left panel (inputs) | Right panel (results + plots)
   - Menu bar: File → Export Plots, Help → About
   - Window icon and professional styling with `ttk` themes

9. **`input_panel.py` — Parameter entry**
   - Labeled entries with units:
     - `a` (m): waveguide width
     - `b` (m): waveguide height
     - `m`: mode index (integer spinner, min=1)
     - `n`: mode index (integer spinner, min=1)
     - `εr`: relative permittivity (float, default=1.0)
     - `μr`: relative permeability (float, default=1.0)
     - `f` (Hz): operating frequency (optional, for β calculation)
   - **"Analyze" button** → triggers computation
   - **"Reset" button** → clears all fields
   - Input validation with error popups

10. **`results_panel.py` — Output display**
    - Cutoff frequency for selected `(m, n)`: large, prominent display
    - Table of first N modes with cutoff frequencies
    - Propagation status: "Propagating ✓" or "Evanescent ✗" based on operating freq
    - Propagation constant β, guide wavelength λ_g

11. **`plot_panel.py` — Tabbed plot viewer**
    - `ttk.Notebook` with tabs:
      - Tab 1: **E-field (x-y plane)**
      - Tab 2: **H-field (x-y plane)**
      - Tab 3: **E-field (y-z plane)**
      - Tab 4: **H-field (y-z plane)**
      - Tab 5: **Cutoff Frequency Chart**
    - Each tab embeds a `FigureCanvasTkAgg`
    - Matplotlib navigation toolbar for zoom/pan
    - "Save Plot" button per tab (PNG, PDF, SVG)

### Phase 4: LaTeX Report (`report/`)
**Priority: HIGH · Estimated: 3–4 hours**

12. **`theory.tex`** — Full derivation
    - Start from Maxwell's equations → Helmholtz equation
    - Boundary conditions for PEC waveguide
    - Separation of variables → modal solutions
    - TM mode specific: Hz = 0, derive all 5 field components
    - Cutoff condition derivation
    - Dispersion relation, phase/group velocity

13. **`implementation.tex`** — Software description
    - Architecture diagram (draw with TikZ or embed PNG)
    - GUI screenshots
    - Code snippets for key algorithms
    - Validation methodology

14. **`results.tex`** — Analysis & plots
    - Cutoff frequency table for sample parameters (e.g., a=2.286cm, b=1.016cm — WR-90)
    - All 4 field plots for TM₁₁, TM₂₁, TM₁₂ modes
    - Discussion of field patterns and physical interpretation
    - Comparison with textbook results

15. **`introduction.tex` & `conclusion.tex`**
    - Motivation, applications (microwave engineering, radar, satellite comm)
    - Summary of results, possible extensions (TE modes, lossy dielectrics)

### Phase 5: Integration & Polish
**Priority: MEDIUM · Estimated: 1–2 hours**

16. **Auto-export plots** from GUI to `report/figures/` for LaTeX inclusion
17. **README.md** with installation instructions and usage guide
18. **Testing** with known waveguide parameters (WR-90, WR-284)
19. **Edge case handling**: evanescent modes, degenerate modes

---

## 6. GUI Mockup Layout

```
┌──────────────────────────────────────────────────────────────────┐
│  Rectangular Waveguide TM Mode Analyzer          [File] [Help]  │
├──────────────────┬───────────────────────────────────────────────┤
│  PARAMETERS      │  ┌─────┬─────┬─────┬─────┬─────┐            │
│                  │  │E x-y│H x-y│E y-z│H y-z│Chart│            │
│  Width a [m]:    │  ├─────┴─────┴─────┴─────┴─────┤            │
│  [___0.02286__]  │  │                               │            │
│                  │  │     ┌───────────────────┐     │            │
│  Height b [m]:   │  │     │                   │     │            │
│  [___0.01016__]  │  │     │   MATPLOTLIB      │     │            │
│                  │  │     │   FIGURE CANVAS    │     │            │
│  Mode m (≥1):    │  │     │                   │     │            │
│  [▲ 1 ▼]        │  │     │   (Quiver/Stream/ │     │            │
│                  │  │     │    Contour Plot)   │     │            │
│  Mode n (≥1):    │  │     │                   │     │            │
│  [▲ 1 ▼]        │  │     └───────────────────┘     │            │
│                  │  │     [Toolbar: 🔍 🏠 💾 ]       │            │
│  εr:             │  ├───────────────────────────────┤            │
│  [___1.0______]  │  │  RESULTS                      │            │
│                  │  │  ──────────────────────────    │            │
│  μr:             │  │  Cutoff freq (TM₁₁): 9.49 GHz│            │
│  [___1.0______]  │  │  Status: Propagating ✓        │            │
│                  │  │  β = 158.2 rad/m               │            │
│  Freq f [GHz]:   │  │  λ_g = 39.7 mm                │            │
│  [___10.0_____]  │  │                               │            │
│                  │  │  Mode Table:                    │            │
│  [  Analyze  ]   │  │  TM₁₁  9.49 GHz              │            │
│  [   Reset   ]   │  │  TM₂₁ 14.76 GHz              │            │
│                  │  │  TM₁₂ 16.15 GHz              │            │
│                  │  │  ...                           │            │
└──────────────────┴───────────────────────────────────────────────┘
```

---

## 7. Key Technical Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Language | **Python 3** | Allowed by project, best for scientific computing |
| GUI Framework | **Tkinter + ttk** | No external deps, professional with themes |
| Plotting | **Matplotlib** | Publication quality, embeds in Tkinter |
| Math | **NumPy** | Vectorized field computation |
| Report | **LaTeX** | Required for top-end presentation |
| Field visualization | **Quiver + Contour + Streamplot** | Shows both direction and magnitude |

---

## 8. Plots to Generate (per mode)

### x-y Plane (transverse cross-section at z = 0)

| Plot | Components | Visualization |
|------|-----------|---------------|
| E-field | E<sub>x</sub>, E<sub>y</sub> (arrows) + E<sub>z</sub> (color fill) | `quiver` + `contourf` |
| H-field | H<sub>x</sub>, H<sub>y</sub> (streamlines) | `streamplot` with color = magnitude |

### y-z Plane (longitudinal section at x = a/2)

| Plot | Components | Visualization |
|------|-----------|---------------|
| E-field | E<sub>y</sub>, E<sub>z</sub> (arrows) | `quiver` with magnitude color |
| H-field | H<sub>x</sub> (arrows, out-of-plane indicated) | `quiver` or annotated contour |

---

## 9. LaTeX Report Structure

```
1. Introduction
   1.1 Background & Motivation
   1.2 Objectives
2. Theory
   2.1 Maxwell's Equations in Source-Free Region
   2.2 Wave Equation & Helmholtz Equation
   2.3 Rectangular Waveguide Boundary Conditions
   2.4 TM Mode Solution (Hz = 0)
   2.5 Field Components Derivation
   2.6 Cutoff Frequency & Dispersion
   2.7 Phase Velocity, Group Velocity, Guide Wavelength
3. Implementation
   3.1 Software Architecture
   3.2 GUI Design
   3.3 Computational Algorithm
4. Results & Discussion
   4.1 Cutoff Frequency Analysis
   4.2 Field Distribution — TM₁₁
   4.3 Field Distribution — TM₂₁, TM₁₂
   4.4 Physical Interpretation
   4.5 Validation with Standard Results
5. Conclusion
   5.1 Summary
   5.2 Future Scope
References
Appendix A: Complete Source Code
```

---

## 10. Timeline (remaining ~25 hours)

| Time Block | Task | Hours |
|------------|------|-------|
| **Now** | Core physics engine (`core/`) | 2 |
| **+2h** | Plotting module (`plotting/`) | 3 |
| **+5h** | GUI application (`gui/`) | 3 |
| **+8h** | Integration & testing | 1 |
| **+9h** | LaTeX report writing | 4 |
| **+13h** | Report figures, polish, edge cases | 2 |
| **Buffer** | Debugging, demo preparation | ~10h |

---

## 11. Validation Strategy

- Compare cutoff frequencies with standard WR-90 waveguide data
- Verify TM₁₁ field patterns match Pozar / Balanis textbook figures
- Check boundary conditions: E<sub>z</sub> = 0 at walls, E tangential = 0 at PEC
- Verify that modes below cutoff show evanescent behavior (imaginary β)

---

## 12. "Top-End" Extras (Differentiators)

- [ ] **3D field visualization** option (using `mpl_toolkits.mplot3d`)
- [ ] **Animation** of wave propagation along z (using `FuncAnimation`)
- [ ] **Mode comparison** — overlay multiple modes on one plot
- [ ] **Dispersion diagram** (β vs ω / f vs β plot)
- [ ] **Export to PDF** button in GUI (saves all plots + summary)
- [ ] **Dark/light theme** toggle in GUI
- [ ] **Real-time update** as parameters change (with slider widgets)

> [!IMPORTANT]
> For TM modes in rectangular waveguide, both **m ≥ 1** and **n ≥ 1** are required (unlike TE modes where one can be zero). The lowest TM mode is TM₁₁. The GUI must enforce this constraint.
