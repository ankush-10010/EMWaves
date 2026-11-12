"""Cutoff frequency analysis and mode tables for rectangular waveguides (TM modes).

This module provides theoretical cutoff wavenumber (kc) and cutoff frequency (fc)
calculations for transverse magnetic (TM_mn) modes in a perfectly conducting
rectangular waveguide filled with a homogeneous, lossless dielectric medium.

Governing Physics:
------------------
For TM modes, the longitudinal magnetic field Hz is zero everywhere (Hz = 0),
and the longitudinal electric field Ez satisfies the Helmholtz equation with
Dirichlet boundary conditions on the conducting walls:
    Ez = 0 at x = 0, x = a and y = 0, y = b.

The modal field distribution requires non-zero integers m >= 1 and n >= 1:
    Ez(x, y, z, t) = E_0 * sin(m*pi*x / a) * sin(n*pi*y / b) * exp(j*(omega*t - beta*z))

Cutoff Wavenumber:
    kc = sqrt((m*pi / a)^2 + (n*pi / b)^2)

Cutoff Frequency:
    fc = (1 / (2*pi*sqrt(mu_0*mu_r*eps_0*eps_r))) * sqrt((m*pi/a)^2 + (n*pi/b)^2)
       = (1 / (2*pi*sqrt(mu*eps))) * kc
       = (c / (2*sqrt(eps_r*mu_r))) * sqrt((m/a)^2 + (n/b)^2)

A mode propagates without attenuation when the operating frequency f > fc.
For TM modes, m = 1, n = 1 (TM_11) is always the dominant mode (lowest fc).
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Tuple

import numpy as np

# Physical constants in SI units
C_0: float = 299792458.0  # Speed of light in vacuum (m/s)
MU_0: float = 4.0 * math.pi * 1e-7  # Permeability of free space (H/m)
EPS_0: float = 1.0 / (MU_0 * C_0**2)  # Permittivity of free space (F/m) ~ 8.8541878176e-12

__all__ = [
    "C_0",
    "MU_0",
    "EPS_0",
    "cutoff_wavenumber",
    "cutoff_frequency",
    "compute_cutoff_table",
    "format_cutoff_table",
    "get_propagating_modes",
    "get_dominant_tm_mode",
]


def cutoff_wavenumber(m: int, n: int, a: float, b: float) -> float:
    """Compute the cutoff wavenumber kc for TM_mn mode.

    Parameters
    ----------
    m : int
        Mode index along the x-dimension (m >= 1).
    n : int
        Mode index along the y-dimension (n >= 1).
    a : float
        Waveguide width in meters (a > 0).
    b : float
        Waveguide height in meters (b > 0).

    Returns
    -------
    float
        Cutoff wavenumber kc in radians per meter (rad/m).

    Raises
    ------
    ValueError
        If waveguide dimensions a or b are non-positive, or if m or n < 1.
    """
    if a <= 0 or b <= 0:
        raise ValueError("Waveguide dimensions 'a' and 'b' must be positive numbers.")
    if m < 1 or n < 1:
        raise ValueError(f"TM mode indices must be integers >= 1 (got m={m}, n={n}).")

    kc = np.sqrt((m * np.pi / a) ** 2 + (n * np.pi / b) ** 2)
    return float(kc)


def cutoff_frequency(
    m: int,
    n: int,
    a: float,
    b: float,
    eps_r: float = 1.0,
    mu_r: float = 1.0,
) -> float:
    """Compute the cutoff frequency fc for TM_mn mode in Hertz.

    Uses the analytical formula:
        fc = (1 / (2*pi*sqrt(mu_0*mu_r*eps_0*eps_r))) * sqrt((m*pi/a)^2 + (n*pi/b)^2)

    Parameters
    ----------
    m : int
        Mode index along x (m >= 1).
    n : int
        Mode index along y (n >= 1).
    a : float
        Waveguide width in meters (a > 0).
    b : float
        Waveguide height in meters (b > 0).
    eps_r : float, optional
        Relative permittivity of dielectric (default is 1.0).
    mu_r : float, optional
        Relative permeability of medium (default is 1.0).

    Returns
    -------
    float
        Cutoff frequency fc in Hertz (Hz).

    Raises
    ------
    ValueError
        If dimensions or constitutive parameters are non-positive, or m, n < 1.
    """
    if eps_r <= 0 or mu_r <= 0:
        raise ValueError("Relative permittivity 'eps_r' and permeability 'mu_r' must be positive.")

    kc = cutoff_wavenumber(m, n, a, b)
    fc = (1.0 / (2.0 * math.pi * math.sqrt(MU_0 * mu_r * EPS_0 * eps_r))) * kc
    return float(fc)


def compute_cutoff_table(
    a: float,
    b: float,
    eps_r: float,
    mu_r: float,
    m_max: int = 5,
    n_max: int = 5,
) -> List[Dict[str, Any]]:
    """Compute cutoff frequency for all TM_mn modes where m=1..m_max, n=1..n_max.

    For TM modes, both m >= 1 and n >= 1. The cutoff frequency is calculated
    using the formula:
        fc = (1 / (2*pi*sqrt(mu_0*mu_r*eps_0*eps_r))) * sqrt((m*pi/a)^2 + (n*pi/b)^2)

    Parameters
    ----------
    a : float
        Waveguide width (x-dimension) in meters.
    b : float
        Waveguide height (y-dimension) in meters.
    eps_r : float
        Relative permittivity of filling medium (dimensionless, > 0).
    mu_r : float
        Relative permeability of filling medium (dimensionless, > 0).
    m_max : int, optional
        Maximum mode index along width (default is 5).
    n_max : int, optional
        Maximum mode index along height (default is 5).

    Returns
    -------
    list of dict
        List of mode dictionaries sorted by ascending cutoff frequency:
        [
            {
                'm': m,
                'n': n,
                'fc': freq_hz,
                'fc_ghz': freq_ghz,
                'kc': wavenumber
            },
            ...
        ]

    Raises
    ------
    ValueError
        If a <= 0, b <= 0, eps_r <= 0, mu_r <= 0, m_max < 1, or n_max < 1.
    """
    if a <= 0 or b <= 0:
        raise ValueError("Waveguide dimensions 'a' and 'b' must be positive numbers.")
    if eps_r <= 0 or mu_r <= 0:
        raise ValueError("Relative permittivity 'eps_r' and permeability 'mu_r' must be positive.")
    if m_max < 1 or n_max < 1:
        raise ValueError("Mode bounds 'm_max' and 'n_max' must be integers >= 1.")

    # Use numpy for mode index grid
    m_arr = np.arange(1, m_max + 1, dtype=int)
    n_arr = np.arange(1, n_max + 1, dtype=int)

    table: List[Dict[str, Any]] = []
    for m in m_arr:
        for n in n_arr:
            m_int = int(m)
            n_int = int(n)
            kc = cutoff_wavenumber(m_int, n_int, a, b)
            fc = cutoff_frequency(m_int, n_int, a, b, eps_r=eps_r, mu_r=mu_r)
            fc_ghz = float(fc / 1e9)
            table.append({
                "m": m_int,
                "n": n_int,
                "fc": fc,
                "fc_ghz": fc_ghz,
                "kc": kc,
            })

    # Sort in ascending order of cutoff frequency fc
    table.sort(key=lambda item: (item["fc"], item["m"], item["n"]))
    return table


def format_cutoff_table(table: List[Dict[str, Any]]) -> str:
    """Format the cutoff frequency table into a nicely formatted string table.

    Parameters
    ----------
    table : list of dict
        Output from `compute_cutoff_table`, containing dictionaries with keys:
        'm', 'n', 'fc', 'fc_ghz', 'kc'.

    Returns
    -------
    str
        Formatted string table suitable for console display or text panels.
    """
    divider = "+---------+-----+-----+----------------+---------------+--------------------+"
    header = "| Mode    |   m |   n |     kc (rad/m) |      fc (GHz) |            fc (Hz) |"

    if not table:
        empty_msg = "(No modes to display)"
        return f"{divider}\n{header}\n{divider}\n| {empty_msg:^{len(divider)-4}} |\n{divider}"

    lines = [divider, header, divider]
    for row in table:
        mode_label = f"TM_{row['m']}{row['n']}"
        m_val = row["m"]
        n_val = row["n"]
        kc_val = row["kc"]
        fc_ghz_val = row["fc_ghz"]
        fc_val = row["fc"]

        line = (
            f"| {mode_label:<7} | {m_val:>3d} | {n_val:>3d} | "
            f"{kc_val:>14.4f} | {fc_ghz_val:>13.4f} | {fc_val:>18.6e} |"
        )
        lines.append(line)

    lines.append(divider)
    return "\n".join(lines)


def get_propagating_modes(
    a: float,
    b: float,
    eps_r: float,
    mu_r: float,
    freq: float,
    m_max: int = 5,
    n_max: int = 5,
) -> List[Dict[str, Any]]:
    """Return only modes where freq > fc (propagating modes).

    Parameters
    ----------
    a : float
        Waveguide width (x-dimension) in meters.
    b : float
        Waveguide height (y-dimension) in meters.
    eps_r : float
        Relative permittivity of filling medium (> 0).
    mu_r : float
        Relative permeability of filling medium (> 0).
    freq : float
        Operating frequency in Hertz (Hz).
    m_max : int, optional
        Maximum mode index along width (default is 5).
    n_max : int, optional
        Maximum mode index along height (default is 5).

    Returns
    -------
    list of dict
        List of propagating mode dictionaries (freq > fc), sorted by
        ascending cutoff frequency.

    Raises
    ------
    ValueError
        If operating frequency is negative, or waveguide dimensions / media
        parameters are invalid.
    """
    if freq < 0:
        raise ValueError(f"Operating frequency 'freq' cannot be negative (got {freq} Hz).")

    table = compute_cutoff_table(a, b, eps_r, mu_r, m_max=m_max, n_max=n_max)
    return [mode for mode in table if freq > mode["fc"]]


def get_dominant_tm_mode(
    a: float,
    b: float,
    eps_r: float = 1.0,
    mu_r: float = 1.0,
) -> Tuple[Tuple[int, int], float]:
    """Return (1, 1) for TM modes and the corresponding cutoff frequency.

    For transverse magnetic (TM) modes in rectangular waveguides, both m >= 1
    and n >= 1 are required, making TM_11 the dominant TM mode with the lowest
    cutoff frequency.

    Parameters
    ----------
    a : float
        Waveguide width in meters (a > 0).
    b : float
        Waveguide height in meters (b > 0).
    eps_r : float, optional
        Relative permittivity of filling medium (default is 1.0).
    mu_r : float, optional
        Relative permeability of filling medium (default is 1.0).

    Returns
    -------
    tuple of ((int, int), float)
        A tuple ((1, 1), fc) where (1, 1) is the dominant TM mode index pair
        and fc is the cutoff frequency in Hertz (Hz).

    Raises
    ------
    ValueError
        If waveguide dimensions or material parameters are non-positive.
    """
    fc = cutoff_frequency(1, 1, a, b, eps_r=eps_r, mu_r=mu_r)
    return (1, 1), fc
