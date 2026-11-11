"""
Rectangular Waveguide TM Mode Analyzer - Core Physics Engine.

This module provides the physics calculations for transverse magnetic (TM) modes
in a rectangular waveguide with perfectly conducting (PEC) walls filled with a
lossless, linear, homogeneous, and isotropic dielectric medium.

For TM modes, the axial magnetic field vanishes everywhere (Hz = 0), and the
axial electric field Ez satisfies the Helmholtz equation:
    (d^2/dx^2 + d^2/dy^2 + kc^2) Ez = 0

Subject to PEC boundary conditions Ez(0, y) = Ez(a, y) = Ez(x, 0) = Ez(x, b) = 0,
the modal solutions are:
    Ez(x, y, z, t) = E0 * sin(m*pi*x / a) * sin(n*pi*y / b) * exp(j*(omega*t - beta*z))
where m >= 1 and n >= 1 are positive integers.
"""

from __future__ import annotations

import cmath
import math
from typing import Any, Dict, List, Optional, Union

import numpy as np

# Physical constants from scipy.constants with fallback to manual constants
try:
    from scipy.constants import c as _c_val, epsilon_0 as _eps0_val, mu_0 as _mu0_val, pi as _pi_val

    c: float = float(_c_val)
    eps_0: float = float(_eps0_val)
    mu_0: float = float(_mu0_val)
    pi: float = float(_pi_val)
except ImportError:
    pi: float = math.pi
    c: float = 299792458.0  # Speed of light in vacuum (m/s)
    mu_0: float = 4.0 * math.pi * 1e-7  # Vacuum magnetic permeability (H/m)
    eps_0: float = 1.0 / (mu_0 * (c**2))  # Vacuum electric permittivity (F/m)

# Standard constant alias
epsilon_0: float = eps_0


class RectangularWaveguide:
    """Rectangular waveguide analyzer for Transverse Magnetic (TM) modes.

    Attributes
    ----------
    a : float
        Waveguide width along x-axis in meters (a > 0).
    b : float
        Waveguide height along y-axis in meters (b > 0).
    eps_r : float
        Relative permittivity of the filling dielectric medium (eps_r >= 1.0).
    mu_r : float
        Relative permeability of the filling dielectric medium (mu_r >= 1.0).
    epsilon : float
        Permittivity of the medium in F/m (epsilon = eps_0 * eps_r).
    mu : float
        Permeability of the medium in H/m (mu = mu_0 * mu_r).
    c_medium : float
        Speed of light in the filling medium in m/s (c_medium = 1 / sqrt(mu * epsilon)).
    eta_medium : float
        Intrinsic plane wave impedance of the medium in Ohms (sqrt(mu / epsilon)).
    """

    def __init__(
        self,
        a: float,
        b: float,
        eps_r: float = 1.0,
        mu_r: float = 1.0,
    ) -> None:
        """Initialize the rectangular waveguide dimensions and medium parameters.

        Parameters
        ----------
        a : float
            Waveguide width in meters (along x-axis, PEC walls at x=0 and x=a).
        b : float
            Waveguide height in meters (along y-axis, PEC walls at y=0 and y=b).
        eps_r : float, optional
            Relative permittivity of dielectric medium (default: 1.0 for vacuum/air).
        mu_r : float, optional
            Relative permeability of dielectric medium (default: 1.0 for non-magnetic).

        Raises
        ------
        ValueError
            If dimensions a or b are non-positive, or if eps_r or mu_r are non-positive.
        """
        if a <= 0.0 or b <= 0.0:
            raise ValueError(
                f"Waveguide dimensions 'a' and 'b' must be strictly positive. Got a={a}, b={b}."
            )
        if eps_r <= 0.0:
            raise ValueError(
                f"Relative permittivity 'eps_r' must be strictly positive. Got eps_r={eps_r}."
            )
        if mu_r <= 0.0:
            raise ValueError(
                f"Relative permeability 'mu_r' must be strictly positive. Got mu_r={mu_r}."
            )

        self.a: float = float(a)
        self.b: float = float(b)
        self.eps_r: float = float(eps_r)
        self.mu_r: float = float(mu_r)

        # Derived physical quantities
        self.epsilon: float = eps_0 * self.eps_r
        self.mu: float = mu_0 * self.mu_r
        self.c_medium: float = 1.0 / math.sqrt(self.mu * self.epsilon)
        self.eta_medium: float = math.sqrt(self.mu / self.epsilon)

    def _validate_mode(self, m: int, n: int) -> tuple[int, int]:
        """Validate and cast mode indices for TM modes.

        For TM modes in a rectangular waveguide, boundary conditions require
        m >= 1 and n >= 1 because if m=0 or n=0, sin(0) = 0 and the field
        identically vanishes everywhere.

        Parameters
        ----------
        m : int
            Mode index along x.
        n : int
            Mode index along y.

        Returns
        -------
        tuple[int, int]
            Validated integer mode indices (m, n).

        Raises
        ------
        TypeError
            If m or n are not integers.
        ValueError
            If m < 1 or n < 1.
        """
        if isinstance(m, float) and m.is_integer():
            m = int(m)
        if isinstance(n, float) and n.is_integer():
            n = int(n)

        if not isinstance(m, (int, np.integer)) or not isinstance(n, (int, np.integer)):
            raise TypeError(
                f"Mode indices (m, n) must be integers, got ({type(m).__name__}, {type(n).__name__})."
            )
        if m < 1 or n < 1:
            raise ValueError(
                f"TM modes require mode indices m >= 1 and n >= 1. Got m={m}, n={n}."
            )
        return int(m), int(n)

    def cutoff_wavenumber(self, m: int, n: int) -> float:
        """Calculate the cutoff wavenumber (kc) for TM_mn mode.

        Formula:
            kc = sqrt((m*pi/a)^2 + (n*pi/b)^2)

        Parameters
        ----------
        m : int
            Mode index along x-axis (m >= 1).
        n : int
            Mode index along y-axis (n >= 1).

        Returns
        -------
        float
            Cutoff wavenumber kc in radians per meter (rad/m).
        """
        m_val, n_val = self._validate_mode(m, n)
        return float(math.sqrt((m_val * pi / self.a) ** 2 + (n_val * pi / self.b) ** 2))

    def cutoff_frequency(self, m: int, n: int) -> float:
        """Calculate the cutoff frequency (fc) for TM_mn mode.

        Formula:
            fc = kc / (2*pi*sqrt(mu*epsilon)) = (c_medium / 2) * sqrt((m/a)^2 + (n/b)^2)

        Parameters
        ----------
        m : int
            Mode index along x-axis (m >= 1).
        n : int
            Mode index along y-axis (n >= 1).

        Returns
        -------
        float
            Cutoff frequency fc in Hertz (Hz).
        """
        kc = self.cutoff_wavenumber(m, n)
        return float(kc / (2.0 * pi * math.sqrt(self.mu * self.epsilon)))

    def wavenumber(self, freq: float) -> float:
        """Calculate the unbounded medium wavenumber (k) at frequency freq.

        Formula:
            k = omega * sqrt(mu * epsilon) = 2*pi*freq / c_medium

        Parameters
        ----------
        freq : float
            Operating frequency in Hertz (Hz).

        Returns
        -------
        float
            Unbounded medium wavenumber in radians per meter (rad/m).
        """
        if freq < 0.0:
            raise ValueError(f"Frequency must be non-negative. Got freq={freq}.")
        omega = 2.0 * pi * freq
        return float(omega * math.sqrt(self.mu * self.epsilon))

    def propagation_constant(
        self, m: int, n: int, freq: Union[float, np.ndarray]
    ) -> Union[float, complex, np.ndarray]:
        """Calculate the propagation constant (beta) for TM_mn mode at frequency freq.

        Formula:
            beta = sqrt(k^2 - kc^2) where k = omega * sqrt(mu * epsilon)

        If freq > fc (above cutoff / propagating):
            beta is real and positive (rad/m).
        If freq == fc (at cutoff):
            beta = 0.0 rad/m.
        If freq < fc (below cutoff / evanescent):
            beta is purely imaginary (1j * alpha, where alpha = sqrt(kc^2 - k^2)).

        Parameters
        ----------
        m : int
            Mode index along x-axis (m >= 1).
        n : int
            Mode index along y-axis (n >= 1).
        freq : float or np.ndarray
            Operating frequency in Hertz (Hz).

        Returns
        -------
        float or complex or np.ndarray
            Propagation constant beta in rad/m. Returns a complex number
            (purely imaginary) when below cutoff.
        """
        kc = self.cutoff_wavenumber(m, n)

        if np.isscalar(freq):
            f_val = float(freq)
            if f_val < 0.0:
                raise ValueError(f"Frequency must be non-negative. Got freq={f_val}.")
            omega = 2.0 * pi * f_val
            k = omega * math.sqrt(self.mu * self.epsilon)
            diff = k**2 - kc**2

            # Guard against floating-point roundoff near cutoff
            if abs(diff) <= 1e-12 * max(k**2, kc**2, 1.0):
                return 0.0
            elif diff > 0.0:
                return float(math.sqrt(diff))
            else:
                return complex(0.0, math.sqrt(-diff))
        else:
            freq_arr = np.asarray(freq, dtype=float)
            if np.any(freq_arr < 0.0):
                raise ValueError("All frequencies must be non-negative.")
            omega_arr = 2.0 * pi * freq_arr
            k_arr = omega_arr * math.sqrt(self.mu * self.epsilon)
            diff_arr = k_arr**2 - kc**2
            # Use scimath.sqrt for clean real/complex branching over array
            return np.lib.scimath.sqrt(diff_arr)

    def wave_impedance(
        self, m: int, n: int, freq: Union[float, np.ndarray]
    ) -> Union[float, complex, np.ndarray]:
        """Calculate the wave impedance (eta_TM) for TM_mn mode at frequency freq.

        Formula:
            eta_TM = beta / (omega * epsilon)

        For propagating modes (f > fc):
            eta_TM = eta_medium * sqrt(1 - (fc/f)^2)  (real and positive, in Ohms)
        For evanescent modes (f < fc):
            eta_TM is purely imaginary (reactive, capacitive impedance).
        At cutoff (f = fc):
            eta_TM = 0.0 Ohms.

        Parameters
        ----------
        m : int
            Mode index along x-axis (m >= 1).
        n : int
            Mode index along y-axis (n >= 1).
        freq : float or np.ndarray
            Operating frequency in Hertz (Hz).

        Returns
        -------
        float or complex or np.ndarray
            TM mode wave impedance in Ohms (Omega).
        """
        if np.isscalar(freq):
            f_val = float(freq)
            if f_val <= 0.0:
                raise ValueError(f"Frequency must be strictly positive. Got freq={f_val}.")
            omega = 2.0 * pi * f_val
            beta = self.propagation_constant(m, n, f_val)
            if beta == 0.0:
                return 0.0
            return beta / (omega * self.epsilon)
        else:
            freq_arr = np.asarray(freq, dtype=float)
            if np.any(freq_arr <= 0.0):
                raise ValueError("All frequencies must be strictly positive.")
            omega_arr = 2.0 * pi * freq_arr
            beta_arr = self.propagation_constant(m, n, freq_arr)
            return beta_arr / (omega_arr * self.epsilon)

    def guide_wavelength(
        self, m: int, n: int, freq: float
    ) -> Optional[float]:
        """Calculate the guide wavelength (lambda_g) for TM_mn mode at frequency freq.

        Formula:
            lambda_g = 2*pi / beta

        The guide wavelength is physically meaningful only for propagating waves
        (where beta is real and positive). For evanescent waves (below or at cutoff),
        this method returns None.

        Parameters
        ----------
        m : int
            Mode index along x-axis (m >= 1).
        n : int
            Mode index along y-axis (n >= 1).
        freq : float
            Operating frequency in Hertz (Hz).

        Returns
        -------
        Optional[float]
            Guide wavelength in meters (m) if propagating, or None if evanescent.
        """
        if not self.is_propagating(m, n, freq):
            return None

        beta = self.propagation_constant(m, n, freq)
        if isinstance(beta, complex) or beta <= 0.0:
            return None

        return float(2.0 * pi / beta)

    def phase_velocity(
        self, m: int, n: int, freq: float
    ) -> Union[float, complex]:
        """Calculate the phase velocity (vp) for TM_mn mode at frequency freq.

        Formula:
            vp = omega / beta

        For propagating modes:
            vp = c_medium / sqrt(1 - (fc/f)^2) > c_medium
        At cutoff:
            vp approaches infinity (returns float('inf')).
        Below cutoff (evanescent):
            returns complex (imaginary phase velocity).

        Parameters
        ----------
        m : int
            Mode index along x-axis (m >= 1).
        n : int
            Mode index along y-axis (n >= 1).
        freq : float
            Operating frequency in Hertz (Hz).

        Returns
        -------
        float or complex
            Phase velocity in meters per second (m/s).
        """
        if freq <= 0.0:
            raise ValueError(f"Frequency must be strictly positive. Got freq={freq}.")
        omega = 2.0 * pi * freq
        beta = self.propagation_constant(m, n, freq)

        if beta == 0.0:
            return float("inf")

        return omega / beta

    def group_velocity(
        self, m: int, n: int, freq: float
    ) -> Union[float, complex]:
        """Calculate the group velocity (vg) for TM_mn mode at frequency freq.

        Formula:
            vg = c_medium^2 / vp

        In the waveguide medium, phase and group velocities satisfy:
            vp * vg = c_medium^2

        For propagating modes:
            vg = c_medium * sqrt(1 - (fc/f)^2) < c_medium
        At cutoff:
            vg = 0.0 m/s.
        Below cutoff (evanescent):
            returns complex (imaginary group velocity).

        Parameters
        ----------
        m : int
            Mode index along x-axis (m >= 1).
        n : int
            Mode index along y-axis (n >= 1).
        freq : float
            Operating frequency in Hertz (Hz).

        Returns
        -------
        float or complex
            Group velocity in meters per second (m/s).
        """
        vp = self.phase_velocity(m, n, freq)
        if math.isinf(vp.real if isinstance(vp, complex) else vp):
            return 0.0
        return (self.c_medium**2) / vp

    def is_propagating(self, m: int, n: int, freq: float) -> bool:
        """Determine whether TM_mn mode is propagating at frequency freq.

        A mode propagates if and only if the operating frequency strictly
        exceeds the mode's cutoff frequency (freq > fc).

        Parameters
        ----------
        m : int
            Mode index along x-axis (m >= 1).
        n : int
            Mode index along y-axis (n >= 1).
        freq : float
            Operating frequency in Hertz (Hz).

        Returns
        -------
        bool
            True if propagating (freq > fc), False if evanescent or at cutoff.
        """
        fc = self.cutoff_frequency(m, n)
        return bool(freq > fc)

    def mode_table(self, m_max: int, n_max: int) -> List[Dict[str, Any]]:
        """Generate a table of TM modes up to (m_max, n_max) sorted by cutoff frequency.

        For TM modes in rectangular waveguides, mode indices must satisfy m >= 1
        and n >= 1. The (0, 0), (m, 0), and (0, n) modes do not exist for TM
        polarization and are therefore omitted.

        Parameters
        ----------
        m_max : int
            Maximum mode index along x-axis (m_max >= 1).
        n_max : int
            Maximum mode index along y-axis (n_max >= 1).

        Returns
        -------
        List[Dict[str, Any]]
            List of dictionaries containing mode information sorted in ascending
            order of cutoff frequency. Each dictionary contains:
            - 'm': mode index along x (int)
            - 'n': mode index along y (int)
            - 'mode': compact mode string, e.g. 'TM11' (str)
            - 'mode_name': subscript mode label, e.g. 'TM_11' (str)
            - 'cutoff_frequency': cutoff frequency in Hz (float)
            - 'fc': alias for cutoff frequency in Hz (float)
            - 'fc_ghz': cutoff frequency in GHz (float)
            - 'cutoff_wavenumber': cutoff wavenumber kc in rad/m (float)
            - 'kc': alias for cutoff wavenumber in rad/m (float)

        Raises
        ------
        ValueError
            If m_max < 1 or n_max < 1.
        """
        if m_max < 1 or n_max < 1:
            raise ValueError(
                f"Mode bounds must satisfy m_max >= 1 and n_max >= 1. Got m_max={m_max}, n_max={n_max}."
            )

        modes: List[Dict[str, Any]] = []
        for m in range(1, m_max + 1):
            for n in range(1, n_max + 1):
                kc = self.cutoff_wavenumber(m, n)
                fc = self.cutoff_frequency(m, n)
                modes.append(
                    {
                        "m": m,
                        "n": n,
                        "mode": f"TM{m}{n}",
                        "mode_name": f"TM_{m}{n}",
                        "cutoff_frequency": fc,
                        "fc": fc,
                        "fc_ghz": fc / 1e9,
                        "cutoff_wavenumber": kc,
                        "kc": kc,
                    }
                )

        modes.sort(key=lambda item: item["cutoff_frequency"])
        return modes

    def __repr__(self) -> str:
        """Return a developer-friendly string representation of the waveguide."""
        return (
            f"RectangularWaveguide(a={self.a}, b={self.b}, "
            f"eps_r={self.eps_r}, mu_r={self.mu_r})"
        )

    def __str__(self) -> str:
        """Return a user-friendly summary of the waveguide."""
        return (
            f"Rectangular Waveguide: {self.a*1e3:.2f} mm x {self.b*1e3:.2f} mm, "
            f"eps_r={self.eps_r}, mu_r={self.mu_r}, c_medium={self.c_medium:.3e} m/s"
        )
