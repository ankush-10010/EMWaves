"""Electromagnetic field computation for TM_mn modes in rectangular waveguide.

Computes all five non-zero field components (Ez, Ex, Ey, Hx, Hy) on 2D meshgrids
for visualization in x-y (transverse) and y-z (longitudinal) cross-sections.

TM Mode Phasor Fields (Hz = 0):
    E_z = E0 sin(kx*x) sin(ky*y) e^{-j*beta*z}
    E_x = (-j*beta*kx / kc^2) E0 cos(kx*x) sin(ky*y) e^{-j*beta*z}
    E_y = (-j*beta*ky / kc^2) E0 sin(kx*x) cos(ky*y) e^{-j*beta*z}
    H_x = (j*omega*eps*ky / kc^2) E0 sin(kx*x) cos(ky*y) e^{-j*beta*z}
    H_y = (-j*omega*eps*kx / kc^2) E0 cos(kx*x) sin(ky*y) e^{-j*beta*z}

where kx = m*pi/a, ky = n*pi/b, kc^2 = kx^2 + ky^2.

Instantaneous fields: Re{ phasor * e^{j*omega*t} }
"""

from __future__ import annotations

import math
from typing import Dict, Optional

import numpy as np

# Physical constants
C_0 = 299792458.0
MU_0 = 4.0 * math.pi * 1e-7
EPS_0 = 1.0 / (MU_0 * C_0**2)


class TMFieldCalculator:
    """Compute TM_mn mode field distributions on 2D grids.

    Parameters
    ----------
    a : float
        Waveguide width (x-direction) in meters.
    b : float
        Waveguide height (y-direction) in meters.
    m : int
        Mode index along x (m >= 1).
    n : int
        Mode index along y (n >= 1).
    eps_r : float
        Relative permittivity of dielectric filling.
    mu_r : float
        Relative permeability of dielectric filling.
    freq : float
        Operating frequency in Hz.
    E0 : float
        Peak amplitude of E_z (V/m). Default 1.0.
    """

    def __init__(
        self,
        a: float,
        b: float,
        m: int,
        n: int,
        eps_r: float,
        mu_r: float,
        freq: float,
        E0: float = 1.0,
    ) -> None:
        if a <= 0 or b <= 0:
            raise ValueError(f"Dimensions must be positive: a={a}, b={b}")
        if m < 1 or n < 1:
            raise ValueError(f"TM modes require m>=1, n>=1: m={m}, n={n}")
        if eps_r <= 0 or mu_r <= 0:
            raise ValueError(f"Material params must be positive: eps_r={eps_r}, mu_r={mu_r}")
        if freq <= 0:
            raise ValueError(f"Frequency must be positive: freq={freq}")

        self.a = float(a)
        self.b = float(b)
        self.m = int(m)
        self.n = int(n)
        self.eps_r = float(eps_r)
        self.mu_r = float(mu_r)
        self.freq = float(freq)
        self.E0 = float(E0)

        # Derived quantities
        self.epsilon = EPS_0 * self.eps_r
        self.mu = MU_0 * self.mu_r
        self.omega = 2.0 * math.pi * self.freq

        # Modal wavenumbers
        self.kx = self.m * math.pi / self.a
        self.ky = self.n * math.pi / self.b
        self.kc2 = self.kx**2 + self.ky**2
        self.kc = math.sqrt(self.kc2)

        # Medium wavenumber
        self.k = self.omega * math.sqrt(self.mu * self.epsilon)

        # Cutoff frequency
        self.fc = self.kc / (2.0 * math.pi * math.sqrt(self.mu * self.epsilon))

        # Propagation constant (complex-valued for generality)
        k2_minus_kc2 = self.k**2 - self.kc2
        if k2_minus_kc2 >= 0:
            self.beta = math.sqrt(k2_minus_kc2)
            self.is_propagating = True
        else:
            self.beta = 0.0
            self.alpha = math.sqrt(-k2_minus_kc2)  # attenuation constant
            self.is_propagating = False

        # Guide wavelength
        if self.is_propagating and self.beta > 0:
            self.lambda_g = 2.0 * math.pi / self.beta
        else:
            self.lambda_g = None

    def compute_xy_fields(
        self,
        nx: int = 50,
        ny: int = 50,
        z: float = 0.0,
        t: float = 0.0,
    ) -> Dict[str, np.ndarray]:
        """Compute field components on the x-y transverse cross-section.

        Parameters
        ----------
        nx, ny : int
            Number of grid points in x and y directions.
        z : float
            Fixed z-coordinate for the cross-section (meters).
        t : float
            Time instant (seconds).

        Returns
        -------
        dict
            Keys: 'Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'x', 'y' (all 2D arrays),
            plus 'x1d', 'y1d' (1D coordinate arrays).
        """
        x1d = np.linspace(0, self.a, nx)
        y1d = np.linspace(0, self.b, ny)
        X, Y = np.meshgrid(x1d, y1d)

        # Spatial mode patterns
        sin_kx_x = np.sin(self.kx * X)
        cos_kx_x = np.cos(self.kx * X)
        sin_ky_y = np.sin(self.ky * Y)
        cos_ky_y = np.cos(self.ky * Y)

        # Phasor propagation factor at (z, t)
        if self.is_propagating:
            phase = np.exp(1j * (self.omega * t - self.beta * z))
        else:
            phase = np.exp(1j * self.omega * t) * np.exp(-self.alpha * z)

        # Complex phasor fields
        Ez_phasor = self.E0 * sin_kx_x * sin_ky_y * phase
        Ex_phasor = (-1j * self.beta * self.kx / self.kc2) * self.E0 * cos_kx_x * sin_ky_y * phase
        Ey_phasor = (-1j * self.beta * self.ky / self.kc2) * self.E0 * sin_kx_x * cos_ky_y * phase
        Hx_phasor = (1j * self.omega * self.epsilon * self.ky / self.kc2) * self.E0 * sin_kx_x * cos_ky_y * phase
        Hy_phasor = (-1j * self.omega * self.epsilon * self.kx / self.kc2) * self.E0 * cos_kx_x * sin_ky_y * phase

        return {
            'Ex': np.real(Ex_phasor),
            'Ey': np.real(Ey_phasor),
            'Ez': np.real(Ez_phasor),
            'Hx': np.real(Hx_phasor),
            'Hy': np.real(Hy_phasor),
            'x': X,
            'y': Y,
            'x1d': x1d,
            'y1d': y1d,
        }

    def compute_yz_fields(
        self,
        ny: int = 40,
        nz: int = 100,
        x: Optional[float] = None,
        t: float = 0.0,
        n_wavelengths: float = 3.0,
    ) -> Dict[str, np.ndarray]:
        """Compute field components on the y-z longitudinal cross-section.

        Parameters
        ----------
        ny, nz : int
            Number of grid points in y and z directions.
        x : float or None
            Fixed x-coordinate. Defaults to a/2 for maximum field.
        t : float
            Time instant (seconds).
        n_wavelengths : float
            Number of guide wavelengths to show along z.

        Returns
        -------
        dict
            Keys: 'Ey', 'Ez', 'Hx', 'Hy', 'y', 'z' (all 2D arrays),
            plus 'y1d', 'z1d'.
        """
        if x is None:
            x = self.a / 2.0

        # Determine z range
        if self.lambda_g is not None:
            z_max = n_wavelengths * self.lambda_g
        else:
            # Evanescent: show a few decay lengths
            z_max = n_wavelengths * 2.0 * math.pi / self.kc if self.kc > 0 else self.a * 5

        y1d = np.linspace(0, self.b, ny)
        z1d = np.linspace(0, z_max, nz)
        Y, Z = np.meshgrid(y1d, z1d)

        # Evaluate mode functions at fixed x
        sin_kx_x = math.sin(self.kx * x)
        cos_kx_x = math.cos(self.kx * x)
        sin_ky_y = np.sin(self.ky * Y)
        cos_ky_y = np.cos(self.ky * Y)

        # Propagation factor
        if self.is_propagating:
            phase = np.exp(1j * (self.omega * t - self.beta * Z))
        else:
            phase = np.exp(1j * self.omega * t) * np.exp(-self.alpha * Z)

        # Phasor fields
        Ez_phasor = self.E0 * sin_kx_x * sin_ky_y * phase
        Ey_phasor = (-1j * self.beta * self.ky / self.kc2) * self.E0 * sin_kx_x * cos_ky_y * phase
        Hx_phasor = (1j * self.omega * self.epsilon * self.ky / self.kc2) * self.E0 * sin_kx_x * cos_ky_y * phase
        Hy_phasor = (-1j * self.omega * self.epsilon * self.kx / self.kc2) * self.E0 * cos_kx_x * sin_ky_y * phase

        return {
            'Ey': np.real(Ey_phasor),
            'Ez': np.real(Ez_phasor),
            'Hx': np.real(Hx_phasor),
            'Hy': np.real(Hy_phasor),
            'y': Y,
            'z': Z,
            'y1d': y1d,
            'z1d': z1d,
        }

    def compute_xz_fields(
        self,
        nx: int = 40,
        nz: int = 100,
        y: Optional[float] = None,
        t: float = 0.0,
        n_wavelengths: float = 3.0,
    ) -> Dict[str, np.ndarray]:
        """Compute field components on the x-z longitudinal cross-section.

        Parameters
        ----------
        nx, nz : int
            Number of grid points in x and z directions.
        y : float or None
            Fixed y-coordinate. Defaults to b/2 for maximum field.
        t : float
            Time instant (seconds).
        n_wavelengths : float
            Number of guide wavelengths to show along z.

        Returns
        -------
        dict
            Keys: 'Ex', 'Ez', 'Hx', 'Hy', 'x', 'z' (all 2D arrays),
            plus 'x1d', 'z1d'.
        """
        if y is None:
            y = self.b / 2.0

        if self.lambda_g is not None:
            z_max = n_wavelengths * self.lambda_g
        else:
            z_max = n_wavelengths * 2.0 * math.pi / self.kc if self.kc > 0 else self.a * 5

        x1d = np.linspace(0, self.a, nx)
        z1d = np.linspace(0, z_max, nz)
        X, Z = np.meshgrid(x1d, z1d)

        sin_kx_x = np.sin(self.kx * X)
        cos_kx_x = np.cos(self.kx * X)
        sin_ky_y = math.sin(self.ky * y)
        cos_ky_y = math.cos(self.ky * y)

        if self.is_propagating:
            phase = np.exp(1j * (self.omega * t - self.beta * Z))
        else:
            phase = np.exp(1j * self.omega * t) * np.exp(-self.alpha * Z)

        Ez_phasor = self.E0 * sin_kx_x * sin_ky_y * phase
        Ex_phasor = (-1j * self.beta * self.kx / self.kc2) * self.E0 * cos_kx_x * sin_ky_y * phase
        Hx_phasor = (1j * self.omega * self.epsilon * self.ky / self.kc2) * self.E0 * sin_kx_x * cos_ky_y * phase
        Hy_phasor = (-1j * self.omega * self.epsilon * self.kx / self.kc2) * self.E0 * cos_kx_x * sin_ky_y * phase

        return {
            'Ex': np.real(Ex_phasor),
            'Ez': np.real(Ez_phasor),
            'Hx': np.real(Hx_phasor),
            'Hy': np.real(Hy_phasor),
            'x': X,
            'z': Z,
            'x1d': x1d,
            'z1d': z1d,
        }
