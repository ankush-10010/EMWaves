# Core physics engine for rectangular waveguide TM mode analysis
from core.cutoff import (
    compute_cutoff_table,
    format_cutoff_table,
    get_dominant_tm_mode,
    get_propagating_modes,
)
from core.waveguide import RectangularWaveguide

__all__ = [
    "RectangularWaveguide",
    "compute_cutoff_table",
    "format_cutoff_table",
    "get_propagating_modes",
    "get_dominant_tm_mode",
]
