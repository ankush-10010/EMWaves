#!/usr/bin/env python3
"""
Rectangular Waveguide TM Mode Analyzer — Entry Point.

Launches the GUI application for analyzing TM_mn modes in a
rectangular waveguide with perfectly conducting walls.

Usage:
    python main.py
"""

import sys
import os

# Ensure project root is on the path
project_root = os.path.dirname(os.path.abspath(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)


def main():
    """Launch the Rectangular Waveguide TM Mode Analyzer GUI."""
    try:
        import matplotlib
        matplotlib.use('TkAgg')
    except ImportError:
        print("ERROR: matplotlib is required. Install with: pip install matplotlib")
        sys.exit(1)

    try:
        import numpy  # noqa: F401
    except ImportError:
        print("ERROR: numpy is required. Install with: pip install numpy")
        sys.exit(1)

    try:
        import tkinter  # noqa: F401
    except ImportError:
        print("ERROR: tkinter is required. Install with your system package manager.")
        sys.exit(1)

    from plotting.utils import apply_style
    apply_style()

    from gui.app import WaveguideApp
    app = WaveguideApp()
    app.mainloop()


if __name__ == '__main__':
    main()
