#!/usr/bin/env python
# coding: utf-8

try:
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("Please install matplotlib")

from typing import List, Optional, Tuple
import math
import os

from .xyzs import MultiPathSet, PathMultiXYZ, XYZ_with_Energy


# ---------------------------
# Helpers
# ---------------------------

_BOHR_PER_ANG = 1.8897259886

def _relative(energies: List[Optional[float]]) -> List[Optional[float]]:
    """Shift energies so the minimum non-None value is zero."""
    vals = [e for e in energies if e is not None]
    if not vals:
        return energies[:]
    baseline = min(vals)
    return [(e - baseline) if e is not None else None for e in energies]

def _to_plot_y(energies: List[Optional[float]]) -> List[float]:
    """Replace None with NaN for matplotlib to break lines."""
    return [e if e is not None else math.nan for e in energies]

def _frame_rmsd(a: XYZ_with_Energy, b: XYZ_with_Energy) -> float:
    """RMSD between two frames (no alignment; assumes same atom order)."""
    if len(a.coords) != len(b.coords):
        raise ValueError(f"Inconsistent atom counts between frames: {len(a.coords)} vs {len(b.coords)}")
    n = len(a.coords)
    acc = 0.0
    for (x1, y1, z1), (x2, y2, z2) in zip(a.coords, b.coords):
        dx = x1 - x2
        dy = y1 - y2
        dz = z1 - z2
        acc += dx*dx + dy*dy + dz*dz
    return math.sqrt(acc / n) if n > 0 else 0.0

def _cumulative_displacement(frames: List[XYZ_with_Energy]) -> Tuple[List[float], str]:
    """
    Compute cumulative displacement (RMSD) along a path:
      x[0] = 0; x[k] = sum_{j=1..k} RMSD(frame[j], frame[j-1])
    Units:
      - If VISNEB_COORD_UNIT in {'ang','angstrom'} (default), convert to Bohr.
      - If VISNEB_COORD_UNIT in {'bohr','au','a0'}, keep as Bohr.
    Returns (x_values, unit_label).
    """
    unit_env = os.environ.get("VISNEB_COORD_UNIT", "ang").strip().lower()
    assume_angstrom = unit_env in {"ang", "angstrom", "a", ""}

    xs: List[float] = [0.0]
    for k in range(1, len(frames)):
        d = _frame_rmsd(frames[k], frames[k-1])
        xs.append(xs[-1] + d)

    if assume_angstrom:
        xs = [v * _BOHR_PER_ANG for v in xs]
        unit = "Bohr"
    else:
        unit = "Bohr"  # already in Bohr

    return xs, unit

def _plot_one_path(ax, xs: List[float], ys: List[float], idx: int, n_paths: int) -> None:
    """Legacy styling: first=black, middle=gray, last=red."""
    if idx == 0:
        color = (0.0, 0.0, 0.0); lw, ms = 1.2, 6.0
    elif idx == n_paths - 1:
        color = (0.6, 0.0, 0.0); lw, ms = 1.5, 8.0
    else:
        color = (0.4, 0.4, 0.4); lw, ms = 1.0, 4.0

    ax.plot(xs, ys, '-', color=color, linewidth=lw)
    ax.plot(xs, ys, '.', color=color, markersize=ms)


# ---------------------------
# Public API (unchanged)
# ---------------------------

def plot_path(multixyz: PathMultiXYZ) -> None:
    """
    Plot all paths in a PathMultiXYZ with legacy look:
      - X: cumulative displacement along each path (RMSD-based), in Bohr
      - Y: relative energy (Ha), per-path minimum set to 0
      - First path black, intermediates gray, last red
    """
    if multixyz.n_paths == 0:
        raise ValueError("No paths to plot.")

    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=150)

    n = multixyz.n_paths
    # paths are lists of XYZ_with_Energy
    for i, path in enumerate(multixyz.paths):
        xs, xunit = _cumulative_displacement(path)
        relE = _relative([p.energy if hasattr(p, "energy") else None for p in path])
        ys = _to_plot_y(relE)
        _plot_one_path(ax, xs, ys, i, n)

    ax.set_xlabel(f"Displacement [{xunit}]", fontsize=12)
    ax.set_ylabel("Energy [Ha]", fontsize=12)
    title = f"Iter.: 0 to {n - 1}" if n > 1 else "Iter.: 0"
    ax.set_title(title)
    # No grid/legend to match legacy
    fig.tight_layout()
    plt.show()


def plot_paths(multipath: MultiPathSet) -> None:
    """
    Plot all paths in a MultiPathSet with legacy look:
      - X: cumulative displacement along each path (RMSD-based), in Bohr
      - Y: relative energy (Ha), per-path minimum set to 0
      - First path black, intermediates gray, last red
    """
    if len(multipath) == 0:
        raise ValueError("No paths to plot.")

    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=150)

    n = multipath.n_paths
    for i, path in enumerate(multipath.paths):
        xs, xunit = _cumulative_displacement(path)
        relE = _relative([p.energy if hasattr(p, "energy") else None for p in path])
        ys = _to_plot_y(relE)
        _plot_one_path(ax, xs, ys, i, n)

    ax.set_xlabel(f"Displacement [{xunit}]", fontsize=12)
    ax.set_ylabel("Energy [Ha]", fontsize=12)
    ax.set_title(f"Iter.: 0 to {n - 1}" if n > 1 else "Iter.: 0")
    fig.tight_layout()
    plt.show()
