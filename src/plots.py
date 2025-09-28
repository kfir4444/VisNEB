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
# Units / config
# ---------------------------
_BOHR_PER_ANG = 1.8897259886

def _assume_angstroms() -> bool:
    unit_env = os.environ.get("VISNEB_COORD_UNIT", "ang").strip().lower()
    return unit_env in {"", "ang", "angstrom", "a"}

def _smooth_enabled() -> bool:
    v = os.environ.get("VISNEB_SMOOTH", "1").strip()
    return v not in {"0", "false", "False", "no", "No"}

# ---------------------------
# Geometry helpers
# ---------------------------
def _frame_rmsd(a: XYZ_with_Energy, b: XYZ_with_Energy) -> float:
    """RMSD between two frames (no alignment; assumes same atom order)."""
    if len(a.coords) != len(b.coords):
        raise ValueError(f"Inconsistent atom counts between frames: {len(a.coords)} vs {len(b.coords)}")
    n = len(a.coords)
    if n == 0:
        return 0.0
    acc = 0.0
    for (x1, y1, z1), (x2, y2, z2) in zip(a.coords, b.coords):
        dx = x1 - x2; dy = y1 - y2; dz = z1 - z2
        acc += dx*dx + dy*dy + dz*dz
    return math.sqrt(acc / n)

def _cumulative_displacement(frames: List[XYZ_with_Energy]) -> Tuple[List[float], str]:
    xs: List[float] = [0.0]
    for k in range(1, len(frames)):
        xs.append(xs[-1] + _frame_rmsd(frames[k], frames[k-1]))
    if _assume_angstroms():
        xs = [v * _BOHR_PER_ANG for v in xs]
    return xs, "Bohr"

# ---------------------------
# Energy helpers
# ---------------------------
def _relative(energies: List[Optional[float]]) -> List[Optional[float]]:
    vals = [e for e in energies if e is not None]
    if not vals:
        return energies[:]
    base = min(vals)
    return [(e - base) if e is not None else None for e in energies]

def _nan_to_none(xs: List[float], ys: List[Optional[float]]) -> Tuple[List[float], List[float]]:
    x2: List[float] = []
    y2: List[float] = []
    for x, y in zip(xs, ys):
        if y is None or (isinstance(y, float) and (math.isnan(y) or math.isinf(y))):
            continue
        x2.append(x); y2.append(float(y))
    return x2, y2

def _unique_increasing(xs: List[float], ys: List[float]) -> Tuple[List[float], List[float]]:
    """Drop duplicate x to ensure strictly increasing for interpolation."""
    if not xs:
        return xs, ys
    xu = [xs[0]]; yu = [ys[0]]
    for i in range(1, len(xs)):
        if xs[i] > xu[-1]:
            xu.append(xs[i]); yu.append(ys[i])
        elif xs[i] < xu[-1]:
            # If out of order, skip (shouldn't happen with cumulative RMSD)
            continue
        else:
            # equal x: keep the latest y (or average)
            yu[-1] = ys[i]
    return xu, yu

# ---------------------------
# Shape-preserving cubic Hermite (PCHIP-like) smoothing (no SciPy)
# ---------------------------
def _pchip_slopes(x: List[float], y: List[float]) -> List[float]:
    n = len(x)
    if n < 2:
        return [0.0] * n
    h = [x[i+1] - x[i] for i in range(n-1)]
    delta = [(y[i+1] - y[i]) / h[i] for i in range(n-1)]

    m = [0.0] * n
    # Endpoints (Fritsch-Butland style)
    if n == 2:
        m[0] = delta[0]; m[1] = delta[0]
        return m

    # m0
    m0 = ((2*h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1]) if (h[0] + h[1]) != 0 else 0.0
    if (m0 * delta[0]) <= 0:
        m0 = 0.0
    elif abs(m0) > 3 * abs(delta[0]):
        m0 = 3 * delta[0]
    m[0] = m0

    # interior
    for i in range(1, n-1):
        if delta[i-1] == 0.0 or delta[i] == 0.0 or (delta[i-1] * delta[i]) < 0.0:
            m[i] = 0.0
        else:
            w1 = 2*h[i] + h[i-1]
            w2 = h[i] + 2*h[i-1]
            m[i] = (w1 + w2) / (w1/delta[i-1] + w2/delta[i])

    # mn-1
    mn = ((2*h[-1] + h[-2]) * delta[-1] - h[-1] * delta[-2]) / (h[-2] + h[-1]) if (h[-2] + h[-1]) != 0 else 0.0
    if (mn * delta[-1]) <= 0:
        mn = 0.0
    elif abs(mn) > 3 * abs(delta[-1]):
        mn = 3 * delta[-1]
    m[-1] = mn

    return m

def _hermite_segment(x0, x1, y0, y1, m0, m1, t):
    """Cubic Hermite basis on [x0,x1] for array t in that interval."""
    h = x1 - x0
    s = (t - x0) / h
    s2 = s * s
    s3 = s2 * s
    h00 =  2*s3 - 3*s2 + 1
    h10 =      s3 - 2*s2 + s
    h01 = -2*s3 + 3*s2
    h11 =      s3 -   s2
    return h00*y0 + h10*h*m0 + h01*y1 + h11*h*m1

def _densify_curve(x: List[float], y: List[float], points_per_seg: int = 20) -> Tuple[List[float], List[float]]:
    """
    Build a smooth (PCHIP-like) curve with ~points_per_seg samples per segment.
    Falls back to straight-line if too few points.
    """
    n = len(x)
    if n <= 2 or not _smooth_enabled():
        # Linear interpolation with denser sampling
        xnew: List[float] = []
        ynew: List[float] = []
        for i in range(n-1):
            x0, x1 = x[i], x[i+1]
            y0, y1 = y[i], y[i+1]
            steps = max(2, points_per_seg)
            for k in range(steps):
                t = k / (steps - 1)
                xt = x0 + t * (x1 - x0)
                yt = y0 + t * (y1 - y0)
                xnew.append(xt); ynew.append(yt)
        return xnew, ynew

    m = _pchip_slopes(x, y)
    xout: List[float] = []
    yout: List[float] = []
    for i in range(n-1):
        x0, x1 = x[i], x[i+1]
        y0, y1 = y[i], y[i+1]
        steps = max(2, points_per_seg)
        # sample including endpoint to ensure continuity
        for k in range(steps):
            t = x0 + (x1 - x0) * (k / (steps - 1))
            xout.append(t)
            yout.append(_hermite_segment(x0, x1, y0, y1, m[i], m[i+1], t))
    return xout, yout

# ---------------------------
# Styling
# ---------------------------
def _plot_one_path(ax, x_nodes: List[float], y_nodes: List[float], idx: int, n_paths: int) -> None:
    """Legacy styling: first=black, middle=gray, last=red."""
    if idx == 0:
        color = (0.0, 0.0, 0.0); lw, ms = 1.2, 6.0
        label = "First Iteration"
    elif idx == n_paths - 1:
        color = (0.6, 0.0, 0.0); lw, ms = 1.5, 8.0
        label = "Converged"
    else:
        color = (0.4, 0.4, 0.4); lw, ms = 1.0, 4.0
        label = ""

    # Smooth densified line (legacy "spline" vibe)
    xs_s, ys_s = _densify_curve(x_nodes, y_nodes, points_per_seg=30)
    ax.plot(xs_s, ys_s, '-', color=color, linewidth=lw, label=label)

    # Original image points
    ax.plot(x_nodes, y_nodes, '.', color=color, markersize=ms)

# ---------------------------
# Public API (unchanged signatures)
# ---------------------------
def plot_path(multixyz: PathMultiXYZ) -> None:
    if multixyz.n_paths == 0:
        raise ValueError("No paths to plot.")

    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=150)

    n = multixyz.n_paths
    for i, path in enumerate(multixyz.paths):
        xs, xunit = _cumulative_displacement(path)
        relE = _relative([p.energy if hasattr(p, "energy") else None for p in path])
        xs2, ys2 = _nan_to_none(xs, relE)
        xs2, ys2 = _unique_increasing(xs2, ys2)
        if len(xs2) < 2:
            continue
        _plot_one_path(ax, xs2, ys2, i, n)

    ax.legend(loc="best")
    ax.set_xlabel(f"Displacement [{xunit}]", fontsize=12)
    ax.set_ylabel("Energy [Ha]", fontsize=12)
    ax.set_title(f"Iter.: 0 to {n - 1}" if n > 1 else "Iter.: 0")
    fig.tight_layout()

def plot_paths(multipath: MultiPathSet) -> None:
    if len(multipath) == 0:
        raise ValueError("No paths to plot.")

    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=150)

    n = multipath.n_paths
    for i, path in enumerate(multipath.paths):
        xs, xunit = _cumulative_displacement(path)
        relE = _relative([p.energy if hasattr(p, "energy") else None for p in path])
        xs2, ys2 = _nan_to_none(xs, relE)
        xs2, ys2 = _unique_increasing(xs2, ys2)
        if len(xs2) < 2:
            continue
        _plot_one_path(ax, xs2, ys2, i, n)
    ax.legend(loc="best")
    ax.set_xlabel(f"Displacement [{xunit}]", fontsize=12)
    ax.set_ylabel("Energy [Ha]", fontsize=12)
    ax.set_title(f"Iter.: 0 to {n - 1}" if n > 1 else "Iter.: 0")
    fig.tight_layout()
