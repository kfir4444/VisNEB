#!/usr/bin/env python
# coding: utf-8

"""
VisNEB CLI (simple): run as
  python visneb.py path  --xyz path.xyz [--nodes-per-path N] [--loose] [--save out.png --dpi 200]
  python visneb.py paths --xyz path.xyz --nodes-per-path N [--loose] [--save out.png --dpi 200]
"""

import argparse
import os
import sys
from pathlib import Path

# --- Make `src` importable without installation ---
ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

try:
    from src.xyzs import PathMultiXYZ, MultiPathSet  # src/visneb/xyzs.py
    from src.plots import plot_path, plot_paths      # src/visneb/plots.py
except ModuleNotFoundError:
    raise ValueError("Couldn't find the files...")

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="visneb",
        description="VisNEB — quick plotting utilities for ORCA NEB multi-XYZ paths."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_common(sp: argparse.ArgumentParser):
        sp.add_argument("--xyz", required=True, type=Path, help="Path to the multi-XYZ file.")
        sp.add_argument("--save", type=Path, default=None, help="Save the figure to a file instead of showing.")
        sp.add_argument("--dpi", type=int, default=150, help="Figure DPI when saving (default: 150)")

    s1 = sub.add_parser("path", help="Plot using PathMultiXYZ (nodes per path optional).")
    add_common(s1)
    s1.add_argument("--nodes-per-path", type=int, default=None,
                    help="If provided, group frames into fixed-size paths; otherwise one path.")
    s1.add_argument("--loose", action="store_true",
                    help="Allow a trailing partial path when --nodes-per-path is set (strict=False).")

    s2 = sub.add_parser("paths", help="Plot using MultiPathSet (fixed nodes per path).")
    add_common(s2)
    s2.add_argument("--nodes-per-path", type=int, required=True, help="Number of images (nodes) per path.")
    s2.add_argument("--loose", action="store_true",
                    help="Allow a trailing partial path (strict=False).")

    return p.parse_args()

def main() -> int:
    args = _parse_args()

    # If saving to file, prefer a non-interactive backend so plt.show() doesn't block on headless machines.
    if args.save:
        os.environ.setdefault("MPLBACKEND", "Agg")

    import matplotlib.pyplot as plt

    try:
        if args.cmd == "path":
            pm = PathMultiXYZ.from_path(
                args.xyz, nodes_per_path=args.nodes_per_path, strict=not args.loose
            )
            plot_path(pm)

        elif args.cmd == "paths":
            mps = MultiPathSet.from_path(
                args.xyz, nodes_per_path=args.nodes_per_path, strict=not args.loose
            )
            plot_paths(mps)

        else:
            raise ValueError(f"Unknown command: {args.cmd}")

        if args.save:
            args.save.parent.mkdir(parents=True, exist_ok=True)
            plt.gcf().savefig(args.save, dpi=args.dpi, bbox_inches="tight")
        else:
            plt.show()

        return 0

    except Exception as e:
        print(f"[VisNEB] Error: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())