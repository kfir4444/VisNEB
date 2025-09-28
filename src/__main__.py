#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import sys
from pathlib import Path

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="visneb",
        description="VisNEB — quick plotting utilities for ORCA NEB multi-XYZ paths."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # Common options
    def add_common(s: argparse.ArgumentParser):
        s.add_argument("--xyz", required=True, type=Path,
                       help="Path to the multi-XYZ file.")
        s.add_argument("--save", type=Path, default=None,
                       help="If given, save the figure to this file (e.g., plot.png) instead of showing it.")
        s.add_argument("--dpi", type=int, default=150,
                       help="Figure DPI (for --save). Default: 150")

    # path: PathMultiXYZ (nodes_per_path optional)
    s_path = sub.add_parser("path", help="Plot paths from a file into PathMultiXYZ (variable nodes ok).")
    add_common(s_path)
    s_path.add_argument("--nodes-per-path", type=int, default=None,
                        help="If provided, group frames into fixed-size paths. If omitted, treats all frames as one path.")
    s_path.add_argument("--loose", action="store_true",
                        help="Allow a trailing partial path when --nodes-per-path is set (strict=False).")

    # paths: MultiPathSet (requires nodes_per_path)
    s_paths = sub.add_parser("paths", help="Plot paths from a file into MultiPathSet (fixed nodes per path).")
    add_common(s_paths)
    s_paths.add_argument("--nodes-per-path", type=int, required=True,
                         help="Number of images (nodes) per path.")
    s_paths.add_argument("--loose", action="store_true",
                         help="Allow a trailing partial path (strict=False).")

    return p.parse_args()

def main() -> int:
    args = _parse_args()

    # If saving, use a non-interactive backend so plt.show() won't block.
    if args.save:
        os.environ.setdefault("MPLBACKEND", "Agg")

    # Imports AFTER backend selection:
    from .xyzs import PathMultiXYZ, MultiPathSet
    from .plots import plot_path, plot_paths
    import matplotlib.pyplot as plt

    try:
        if args.cmd == "path":
            pm = PathMultiXYZ.from_path(
                args.xyz,
                nodes_per_path=args.nodes_per_path,
                strict=not args.loose,
            )
            plot_path(pm)

        elif args.cmd == "paths":
            mps = MultiPathSet.from_path(
                args.xyz,
                nodes_per_path=args.nodes_per_path,
                strict=not args.loose,
            )
            plot_paths(mps)

        else:
            raise ValueError(f"Unknown command: {args.cmd}")

        if args.save:
            # Save the most recently created figure
            out = Path(args.save)
            out.parent.mkdir(parents=True, exist_ok=True)
            plt.gcf().savefig(out, dpi=args.dpi, bbox_inches="tight")
        return 0

    except Exception as e:
        print(f"[VisNEB] Error: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    raise SystemExit(main())
