#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations
from pathlib import Path
from typing import List, Tuple, Optional, Iterable
import re

class XYZ:
    def __init__(self, atoms: List[str], coords: List[Tuple[float, float, float]], comment: str = ""):
        self.atoms = atoms
        self.coords = coords
        self.comment = comment

    @classmethod
    def from_path(cls, path: str | Path) -> "XYZ":
        path = Path(path)
        with path.open("r") as f:
            lines = [line.rstrip("\n") for line in f]
        return cls._from_block(lines)

    @classmethod
    def _from_block(cls, lines: List[str]) -> "XYZ":
        # lines is a whole XYZ file or a single block [n, comment, coords...]
        if not lines:
            raise ValueError("Empty XYZ block")
        # If it's a whole file, keep only the first block
        n_atoms = int(lines[0].strip())
        comment = lines[1].strip() if len(lines) > 1 else ""
        atom_lines = lines[2:2 + n_atoms]
        if len(atom_lines) < n_atoms:
            raise ValueError("Incomplete XYZ atom lines")

        atoms, coords = [], []
        for ln in atom_lines:
            parts = ln.split()
            if len(parts) < 4:
                raise ValueError(f"Invalid atom line: {ln}")
            atom = parts[0]
            x, y, z = map(float, parts[1:4])
            atoms.append(atom)
            coords.append((x, y, z))
        return cls(atoms, coords, comment)

    def to_string(self) -> str:
        lines = [str(len(self.atoms)), self.comment]
        for atom, (x, y, z) in zip(self.atoms, self.coords):
            lines.append(f"{atom:2s} {x:15.8f} {y:15.8f} {z:15.8f}")
        return "\n".join(lines)

    def __repr__(self):
        return f"<XYZ {len(self.atoms)} atoms, comment='{self.comment}'>"


class XYZ_with_Energy(XYZ):
    def __init__(self, atoms, coords, comment=""):
        super().__init__(atoms, coords, comment)
        self.energy: Optional[float] = self._parse_energy(comment)

    @staticmethod
    def _parse_energy(comment: str) -> Optional[float]:
        # Match 'E <number>' where number can be float or scientific
        m = re.search(r"\bE\s+(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", comment)
        return float(m.group(1)) if m else None

    @classmethod
    def _from_block(cls, lines: List[str]) -> "XYZ_with_Energy":
        # Reuse XYZ parsing then attach energy
        base = XYZ._from_block(lines)
        return cls(base.atoms, base.coords, base.comment)

    def __repr__(self):
        return f"<XYZ_with_Energy {len(self.atoms)} atoms, E={self.energy}>"


class MultiXYZ:
    def __init__(self, xyz_list: List[XYZ]):
        self.xyz_list = xyz_list

    @staticmethod
    def _iter_blocks(all_lines: List[str]) -> Iterable[List[str]]:
        i, n = 0, len(all_lines)
        while i < n:
            # skip blank lines
            while i < n and not all_lines[i].strip():
                i += 1
            if i >= n:
                break
            try:
                n_atoms = int(all_lines[i].strip())
            except ValueError:
                raise ValueError(f"Invalid atom count at line {i+1}")
            block = all_lines[i:i + 2 + n_atoms]
            if len(block) < 2 + n_atoms:
                raise ValueError("Truncated multi-XYZ block")
            yield block
            i += 2 + n_atoms

    @classmethod
    def from_path(cls, path: str | Path, frame_cls: type[XYZ] = XYZ) -> "MultiXYZ":
        """Load a multi-XYZ file into a list of frames of type `frame_cls`."""
        path = Path(path)
        with path.open("r") as f:
            lines = [line.rstrip("\n") for line in f]
        frames: List[XYZ] = []
        for block in cls._iter_blocks(lines):
            frames.append(frame_cls._from_block(block))
        return cls(frames)

    def __len__(self):
        return len(self.xyz_list)

    def __getitem__(self, idx):
        return self.xyz_list[idx]

    def __iter__(self):
        return iter(self.xyz_list)

    def __repr__(self):
        return f"<MultiXYZ with {len(self)} frames>"


class PathMultiXYZ(MultiXYZ):
    """
    A MultiXYZ that groups frames into paths.
    If nodes_per_path is given, it partitions accordingly.
    If nodes_per_path is None, the entire set of frames is treated as one path.
    Frames are expected to be XYZ_with_Energy so energies are available.
    """
    def __init__(self, xyz_list: List[XYZ_with_Energy], nodes_per_path: Optional[int] = None, *, strict: bool = True):
        super().__init__(xyz_list)
        self.nodes_per_path = nodes_per_path
        self.strict = bool(strict)

        if nodes_per_path is None:
            # single path containing all frames
            self._paths = [self.xyz_list]
        else:
            if nodes_per_path <= 0:
                raise ValueError("nodes_per_path must be positive")
            total = len(self.xyz_list)
            if self.strict and total % nodes_per_path != 0:
                raise ValueError(
                    f"Frame count ({total}) is not divisible by nodes_per_path ({nodes_per_path}). "
                    f"Use strict=False to keep a trailing partial path."
                )
            self._paths = []
            for i in range(0, total, nodes_per_path):
                chunk = self.xyz_list[i:i + nodes_per_path]
                if not chunk:
                    continue
                if self.strict and len(chunk) != nodes_per_path:
                    continue
                self._paths.append(chunk)

        # Cache energies per path
        self._energies_per_path: List[List[Optional[float]]] = [
            [frame.energy if hasattr(frame, "energy") else None for frame in path]
            for path in self._paths
        ]

    @classmethod
    def from_path(cls, path: str | Path, nodes_per_path: Optional[int] = None, *, strict: bool = True) -> "PathMultiXYZ":
        """Load a multi-XYZ file, parse frames as XYZ_with_Energy, and group into paths."""
        multi = MultiXYZ.from_path(path, frame_cls=XYZ_with_Energy)
        return cls(xyz_list=list(multi), nodes_per_path=nodes_per_path, strict=strict)

    # Accessors
    @property
    def n_paths(self) -> int:
        return len(self._paths)

    def get_path(self, i: int) -> List[XYZ_with_Energy]:
        return self._paths[i]

    def get_path_energies(self, i: int) -> List[Optional[float]]:
        return self._energies_per_path[i]

    @property
    def paths(self) -> List[List[XYZ_with_Energy]]:
        return self._paths

    @property
    def energies(self) -> List[List[Optional[float]]]:
        return self._energies_per_path

    def __repr__(self):
        node_info = f"{self.nodes_per_path} nodes" if self.nodes_per_path else "variable nodes"
        return f"<PathMultiXYZ {self.n_paths} paths ({node_info}, total {len(self)} frames)>"


class MultiPathSet:
    """
    Represents a file containing multiple paths (e.g., many NEB paths) laid out
    one after the other. Requires `nodes_per_path` to partition frames.

    Example:
      total_frames = 868, nodes_per_path = 14  -> n_paths = 62
    """
    def __init__(self, paths: List[List[XYZ_with_Energy]], nodes_per_path: int):
        if nodes_per_path <= 0:
            raise ValueError("nodes_per_path must be a positive integer")
        if not paths:
            raise ValueError("No paths were constructed (empty input).")
        self._paths = paths
        self.nodes_per_path = int(nodes_per_path)
        # Cache energies in same shape
        self._energies = [
            [frame.energy if hasattr(frame, "energy") else None for frame in path]
            for path in self._paths
        ]

    @classmethod
    def from_path(cls, path: str | Path, nodes_per_path: int, *, strict: bool = True) -> "MultiPathSet":
        """
        Load a multi-XYZ file, split into paths of equal size `nodes_per_path`.

        - strict=True: require total_frames % nodes_per_path == 0 (no partial tail).
        - strict=False: allow a trailing partial path at the end.
        """
        multi = MultiXYZ.from_path(path, frame_cls=XYZ_with_Energy)
        frames: List[XYZ_with_Energy] = list(multi)  # type: ignore[assignment]
        total = len(frames)

        if strict and total % nodes_per_path != 0:
            raise ValueError(
                f"Frame count ({total}) is not divisible by nodes_per_path ({nodes_per_path}). "
                f"Use strict=False to keep a trailing partial path."
            )

        paths: List[List[XYZ_with_Energy]] = []
        for i in range(0, total, nodes_per_path):
            chunk = frames[i:i + nodes_per_path]
            if not chunk:
                continue
            if strict and len(chunk) != nodes_per_path:
                # Shouldn't occur due to the earlier check
                continue
            paths.append(chunk)

        return cls(paths=paths, nodes_per_path=nodes_per_path)

    @property
    def n_paths(self) -> int:
        return len(self._paths)

    @property
    def total_frames(self) -> int:
        return sum(len(p) for p in self._paths)

    @property
    def paths(self) -> List[List[XYZ_with_Energy]]:
        return self._paths

    @property
    def energies(self) -> List[List[Optional[float]]]:
        return self._energies

    def get_path(self, idx: int) -> List[XYZ_with_Energy]:
        return self._paths[idx]

    def get_path_energies(self, idx: int) -> List[Optional[float]]:
        return self._energies[idx]

    def path_min_energy(self, idx: int) -> Optional[float]:
        es = [e for e in self._energies[idx] if e is not None]
        return min(es) if es else None

    def __len__(self):
        return self.n_paths

    def __iter__(self):
        return iter(self._paths)

    def __repr__(self):
        return f"<MultiPathSet {self.n_paths} paths × {self.nodes_per_path} nodes (total {self.total_frames} frames)>"
