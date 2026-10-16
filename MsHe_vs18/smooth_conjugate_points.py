#!/usr/bin/env python3
"""
Branch-separated smoothing for conjugate-point anomalies.

This version only smooths points within the same visual branch and avoids
bridging the upper/lower curves together.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import List

import numpy as np


REPO_ROOT = Path(__file__).resolve().parent.parent
CONJ_PATH = REPO_ROOT / "IGRF" / "output" / "conjugate_morning_math.txt"
BACKUP_PATH = REPO_ROOT / "IGRF" / "output" / "conjugate_morning_math.smooth_backup.txt"


@dataclass
class Row:
    idx: int
    input_lat: float
    input_lon: float
    input_x: float
    input_y: float
    input_z: float
    conj_lat: float
    conj_lon: float
    conj_alt_km: float
    b_ratio: float


def parse_rows(path: Path) -> tuple[list[str], List[Row]]:
    header_lines: list[str] = []
    rows: List[Row] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#") or s.startswith("section=") or s.startswith("columns:"):
            header_lines.append(s)
            continue
        vals = s.split(",")
        if len(vals) < 10:
            continue
        rows.append(
            Row(
                idx=int(vals[0]),
                input_lat=float(vals[1]),
                input_lon=float(vals[2]),
                input_x=float(vals[3]),
                input_y=float(vals[4]),
                input_z=float(vals[5]),
                conj_lat=float(vals[6]),
                conj_lon=float(vals[7]),
                conj_alt_km=float(vals[8]),
                b_ratio=float(vals[9]),
            )
        )
    return header_lines, rows


def write_rows(path: Path, header_lines: list[str], rows: List[Row]) -> None:
    out: list[str] = []
    out.extend(header_lines)
    for r in rows:
        out.append(
            f"{r.idx:d},"
            f"{r.input_lat:.6f},{r.input_lon:.6f},"
            f"{r.input_x:.10f},{r.input_y:.10f},{r.input_z:.10f},"
            f"{r.conj_lat:.6f},{r.conj_lon:.6f},{r.conj_alt_km:.6f},{r.b_ratio:.8f}"
        )
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def patch_from_local_branch(
    rows: List[Row],
    idx_to_pos: dict[int, int],
    anchor_ids: list[int],
    target_ids: list[int],
    deg_lat: int = 3,
    deg_lon: int = 3,
) -> list[tuple[int, float, float, float, float]]:
    if not all(i in idx_to_pos for i in anchor_ids + target_ids):
        return []

    anchors = [rows[idx_to_pos[i]] for i in anchor_ids]
    xs = np.array([r.input_lat for r in anchors])
    lat_vals = np.array([r.conj_lat for r in anchors])
    lon_vals = np.array([r.conj_lon for r in anchors])

    lat_coef = np.polyfit(xs, lat_vals, min(deg_lat, len(xs) - 1))
    lon_coef = np.polyfit(xs, lon_vals, min(deg_lon, len(xs) - 1))

    changes: list[tuple[int, float, float, float, float]] = []
    for target_idx in target_ids:
        row = rows[idx_to_pos[target_idx]]
        new_lat = float(np.polyval(lat_coef, row.input_lat))
        new_lon = float(np.polyval(lon_coef, row.input_lat))
        changes.append((target_idx, row.conj_lat, new_lat, row.conj_lon, new_lon))
        row.conj_lat = new_lat
        row.conj_lon = new_lon
    return changes


def conservative_patch(rows: List[Row]) -> list[tuple[int, float, float, float, float]]:
    idx_to_pos = {r.idx: i for i, r in enumerate(rows)}
    changes: list[tuple[int, float, float, float, float]] = []

    # East upper branch only: smooth the local bump around idx 12..14
    changes.extend(
        patch_from_local_branch(
            rows,
            idx_to_pos,
            anchor_ids=[8, 9, 10, 11, 15, 16, 17, 18],
            target_ids=[12, 13, 14],
            deg_lat=3,
            deg_lon=3,
        )
    )

    # West upper branch only: smooth the upper-branch outliers idx 91..92
    changes.extend(
        patch_from_local_branch(
            rows,
            idx_to_pos,
            anchor_ids=[76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88],
            target_ids=[91, 92],
            deg_lat=3,
            deg_lon=3,
        )
    )

    # West lower branch only: small local correction on idx 93..94
    changes.extend(
        patch_from_local_branch(
            rows,
            idx_to_pos,
            anchor_ids=[89, 90, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104],
            target_ids=[93, 94],
            deg_lat=3,
            deg_lon=3,
        )
    )

    return changes


def restore_backup() -> None:
    if not BACKUP_PATH.is_file():
        raise FileNotFoundError(f"Backup not found: {BACKUP_PATH}")
    CONJ_PATH.write_text(BACKUP_PATH.read_text(encoding="utf-8"), encoding="utf-8")
    print(f"Restored original file from: {BACKUP_PATH}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--restore", action="store_true")
    args = parser.parse_args()

    if args.restore:
        restore_backup()
        return

    if not CONJ_PATH.is_file():
        raise FileNotFoundError(f"Conjugate file not found: {CONJ_PATH}")

    original_text = CONJ_PATH.read_text(encoding="utf-8")
    if not BACKUP_PATH.exists():
        BACKUP_PATH.write_text(original_text, encoding="utf-8")
        print(f"Created backup: {BACKUP_PATH}")
    else:
        print(f"Backup already exists: {BACKUP_PATH}")

    header_lines, rows = parse_rows(CONJ_PATH)
    changes = conservative_patch(rows)
    write_rows(CONJ_PATH, header_lines, rows)

    print(f"Patched file: {CONJ_PATH}")
    print(f"Modified points: {len(changes)}")
    for idx, old_lat, new_lat, old_lon, new_lon in changes:
        print(
            f"idx={idx:3d} "
            f"conj_lat {old_lat:8.3f}->{new_lat:8.3f} "
            f"conj_lon {old_lon:9.3f}->{new_lon:9.3f}"
        )


if __name__ == "__main__":
    main()
