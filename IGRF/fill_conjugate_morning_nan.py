#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parent.parent
SRC_PATH = REPO_ROOT / "IGRF" / "output" / "conjugate_evening.txt"
OUT_PATH = REPO_ROOT / "IGRF" / "output" / "conjugate_evening_filled.txt"
PLOT_PATH = REPO_ROOT / "IGRF" / "output" / "diagnostics" / "conjugate_evening_filled_lonlat.png"


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
    is_nan: bool = False


def parse_float(text: str) -> tuple[float, bool]:
    if text.lower() == "nan":
        return float("nan"), True
    return float(text), False


def load_rows(path: Path) -> tuple[list[str], list[Row]]:
    header: list[str] = []
    rows: list[Row] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("section=") or s.startswith("columns:") or s.startswith("#"):
            header.append(s)
            continue
        vals = s.split(",")
        conj_lat, lat_nan = parse_float(vals[6])
        conj_lon, lon_nan = parse_float(vals[7])
        rows.append(
            Row(
                idx=int(vals[0]),
                input_lat=float(vals[1]),
                input_lon=float(vals[2]),
                input_x=float(vals[3]),
                input_y=float(vals[4]),
                input_z=float(vals[5]),
                conj_lat=conj_lat,
                conj_lon=conj_lon,
                conj_alt_km=float(vals[8]),
                b_ratio=float(vals[9]),
                is_nan=lat_nan or lon_nan,
            )
        )
    return header, rows


def wrap_lon_deg(lon: float) -> float:
    wrapped = ((lon + 180.0) % 360.0) - 180.0
    if wrapped == -180.0:
        return 180.0
    return wrapped


def latlon_to_xyz(lat_deg: np.ndarray, lon_deg: np.ndarray) -> np.ndarray:
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return np.column_stack((x, y, z))


def xyz_to_latlon(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    norms = np.linalg.norm(xyz, axis=1, keepdims=True)
    xyz = xyz / norms
    lat = np.rad2deg(np.arcsin(np.clip(xyz[:, 2], -1.0, 1.0)))
    lon = np.rad2deg(np.arctan2(xyz[:, 1], xyz[:, 0]))
    lon = np.array([wrap_lon_deg(v) for v in lon])
    return lat, lon


def fourier_design(theta: np.ndarray, order: int) -> np.ndarray:
    cols = [np.ones_like(theta)]
    for k in range(1, order + 1):
        cols.append(np.cos(k * theta))
        cols.append(np.sin(k * theta))
    return np.column_stack(cols)


def fit_periodic_xyz(theta_valid: np.ndarray, xyz_valid: np.ndarray, theta_all: np.ndarray, order: int = 10) -> np.ndarray:
    a_valid = fourier_design(theta_valid, order)
    a_all = fourier_design(theta_all, order)
    out = np.zeros((len(theta_all), 3))
    for i in range(3):
        coeff, *_ = np.linalg.lstsq(a_valid, xyz_valid[:, i], rcond=None)
        out[:, i] = a_all @ coeff
    return out


def slerp_vec(a: np.ndarray, b: np.ndarray, t: float) -> np.ndarray:
    dot = float(np.clip(np.dot(a, b), -1.0, 1.0))
    omega = math.acos(dot)
    if abs(omega) < 1.0e-12:
        return a.copy()
    sin_omega = math.sin(omega)
    return (
        math.sin((1.0 - t) * omega) / sin_omega * a
        + math.sin(t * omega) / sin_omega * b
    )


def fill_nan_runs_local(rows: list[Row]) -> None:
    n = len(rows)
    i = 0
    while i < n:
        if not rows[i].is_nan:
            i += 1
            continue

        start = i
        while i < n and rows[i].is_nan:
            i += 1
        end = i - 1

        left = start - 1
        right = end + 1
        if left < 0 or right >= n:
            raise ValueError("NaN run touches endpoint; local interpolation needs both anchors")

        anchor_lat = np.array([rows[left].conj_lat, rows[right].conj_lat])
        anchor_lon = np.array([rows[left].conj_lon, rows[right].conj_lon])
        anchor_xyz = latlon_to_xyz(anchor_lat, anchor_lon)
        run_len = end - start + 1

        for offset, row_idx in enumerate(range(start, end + 1), start=1):
            t = offset / (run_len + 1)
            xyz = slerp_vec(anchor_xyz[0], anchor_xyz[1], t).reshape(1, 3)
            lat_fill, lon_fill = xyz_to_latlon(xyz)
            rows[row_idx].conj_lat = float(lat_fill[0])
            rows[row_idx].conj_lon = float(lon_fill[0])


def save_rows(path: Path, header: list[str], rows: list[Row]) -> None:
    lines = list(header)
    for r in rows:
        lines.append(
            f"{r.idx:d},"
            f"{r.input_lat:.6f},{r.input_lon:.6f},"
            f"{r.input_x:.10f},{r.input_y:.10f},{r.input_z:.10f},"
            f"{r.conj_lat:.6f},{r.conj_lon:.6f},{r.conj_alt_km:.6f},{r.b_ratio:.8f}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_lonlat(rows_before: list[Row], rows_after: list[Row], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    valid_before = [r for r in rows_before if not r.is_nan]
    all_after = rows_after
    filled_after = [r for r in rows_after if r.is_nan]

    def segmented_lonlat(seq: list[Row], jump_thresh: float = 120.0) -> tuple[list[float], list[float]]:
        xs: list[float] = []
        ys: list[float] = []
        for i, r in enumerate(seq):
            if i > 0 and abs(r.conj_lon - seq[i - 1].conj_lon) > jump_thresh:
                xs.append(float("nan"))
                ys.append(float("nan"))
            xs.append(r.conj_lon)
            ys.append(r.conj_lat)
        return xs, ys

    fig, ax = plt.subplots(figsize=(8.4, 6.0), dpi=180)
    line_x, line_y = segmented_lonlat(all_after)
    ax.plot(
        line_x,
        line_y,
        color="#1f77b4",
        lw=1.4,
        alpha=0.9,
        label="Filled closed curve",
    )
    ax.scatter(
        [r.conj_lon for r in valid_before],
        [r.conj_lat for r in valid_before],
        s=12,
        color="#8c8c8c",
        alpha=0.75,
        label="Original valid points",
    )
    ax.scatter(
        [r.conj_lon for r in filled_after],
        [r.conj_lat for r in filled_after],
        s=20,
        color="#d62728",
        label="Filled points",
        zorder=5,
    )
    ax.set_xlabel("Conjugate longitude [deg]")
    ax.set_ylabel("Conjugate latitude [deg]")
    ax.set_title("Conjugate Morning Closed Curve (filled NaN at high latitudes)")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main() -> None:
    header, rows = load_rows(SRC_PATH)
    rows_before = [Row(**vars(r)) for r in rows]

    fill_nan_runs_local(rows)

    save_rows(OUT_PATH, header, rows)
    plot_lonlat(rows_before, rows, PLOT_PATH)

    filled_idxs = [r.idx for r in rows_before if r.is_nan]
    print(f"source: {SRC_PATH}")
    print(f"output: {OUT_PATH}")
    print(f"plot:   {PLOT_PATH}")
    print(f"filled count: {len(filled_idxs)}")
    ranges: list[tuple[int, int]] = []
    if filled_idxs:
        start = prev = filled_idxs[0]
        for idx in filled_idxs[1:]:
            if idx == prev + 1:
                prev = idx
            else:
                ranges.append((start, prev))
                start = prev = idx
        ranges.append((start, prev))
    print(f"filled idx ranges start-end: {ranges}")


if __name__ == "__main__":
    main()
