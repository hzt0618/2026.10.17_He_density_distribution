from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
RUN_DATE = "20260320"
OUTPUT_DIR = ROOT / "output" / RUN_DATE
RUN_SPECS = [
    {
        "name": "morning",
        "src_name": "conjugate_morning.txt",
        "filled_name": "conjugate_morning_filled.txt",
        "plot_name": "conjugate_morning_filled_lonlat.png",
        "title": "Conjugate Morning Filled Points in Lon-Lat Map",
    },
    {
        "name": "evening",
        "src_name": "conjugate_evening.txt",
        "filled_name": "conjugate_evening_filled.txt",
        "plot_name": "conjugate_evening_filled_lonlat.png",
        "title": "Conjugate Evening Filled Points in Lon-Lat Map",
    },
]


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
        if len(vals) < 10:
            continue

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


def fill_nan_runs_local(rows: list[Row]) -> int:
    filled_count = 0
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
            filled_count += 1

    return filled_count


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


def plot_lonlat(rows_before: list[Row], rows_after: list[Row], out_path: Path, title: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    valid_before = [r for r in rows_before if not r.is_nan]
    filled_after = [r for r in rows_after if r.is_nan]
    line_x, line_y = segmented_lonlat(rows_after)

    fig, ax = plt.subplots(figsize=(11, 5.8), dpi=180)
    ax.set_facecolor("#f6fbff")
    ax.grid(True, color="#c9d8e6", linewidth=0.7, alpha=0.7)

    ax.plot(
        line_x,
        line_y,
        color="#2b6cb0",
        lw=1.3,
        alpha=0.9,
        label="Filled curve",
        zorder=2,
    )
    ax.scatter(
        [r.conj_lon for r in valid_before],
        [r.conj_lat for r in valid_before],
        s=14,
        color="#7f8c8d",
        alpha=0.8,
        label="Original valid points",
        zorder=3,
    )
    if filled_after:
        ax.scatter(
            [r.conj_lon for r in filled_after],
            [r.conj_lat for r in filled_after],
            s=24,
            color="#d62828",
            label="Filled points",
            zorder=4,
        )

    ax.scatter([rows_after[0].conj_lon], [rows_after[0].conj_lat], color="#f4a261", s=42, zorder=5, label="Start")
    ax.scatter([rows_after[-1].conj_lon], [rows_after[-1].conj_lat], color="#2a9d8f", s=42, zorder=5, label="End")

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(range(-180, 181, 30))
    ax.set_yticks(range(-90, 91, 15))
    ax.set_xlabel("Conjugate longitude [deg]")
    ax.set_ylabel("Conjugate latitude [deg]")
    ax.set_title(title)
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def process_one(spec: dict[str, Path | str]) -> None:
    src = OUTPUT_DIR / str(spec["src_name"])
    filled = OUTPUT_DIR / str(spec["filled_name"])
    plot = OUTPUT_DIR / str(spec["plot_name"])
    title = str(spec["title"])

    header, rows = load_rows(src)
    rows_before = [Row(**vars(r)) for r in rows]

    filled_count = fill_nan_runs_local(rows)
    save_rows(filled, header, rows)
    plot_lonlat(rows_before, rows, plot, title)

    print(f"source: {src}")
    print(f"filled: {filled}")
    print(f"plot:   {plot}")
    print(f"filled count: {filled_count}")


def main() -> None:
    for spec in RUN_SPECS:
        process_one(spec)


if __name__ == "__main__":
    main()
