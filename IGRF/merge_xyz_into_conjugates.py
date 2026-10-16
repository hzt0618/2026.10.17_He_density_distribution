from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    terminator, morning, evening = default_paths()
    parser = argparse.ArgumentParser(
        description="Inject x/y/z columns from a terminator file into conjugate output files."
    )
    parser.add_argument("--terminator", type=Path, default=terminator, help="Terminator source file")
    parser.add_argument("--morning", type=Path, default=morning, help="Morning conjugate file")
    parser.add_argument("--evening", type=Path, default=evening, help="Evening conjugate file")
    parser.add_argument("--no-backup", action="store_true", help="Do not create .bak files before overwriting")
    return parser.parse_args()


def default_paths() -> tuple[Path, Path, Path]:
    script_dir = Path(__file__).resolve().parent
    code_dir = script_dir.parent
    terminator = code_dir / "terminators" / "terminator_20251009_19.50_109.00.txt"
    morning = script_dir / "output" / "MAIN" / "conjugate_morning.txt"
    evening = script_dir / "output" / "MAIN" / "conjugate_evening.txt"
    return terminator, morning, evening


def is_close(a: float, b: float, tol: float = 1e-6) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tol


def parse_terminator_sections(path: Path) -> dict[str, dict[int, tuple[float, float, float, float, float]]]:
    sections: dict[str, dict[int, tuple[float, float, float, float, float]]] = {}
    current_name: str | None = None

    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("["):
            current_name = line
            sections[current_name] = {}
            continue
        if current_name is None:
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) != 6:
            continue
        idx = int(float(parts[0]))
        lon = float(parts[1])
        lat = float(parts[2])
        x = float(parts[3])
        y = float(parts[4])
        z = float(parts[5])
        sections[current_name][idx] = (lat, lon, x, y, z)

    if len(sections) < 2:
        raise ValueError(f"Expected at least two sections in {path}")
    return sections


def read_conjugate(path: Path) -> tuple[list[str], list[dict[str, float]]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    header: list[str] = []
    rows: list[dict[str, float]] = []

    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("section="):
            header.append(line)
            continue
        if line.startswith("columns:"):
            continue

        parts = [p.strip() for p in line.split(",")]
        if len(parts) not in (7, 10):
            raise ValueError(f"Unexpected column count {len(parts)} in {path}: {line}")

        if len(parts) == 7:
            idx, input_lat, input_lon, conj_lat, conj_lon, conj_alt_km, b_ratio = parts
            row = {
                "idx": int(float(idx)),
                "input_lat": float(input_lat),
                "input_lon": float(input_lon),
                "input_x": math.nan,
                "input_y": math.nan,
                "input_z": math.nan,
                "conj_lat": float(conj_lat),
                "conj_lon": float(conj_lon),
                "conj_alt_km": float(conj_alt_km),
                "b_ratio": float(b_ratio),
            }
        else:
            (
                idx,
                input_lat,
                input_lon,
                input_x,
                input_y,
                input_z,
                conj_lat,
                conj_lon,
                conj_alt_km,
                b_ratio,
            ) = parts
            row = {
                "idx": int(float(idx)),
                "input_lat": float(input_lat),
                "input_lon": float(input_lon),
                "input_x": float(input_x),
                "input_y": float(input_y),
                "input_z": float(input_z),
                "conj_lat": float(conj_lat),
                "conj_lon": float(conj_lon),
                "conj_alt_km": float(conj_alt_km),
                "b_ratio": float(b_ratio),
            }
        rows.append(row)

    if not header:
        raise ValueError(f"Missing section header in {path}")
    return header, rows


def merge_xyz(
    conj_path: Path,
    terminator_points: dict[int, tuple[float, float, float, float, float]],
    backup: bool = True,
) -> None:
    header, rows = read_conjugate(conj_path)
    output_lines = [header[0], "columns: idx,input_lat,input_lon,input_x,input_y,input_z,conj_lat,conj_lon,conj_alt_km,B_ratio"]

    for row in rows:
        idx = int(row["idx"])
        if idx not in terminator_points:
            raise KeyError(f"Point idx={idx} not found in terminator source for {conj_path}")

        src_lat, src_lon, src_x, src_y, src_z = terminator_points[idx]
        if not is_close(row["input_lat"], src_lat) or not is_close(row["input_lon"], src_lon):
            raise ValueError(
                f"Lat/lon mismatch for idx={idx} in {conj_path}: "
                f"conjugate=({row['input_lat']:.6f},{row['input_lon']:.6f}) "
                f"terminator=({src_lat:.6f},{src_lon:.6f})"
            )

        output_lines.append(
            f"{idx},"
            f"{row['input_lat']:.6f},{row['input_lon']:.6f},"
            f"{src_x:.10f},{src_y:.10f},{src_z:.10f},"
            f"{row['conj_lat']:.6f},{row['conj_lon']:.6f},"
            f"{row['conj_alt_km']:.6f},{row['b_ratio']:.8f}"
        )

    if backup:
        backup_path = conj_path.with_suffix(conj_path.suffix + ".bak")
        shutil.copy2(conj_path, backup_path)

    conj_path.write_text("\n".join(output_lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    terminator_path = args.terminator.resolve()
    morning_path = args.morning.resolve()
    evening_path = args.evening.resolve()

    sections = parse_terminator_sections(terminator_path)
    section_names = list(sections.keys())
    if len(section_names) < 2:
        raise ValueError(f"Need two terminator sections in {terminator_path}")

    merge_xyz(morning_path, sections[section_names[0]], backup=not args.no_backup)
    merge_xyz(evening_path, sections[section_names[1]], backup=not args.no_backup)

    print(f"Updated: {morning_path}")
    print(f"Updated: {evening_path}")
    print(f"Source:  {terminator_path}")


if __name__ == "__main__":
    main()
