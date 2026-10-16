from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
DEFAULT_FILE_A = BASE_DIR / "glow model" / "output" / "Lat80.00_Lon143.16_20260119" / "019" / "Lat80.00_Lon143.16_20260119_vs276.txt"
DEFAULT_FILE_B = BASE_DIR / "glow model" / "output" / "high" / "Lat80.00_Lon143.16_20260119" / "019" / "Lat80.00_Lon143.16_20260119_vs276.txt"
DEFAULT_BLOCK_STATS = BASE_DIR / "compare_vs276_block_stats.csv"
DEFAULT_PROFILE_STATS = BASE_DIR / "compare_vs276_profile_stats.csv"


@dataclass
class GlowOutput:
    header: list[float]
    jmax: int
    nbins: int
    uflx_rows: list[list[float]]
    dflx_rows: list[list[float]]
    ener: list[float]
    delta_e: list[float]
    heion: list[float]
    ze: list[float]
    electem: list[float]


def parse_float_row(line: str) -> list[float]:
    return [float(value) for value in line.split()]


def load_glow_output(file_path: Path) -> GlowOutput:
    lines = [line.strip() for line in file_path.read_text(encoding="utf-8", errors="ignore").splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"Empty file: {file_path}")

    header = parse_float_row(lines[0])
    if len(header) < 14:
        raise ValueError(f"Unexpected header format in {file_path}")

    jmax = int(round(header[0]))
    nbins = int(round(header[1]))
    expected_lines = 1 + jmax + jmax + 5
    if len(lines) != expected_lines:
        raise ValueError(
            f"Unexpected line count in {file_path}: got {len(lines)}, expected {expected_lines}"
        )

    uflx_rows = [parse_float_row(line) for line in lines[1 : 1 + jmax]]
    dflx_rows = [parse_float_row(line) for line in lines[1 + jmax : 1 + 2 * jmax]]
    ener = parse_float_row(lines[1 + 2 * jmax])
    delta_e = parse_float_row(lines[2 + 2 * jmax])
    heion = parse_float_row(lines[3 + 2 * jmax])
    ze = parse_float_row(lines[4 + 2 * jmax])
    electem = parse_float_row(lines[5 + 2 * jmax])

    for index, row in enumerate(uflx_rows, start=1):
        if len(row) != nbins + 1:
            raise ValueError(f"Unexpected uflx column count at row {index} in {file_path}")
    for index, row in enumerate(dflx_rows, start=1):
        if len(row) != nbins + 1:
            raise ValueError(f"Unexpected dflx column count at row {index} in {file_path}")
    for name, row, expected in (
        ("ener", ener, nbins),
        ("del", delta_e, nbins),
        ("heion", heion, jmax),
        ("ze", ze, jmax),
        ("electem", electem, jmax),
    ):
        if len(row) != expected:
            raise ValueError(f"Unexpected {name} length in {file_path}: got {len(row)}, expected {expected}")

    return GlowOutput(
        header=header,
        jmax=jmax,
        nbins=nbins,
        uflx_rows=uflx_rows,
        dflx_rows=dflx_rows,
        ener=ener,
        delta_e=delta_e,
        heion=heion,
        ze=ze,
        electem=electem,
    )


def max_abs_diff(values_a: list[float], values_b: list[float]) -> float:
    return max(abs(a - b) for a, b in zip(values_a, values_b)) if values_a else 0.0


def compare_vector(values_a: list[float], values_b: list[float]) -> tuple[float, float, int]:
    total_abs = 0.0
    max_abs = 0.0
    max_index = 1
    for index, (a_val, b_val) in enumerate(zip(values_a, values_b), start=1):
        abs_diff = abs(a_val - b_val)
        total_abs += abs_diff
        if abs_diff > max_abs:
            max_abs = abs_diff
            max_index = index
    mean_abs = total_abs / len(values_a) if values_a else 0.0
    return mean_abs, max_abs, max_index


def compare_flux_block(rows_a: list[list[float]], rows_b: list[list[float]], block_name: str) -> tuple[dict[str, float | str], list[dict[str, float | str]]]:
    block_total_abs = 0.0
    block_count = 0
    block_max_abs = 0.0
    block_max_row = 1
    block_max_col = 1
    profile_stats: list[dict[str, float | str]] = []

    for row_index, (row_a, row_b) in enumerate(zip(rows_a, rows_b), start=1):
        z_a = row_a[0]
        z_b = row_b[0]
        if abs(z_a - z_b) > 1e-6:
            raise ValueError(f"{block_name} altitude mismatch at row {row_index}: {z_a} vs {z_b}")

        flux_a = row_a[1:]
        flux_b = row_b[1:]
        mean_abs, row_max_abs, row_max_col = compare_vector(flux_a, flux_b)
        profile_stats.append(
            {
                "block": block_name,
                "row_index": row_index,
                "altitude_km": z_a,
                "mean_abs_diff": mean_abs,
                "max_abs_diff": row_max_abs,
                "max_abs_col": row_max_col,
            }
        )

        block_total_abs += mean_abs * len(flux_a)
        block_count += len(flux_a)
        if row_max_abs > block_max_abs:
            block_max_abs = row_max_abs
            block_max_row = row_index
            block_max_col = row_max_col

    summary = {
        "block": block_name,
        "elements": block_count,
        "mean_abs_diff": block_total_abs / block_count if block_count else 0.0,
        "max_abs_diff": block_max_abs,
        "max_row": block_max_row,
        "max_col": block_max_col,
    }
    return summary, profile_stats


def compare_outputs(output_a: GlowOutput, output_b: GlowOutput) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]]]:
    if output_a.jmax != output_b.jmax or output_a.nbins != output_b.nbins:
        raise ValueError(
            f"Shape mismatch: A(jmax={output_a.jmax}, nbins={output_a.nbins}) vs "
            f"B(jmax={output_b.jmax}, nbins={output_b.nbins})"
        )

    block_stats: list[dict[str, float | str]] = []
    profile_stats: list[dict[str, float | str]] = []

    header_labels = [
        "jmax", "nbins", "sza", "idate", "ut", "glat", "glong",
        "f107a", "f107", "f107p", "ap", "ef", "ec", "iii",
    ]
    header_diffs = [abs(a - b) for a, b in zip(output_a.header, output_b.header)]
    max_header_index = max(range(len(header_diffs)), key=header_diffs.__getitem__)
    block_stats.append(
        {
            "block": "header",
            "elements": len(header_diffs),
            "mean_abs_diff": sum(header_diffs) / len(header_diffs),
            "max_abs_diff": header_diffs[max_header_index],
            "max_row": 1,
            "max_col": max_header_index + 1,
            "label": header_labels[max_header_index] if max_header_index < len(header_labels) else f"col_{max_header_index+1}",
        }
    )

    for block_name, rows_a, rows_b in (
        ("uflx", output_a.uflx_rows, output_b.uflx_rows),
        ("dflx", output_a.dflx_rows, output_b.dflx_rows),
    ):
        block_summary, block_profiles = compare_flux_block(rows_a, rows_b, block_name)
        block_stats.append(block_summary)
        profile_stats.extend(block_profiles)

    for block_name, values_a, values_b, axis_kind in (
        ("ener", output_a.ener, output_b.ener, "energy_bin"),
        ("del", output_a.delta_e, output_b.delta_e, "energy_bin"),
        ("heion", output_a.heion, output_b.heion, "altitude_index"),
        ("ze", output_a.ze, output_b.ze, "altitude_index"),
        ("electem", output_a.electem, output_b.electem, "altitude_index"),
    ):
        mean_abs, block_max_abs, block_max_index = compare_vector(values_a, values_b)
        block_stats.append(
            {
                "block": block_name,
                "elements": len(values_a),
                "mean_abs_diff": mean_abs,
                "max_abs_diff": block_max_abs,
                "max_row": block_max_index,
                "max_col": 1,
                "label": axis_kind,
            }
        )

    return block_stats, profile_stats


def write_csv(file_path: Path, rows: list[dict[str, float | str]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)

    with file_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def print_header_comparison(output_a: GlowOutput, output_b: GlowOutput) -> None:
    labels = [
        "jmax", "nbins", "sza", "idate", "ut", "glat", "glong",
        "f107a", "f107", "f107p", "ap", "ef", "ec", "iii",
    ]
    print("=== Header Compare ===")
    for index, (value_a, value_b) in enumerate(zip(output_a.header, output_b.header)):
        label = labels[index] if index < len(labels) else f"col_{index + 1}"
        print(f"{label:6s}: A={value_a:.6e}  B={value_b:.6e}  diff={abs(value_a - value_b):.6e}")


def print_block_summary(block_stats: list[dict[str, float | str]]) -> None:
    print("\n=== Block Summary ===")
    for row in block_stats:
        label_suffix = f", label={row['label']}" if "label" in row else ""
        print(
            f"{str(row['block']):8s} "
            f"elements={int(row['elements']) if 'elements' in row else 0:<6d} "
            f"mean_abs_diff={float(row['mean_abs_diff']):.6e} "
            f"max_abs_diff={float(row['max_abs_diff']):.6e} "
            f"max_row={int(row['max_row']) if 'max_row' in row else 0} "
            f"max_col={int(row['max_col']) if 'max_col' in row else 0}"
            f"{label_suffix}"
        )


def main() -> None:
    file_a = DEFAULT_FILE_A
    file_b = DEFAULT_FILE_B

    output_a = load_glow_output(file_a)
    output_b = load_glow_output(file_b)
    block_stats, profile_stats = compare_outputs(output_a, output_b)

    write_csv(DEFAULT_BLOCK_STATS, block_stats)
    write_csv(DEFAULT_PROFILE_STATS, profile_stats)

    print(f"File A: {file_a}")
    print(f"File B: {file_b}")
    print_header_comparison(output_a, output_b)
    print_block_summary(block_stats)
    print(f"\nBlock stats CSV   : {DEFAULT_BLOCK_STATS}")
    print(f"Profile stats CSV : {DEFAULT_PROFILE_STATS}")


if __name__ == "__main__":
    main()
