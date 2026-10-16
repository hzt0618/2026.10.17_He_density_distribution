from ZhangPaxtonModel import model_energy_from_geo
from datetime import datetime, timedelta, timezone
from pathlib import Path
import pandas as pd
import spaceweather as sw
import sys

SCRIPT_DIR = Path(__file__).resolve().parent
CODE_DIR = SCRIPT_DIR.parent.parent
IGRF_OUTPUT_DIR = CODE_DIR / "IGRF" / "output"
RUN_DATE = "2025-03-20"
RUN_DIR_NAME = RUN_DATE.replace("-", "")
IGRF_RUN_DIR = IGRF_OUTPUT_DIR / RUN_DIR_NAME

# =========================
# User-editable config
# =========================
CONFIG = {
    # Date in UTC for daily run (YYYY-MM-DD). CLI arg still overrides this.
    "date_utc": RUN_DATE,
    # Kp control
    "use_spaceweather_kp": True,
    "update_spaceweather_data": False,
    "kp_fallback": 3.0,
    # Point source
    "points_files": [
        IGRF_RUN_DIR / "conjugate_morning_filled.txt",
        IGRF_RUN_DIR / "conjugate_evening_filled.txt",
    ],
    "lat_column": "input_lat",
    "lon_column": "input_lon",
    "alt_km": 1000.0,
    "time_step_minutes": 10,
    # Output
    "output_dir": SCRIPT_DIR / "output",
    "output_prefix": "ef_ec",
    "output_ext": ".txt",
}


def infer_run_label(points_file):
    """Infer a short morning/evening style label from the points filename."""
    name = Path(points_file).stem.lower()
    if "morning" in name:
        return "morning"
    if "evening" in name:
        return "evening"
    return "points"


def main():
    # Keep this single-point demo unchanged.
    lat = 70.0
    lon = -50.0
    alt_km = 600.0
    time_utc = "2024-03-14"
    kp = 5.3

    result = model_energy_from_geo(
        lat=lat,
        lon=lon,
        alt_km=alt_km,
        time_utc=time_utc,
        kp=kp,
    )

    Q = result["energy_flux_erg_s_cm2"]
    E0 = result["mean_energy_keV"]
    ef = Q
    ec = E0 * 1000

    print("Input:")
    print(f"lat={lat}, lon={lon}, alt_km={alt_km}, time_utc={time_utc}, kp={kp}")

    print("\nOutput:")
    print(f"MLAT = {result['mlat']}")
    print(f"MLON = {result['mlon']}")
    print(f"MLT  = {result['mlt']}")
    print(f"Mean energy (keV)        = {E0}")
    print(f"Energy flux (erg/s/cm^2) = {Q}")
    print(f"ef (= Q)                 = {ef}")
    print(f"ec = E0 * 1000           = {ec}")


def _normalize_time_to_utc_naive(time_utc):
    """Convert ISO string/datetime to naive UTC timestamp for index matching."""
    if isinstance(time_utc, str):
        dt = datetime.fromisoformat(time_utc.replace("Z", "+00:00"))
    elif isinstance(time_utc, datetime):
        dt = time_utc
    else:
        raise TypeError("time_utc must be datetime or ISO-8601 string.")

    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    else:
        dt = dt.astimezone(timezone.utc)

    return pd.Timestamp(dt).tz_convert(None)


def _kp_bin_start_utc(time_utc):
    """Return the UTC 3-hour Kp bin start for a given time."""
    target_utc = _normalize_time_to_utc_naive(time_utc)
    bin_hour = (target_utc.hour // 3) * 3
    return target_utc.replace(hour=bin_hour, minute=0, second=0, microsecond=0)


def get_kp_for_ut(time_utc, update=False):
    """Get Kp for the explicit UTC 3-hour bin containing the given time."""
    df = sw.ap_kp_3h(update=update)
    bin_start_utc = _kp_bin_start_utc(time_utc)
    bin_midpoint_utc = bin_start_utc + pd.Timedelta(hours=1, minutes=30)
    valid = df[(df["Kp"] >= 0) & (df.index == bin_midpoint_utc)]
    if valid.empty:
        raise ValueError(f"No valid Kp found for 3-hour bin starting at {bin_start_utc.isoformat()}.")
    return float(valid["Kp"].iloc[0]), bin_start_utc


def load_points(points_file, lat_column="input_lat", lon_column="input_lon"):
    """Load geographic points from the conjugate file."""
    points_file = Path(points_file)
    df = pd.read_csv(points_file, skiprows=2, header=None)
    df.columns = [
        "idx",
        "input_lat",
        "input_lon",
        "input_x",
        "input_y",
        "input_z",
        "conj_lat",
        "conj_lon",
        "conj_alt_km",
        "B_ratio",
    ]

    missing = [col for col in (lat_column, lon_column) if col not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in points file: {missing}")

    points = df[["idx", lat_column, lon_column]].copy()
    points = points.rename(columns={lat_column: "lat", lon_column: "lon"})
    return points


def run_points_to_txt(time_utc, kp, points_df, alt, out_path):
    """Compute ef/ec for file-loaded geographic points at a single time, then save to txt."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    alt = float(alt)
    total = len(points_df)
    ok_count = 0
    fail_count = 0

    with out_path.open("w", encoding="utf-8") as f:
        f.write(f"# time_utc={time_utc}, kp={kp}, alt_km={alt}\n")
        f.write("idx,lat,lon,alt_km,mlat,mlon,mlt,ef,ec,status\n")

        for point in points_df.itertuples(index=False):
            try:
                result = model_energy_from_geo(
                    lat=float(point.lat),
                    lon=float(point.lon),
                    alt_km=alt,
                    time_utc=time_utc,
                    kp=kp,
                )
                q = result["energy_flux_erg_s_cm2"]
                e0 = result["mean_energy_keV"]
                ef = q
                ec = e0 * 1000
                f.write(
                    f"{int(point.idx)},{float(point.lat):.6f},{float(point.lon):.6f},{alt:.1f},"
                    f"{result['mlat']:.6f},{result['mlon']:.6f},{result['mlt']:.6f},"
                    f"{ef:.10f},{ec:.6f},ok\n"
                )
                ok_count += 1
            except Exception as exc:
                f.write(
                    f"{int(point.idx)},{float(point.lat):.6f},{float(point.lon):.6f},{alt:.1f},"
                    f"nan,nan,nan,nan,nan,fail:{type(exc).__name__}\n"
                )
                fail_count += 1

    print(f"Saved {total} rows to: {out_path}")
    print(f"Success: {ok_count}, Failed: {fail_count}")


def run_day_points_to_txt(
    date_utc,
    points_df,
    alt,
    out_path,
    use_spaceweather_kp=True,
    update_spaceweather_data=False,
    kp_fallback=3.0,
    time_step_minutes=10,
):
    """Run a full UTC-day grid for file-loaded points and save all results into one txt file."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(date_utc, str):
        day = datetime.fromisoformat(date_utc).date()
    elif isinstance(date_utc, datetime):
        day = date_utc.date()
    else:
        raise TypeError("date_utc must be datetime or ISO-8601 date string like YYYY-MM-DD.")

    alt = float(alt)
    time_step_minutes = int(time_step_minutes)
    if time_step_minutes <= 0 or 1440 % time_step_minutes != 0:
        raise ValueError("time_step_minutes must be a positive divisor of 1440.")

    steps_per_day = 1440 // time_step_minutes
    total_rows = steps_per_day * len(points_df)
    ok_count = 0
    fail_count = 0

    with out_path.open("w", encoding="utf-8") as f:
        f.write(f"# date_utc={day.isoformat()}, alt_km={alt}, points={len(points_df)}\n")
        f.write("time_utc,kp,idx,lat,lon,alt_km,mlat,mlon,mlt,ef,ec,status\n")

        day_start = datetime.combine(day, datetime.min.time())
        for step in range(steps_per_day):
            time_dt = day_start + timedelta(minutes=step * time_step_minutes)
            time_utc = time_dt.strftime("%Y-%m-%dT%H:%M:%S")
            if use_spaceweather_kp:
                try:
                    kp, _ = get_kp_for_ut(time_utc, update=update_spaceweather_data)
                except Exception:
                    kp = float(kp_fallback)
            else:
                kp = float(kp_fallback)

            for point in points_df.itertuples(index=False):
                try:
                    result = model_energy_from_geo(
                        lat=float(point.lat),
                        lon=float(point.lon),
                        alt_km=alt,
                        time_utc=time_utc,
                        kp=kp,
                    )
                    q = result["energy_flux_erg_s_cm2"]
                    e0 = result["mean_energy_keV"]
                    ef = q
                    ec = e0 * 1000
                    f.write(
                        f"{time_utc},{kp:.3f},{int(point.idx)},{float(point.lat):.6f},{float(point.lon):.6f},{alt:.1f},"
                        f"{result['mlat']:.6f},{result['mlon']:.6f},{result['mlt']:.6f},"
                        f"{ef:.10f},{ec:.6f},ok\n"
                    )
                    ok_count += 1
                except Exception as exc:
                    f.write(
                        f"{time_utc},{kp:.3f},{int(point.idx)},{float(point.lat):.6f},{float(point.lon):.6f},{alt:.1f},"
                        f"nan,nan,nan,nan,nan,fail:{type(exc).__name__}\n"
                    )
                    fail_count += 1

    print(f"Saved {total_rows} rows to: {out_path}")
    print(f"Success: {ok_count}, Failed: {fail_count}")


if __name__ == "__main__":
    # date input: first arg, e.g. python main_term.py 2025-10-09
    shared_date_utc = sys.argv[1] if len(sys.argv) > 1 else CONFIG["date_utc"]
    shared_kp_fallback = float(CONFIG["kp_fallback"])
    use_spaceweather_kp = bool(CONFIG["use_spaceweather_kp"])
    update_spaceweather_data = bool(CONFIG["update_spaceweather_data"])
    time_step_minutes = int(CONFIG["time_step_minutes"])
    alt = float(CONFIG["alt_km"])

    date_stamp = datetime.fromisoformat(shared_date_utc).strftime("%Y%m%d")
    output_dir = Path(CONFIG["output_dir"]).expanduser().resolve()

    for points_file_cfg in CONFIG["points_files"]:
        points_file = Path(points_file_cfg).expanduser().resolve()
        points_df = load_points(
            points_file=points_file,
            lat_column=str(CONFIG["lat_column"]),
            lon_column=str(CONFIG["lon_column"]),
        )

        run_label = infer_run_label(points_file)
        output_file = output_dir / f"{CONFIG['output_prefix']}_{date_stamp}_{run_label}{CONFIG['output_ext']}"
        run_day_points_to_txt(
            date_utc=shared_date_utc,
            points_df=points_df,
            alt=alt,
            out_path=output_file,
            use_spaceweather_kp=use_spaceweather_kp,
            update_spaceweather_data=update_spaceweather_data,
            kp_fallback=shared_kp_fallback,
            time_step_minutes=time_step_minutes,
        )

    print("done")
