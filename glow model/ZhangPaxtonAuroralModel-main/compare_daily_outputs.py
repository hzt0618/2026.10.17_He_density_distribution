from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "output"
DATE_TAGS = ["20260117", "20260118", "20260119", "20260120", "20260121"]


def load_day_table(date_tag: str) -> pd.DataFrame:
    file_path = OUTPUT_DIR / f"{date_tag}_24h.txt"
    if not file_path.is_file():
        raise FileNotFoundError(f"Missing file: {file_path}")

    df = pd.read_csv(file_path, comment="#")
    if "time_utc" not in df.columns:
        raise ValueError(f"Missing time_utc column in: {file_path}")

    df = df[df["status"] == "ok"].copy()
    df["time_utc"] = pd.to_datetime(df["time_utc"], utc=True)
    df["ut_hour"] = (
        df["time_utc"].dt.hour
        + df["time_utc"].dt.minute / 60.0
        + df["time_utc"].dt.second / 3600.0
    )
    df["date_tag"] = date_tag
    return df


def main() -> None:
    frames = [load_day_table(date_tag) for date_tag in DATE_TAGS]
    combined = pd.concat(frames, ignore_index=True).sort_values("time_utc").reset_index(drop=True)
    combined["day_index"] = (combined["time_utc"] - combined["time_utc"].min()).dt.total_seconds() / 86400.0
    alt_km = combined["alt_km"].iloc[0]

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    series_specs = [
        ("kp", "Kp"),
        ("ef", "Energy Flux (erg/s/cm^2)"),
        ("ec", "Mean Energy (eV)"),
    ]

    for ax, (column, ylabel) in zip(axes, series_specs):
        ax.plot(combined["day_index"], combined[column], marker="o", linewidth=1.5, markersize=3, color="tab:blue")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)

    day_ticks = list(range(len(DATE_TAGS)))
    for ax in axes:
        for x in day_ticks:
            ax.axvline(x=x, color="0.7", linewidth=0.8, linestyle="--")

    axes[-1].set_xlabel("UT Hour")
    axes[-1].set_xlim(combined["day_index"].min(), combined["day_index"].max())
    axes[-1].set_xticks(day_ticks)
    axes[-1].set_xticklabels(DATE_TAGS)
    fig.suptitle(
        f"Time-stitched series for 20260117-20260121 at lon=109, lat=80, alt={alt_km:g} km"
    )
    fig.tight_layout()

    out_png = OUTPUT_DIR / "compare_20260117_20260121.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Saved figure: {out_png}")


if __name__ == "__main__":
    main()
