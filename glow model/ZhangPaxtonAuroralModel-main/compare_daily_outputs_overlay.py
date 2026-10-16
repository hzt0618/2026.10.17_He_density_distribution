from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "output"
DATE_TAGS = ["20260119"]


def load_day_table(date_tag: str) -> pd.DataFrame:
    file_path = OUTPUT_DIR / f"{date_tag}_24h.txt"
    if not file_path.is_file():
        raise FileNotFoundError(f"Missing file: {file_path}")

    df = pd.read_csv(file_path, comment="#")
    df = df[df["status"] == "ok"].copy()
    df["time_utc"] = pd.to_datetime(df["time_utc"], utc=True)
    df["ut_hour"] = (
        df["time_utc"].dt.hour
        + df["time_utc"].dt.minute / 60.0
        + df["time_utc"].dt.second / 3600.0
    )
    df["date_tag"] = date_tag
    return df.sort_values("ut_hour").reset_index(drop=True)


def main() -> None:
    frames = [load_day_table(date_tag) for date_tag in DATE_TAGS]
    date_label = DATE_TAGS[0]
    alt_km = frames[0]["alt_km"].iloc[0]

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    series_specs = [
        ("kp", "Kp"),
        ("ef", "Energy Flux (erg/s/cm^2)"),
        ("ec", "Mean Energy (eV)"),
    ]

    for df in frames:
        label = df["date_tag"].iat[0]
        for ax, (column, ylabel) in zip(axes, series_specs):
            ax.plot(
                df["ut_hour"],
                df[column],
                marker="o",
                linewidth=1.5,
                markersize=3,
                label=label,
            )
            ax.set_ylabel(ylabel)
            ax.grid(True, alpha=0.3)

    for ax in axes:
        ax.legend(loc="best")

    axes[-1].set_xlabel("UT Hour")
    axes[-1].set_xlim(0, 23)
    axes[-1].set_xticks(range(0, 24, 2))
    fig.suptitle(f"Single-day series for {date_label} at lon=109, lat=80, alt={alt_km:g} km")
    fig.tight_layout()

    out_png = OUTPUT_DIR / f"compare_{date_label}.png"
    fig.savefig(out_png, dpi=200)
    #plt.close(fig)
    print(f"Saved figure: {out_png}")


if __name__ == "__main__":
    main()
