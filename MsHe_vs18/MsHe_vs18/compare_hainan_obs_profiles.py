from __future__ import annotations

import csv
import argparse
import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerTuple
from scipy.io import loadmat

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "axes.spines.right": True,
        "axes.spines.top": True,
        "axes.linewidth": 0.75,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
    }
)


def set_full_box(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.75)


@dataclass
class ObsProfile:
    path: Path
    local_time: datetime
    ut_hour: float
    altitude_km: np.ndarray
    density_cm3: np.ndarray
    err_pct: np.ndarray


def parse_obs_file(path: Path) -> ObsProfile:
    local_time = None
    alts = []
    dens = []
    errs = []
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#DataTime:"):
                local_time = datetime.fromisoformat(line.split(":", 1)[1].strip())
            elif line.startswith("#") or line.startswith("Elevation") or not line.strip():
                continue
            else:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        alts.append(float(parts[0]))
                        dens.append(float(parts[1]))
                        errs.append(float(parts[2]))
                    except ValueError:
                        pass
    if local_time is None:
        raise ValueError(f"Missing #DataTime in {path}")
    ut_time = local_time - timedelta(hours=8)
    ut_hour = ut_time.hour + ut_time.minute / 60 + ut_time.second / 3600
    return ObsProfile(
        path=path,
        local_time=local_time,
        ut_hour=ut_hour,
        altitude_km=np.array(alts, dtype=float),
        density_cm3=np.array(dens, dtype=float),
        err_pct=np.array(errs, dtype=float),
    )


def nearest_index(values: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(values - target)))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare Hainan observed and modeled He* profiles.")
    base_default = Path(r"D:\1CH2\1study\DXK\zhanshi\spaceborne\code")
    parser.add_argument("--base", type=Path, default=base_default)
    parser.add_argument("--model-dir", type=Path, default=base_default / "MsHe_vs18" / "output" / "hainna_20251009_sth_05")
    parser.add_argument("--obs-dir", type=Path, default=base_default / "hainan_data" / "20251009" / "ODAZH_ALAL01_NDMH_L2A")
    parser.add_argument("--out-dir", type=Path, default=base_default / "MsHe_vs18" / "output" / "compare_hainan_20251009_model_vs_obs")
    return parser.parse_args()


def find_single_file(model_dir: Path, suffix: str) -> Path:
    matches = sorted(model_dir.glob(f"*{suffix}"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected exactly one '*{suffix}' in {model_dir}, found {len(matches)}")
    return matches[0]


def main() -> None:
    args = parse_args()
    model_dir = args.model_dir
    obs_dir = args.obs_dir
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    mat_path = find_single_file(model_dir, "_nhe_full.mat")
    results_csv = find_single_file(model_dir, "_results.csv")
    model_mat = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    z_km = np.array(model_mat["z_km"], dtype=float).ravel()
    ut_h = np.array(model_mat["ut_h"], dtype=float).ravel()
    n_he_matrix = np.array(model_mat["n_he_matrix"], dtype=float)
    results_rows = list(csv.DictReader(results_csv.open("r", encoding="utf-8-sig")))

    obs_profiles = [parse_obs_file(p) for p in sorted(obs_dir.glob("*.TXT"))]

    summary_rows = []
    for obs in obs_profiles:
        mi = nearest_index(ut_h, obs.ut_hour)
        model_ut = float(ut_h[mi])
        model_vs = int(results_rows[mi]["vs_idx"])
        model_peak = float(results_rows[mi]["n_he_max_cm3"])
        model_peak_alt = float(results_rows[mi]["n_he_max_alt_km"])

        valid_obs = np.isfinite(obs.density_cm3)
        valid_obs_hi = valid_obs & (obs.altitude_km >= 400)
        obs_peak = float(np.max(obs.density_cm3[valid_obs_hi])) if np.any(valid_obs_hi) else math.nan
        obs_peak_alt = float(obs.altitude_km[valid_obs_hi][np.argmax(obs.density_cm3[valid_obs_hi])]) if np.any(valid_obs_hi) else math.nan

        fig, ax = plt.subplots(figsize=(3.46, 2.72))
        model_profile = n_he_matrix[mi, :]
        valid_model = np.isfinite(model_profile) & (model_profile > 0) & (z_km > 200)
        model_color = "black"
        obs_color = "#d62728"
        model_line, = ax.plot(
            model_profile[valid_model],
            z_km[valid_model],
            "-",
            lw=1.8,
            color=model_color,
            label="Model",
        )

        obs_keep = valid_obs & (obs.altitude_km > 200)
        obs_den = obs.density_cm3[obs_keep]
        obs_alt = obs.altitude_km[obs_keep]
        obs_err = obs.err_pct[obs_keep]
        err_frac = np.clip(np.abs(obs_err) / 100.0, 0.0, 5.0)
        xerr = np.abs(obs_den) * err_frac
        err_keep = obs_den > 0
        obs_band = None
        if np.any(err_keep):
            obs_band = ax.fill_betweenx(
                obs_alt[err_keep],
                np.maximum(obs_den[err_keep] - xerr[err_keep], 0.0),
                obs_den[err_keep] + xerr[err_keep],
                color=obs_color,
                alpha=0.20,
                linewidth=0,
            )
        obs_line, = ax.plot(
            obs_den,
            obs_alt,
            "o-",
            ms=2.6,
            lw=1.8,
            color=obs_color,
            markeredgewidth=0,
            label="Hainan",
        )

        x_candidates = []
        model_vis = valid_model & (z_km <= 1000)
        if np.any(model_vis):
            x_candidates.extend(model_profile[model_vis].tolist())
        obs_vis = (obs_alt >= 200) & (obs_alt <= 1000)
        if np.any(obs_vis):
            x_candidates.extend(obs_den[obs_vis].tolist())
            err_vis = obs_vis & (obs_den > 0)
            if np.any(err_vis):
                x_candidates.extend(np.maximum(obs_den[err_vis] - xerr[err_vis], 0.0).tolist())
                x_candidates.extend((obs_den[err_vis] + xerr[err_vis]).tolist())
        x_candidates = np.array([x for x in x_candidates if np.isfinite(x)], dtype=float)
        if x_candidates.size > 0:
            x_lo = np.nanmin(x_candidates)
            x_hi = np.nanmax(x_candidates)
            span = x_hi - x_lo
            if span <= 0:
                span = max(abs(x_hi), 0.1)
            pad = max(0.08 * span, 0.05)
            x_lo = x_lo - pad
            x_hi = x_hi + pad
            if np.nanmin(x_candidates) >= 0:
                x_lo = 0.0
            ax.set_xlim(x_lo, x_hi)

        ax.set_xlabel(r"$n_{\mathrm{He}}$ [cm$^{-3}$]")
        ax.set_ylabel("Altitude [km]")
        ax.grid(False)
        obs_handle = (obs_line, obs_band) if obs_band is not None else obs_line
        ax.legend(
            [obs_handle, model_line],
            ["Hainan", "Model"],
            loc="lower right",
            frameon=False,
            handlelength=1.7,
            borderpad=0.3,
            labelspacing=0.3,
            handler_map={tuple: HandlerTuple(ndivide=None)},
        )
        ax.set_ylim(200, 1000)
        ax.tick_params(direction="out")
        set_full_box(ax)

        stem = obs.path.stem
        png_path = out_dir / f"{stem}_compare.png"
        pdf_path = out_dir / f"{stem}_compare.pdf"
        fig.tight_layout(pad=0.3)
        fig.savefig(png_path, dpi=600)
        fig.savefig(pdf_path)
        plt.close(fig)

        summary_rows.append(
            {
                "obs_file": obs.path.name,
                "obs_local_time": obs.local_time.isoformat(sep="T"),
                "obs_ut_h": f"{obs.ut_hour:.6f}",
                "matched_model_vs": model_vs,
                "matched_model_ut_h": f"{model_ut:.6f}",
                "model_peak_cm3": f"{model_peak:.6f}",
                "model_peak_alt_km": f"{model_peak_alt:.1f}",
                "obs_peak_gt400_cm3": "" if math.isnan(obs_peak) else f"{obs_peak:.6f}",
                "obs_peak_gt400_alt_km": "" if math.isnan(obs_peak_alt) else f"{obs_peak_alt:.1f}",
                "plot_png": png_path.name,
                "plot_pdf": pdf_path.name,
            }
        )

    summary_csv = out_dir / "comparison_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    overview_txt = out_dir / "comparison_overview.txt"
    with overview_txt.open("w", encoding="utf-8") as f:
        f.write("Assumption: observation DataTime is local Beijing time (UTC+8); model uses UT.\n")
        f.write("Each observation profile is paired with the nearest model UT sample.\n")
        f.write("Observation peak statistics in the CSV use positive densities above 400 km only.\n")


if __name__ == "__main__":
    main()
