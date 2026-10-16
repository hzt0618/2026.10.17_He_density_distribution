from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


FOLDER = Path(
    r"D:\1CH2\1study\DXK\zhanshi\spaceborne\code\MsHe_vs18\output\nor_20230225-28_sth_05"
)
TAGS = [
    "nor_20230225_sth_05",
    "nor_20230226_sth_05",
    "nor_20230227_sth_05",
    "nor_20230228_sth_05",
]


def main() -> None:
    all_ut = []
    all_heat = []
    z_ref = None

    for day_index, tag in enumerate(TAGS):
        mat_path = FOLDER / f"{tag}_heatmap_nhe.mat"
        data = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
        ut_h = np.array(data["ut_h_sorted"], dtype=float).ravel()
        z_km = np.array(data["z_km_grid"], dtype=float).ravel()
        n_he = np.array(data["n_he_heatmap"], dtype=float)
        n_he[~np.isfinite(n_he)] = np.nan

        if z_ref is None:
            z_ref = z_km
        elif z_ref.shape != z_km.shape or not np.allclose(z_ref, z_km):
            raise RuntimeError(f"Altitude grid mismatch in {mat_path}")

        all_ut.append(ut_h + day_index * 24.0)
        all_heat.append(n_he)

    x = np.concatenate(all_ut)
    heat = np.vstack(all_heat)
    order = np.argsort(x)
    x = x[order]
    heat = heat[order, :]

    valid = heat[np.isfinite(heat) & (heat > 0)]
    if valid.size == 0:
        raise RuntimeError("No finite positive n_he values found.")
    clim_low, clim_high = np.percentile(valid, [2, 98])
    if clim_high <= clim_low:
        clim_low = float(np.nanmin(valid))
        clim_high = float(np.nanmax(valid))
    if clim_high <= clim_low:
        clim_high = clim_low + np.finfo(float).eps

    fig, ax = plt.subplots(figsize=(16, 5.6), facecolor="white")
    x_grid, z_grid = np.meshgrid(x, z_ref)
    contour = ax.contourf(
        x_grid,
        z_grid,
        heat.T,
        levels=30,
        cmap="turbo",
        vmin=clim_low,
        vmax=clim_high,
    )

    for boundary in [24, 48, 72]:
        ax.axvline(boundary, color="k", linestyle="--", linewidth=1.0, alpha=0.65)

    ax.set_xlim(0, 96)
    ax.set_ylim(float(np.nanmin(z_ref)), float(np.nanmax(z_ref)))
    ax.set_xlabel("Continuous UT [hour from 2023-02-25 00:00]")
    ax.set_ylabel("Altitude [km]")
    ax.set_title("He(2^3S) density | nor sth 05 | 2023-02-25 to 2023-02-28")
    ax.set_xticks([0, 12, 24, 36, 48, 60, 72, 84, 96])
    ax.set_xticklabels([
        "02-25 00",
        "02-25 12",
        "02-26 00",
        "02-26 12",
        "02-27 00",
        "02-27 12",
        "02-28 00",
        "02-28 12",
        "03-01 00",
    ])
    ax.grid(True, alpha=0.25)

    cbar = fig.colorbar(contour, ax=ax)
    cbar.set_label(r"$n_{he}$ [cm$^{-3}$]")

    fig.tight_layout()
    out_png = FOLDER / "nor_20230225_20230228_sth_05_fig1_heatmap_continuous.png"
    out_pdf = FOLDER / "nor_20230225_20230228_sth_05_fig1_heatmap_continuous.pdf"
    fig.savefig(out_png, dpi=220)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(out_png)
    print(out_pdf)


if __name__ == "__main__":
    main()
