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
]


def main() -> None:
    panels = []
    valid_values = []

    for tag in TAGS:
        mat_path = FOLDER / f"{tag}_heatmap_nhe.mat"
        data = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
        ut_h = np.array(data["ut_h_sorted"], dtype=float).ravel()
        z_km = np.array(data["z_km_grid"], dtype=float).ravel()
        n_he = np.array(data["n_he_heatmap"], dtype=float)
        n_he[~np.isfinite(n_he)] = np.nan
        panels.append((tag, ut_h, z_km, n_he))

        vals = n_he[np.isfinite(n_he) & (n_he > 0)]
        if vals.size:
            valid_values.append(vals)

    if not valid_values:
        raise RuntimeError("No finite positive n_he values found.")

    all_valid = np.concatenate(valid_values)
    clim_low, clim_high = np.percentile(all_valid, [2, 98])
    if clim_high <= clim_low:
        clim_low = float(np.nanmin(all_valid))
        clim_high = float(np.nanmax(all_valid))
    if clim_high <= clim_low:
        clim_high = clim_low + np.finfo(float).eps

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(12.5, 10.5),
        sharex=True,
        sharey=True,
        constrained_layout=True,
        facecolor="white",
    )

    contour = None
    for ax, (tag, ut_h, z_km, n_he) in zip(axes, panels):
        date_label = tag.split("_")[1]
        x_grid, z_grid = np.meshgrid(ut_h, z_km)
        contour = ax.contourf(
            x_grid,
            z_grid,
            n_he.T,
            levels=30,
            cmap="turbo",
            vmin=clim_low,
            vmax=clim_high,
        )
        ax.set_title(f"{date_label} | {tag}", fontsize=12)
        ax.set_ylabel("Altitude [km]")
        ax.grid(True, alpha=0.25)
        ax.set_xlim(0, 24)

    axes[-1].set_xlabel("UT [hour]")
    fig.suptitle("He(2^3S) density heatmaps | nor sth 05 | 2023-02-25 to 2023-02-27", fontsize=14)

    if contour is not None:
        cbar = fig.colorbar(contour, ax=axes, location="right", shrink=0.96)
        cbar.set_label(r"$n_{he}$ [cm$^{-3}$]")

    out_png = FOLDER / "nor_20230225_20230227_sth_05_fig1_heatmap_combined.png"
    out_pdf = FOLDER / "nor_20230225_20230227_sth_05_fig1_heatmap_combined.pdf"
    fig.savefig(out_png, dpi=220)
    fig.savefig(out_pdf)
    plt.close(fig)
    print(out_png)
    print(out_pdf)


if __name__ == "__main__":
    main()
