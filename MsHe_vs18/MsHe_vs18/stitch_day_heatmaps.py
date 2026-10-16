from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


def main() -> None:
    folder = Path(r"D:\1CH2\1study\DXK\zhanshi\spaceborne\code\MsHe_vs18\output\Lat78.10_Lon16")
    tags = ["20230225", "20230226", "20230227", "20230228"]

    all_ut = []
    all_heat = []
    z_km = None
    day_offsets = []

    for day_index, tag in enumerate(tags):
        mat = loadmat(folder / f"heatmap_nhe_{tag}.mat", squeeze_me=True, struct_as_record=False)
        ut = np.array(mat["ut_h_sorted"], dtype=float).ravel()
        heat = np.array(mat["n_he_heatmap"], dtype=float)
        z = np.array(mat["z_km_grid"], dtype=float).ravel()
        if z_km is None:
            z_km = z
        day_offsets.append(day_index * 24.0)
        all_ut.append(ut + day_index * 24.0)
        all_heat.append(heat)

    x = np.concatenate(all_ut)
    h = np.vstack(all_heat)

    fig, ax = plt.subplots(figsize=(15, 5.5), facecolor="white")
    X, Z = np.meshgrid(x, z_km)
    cf = ax.contourf(X, Z, h.T, 30, cmap="turbo", linewidths=0)

    valid = h[np.isfinite(h) & (h > 0)]
    if valid.size > 0:
        c_lo = np.percentile(valid, 2)
        c_hi = np.percentile(valid, 98)
        if c_hi <= c_lo:
            c_lo = float(np.min(valid))
            c_hi = float(np.max(valid))
        cf.set_clim(c_lo, c_hi)

    for boundary in [24, 48, 72]:
        ax.axvline(boundary, color="k", linestyle="--", linewidth=1.0, alpha=0.7)

    ax.set_xlim(float(np.min(x)), float(np.max(x)))
    ax.set_ylim(float(np.min(z_km)), float(np.max(z_km)))
    ax.set_xlabel("Continuous UT [hour from 2023-02-25 00:00]")
    ax.set_ylabel("Altitude [km]")
    ax.set_title("He(2^3S) density | Lat=78.10, Lon=16.00 | 2023-02-25 to 2023-02-28")
    ax.set_xticks([0, 12, 24, 36, 48, 60, 72, 84, 96])
    ax.set_xticklabels(["02-25 00", "02-25 12", "02-26 00", "02-26 12", "02-27 00", "02-27 12", "02-28 00", "02-28 12", "03-01 00"])
    ax.grid(True, alpha=0.25)

    cb = fig.colorbar(cf, ax=ax)
    cb.set_label("n_he [cm$^{-3}$]")

    plt.tight_layout()
    out_png = folder / "fig1_heatmap_20230225_20230228_continuous.png"
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
