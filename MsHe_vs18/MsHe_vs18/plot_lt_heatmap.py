from __future__ import annotations

import argparse
import os
from pathlib import Path
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy.io import loadmat


plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "font.size": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.75,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
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


def load_mat(path: Path) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return loadmat(path, squeeze_me=True, struct_as_record=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot an LT-altitude He(2^3S) heatmap from a merged LT output folder.")
    parser.add_argument("--run-tag", default="hainna_20251009_LT_0_24_sth_05")
    args = parser.parse_args()

    base = Path(r"D:\1CH2\1study\DXK\zhanshi\spaceborne\code\MsHe_vs18\output")
    run_tag = args.run_tag
    out_dir = base / run_tag
    data_path = out_dir / f"{run_tag}_heatmap_nhe.mat"

    data = load_mat(data_path)
    lt_h = np.asarray(data["lt_h_sorted"], dtype=float).ravel()
    z_km = np.asarray(data["z_km_grid"], dtype=float).ravel()
    n_he = np.asarray(data["n_he_heatmap"], dtype=float)

    order = np.argsort(lt_h)
    lt_h = lt_h[order]
    n_he = n_he[order, :]
    n_he[~np.isfinite(n_he)] = np.nan

    fig, ax = plt.subplots(figsize=(3.46, 2.35))
    lt_grid, z_grid = np.meshgrid(lt_h, z_km)
    im = ax.contourf(lt_grid, z_grid, n_he.T, 30, cmap="turbo")

    valid = n_he[np.isfinite(n_he) & (n_he > 0)]
    if valid.size:
        lo, hi = np.nanpercentile(valid, [2, 98])
        if hi <= lo:
            lo, hi = np.nanmin(valid), np.nanmax(valid)
        if hi <= lo:
            hi = lo + np.finfo(float).eps
        im.set_clim(lo, hi)

    ax.set_xlim(0, 24)
    ax.set_ylim(float(np.nanmin(z_km)), float(np.nanmax(z_km)))
    ax.set_xticks(np.arange(0, 25, 4))
    ax.set_xlabel("Local time [h]")
    ax.set_ylabel("Altitude [km]")
    ax.tick_params(direction="out")
    ax.grid(False)

    cb = fig.colorbar(im, ax=ax, pad=0.025)
    cb.set_label(r"$n_{\mathrm{He}}$ [cm$^{-3}$]")
    cb.locator = MaxNLocator(nbins=5)
    cb.update_ticks()
    cb.ax.tick_params(labelsize=8, width=0.7, length=2.5, direction="out")
    cb.outline.set_linewidth(0.7)

    png_path = out_dir / f"{run_tag}_fig1_heatmap_LT.png"
    pdf_path = out_dir / f"{run_tag}_fig1_heatmap_LT.pdf"
    saved_pdf_path = pdf_path
    try:
        with saved_pdf_path.open("wb") as fh:
            fig.savefig(fh, format="pdf", bbox_inches="tight", pad_inches=0.02)
    except PermissionError as exc:
        saved_pdf_path = out_dir / f"{run_tag}_fig1_heatmap_LT_updated_{os.getpid()}.pdf"
        print(f"Skipped locked PDF: {exc}")
        with saved_pdf_path.open("wb") as fh:
            fig.savefig(fh, format="pdf", bbox_inches="tight", pad_inches=0.02)
    try:
        fig.savefig(png_path, dpi=600, bbox_inches="tight", pad_inches=0.02)
    except PermissionError as exc:
        print(f"Skipped locked PNG: {exc}")
    plt.close(fig)

    print(png_path)
    print(saved_pdf_path)


if __name__ == "__main__":
    main()
