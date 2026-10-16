#!/usr/bin/env python
"""Update and export space weather indices to local files.

Default output:
  - D:/1CH2/1study/DXK/zhanshi/星载/15_程序集合/MsHe_vs18/f107/data/sw_daily_latest.csv
  - D:/1CH2/1study/DXK/zhanshi/星载/15_程序集合/MsHe_vs18/f107/data/sw_ap_kp_3h_latest.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import requests
import spaceweather as sw
import urllib3


def _install_relaxed_ssl_download() -> None:
    """Work around expired upstream certificates on CelesTrak."""
    try:
        import spaceweather.celestrak as celestrak
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(f"Unable to load spaceweather.celestrak: {exc}") from exc

    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

    def _dl_file_noverify(swpath, url):
        with requests.get(url, stream=True, verify=False, timeout=120) as response:
            response.raise_for_status()
            with open(swpath, "wb") as fd:
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        fd.write(chunk)

    celestrak._dl_file = _dl_file_noverify


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download/update latest space weather data and save CSV files."
    )
    parser.add_argument(
        "--out-dir",
        default=r"D:\1CH2\1study\DXK\zhanshi\星载\15_程序集合\MsHe_vs18\f107\data",
        help="Output directory for CSV files.",
    )
    parser.add_argument(
        "--daily-name",
        default="sw_daily_latest.csv",
        help="Filename for daily indices CSV.",
    )
    parser.add_argument(
        "--kp3h-name",
        default="sw_ap_kp_3h_latest.csv",
        help="Filename for 3-hour Ap/Kp CSV.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    daily_path = out_dir / args.daily_name
    kp3h_path = out_dir / args.kp3h_name

    try:
        _install_relaxed_ssl_download()

        print("Updating space weather data from remote source...")
        sw.update_data()

        print("Loading daily indices (Ap/Kp/F10.7)...")
        daily_df = sw.sw_daily()

        print("Loading 3-hourly Ap/Kp indices...")
        kp3h_df = sw.ap_kp_3h()

        daily_df.to_csv(daily_path, encoding="utf-8")
        kp3h_df.to_csv(kp3h_path, encoding="utf-8")

        print(f"Saved: {daily_path}")
        print(f"Saved: {kp3h_path}")
        if len(daily_df) > 0:
            print("Latest daily row:")
            print(daily_df.tail(1).to_string())
        if len(kp3h_df) > 0:
            print("Latest 3h row:")
            print(kp3h_df.tail(1).to_string())
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to update/export space weather data: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
