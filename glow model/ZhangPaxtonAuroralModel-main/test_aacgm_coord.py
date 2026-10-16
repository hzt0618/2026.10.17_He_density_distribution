from datetime import datetime, timezone

import aacgmv2


# User-editable test input
LAT = 80.0
LON = 109.0
ALT_KM = 100.0
TIME_UTC = datetime(2026, 1, 20, 0, 0, 0, tzinfo=timezone.utc)
METHOD = "ALLOWTRACE"


def main() -> None:
    mlat, mlon, mlt = aacgmv2.get_aacgm_coord(
        LAT,
        LON,
        ALT_KM,
        TIME_UTC,
        method=METHOD,
    )

    print("Input:")
    print(f"  lat      = {LAT}")
    print(f"  lon      = {LON}")
    print(f"  alt_km   = {ALT_KM}")
    print(f"  time_utc = {TIME_UTC.isoformat()}")
    print(f"  method   = {METHOD}")

    print("\nOutput:")
    print(f"  mlat = {mlat}")
    print(f"  mlon = {mlon}")
    print(f"  mlt  = {mlt}")


if __name__ == "__main__":
    main()
