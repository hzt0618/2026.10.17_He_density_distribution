from pathlib import Path

import pandas as pd


def extract_daily_fields(csv_path: Path) -> dict[str, object]:
    df = pd.read_csv(csv_path)

    ap_columns = ["Ap0", "Ap3", "Ap6", "Ap9", "Ap12", "Ap15", "Ap18", "Ap21"]

    year = df["year"].to_numpy()
    month = df["month"].to_numpy()
    day = df["day"].to_numpy()
    ap_8 = df[ap_columns].to_numpy()
    f107_obs = df["f107_obs"].to_numpy()
    f107_81ctr_obs = df["f107_81ctr_obs"].to_numpy()
    f107_81lst_obs = df["f107_81lst_obs"].to_numpy()

    return {
        "year": year,
        "month": month,
        "day": day,
        "ap_8": ap_8,
        "f107_obs": f107_obs,
        "f107_81ctr_obs": f107_81ctr_obs,
        "f107_81lst_obs": f107_81lst_obs,
    }


def main() -> int:
    script_dir = Path(__file__).resolve().parent
    candidate_paths = [
        script_dir.parent.parent / "sw_daily_latest.csv",
        script_dir / "data" / "sw_daily_latest.csv",
    ]
    csv_path = next((p for p in candidate_paths if p.exists()), None)
    if csv_path is None:
        raise FileNotFoundError(
            "CSV not found. Checked:\n"
            + "\n".join(str(p) for p in candidate_paths)
        )

    fields = extract_daily_fields(csv_path)

    print(f"Loaded: {csv_path}")
    print(f"rows: {len(fields['year'])}")
    print(
        "variables: year, month, day, ap_8, f107_obs, "
        "f107_81ctr_obs, f107_81lst_obs"
    )
    print(f"year[:3] = {fields['year'][:3]}")
    print(f"ap_8[:2] =\n{fields['ap_8'][:2]}")
    print(f"f107_obs[:3] = {fields['f107_obs'][:3]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
