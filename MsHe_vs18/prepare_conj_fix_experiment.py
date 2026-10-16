#!/usr/bin/env python3
"""
Prepare or restore a minimal conjugate-point experiment for the 2025-10-09
terminator run near 109E.

Default action:
    - back up IGRF/output/conjugate_morning_math.txt once
    - replace the 5.5N row with a smooth test conjugate solution

Restore action:
    python prepare_conj_fix_experiment.py --restore
"""

from __future__ import annotations

import argparse
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
CONJ_PATH = REPO_ROOT / "IGRF" / "output" / "conjugate_morning_math.txt"
BACKUP_PATH = REPO_ROOT / "IGRF" / "output" / "conjugate_morning_math.bak.codex.txt"

TARGET_PREFIX = "8,5.500711,107.470826,"
ORIGINAL_LINE = (
    "8,5.500711,107.470826,-0.2988376342,0.9494773788,0.0958581007,"
    "30.344105,108.226936,0.000000,2.67007525"
)
PATCHED_LINE = (
    "8,5.500711,107.470826,-0.2988376342,0.9494773788,0.0958581007,"
    "-19.050000,108.365000,0.000000,2.67007525"
)


def restore_backup() -> None:
    if not BACKUP_PATH.is_file():
        raise FileNotFoundError(f"Backup not found: {BACKUP_PATH}")
    CONJ_PATH.write_text(BACKUP_PATH.read_text(encoding="utf-8"), encoding="utf-8")
    print(f"Restored original file from: {BACKUP_PATH}")


def apply_patch() -> None:
    if not CONJ_PATH.is_file():
        raise FileNotFoundError(f"Conjugate file not found: {CONJ_PATH}")

    text = CONJ_PATH.read_text(encoding="utf-8")
    if not BACKUP_PATH.exists():
        BACKUP_PATH.write_text(text, encoding="utf-8")
        print(f"Created backup: {BACKUP_PATH}")
    else:
        print(f"Backup already exists: {BACKUP_PATH}")

    lines = text.splitlines()
    replaced = False
    for i, line in enumerate(lines):
        if line.startswith(TARGET_PREFIX):
            lines[i] = PATCHED_LINE
            replaced = True
            break

    if not replaced:
        raise RuntimeError(f"Target row not found in: {CONJ_PATH}")

    CONJ_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Patched row in: {CONJ_PATH}")
    print("Old row:")
    print(f"  {ORIGINAL_LINE}")
    print("New row:")
    print(f"  {PATCHED_LINE}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--restore",
        action="store_true",
        help="restore conjugate_morning_math.txt from the backup created by this script",
    )
    args = parser.parse_args()

    if args.restore:
        restore_backup()
    else:
        apply_patch()


if __name__ == "__main__":
    main()
