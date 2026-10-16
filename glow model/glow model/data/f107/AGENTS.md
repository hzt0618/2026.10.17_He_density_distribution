# Repository Guidelines

## Project Structure & Module Organization
- `update_spaceweather.py`: main entry point that downloads and exports latest space weather datasets.
- `data/`: generated CSV outputs (`sw_daily_latest.csv`, `sw_ap_kp_3h_latest.csv`). Treat these as runtime artifacts unless a change explicitly updates baseline data.
- No dedicated test package exists yet; add tests under a top-level `tests/` directory when introducing non-trivial logic.

## Build, Test, and Development Commands
- `python -m venv .venv` then `.\.venv\Scripts\Activate.ps1`: create and activate a local virtual environment (Windows PowerShell).
- `python -m pip install --upgrade pip spaceweather pandas`: install runtime dependencies used by the script.
- `python update_spaceweather.py`: fetch remote data and write CSV files to `data/`.
- `python update_spaceweather.py --out-dir .\data --daily-name daily.csv --kp3h-name kp3h.csv`: example custom output run.

## Coding Style & Naming Conventions
- Follow PEP 8 with 4-space indentation and clear, small functions.
- Use type hints for public functions (`def main() -> int`) and `pathlib.Path` for filesystem paths.
- Prefer snake_case for variables/functions and descriptive CLI option names (`--kp3h-name`).
- Keep logging/user output concise and action-focused.

## Testing Guidelines
- Preferred framework: `pytest`.
- Place tests in `tests/` and name files `test_*.py` (for example, `tests/test_update_spaceweather.py`).
- Prioritize argument parsing, output path handling, and failure-path behavior (network/API errors).
- Run tests with `python -m pytest -q`.

## Commit & Pull Request Guidelines
- Use short, imperative commit subjects (for example, `Add CLI option for output directory`).
- Keep one logical change per commit; include data file updates only when intentional and described.
- PRs should include: purpose, behavior changes, commands run, and sample output paths.
- Link related issue/task IDs when available and include before/after notes for CLI behavior.

## Security & Configuration Tips
- Do not hardcode machine-specific absolute paths in new code; prefer CLI flags and relative defaults.
- Validate output directories before writing, and avoid committing secrets or tokens (none are required for current flow).
