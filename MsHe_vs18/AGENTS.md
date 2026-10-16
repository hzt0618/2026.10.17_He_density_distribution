# Repository Guidelines

## Project Structure & Module Organization
- Main code lives in `MsHe_vs18/` and is MATLAB-first (`*.m`).
- Primary entry scripts are `MsHe.m`, `MsHe_all.m`, and versioned variants such as `MsHe_vs18.m`.
- Physics helper modules are colocated with scripts (for example `he_produc*.m`, `he_loss.m`, `g_factor*.m`, `ExtensionHeight*.m`).
- Data assets (`*.mat`, `*.txt`) are stored in the same directory and are consumed directly by scripts.
- There is no dedicated `tests/` folder yet; validation is currently script-driven.

## Build, Test, and Development Commands
- Run the main model (batch):
  ```powershell
  matlab -batch "cd('MsHe_vs18'); MsHe"
  ```
- Run yearly/all-scenario workflow:
  ```powershell
  matlab -batch "cd('MsHe_vs18'); MsHe_all"
  ```
- Run a specific versioned experiment:
  ```powershell
  matlab -batch "cd('MsHe_vs18'); MsHe_vs18"
  ```
- Quick syntax check for all scripts:
  ```powershell
  matlab -batch "cd('MsHe_vs18'); checkcode('*.m')"
  ```

## Coding Style & Naming Conventions
- Use 4-space indentation and keep one statement per line.
- Keep function names lowercase with underscores for physics utilities (for example `read_glowdata`, `he_dif_e`).
- Keep scenario scripts versioned with `vs` suffixes (`MsHe_vs17_8.m`, `he_denvs7.m`) instead of overwriting baselines.
- Prefer explicit variable names for geophysical quantities (`glat`, `glong`, `doy`, `eflx`).

## Testing Guidelines
- No formal `matlab.unittest` suite exists yet; use reproducible smoke runs.
- For every change, run at least one baseline script (`MsHe` or `MsHe_vs18`) and one affected variant.
- Compare key outputs (profiles/figures or saved arrays) against prior results before merging.
- If adding tests, place them under `MsHe_vs18/tests/` using `test_*.m` naming.

## Commit & Pull Request Guidelines
- Current repository workflow is version-oriented; keep commits scoped to one model or helper change.
- Use imperative commit subjects, e.g., `model: update he_produc2 recombination term`.
- In PRs include: purpose, changed scripts/data files, run command(s), and before/after evidence (figure snapshot or numeric delta summary).
- Link related issue/task IDs when available and call out any data dependency changes explicitly.

## Data & Configuration Notes
- Many scripts reference local absolute paths; prefer converting new code to repo-relative paths.
- Do not commit generated temporary outputs unless they are intentional reference artifacts.
