# Repository Guidelines

## Project Structure & Module Organization
- Core source lives in `glow model/` and is organized as Intel Fortran units (`*.f90`, `*.f`).
- Entry point is `glow model/main.f90` (program `test`) calling `run_glowbasic_Heion_args` from `glowbasic_Heion_api.f90`.
- Physics/model modules are split by concern (for example `mzgrid*.f90`, `ephoto.f90`, `ssflux.f90`, `geomag.f90`).
- Runtime inputs are under `glow model/data/` (`iri90/`, solar flux tables, NOEM data).
- Generated artifacts are in `glow model/x64/` and large simulation outputs in folders like `Lat19.50_Lon109.00_20240106/`.

## Build, Test, and Development Commands
- This is a Visual Studio + Intel Fortran project (`glow model.sln`, `glow model.vfproj`).
- Build (Developer Command Prompt):
  - `msbuild "glow model\\glow model.sln" /p:Configuration=Release /p:Platform=x64`
- Build from IDE:
  - `devenv "glow model\\glow model.sln" /Build "Release|x64"`
- Run executable:
  - `."glow model\\x64\\Release\\glow model.exe"`
- Use x64 by default; `vfproj` enables Intel `ifx` for x64 configs.

## Coding Style & Naming Conventions
- Follow existing Fortran style: 2-space indentation, `implicit none`, one module/program per file where practical.
- Keep filenames descriptive and lowercase with underscores (for example `glow_config_module.f90`).
- Prefer explicit `intent(in|out|inout)` and grouped `use ..., only:` imports.
- Keep path handling centralized in helpers like `write_path.f90`; avoid hard-coded absolute paths.

## Testing Guidelines
- No automated unit-test framework is currently configured.
- Validate by running representative cases and checking generated `output_vs*.txt` structure/content (see `data_readme.txt` for expected layout).
- For model changes, compare key fields (`sza`, `uflx`, `dflx`, `heion`, `ze`) against a baseline case.

## Commit & Pull Request Guidelines
- Local `.git` metadata is not present in this snapshot; use a consistent convention going forward.
- Commit format: `type(scope): short summary` (for example `fix(run_glow4): correct UT wrap logic`).
- Keep commits focused (source vs. generated outputs separated).
- PRs should include: purpose, changed files/modules, validation case(s), and before/after sample numbers or plots when physics outputs change.

## Data & Artifact Hygiene
- Do not commit build products from `x64/` or transient IDE state from `.vs/`.
- Large generated output directories should be added only when required for reproducible scientific results.
