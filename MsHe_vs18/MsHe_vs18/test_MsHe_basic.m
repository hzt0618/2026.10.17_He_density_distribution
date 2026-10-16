% test_MsHe_basic.m
% Basic smoke test for MsHe.m

clear;
clc;

script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..', '..');

glow_file = fullfile(repo_root, 'glow model', 'glow model', 'output', ...
    'Lat-60.00_Lon109.00_20251009', '282', 'Lat-60.00_Lon109.00_20251009_vs025.txt');

if ~isfile(glow_file)
    error('Test glow file not found: %s', glow_file);
end

r = MsHe(glow_file, struct('make_plot', false));

assert(isfield(r, 'profile') && isfield(r.profile, 'n_He_cm3'), 'Missing n_he output');
assert(numel(r.profile.n_He_cm3) == r.grid.nz, 'n_he length mismatch');
assert(all(isfinite(r.profile.n_He_cm3)), 'n_he contains non-finite values');

fprintf('test_MsHe_basic passed. max n_he = %.6e at %.1f km\n', ...
    r.profile.n_he_max_cm3, r.profile.n_he_max_alt_km);

