% test_MsHe.m
% Smoke test for MsHe.m

clear;
clc;

script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..', '..');

% Preferred test input (user-provided path)
glow_file = 'D:\\1CH2\\1study\\DXK\\zhanshi\\鏄熻浇\\code\\glow model\\glow model\\output\\HaiNan20240106\\6Module\\output_vs50.txt';

% Fallback test input
if ~isfile(glow_file)
    glow_file = fullfile(script_dir, '..', 'glow_data', 'HaiNan20240106', '006', 'output_vs050.txt');
end

if ~isfile(glow_file)
    error('No valid test input found. Checked preferred and fallback glow files.');
end

fprintf('Testing MsHe with file:\n%s\n\n', glow_file);

opts = struct('make_plot', true);
r = MsHe(glow_file, opts);

% Structural checks
required_top_fields = {'input', 'meta', 'grid', 'raw', 'profile'};
for k = 1:numel(required_top_fields)
    assert(isfield(r, required_top_fields{k}), 'Missing field: r.%s', required_top_fields{k});
end

required_profile_fields = {'n_He_cm3', 'qe', 'qr', 'q', 'Loss', 'Dif', 'Hn_cm'};
for k = 1:numel(required_profile_fields)
    assert(isfield(r.profile, required_profile_fields{k}), ...
        'Missing field: r.profile.%s', required_profile_fields{k});
end

% Numeric checks
z_km = r.grid.z_km(:);
n_he = r.profile.n_He_cm3(:);
assert(numel(z_km) == r.grid.nz, 'Grid size mismatch: z_km vs grid.nz');
assert(numel(n_he) == r.grid.nz, 'Profile size mismatch: n_he vs grid.nz');
assert(all(isfinite(n_he)), 'n_he has non-finite values');

fprintf('MsHe smoke test passed.\n');
fprintf('nz = %d, z range = [%.1f, %.1f] km\n', r.grid.nz, min(z_km), max(z_km));
fprintf('n_he max = %.6e cm^-3 at %.1f km\n', ...
    r.profile.n_he_max_cm3, r.profile.n_he_max_alt_km);

