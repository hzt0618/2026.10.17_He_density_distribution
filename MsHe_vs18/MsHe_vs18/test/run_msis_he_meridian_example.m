%% Run NRLMSISE-00 He density along one longitude meridian.
% Edit the parameters in this block, then run this script.

clear; clc;

% ===================== USER PARAMETERS =====================
date_input = '20251221';        % 'yyyyMMdd' or 'yyyy-MM-dd'
longitude_deg = 109.0;           % degrees east, west is negative
ut_hour = 23.0;                 % UT hour of day

lat_deg = -90:1:90;             % full meridian latitude grid
alt_km = 100:10:1000;           % altitude grid in km
make_plot = true;               % save PNG heatmap

% Leave empty to save outputs in this test folder.
out_dir = '';

% Leave empty to use the repository sw_daily_latest.csv automatically.
sw_csv = '';
% ===========================================================

script_dir = fileparts(mfilename('fullpath'));
if isempty(out_dir)
    out_dir = script_dir;
end

opts = struct();
opts.ut_s = ut_hour * 3600;
opts.lat_deg = lat_deg;
opts.alt_km = alt_km;
opts.out_dir = out_dir;
opts.sw_csv = sw_csv;
opts.make_plot = make_plot;

out = run_msis_he_meridian(date_input, longitude_deg, opts);

fprintf('\nOutput files:\n');
fprintf('CSV: %s\n', out.csv_path);
fprintf('MAT: %s\n', out.mat_path);
if make_plot
    fprintf('PNG: %s\n', out.png_path);
end
