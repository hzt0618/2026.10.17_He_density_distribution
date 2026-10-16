function out = run_msis_he_meridian(date_input, longitude_deg, opts)
%RUN_MSIS_HE_MERIDIAN Get MSIS neutral He density along a full meridian.
%
%   out = run_msis_he_meridian(date_input, longitude_deg)
%   out = run_msis_he_meridian(date_input, longitude_deg, opts)
%
% Inputs
%   date_input    : datetime, 'yyyyMMdd', 'yyyy-MM-dd', or datenum-compatible date
%   longitude_deg : longitude in degrees east
%   opts          : optional struct
%       opts.ut_s       : UT seconds of day, default 12*3600
%       opts.lat_deg    : latitude grid, default -90:1:90
%       opts.alt_km     : altitude grid, default 100:10:1000
%       opts.out_dir    : output folder, default this test folder
%       opts.sw_csv     : solar/geomagnetic CSV path
%       opts.make_plot  : save heatmap PNG, default true
%
% Outputs
%   Files are saved under opts.out_dir:
%       msis_he_meridian_*.csv : long table (lat, alt, He density)
%       msis_he_meridian_*.mat : MATLAB data
%       msis_he_meridian_*.png : latitude-altitude heatmap
%
% Notes
%   atmosnrlmsise00 returns He number density in m^-3. This script stores
%   He density in cm^-3 to match the MsHe code convention.

if nargin < 1 || isempty(date_input)
    error('date_input is required, e.g. ''20230225''.');
end
if nargin < 2 || isempty(longitude_deg)
    error('longitude_deg is required, e.g. 16.0.');
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..', '..', '..');
default_sw_csv = fullfile(repo_root, 'glow model', 'glow model', ...
    'data', 'f107', 'data', 'sw_daily_latest.csv');

opts = apply_defaults(opts, struct( ...
    'ut_s', 12 * 3600, ...
    'lat_deg', -90:1:90, ...
    'alt_km', 100:10:1000, ...
    'out_dir', script_dir, ...
    'sw_csv', default_sw_csv, ...
    'make_plot', true));
if isempty(opts.sw_csv)
    opts.sw_csv = default_sw_csv;
end

dt = normalize_date(date_input);
date_tag = char(string(dt, 'yyyyMMdd'));
doy = day(dt, 'dayofyear');
yr = year(dt);

if ~isfolder(opts.out_dir)
    mkdir(opts.out_dir);
end

sw = get_sw_daily(opts.sw_csv, dt);
f107 = sw.f107_obs;
f107a = sw.f107_81ctr_obs;
ap = sw.Apavg;

lat_deg = opts.lat_deg(:);
alt_km = opts.alt_km(:).';
nlat = numel(lat_deg);
nalt = numel(alt_km);

he_cm3 = nan(nlat, nalt);
temp_k = nan(nlat, nalt);

for i = 1:nlat
    % atmosnrlmsise00 altitude input is meters.
    [T, DEN] = atmosnrlmsise00(alt_km * 1000, lat_deg(i), longitude_deg, ...
        yr, doy, opts.ut_s, f107a, f107, ap);
    temp_k(i, :) = reshape(T(:, 2), 1, []);
    he_cm3(i, :) = reshape(DEN(:, 1), 1, []) * 1e-6;
end

[LAT, ALT] = ndgrid(lat_deg, alt_km);
result_table = table(LAT(:), ALT(:), he_cm3(:), temp_k(:), ...
    'VariableNames', {'lat_deg', 'alt_km', 'he_cm3', 'neutral_temp_K'});

lon_tag = format_lon_tag(longitude_deg);
ut_tag = sprintf('UT%05.0fs', opts.ut_s);
base_name = sprintf('msis_he_meridian_%s_lon%s_%s', date_tag, lon_tag, ut_tag);
csv_path = fullfile(opts.out_dir, [base_name '.csv']);
mat_path = fullfile(opts.out_dir, [base_name '.mat']);
png_path = fullfile(opts.out_dir, [base_name '.png']);

writetable(result_table, csv_path);

out = struct();
out.date = dt;
out.date_tag = date_tag;
out.longitude_deg = longitude_deg;
out.ut_s = opts.ut_s;
out.lat_deg = lat_deg;
out.alt_km = alt_km;
out.he_cm3 = he_cm3;
out.neutral_temp_K = temp_k;
out.sw = sw;
out.csv_path = csv_path;
out.mat_path = mat_path;
out.png_path = '';

save(mat_path, 'out', 'lat_deg', 'alt_km', 'he_cm3', 'temp_k', ...
    'longitude_deg', 'dt', 'doy', 'yr', 'f107', 'f107a', 'ap');

if opts.make_plot
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [120, 120, 1000, 560]);
    [LATP, ALTP] = meshgrid(lat_deg, alt_km);
    plot_data = he_cm3.';
    contourf(LATP, ALTP, plot_data, 30, 'LineStyle', 'none');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
    xlabel('Latitude [deg]');
    ylabel('Altitude [km]');
    title(sprintf('NRLMSISE-00 He density | %s | Lon=%.2f deg | UT=%.2f h', ...
        char(string(dt, 'yyyy-MM-dd')), longitude_deg, opts.ut_s / 3600), ...
        'Interpreter', 'none');
    cb = colorbar;
    cb.Label.String = 'He density [cm^{-3}]';
    colormap(turbo);
    grid on;
    box on;

    valid = plot_data(isfinite(plot_data) & plot_data > 0);
    if ~isempty(valid)
        lo = prctile(valid, 2);
        hi = prctile(valid, 98);
        if hi > lo
            clim([lo, hi]);
        end
    end

    exportgraphics(fig, png_path, 'Resolution', 220);
    close(fig);
    out.png_path = png_path;
end

fprintf('Saved CSV: %s\n', csv_path);
fprintf('Saved MAT: %s\n', mat_path);
if opts.make_plot
    fprintf('Saved PNG: %s\n', png_path);
end
fprintf('F10.7=%.3g, F10.7a=%.3g, Ap=%.3g\n', f107, f107a, ap);
end

function opts = apply_defaults(opts, defaults)
names = fieldnames(defaults);
for i = 1:numel(names)
    name = names{i};
    if ~isfield(opts, name) || isempty(opts.(name))
        opts.(name) = defaults.(name);
    end
end
end

function dt = normalize_date(date_input)
if isa(date_input, 'datetime')
    dt = date_input;
elseif isnumeric(date_input)
    if date_input > 1e7
        dt = datetime(num2str(date_input), 'InputFormat', 'yyyyMMdd');
    else
        dt = datetime(date_input, 'ConvertFrom', 'datenum');
    end
else
    s = char(string(date_input));
    if ~isempty(regexp(s, '^\d{8}$', 'once'))
        dt = datetime(s, 'InputFormat', 'yyyyMMdd');
    else
        dt = datetime(s, 'InputFormat', 'yyyy-MM-dd');
    end
end
dt = dateshift(dt, 'start', 'day');
end

function sw = get_sw_daily(sw_csv, dt)
if ~isfile(sw_csv)
    error('Solar/geomagnetic CSV not found: %s', sw_csv);
end

tbl = readtable(sw_csv, 'VariableNamingRule', 'preserve');
date_col = tbl.Properties.VariableNames{1};
date_vec = datetime(string(tbl.(date_col)), 'InputFormat', 'yyyy/M/d');
idx = find(date_vec == dt, 1, 'first');
if isempty(idx)
    error('Date %s not found in %s', char(string(dt, 'yyyy-MM-dd')), sw_csv);
end

sw = struct();
sw.csv_path = sw_csv;
sw.date = dt;
sw.f107_obs = numeric_field(tbl, 'f107_obs', idx);
sw.f107_81ctr_obs = numeric_field(tbl, 'f107_81ctr_obs', idx);
sw.Apavg = numeric_field(tbl, 'Apavg', idx);

if ~isfinite(sw.f107_obs) || sw.f107_obs <= 0
    error('Invalid f107_obs for %s in %s', char(string(dt, 'yyyy-MM-dd')), sw_csv);
end
if ~isfinite(sw.f107_81ctr_obs) || sw.f107_81ctr_obs <= 0
    error('Invalid f107_81ctr_obs for %s in %s', char(string(dt, 'yyyy-MM-dd')), sw_csv);
end
if ~isfinite(sw.Apavg) || sw.Apavg < 0
    error('Invalid Apavg for %s in %s', char(string(dt, 'yyyy-MM-dd')), sw_csv);
end
end

function x = numeric_field(tbl, name, idx)
if ~ismember(name, tbl.Properties.VariableNames)
    error('Column not found in solar/geomagnetic table: %s', name);
end
v = tbl.(name);
if isnumeric(v)
    x = v(idx);
else
    x = str2double(string(v(idx)));
end
end

function s = format_lon_tag(lon)
prefix = 'E';
if lon < 0
    prefix = 'W';
end
s = sprintf('%s%07.2f', prefix, abs(lon));
s = strrep(s, '.', 'p');
end
