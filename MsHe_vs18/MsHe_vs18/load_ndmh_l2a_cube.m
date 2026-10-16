function out = load_ndmh_l2a_cube(data_dir, save_mat_path, model_heatmap_path)
%LOAD_NDMH_L2A_CUBE Read NDMH_L2A TXT files and build time-height cubes.
%   out = load_ndmh_l2a_cube(data_dir)
%   out = load_ndmh_l2a_cube(data_dir, save_mat_path)
%   out = load_ndmh_l2a_cube(data_dir, save_mat_path, model_heatmap_path)
%
% Outputs:
%   out.time_vec        : [nt x 1] datetime, sorted by DataTime
%   out.z_km            : [nz x 1] altitude grid (km)
%   out.den_3d          : [nt x nz x 1] density (cm^-3), missing->NaN
%   out.denerr_3d       : [nt x nz x 1] density error (%), missing->NaN
%   out.he_cube         : [nt x nz x 2], (:,:,1)=density, (:,:,2)=error
%   out.quality_flag    : [nt x 1] 1/2/3
%   out.file_list       : [nt x 1] source file names
%   out.valid_mask      : [nt x nz] density>0 & finite error & quality~=3
%   out.model_compare   : model-vs-observation matching info
%
% Missing values in source files:
%   Elevation missing: -99999999
%   Density missing  : -999999999999
%   DenEr missing    : -9999999999

if nargin < 1 || isempty(data_dir)
    script_dir = fileparts(mfilename('fullpath'));
    data_dir = fullfile(script_dir, '..', '..', 'hainan_data', '20251009', ...
        'ODAZH_ALAL01_NDMH_L2A');
end

if nargin < 2
    save_mat_path = '';
end

if nargin < 3
    model_heatmap_path = '';
end

if ~isfolder(data_dir)
    error('Data folder not found: %s', data_dir);
end

files = dir(fullfile(data_dir, '*.TXT'));
if isempty(files)
    files = dir(fullfile(data_dir, '*.txt'));
end
if isempty(files)
    error('No TXT files found under: %s', data_dir);
end

all_names = string({files.name})';
all_times = NaT(numel(files), 1);
all_qf = nan(numel(files), 1);
all_z = cell(numel(files), 1);
all_den = cell(numel(files), 1);
all_err = cell(numel(files), 1);

for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    [t_i, qf_i, z_i, den_i, err_i] = parse_one_file(fpath);
    all_times(i) = t_i;
    all_qf(i) = qf_i;
    all_z{i} = z_i;
    all_den{i} = den_i;
    all_err{i} = err_i;
end

% Sort by time. Unknown times are pushed to the end.
invalid_time = isnat(all_times);
safe_time = all_times;
safe_time(invalid_time) = datetime(9999, 1, 1);
[~, ord] = sort(safe_time);

all_names = all_names(ord);
all_times = all_times(ord);
all_qf = all_qf(ord);
all_z = all_z(ord);
all_den = all_den(ord);
all_err = all_err(ord);

% Build a common altitude grid from union of all files.
z_union = [];
for i = 1:numel(all_z)
    z_union = union(z_union, all_z{i});
end
z_km = z_union(:);
z_km = z_km(z_km >= 200);
nt = numel(all_names);
nz = numel(z_km);

den_2d = nan(nt, nz);
denerr_2d = nan(nt, nz);

for i = 1:nt
    [lia, locb] = ismember(all_z{i}, z_km);
    idx = locb(lia);
    den_2d(i, idx) = all_den{i}(lia);
    denerr_2d(i, idx) = all_err{i}(lia);
end

den_3d = reshape(den_2d, [nt, nz, 1]);
denerr_3d = reshape(denerr_2d, [nt, nz, 1]);
he_cube = cat(3, den_2d, denerr_2d);

quality_flag = all_qf(:);
valid_mask = (den_2d > 0) & isfinite(denerr_2d) & (quality_flag ~= 3);

out = struct();
out.time_vec = all_times(:);
out.z_km = z_km;
out.den_3d = den_3d;
out.denerr_3d = denerr_3d;
out.he_cube = he_cube;
out.quality_flag = quality_flag;
out.file_list = all_names(:);
out.valid_mask = valid_mask;

[peak_den, peak_idx] = max(den_2d, [], 2, 'omitnan');
peak_z_km = nan(nt, 1);
for i = 1:nt
    if isfinite(peak_den(i)) && isfinite(peak_idx(i)) && peak_idx(i) >= 1 && peak_idx(i) <= nz
        peak_z_km(i) = z_km(peak_idx(i));
    end
end

out.figure = plot_cubes(all_times, z_km, den_2d, denerr_2d, quality_flag, peak_den, peak_z_km);
out.peak_den_cm3 = peak_den;
out.peak_z_km = peak_z_km;

model_cmp = read_model_heatmap(model_heatmap_path, all_times, z_km);
out.model_compare = model_cmp;
out.figure_each_time = plot_each_time_profiles(all_times, z_km, den_2d, denerr_2d, quality_flag, model_cmp);

if ~isempty(save_mat_path)
    time_vec = out.time_vec; %#ok<NASGU>
    den_3d = out.den_3d; %#ok<NASGU>
    denerr_3d = out.denerr_3d; %#ok<NASGU>
    he_cube = out.he_cube; %#ok<NASGU>
    quality_flag = out.quality_flag; %#ok<NASGU>
    file_list = out.file_list; %#ok<NASGU>
    valid_mask = out.valid_mask; %#ok<NASGU>
    save(save_mat_path, 'time_vec', 'z_km', 'den_3d', 'denerr_3d', ...
        'he_cube', 'quality_flag', 'file_list', 'valid_mask');
end

fprintf('Loaded %d files from %s\n', nt, data_dir);
fprintf('Cube size: den_3d = [%d x %d x 1], denerr_3d = [%d x %d x 1]\n', ...
    size(den_3d, 1), size(den_3d, 2), size(denerr_3d, 1), size(denerr_3d, 2));
if model_cmp.available
    fprintf('Model compare loaded: %s\n', model_cmp.path);
else
    fprintf('Model compare not loaded.\n');
end

end

function fig_list = plot_each_time_profiles(time_vec, z_km, den_2d, denerr_2d, quality_flag, model_cmp)
if nargin < 6
    model_cmp = struct('available', false);
end

nt = size(den_2d, 1);
fig_list = gobjects(nt, 1);

for i = 1:nt
    t_label = 'NaT';
    if ~isnat(time_vec(i))
        t_label = char(string(time_vec(i), 'yyyy-MM-dd HH:mm:ss'));
    end

    fig_list(i) = figure('Color', 'w', 'Position', [180, 120, 760, 560], ...
        'Name', sprintf('Profile_%02d_%s', i, t_label), 'NumberTitle', 'off');

    den_i = den_2d(i, :).';
    err_pct_i = denerr_2d(i, :).';
    abs_err_i = abs(den_i) .* err_pct_i / 100;

    % Keep invalid/negative density points; only missing values remain NaN.
    h_obs = plot(den_i, z_km, '-o', 'LineWidth', 1.1, 'MarkerSize', 3, ...
        'Color', [0.1, 0.35, 0.85]);
    hold on;
    valid_err = isfinite(den_i) & isfinite(z_km) & isfinite(abs_err_i) & (abs_err_i >= 0);
    cap_half_km = 6;
    h_err = plot(nan, nan, '-', 'Color', [0.2, 0.2, 0.2], 'LineWidth', 1.0);
    for k = find(valid_err).'
        x1 = den_i(k) - abs_err_i(k);
        x2 = den_i(k) + abs_err_i(k);
        y0 = z_km(k);
        plot([x1, x2], [y0, y0], '-', 'Color', [0.2, 0.2, 0.2], 'LineWidth', 1.0);
        plot([x1, x1], [y0 - cap_half_km, y0 + cap_half_km], '-', ...
            'Color', [0.2, 0.2, 0.2], 'LineWidth', 1.0);
        plot([x2, x2], [y0 - cap_half_km, y0 + cap_half_km], '-', ...
            'Color', [0.2, 0.2, 0.2], 'LineWidth', 1.0);
    end

    has_model = model_cmp.available && i <= numel(model_cmp.model_den_interp) && ...
        any(isfinite(model_cmp.model_den_interp{i}));
    if has_model
        mdl = model_cmp.model_den_interp{i};
        dtm = model_cmp.match_abs_dt_h(i);
        h_mdl = plot(mdl, z_km, '-', 'LineWidth', 1.7, 'Color', [0.9, 0.2, 0.1]);
        legend([h_obs, h_err, h_mdl], {'Obs density', 'Obs error bar', ...
            sprintf('Model n_{he} (|\\Delta t|=%.3f h)', dtm)}, 'Location', 'best');
    else
        legend([h_obs, h_err], {'Obs density', 'Obs error bar'}, 'Location', 'best');
    end

    hold off;
    xlabel('Density [cm^{-3}]');
    ylabel('Altitude [km]');
    title(sprintf('Density with Error Bars | %s | QF=%g', t_label, quality_flag(i)));
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
    grid on;
    box on;
end
end

function model_cmp = read_model_heatmap(model_heatmap_path, obs_time_vec, obs_z_km)
model_cmp = struct();
model_cmp.available = false;
model_cmp.path = '';
model_cmp.model_ut_h = [];
model_cmp.model_z_km = [];
model_cmp.model_n_he = [];
model_cmp.match_model_idx = [];
model_cmp.match_abs_dt_h = [];
model_cmp.model_den_interp = {};

if nargin < 1 || isempty(model_heatmap_path)
    if any(~isnat(obs_time_vec))
        d0 = datestr(obs_time_vec(find(~isnat(obs_time_vec), 1, 'first')), 'yyyymmdd');
    else
        d0 = datestr(datetime('today'), 'yyyymmdd');
    end
    script_dir = fileparts(mfilename('fullpath'));
    model_heatmap_path = fullfile(script_dir, '..', 'output', sprintf('heatmap_nhe_%s.mat', d0));
end
model_cmp.path = model_heatmap_path;

if ~isfile(model_heatmap_path)
    warning('Model heatmap mat not found: %s', model_heatmap_path);
    return;
end

S = load(model_heatmap_path);
if ~isfield(S, 'ut_h_sorted') || ~isfield(S, 'z_km_grid') || ~isfield(S, 'n_he_heatmap')
    warning('Model mat missing required vars: ut_h_sorted/z_km_grid/n_he_heatmap');
    return;
end

ut_h = S.ut_h_sorted(:);
z_mdl = S.z_km_grid(:);
nhe = S.n_he_heatmap;
if size(nhe, 1) ~= numel(ut_h) || size(nhe, 2) ~= numel(z_mdl)
    warning('Model heatmap shape mismatch in %s', model_heatmap_path);
    return;
end

obs_h = hour(obs_time_vec) + minute(obs_time_vec) ./ 60 + second(obs_time_vec) ./ 3600;
match_idx = nan(numel(obs_h), 1);
match_dt = nan(numel(obs_h), 1);
model_interp = cell(numel(obs_h), 1);

for i = 1:numel(obs_h)
    if isnan(obs_h(i))
        continue;
    end
    [dmin, j] = min(abs(ut_h - obs_h(i)));
    match_idx(i) = j;
    match_dt(i) = dmin;

    mdl_profile = nhe(j, :).';
    model_interp{i} = interp1(z_mdl, mdl_profile, obs_z_km, 'linear', NaN);
end

model_cmp.available = true;
model_cmp.model_ut_h = ut_h;
model_cmp.model_z_km = z_mdl;
model_cmp.model_n_he = nhe;
model_cmp.match_model_idx = match_idx;
model_cmp.match_abs_dt_h = match_dt;
model_cmp.model_den_interp = model_interp;
end

function fig = plot_cubes(time_vec, z_km, den_2d, denerr_2d, quality_flag, peak_den, peak_z_km)
time_num = datenum(time_vec);
valid_t = isfinite(time_num);
valid_z = isfinite(z_km);

fig = struct();

if any(valid_t) && any(valid_z)
    % Density heatmap: time x altitude
    fig.den = figure('Color', 'w', 'Position', [100, 100, 1100, 450]);
    [TT, ZZ] = meshgrid(time_num(valid_t), z_km(valid_z));
    data_den = den_2d(valid_t, valid_z);
    data_den(data_den <= 0) = NaN;
    contourf(TT, ZZ, data_den.', 30, 'LineStyle', 'none');
    datetick('x', 'HH:MM', 'keeplimits');
    xlabel('Time (UT)');
    ylabel('Altitude [km]');
    title('Metastable He Number Density');
    cb = colorbar;
    cb.Label.String = 'Density [cm^{-3}]';
    colormap(turbo);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
    grid on;

    % Error heatmap: time x altitude
    fig.err = figure('Color', 'w', 'Position', [120, 120, 1100, 450]);
    [TT2, ZZ2] = meshgrid(time_num(valid_t), z_km(valid_z));
    data_err = denerr_2d(valid_t, valid_z);
    contourf(TT2, ZZ2, data_err.', 30, 'LineStyle', 'none');
    datetick('x', 'HH:MM', 'keeplimits');
    xlabel('Time (UT)');
    ylabel('Altitude [km]');
    title('Metastable He Density Relative Error');
    cb2 = colorbar;
    cb2.Label.String = 'DenEr [%]';
    colormap(turbo);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
    grid on;
end

% Peak metrics vs time
fig.peak = figure('Color', 'w', 'Position', [140, 140, 1000, 450]);
yyaxis left;
plot(time_vec, peak_den, '-o', 'LineWidth', 1.2, 'MarkerSize', 4);
set(gca, 'YScale', 'log');
ylabel('Peak Density [cm^{-3}]');

yyaxis right;
plot(time_vec, peak_z_km, '-s', 'LineWidth', 1.2, 'MarkerSize', 4);
ylabel('Peak Altitude [km]');

xlabel('Time');
title('Peak He Density and Peak Altitude');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
grid on;
box on;

% Overlay quality flag as text marks for quick screening.
bad_idx = find(quality_flag == 3);
if ~isempty(bad_idx)
    yyaxis right;
    hold on;
    plot(time_vec(bad_idx), peak_z_km(bad_idx), 'rx', 'MarkerSize', 8, 'LineWidth', 1.2);
    hold off;
end
end

function [data_time, quality_flag, z_km, den, err] = parse_one_file(fpath)
txt = fileread(fpath);
lines = regexp(txt, '\r\n|\n|\r', 'split');

data_time = NaT;
quality_flag = NaN;

for i = 1:numel(lines)
    ln = strtrim(lines{i});
    if startsWith(ln, '#DataTime:')
        s = strtrim(extractAfter(string(ln), '#DataTime:'));
        try
            data_time = datetime(s, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
        catch
            data_time = NaT;
        end
    elseif startsWith(ln, '#QualityFlag:')
        t = regexp(ln, '^\#QualityFlag:\s*(\d+)', 'tokens', 'once');
        if ~isempty(t)
            quality_flag = str2double(t{1});
        end
    end
end

start_idx = find(contains(lines, 'Elevation') & contains(lines, 'Density') & contains(lines, 'DenEr'), 1, 'first');
if isempty(start_idx)
    error('Cannot find data header in file: %s', fpath);
end

z_km = [];
den = [];
err = [];

for i = (start_idx + 1):numel(lines)
    ln = strtrim(lines{i});
    if isempty(ln) || startsWith(ln, '#')
        continue;
    end
    nums = sscanf(ln, '%f');
    if numel(nums) < 3
        continue;
    end
    z_i = nums(1);
    den_i = nums(2);
    err_i = nums(3);

    if z_i <= -9e7
        z_i = NaN;
    end
    if den_i <= -9e11
        den_i = NaN;
    end
    if err_i <= -9e9
        err_i = NaN;
    end

    z_km(end + 1, 1) = z_i; %#ok<AGROW>
    den(end + 1, 1) = den_i; %#ok<AGROW>
    err(end + 1, 1) = err_i; %#ok<AGROW>
end

valid_row = isfinite(z_km);
z_km = z_km(valid_row);
den = den(valid_row);
err = err(valid_row);
end
