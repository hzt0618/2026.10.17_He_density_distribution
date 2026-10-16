% test_mshe_param_heatmaps_day.m
% Build parameter heatmaps for one day using the same data source logic as main_day.m.

clear;
clc;

script_dir = fileparts(mfilename('fullpath'));
output_base = fullfile(script_dir, '..', '..', 'glow model', 'glow model', 'output', 'north');

% Only edit date_tag for daily runs; doy_folder is computed automatically.
date_tag = '20260119';
dt = datetime(date_tag, 'InputFormat', 'yyyyMMdd');
doy_folder = sprintf('%03d', day(dt, 'dayofyear'));

% Auto-discover the date folder: Lat*_Lon*_<date_tag>
date_dirs = dir(fullfile(output_base, '**', sprintf('Lat*_Lon*_%s', date_tag)));
date_dirs = date_dirs([date_dirs.isdir]);
if isempty(date_dirs)
    error('No data folder found for date %s under: %s', date_tag, output_base);
end
if numel(date_dirs) > 1
    error('Multiple data folders found for date %s; please keep only one match under: %s', date_tag, output_base);
end

output_root = fullfile(date_dirs(1).folder, date_dirs(1).name);
site_tag = date_dirs(1).name;
tok_site = regexp(site_tag, '^Lat([-+]?\d*\.?\d+)_Lon([-+]?\d*\.?\d+)_(\d{8})$', 'tokens', 'once');
if isempty(tok_site)
    error('Folder name does not match expected pattern: %s', site_tag);
end
fixed_lat = str2double(tok_site{1});
lon_tag = tok_site{2};
source_tag = string(split(output_base, filesep));
source_tag = char(source_tag(end));

% 5-minute cadence for a full day
vs_list = 1:(24 * 60 / 5);

opts = struct('make_plot', false);

n_vs = numel(vs_list);
ok = false(n_vs, 1);
ut_h = nan(n_vs, 1);

z_km = [];

% Core controlling terms from MsHe pipeline
qe_matrix = [];
qe_per_he_matrix = [];
qr_matrix = [];
q_matrix = [];
loss_matrix = [];
dif_matrix = [];
hn_matrix = [];
he_matrix = [];
nhe_matrix = [];

lat_tag = sprintf('%.2f', fixed_lat);

for i = 1:n_vs
    vs_idx = vs_list(i);
    file_name = sprintf('%s_vs%03d.txt', site_tag, vs_idx);
    glow_file = fullfile(output_root, doy_folder, file_name);

    if ~isfile(glow_file)
        fprintf('Missing file (skip): %s\n', glow_file);
        continue;
    end

    try
        r = MsHe(glow_file, opts);

        if isempty(z_km)
            z_km = r.grid.z_km(:);
            nz = numel(z_km);

            qe_matrix = nan(n_vs, nz);
            qe_per_he_matrix = nan(n_vs, nz);
            qr_matrix = nan(n_vs, nz);
            q_matrix = nan(n_vs, nz);
            loss_matrix = nan(n_vs, nz);
            dif_matrix = nan(n_vs, nz);
            hn_matrix = nan(n_vs, nz);
            he_matrix = nan(n_vs, nz);
            nhe_matrix = nan(n_vs, nz);
        end

        ut_h(i) = r.meta.ut_s / 3600;

        qe_matrix(i, :) = r.profile.qe(:).';
        qe_per_he = nan(size(r.profile.qe(:)));
        he_valid = isfinite(r.raw.He(:)) & (r.raw.He(:) > 0);
        qe_per_he(he_valid) = r.profile.qe(he_valid) ./ r.raw.He(he_valid);
        qe_per_he_matrix(i, :) = qe_per_he(:).';
        qr_matrix(i, :) = r.profile.qr(:).';
        q_matrix(i, :) = r.profile.q(:).';
        loss_matrix(i, :) = r.profile.Loss(:).';
        dif_matrix(i, :) = r.profile.Dif(:).';
        hn_matrix(i, :) = r.profile.Hn_cm(:).';
        he_matrix(i, :) = r.raw.He(:).';
        nhe_matrix(i, :) = r.profile.n_He_cm3(:).';

        ok(i) = true;
        fprintf('vs=%03d ok, UT=%.3f h\n', vs_idx, ut_h(i));
    catch ME
        fprintf('vs=%03d failed: %s\n', vs_idx, ME.message);
    end
end

fprintf('\nFinished: %d/%d succeeded.\n', nnz(ok), n_vs);

if isempty(z_km)
    error('No valid files were processed.');
end

valid_idx = find(ok);
ut_valid = ut_h(valid_idx);

[ut_sorted, ord] = sort(ut_valid);

qe_plot = qe_matrix(valid_idx, :);
qe_per_he_plot = qe_per_he_matrix(valid_idx, :);
qr_plot = qr_matrix(valid_idx, :);
q_plot = q_matrix(valid_idx, :);
loss_plot = loss_matrix(valid_idx, :);
dif_plot = dif_matrix(valid_idx, :);
hn_plot = hn_matrix(valid_idx, :);
he_plot = he_matrix(valid_idx, :);
nhe_plot = nhe_matrix(valid_idx, :);

qe_plot = qe_plot(ord, :);
qe_per_he_plot = qe_per_he_plot(ord, :);
qr_plot = qr_plot(ord, :);
q_plot = q_plot(ord, :);
loss_plot = loss_plot(ord, :);
dif_plot = dif_plot(ord, :);
hn_plot = hn_plot(ord, :);
he_plot = he_plot(ord, :);
nhe_plot = nhe_plot(ord, :);

out_root = fullfile(script_dir, '..', 'output', source_tag);
out_dir = fullfile(out_root, date_tag);
if ~isfolder(out_dir)
    mkdir(out_dir);
end

% Save full parameter data
param_mat = fullfile(out_dir, sprintf('mshe_params_full_%s.mat', date_tag));
save(param_mat, 'date_tag', 'doy_folder', 'fixed_lat', 'lon_tag', ...
    'vs_list', 'ok', 'ut_h', 'ut_sorted', 'z_km', ...
    'qe_matrix', 'qe_per_he_matrix', 'qr_matrix', 'q_matrix', 'loss_matrix', 'dif_matrix', 'hn_matrix', 'he_matrix', 'nhe_matrix', ...
    'qe_plot', 'qe_per_he_plot', 'qr_plot', 'q_plot', 'loss_plot', 'dif_plot', 'hn_plot', 'he_plot', 'nhe_plot');
fprintf('Saved data: %s\n', param_mat);

plot_param_heatmap(ut_sorted, z_km, qe_plot, ...
    'qe', 'cm^{-3} s^{-1}', fixed_lat, date_tag, out_dir, false, 300);
plot_param_heatmap(ut_sorted, z_km, qe_per_he_plot, ...
    'qe/He', 's^{-1}', fixed_lat, date_tag, out_dir, false, 300);
plot_param_heatmap(ut_sorted, z_km, qr_plot, ...
    'qr', 'cm^{-3} s^{-1}', fixed_lat, date_tag, out_dir, false, 300);
plot_param_heatmap(ut_sorted, z_km, q_plot, ...
    'q=qe+qr', 'cm^{-3} s^{-1}', fixed_lat, date_tag, out_dir, false, 300);
% plot_param_heatmap(ut_sorted, z_km, loss_plot, ...
%     'Loss', 's^{-1}', fixed_lat, date_tag, out_dir, true, NaN);
% plot_param_heatmap(ut_sorted, z_km, dif_plot, ...
%     'Dif', 'cm^{2} s^{-1}', fixed_lat, date_tag, out_dir, true, NaN);
% plot_param_heatmap(ut_sorted, z_km, hn_plot, ...
%     'Hn', 'cm', fixed_lat, date_tag, out_dir, true, NaN);

plot_param_heatmap(ut_sorted, z_km, he_plot, ...
    'He', 'cm^{-3}', fixed_lat, date_tag, out_dir, true, 200);

% Final result heatmap: metastable helium density (n_he)
% Color limits are determined only from z >= 300 km values.
plot_param_heatmap(ut_sorted, z_km, nhe_plot, ...
    'n_he', 'cm^{-3}', fixed_lat, date_tag, out_dir, false, 300);

fprintf('All heatmaps saved in: %s\n', out_dir);


function plot_param_heatmap(ut_sorted, z_km, value_matrix, tag, unit_txt, fixed_lat, date_tag, out_dir, use_log_scale, clim_z_min_km)
    fig = figure('Color', 'w', 'Position', [120, 120, 1100, 450]);

    plot_matrix = value_matrix;
    color_label = sprintf('%s [%s]', tag, unit_txt);
    title_suffix = tag;
    if strcmp(tag, 'n_he')
        title_suffix = 'n_{metahe}';
    end

    if use_log_scale
        valid = isfinite(plot_matrix) & (plot_matrix > 0);
        tmp = nan(size(plot_matrix));
        tmp(valid) = log10(plot_matrix(valid));
        plot_matrix = tmp;
        color_label = sprintf('log10(%s [%s])', tag, unit_txt);
        title_suffix = sprintf('log10(%s)', tag);
    end

    [UT, Z] = meshgrid(ut_sorted, z_km);
    contourf(UT, Z, plot_matrix.', 30, 'LineStyle', 'none');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
    xlabel('UT [hour]');
    ylabel('Altitude [km]');
    title(sprintf('%s heatmap at lat=%.2f, date %s', title_suffix, fixed_lat, date_tag));
    cb = colorbar;
    cb.Label.String = color_label;
    colormap(turbo);

    if strcmp(tag, 'He')
        ylim([200, max(z_km)]);
    end

    if ~isfinite(clim_z_min_km)
        vals = plot_matrix(isfinite(plot_matrix));
    else
        zmask = z_km >= clim_z_min_km;
        vals = plot_matrix(:, zmask);
        vals = vals(isfinite(vals));
    end

    if isempty(vals)
        caxis([0, 1]);
    else
        if strcmp(tag, 'n_he') || strcmp(tag, 'qe') || strcmp(tag, 'qe/He') || strcmp(tag, 'qr') || strcmp(tag, 'q=qe+qr')
            c_lo = prctile(vals, 12);
            c_hi = prctile(vals, 99.5);
        elseif strcmp(tag, 'He')
            % Use robust percentiles so the He heatmap is not dominated by a few extreme values.
            c_lo = prctile(vals, 5);
            c_hi = prctile(vals, 99.5);
        else
            c_lo = prctile(vals, 2);
            c_hi = prctile(vals, 98);
        end
        if c_hi <= c_lo
            c_lo = min(vals);
            c_hi = max(vals);
        end
        if c_hi <= c_lo
            c_hi = c_lo + eps(c_lo + 1);
        end
        caxis([c_lo, c_hi]);
    end

    grid on;

    safe_tag = regexprep(tag, '[^a-zA-Z0-9_]', '_');
    png_path = fullfile(out_dir, sprintf('heatmap_%s_%s.png', safe_tag, date_tag));
    fig_path = fullfile(out_dir, sprintf('heatmap_%s_%s.fig', safe_tag, date_tag));
    saveas(fig, png_path);
    savefig(fig, fig_path);
    fprintf('Saved figure: %s\n', png_path);
    fprintf('Saved figure: %s\n', fig_path);
end















