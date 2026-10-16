% main_day.m
% Batch run MsHe at fixed latitude, scanning UT via vs files.

clear;
clc;
 
script_dir = fileparts(mfilename('fullpath'));
run_tag = 'hainna_20251010_sth_05';
output_base = fullfile(script_dir, '..', '..', 'glow model', 'glow model', 'output', run_tag);

% Only edit date_tag for daily runs; doy_folder is computed automatically.
date_tag = '20251010';
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
save_root = fullfile(script_dir, '..', 'output', run_tag);
site_tag = date_dirs(1).name;
tok_site = regexp(site_tag, '^Lat([-+]?\d*\.?\d+)_Lon([-+]?\d*\.?\d+)_(\d{8})$', 'tokens', 'once');
if isempty(tok_site)
    error('Folder name does not match expected pattern: %s', site_tag);
end
fixed_lat = str2double(tok_site{1});
lon_tag = tok_site{2};

if ~isfolder(save_root)
    mkdir(save_root);
end

% 5-minute cadence for a full day: 24*60/5 = 288 samples.
vs_list = 1:(24 * 60 / 5);

opts = struct('make_plot', false);
nut = numel(vs_list);
max_n_he = nan(nut, 1);
max_alt_km = nan(nut, 1);
ok = false(nut, 1);
ut_s = nan(nut, 1);
ut_h = nan(nut, 1);
z_km = [];
n_he_matrix = [];

lat_tag = sprintf('%.2f', fixed_lat);

for i = 1:nut
    vs_idx = vs_list(i);
    file_name = sprintf('%s_vs%03d.txt', site_tag, vs_idx);
    glow_file = fullfile(output_root, doy_folder, file_name);

    if ~isfile(glow_file)
        fprintf('Missing file (skip): %s\n', glow_file);
        continue;
    end

    try
        r = MsHe(glow_file, opts);
        required_top_fields = {'input', 'meta', 'grid', 'raw', 'profile'};
        for k = 1:numel(required_top_fields)
            assert(isfield(r, required_top_fields{k}), 'Missing field: r.%s', required_top_fields{k});
        end
        assert(isfield(r.profile, 'n_He_cm3'), 'Missing field: r.profile.n_He_cm3');
        assert(numel(r.profile.n_He_cm3) == r.grid.nz, 'n_he length does not match grid size');
        assert(all(isfinite(r.profile.n_He_cm3)), 'n_he contains non-finite values');

        max_n_he(i) = r.profile.n_he_max_cm3;
        max_alt_km(i) = r.profile.n_he_max_alt_km;
        ut_s(i) = r.meta.ut_s;
        ut_h(i) = r.meta.ut_s / 3600;
        ok(i) = true;

        if isempty(z_km)
            z_km = r.grid.z_km(:);
            n_he_matrix = nan(nut, numel(z_km));
        end
        n_he_matrix(i, :) = r.profile.n_He_cm3(:).';

        fprintf('vs=%03d: UT=%7.1f s (%6.3f h), max_n_he=% .6e cm^-3, alt=%6.1f km\n', ...
            vs_idx, ut_s(i), ut_h(i), max_n_he(i), max_alt_km(i));
    catch ME
        fprintf('vs=%03d failed: %s\n', vs_idx, ME.message);
    end
end

results = table(vs_list(:), ok, ut_s, ut_h, max_n_he, max_alt_km, ...
    'VariableNames', {'vs_idx', 'ok', 'ut_s', 'ut_h', 'n_he_max_cm3', 'n_he_max_alt_km'});

fprintf('\nFinished: %d/%d succeeded.\n', nnz(ok), nut);
disp(results);

if ~isempty(z_km)
    % Wide output: one row per UT sample (vs index), one column per altitude.
    z_names = compose('z_%dkm', round(z_km));
    n_he_table = array2table(n_he_matrix, 'VariableNames', z_names);
    n_he_table = addvars(n_he_table, vs_list(:), ok, ut_s, ut_h, 'Before', 1, ...
        'NewVariableNames', {'vs_idx', 'ok', 'ut_s', 'ut_h'});

    % Long output: (vs_idx, ut_h, z_km, n_He_cm3)
    [vs_grid, z_grid] = ndgrid(vs_list(:), z_km(:));
    [ut_h_grid, ~] = ndgrid(ut_h, z_km(:));
    n_he_profile = table(vs_grid(:), ut_h_grid(:), z_grid(:), n_he_matrix(:), ...
        'VariableNames', {'vs_idx', 'ut_h', 'z_km', 'n_He_cm3'});

    fprintf('\nCreated full-profile outputs:\n');
    fprintf('1) n_he_table   (wide: UT sample x altitude)\n');
    fprintf('2) n_he_profile (long: vs_idx, ut_h, z_km, n_He_cm3)\n');

    valid_idx = find(ok);
    ut_valid = ut_h(valid_idx);
    n_he_valid = n_he_matrix(valid_idx, :);

    % Sort by UT for smooth plotting.
    [ut_valid_sorted, ord] = sort(ut_valid);
    n_he_valid_sorted = n_he_valid(ord, :);

    % Plot 1: UT-altitude heatmap of n_he.
    n_he_plot = n_he_valid_sorted;
    n_he_plot(~isfinite(n_he_plot)) = nan;

    % Match the Optica single-column final size so text stays near 8--9 pt.
    fig_font = 'Times New Roman';
    fig_label_font_pt = 9;
    fig_tick_font_pt = 8;
    fig_line_width = 0.6;
    fig1_width_in = 3.35;
    fig1_height_in = 1.55;
    fig1 = figure('Color', 'w', ...
        'Units', 'inches', ...
        'Position', [1, 1, fig1_width_in, fig1_height_in], ...
        'PaperUnits', 'inches', ...
        'PaperPosition', [0, 0, fig1_width_in, fig1_height_in], ...
        'PaperSize', [fig1_width_in, fig1_height_in]);
    [UT, Z] = meshgrid(ut_valid_sorted, z_km);
    contourf(UT, Z, n_he_plot.', 30, 'LineStyle', 'none');
    ax1 = gca;
    set(ax1, 'FontName', fig_font, 'FontSize', fig_tick_font_pt, ...
        'LineWidth', fig_line_width, 'TickDir', 'out');
    xlabel('UT [hour]', 'FontName', fig_font, 'FontSize', fig_label_font_pt);
    ylabel('Altitude [km]', 'FontName', fig_font, 'FontSize', fig_label_font_pt);
    cb = colorbar;
    cb.Label.String = 'n_{he} [cm^{-3}]';
    cb.FontName = fig_font;
    cb.FontSize = fig_tick_font_pt;
    cb.Label.FontName = fig_font;
    cb.Label.FontSize = fig_label_font_pt;
    colormap(turbo);

    % Auto color range from robust percentiles of valid density values.
    valid_vals = n_he_plot(isfinite(n_he_plot) & n_he_plot > 0);
    if isempty(valid_vals)
        caxis([0, 1]);
    else
        clim_low = prctile(valid_vals, 2);
        clim_high = prctile(valid_vals, 98);
        if clim_high <= clim_low
            clim_low = min(valid_vals);
            clim_high = max(valid_vals);
        end
        if clim_high <= clim_low
            clim_high = clim_low + eps(clim_low + 1);
        end
        caxis([clim_low, clim_high]);
    end

    grid off;

    out_dir = save_root;

    % Save complete outputs (summary + full density-height data).
    results_csv = fullfile(out_dir, sprintf('%s_results.csv', run_tag));
    nhe_table_csv = fullfile(out_dir, sprintf('%s_nhe_table.csv', run_tag));
    nhe_profile_csv = fullfile(out_dir, sprintf('%s_nhe_profile.csv', run_tag));
    full_data_mat = fullfile(out_dir, sprintf('%s_nhe_full.mat', run_tag));
    writetable(results, results_csv);
    writetable(n_he_table, nhe_table_csv);
    writetable(n_he_profile, nhe_profile_csv);
    save(full_data_mat, 'results', 'n_he_table', 'n_he_profile', ...
        'z_km', 'n_he_matrix', 'ut_h', 'ut_s', 'vs_list', ...
        'date_tag', 'doy_folder', 'lon_tag', 'fixed_lat', 'run_tag');
    fprintf('Saved data: %s\n', results_csv);
    fprintf('Saved data: %s\n', nhe_table_csv);
    fprintf('Saved data: %s\n', nhe_profile_csv);
    fprintf('Saved data: %s\n', full_data_mat);

    % Save Plot 1 source data for reuse (UT-altitude heatmap).
    heatmap_mat = fullfile(out_dir, sprintf('%s_heatmap_nhe.mat', run_tag));
    ut_h_sorted = ut_valid_sorted; %#ok<NASGU>
    z_km_grid = z_km; %#ok<NASGU>
    n_he_heatmap = n_he_plot; %#ok<NASGU>
    save(heatmap_mat, 'UT', 'Z', 'n_he_heatmap', 'ut_h_sorted', 'z_km_grid', ...
        'doy_folder', 'date_tag', 'lon_tag', 'fixed_lat', 'run_tag');
    fprintf('Saved Plot 1 data: %s\n', heatmap_mat);

    fig1_png = fullfile(out_dir, sprintf('%s_fig1_heatmap.png', run_tag));
    fig1_pdf = fullfile(out_dir, sprintf('%s_fig1_heatmap.pdf', run_tag));
    fig1_fig = fullfile(out_dir, sprintf('%s_fig1_heatmap.fig', run_tag));
    exportgraphics(fig1, fig1_png, 'Resolution', 600);
    exportgraphics(fig1, fig1_pdf, 'ContentType', 'vector');
    savefig(fig1, fig1_fig);
    fprintf('Saved figure: %s\n', fig1_png);
    fprintf('Saved figure: %s\n', fig1_pdf);
    fprintf('Saved figure: %s\n', fig1_fig);

    % Plot 2: selected UT profiles.
    selected_ut = [min(ut_valid_sorted), mean(ut_valid_sorted), max(ut_valid_sorted)];
    fig2 = figure('Color', 'w', 'Position', [150, 150, 600, 700]);
    hold on;
    for j = 1:numel(selected_ut)
        [~, near_idx] = min(abs(ut_valid_sorted - selected_ut(j)));
        plot(n_he_valid_sorted(near_idx, :), z_km, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('UT=%.3f h', ut_valid_sorted(near_idx)));
    end
    set(gca, 'XScale', 'log', 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
    xlabel('n_{he} [cm^{-3}]');
    ylabel('Altitude [km]');
    title(sprintf('He(2^3S) vertical profiles at lat=%.2f', fixed_lat));
    legend('Location', 'best');
    grid on;
    box on;
    hold off;

    fig2_png = fullfile(out_dir, sprintf('%s_fig2_profiles.png', run_tag));
    fig2_fig = fullfile(out_dir, sprintf('%s_fig2_profiles.fig', run_tag));
    saveas(fig2, fig2_png);
    savefig(fig2, fig2_fig);
    fprintf('Saved figure: %s\n', fig2_png);
    fprintf('Saved figure: %s\n', fig2_fig);
end



