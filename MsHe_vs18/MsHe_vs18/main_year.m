% main_year.m
% Batch run MsHe for a full year at fixed latitude, using only vs001 each day.

clear;
clc;

script_dir = fileparts(mfilename('fullpath'));
output_root = fullfile(script_dir, '..', '..', 'glow model', 'glow model', 'output', 'hainan_2026_sth_05');
[~, input_data_name] = fileparts(output_root);
save_root = fullfile(script_dir, '..', 'output', input_data_name);

target_year = 2026;
lon_tag = '109.00';
fixed_lat = 19.5;
vs_idx = 1;

opts = struct('make_plot', false);

assert(isfolder(output_root), 'Output folder does not exist: %s', output_root);
if ~isfolder(save_root)
    mkdir(save_root);
end

t0 = datetime(target_year, 1, 1);
t1 = datetime(target_year, 12, 31);
all_dates = (t0:t1).';
nday = numel(all_dates);

doy_list = day(all_dates, 'dayofyear');
date_tags = string(datestr(all_dates, 'yyyymmdd'));

max_n_he = nan(nday, 1);
max_alt_km = nan(nday, 1);
ok = false(nday, 1);
ut_s = nan(nday, 1);
ut_h = nan(nday, 1);
z_km = [];
n_he_matrix = [];

lat_tag = sprintf('%.2f', fixed_lat);

for i = 1:nday
    doy_folder = sprintf('%03d', doy_list(i));
    date_tag = char(date_tags(i));
    site_tag = sprintf('Lat%s_Lon%s_%s', lat_tag, lon_tag, date_tag);
    file_name = sprintf('%s_vs%03d.txt', site_tag, vs_idx);
    glow_file = fullfile(output_root, site_tag, doy_folder, file_name);

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
            n_he_matrix = nan(nday, numel(z_km));
        end
        n_he_matrix(i, :) = r.profile.n_He_cm3(:).';

        fprintf('DoY=%03d, date=%s, vs=%03d: UT=%7.1f s (%6.3f h), max_n_he=% .6e cm^-3, alt=%6.1f km\n', ...
            doy_list(i), date_tag, vs_idx, ut_s(i), ut_h(i), max_n_he(i), max_alt_km(i));
    catch ME
        fprintf('DoY=%03d, date=%s, vs=%03d failed: %s\n', doy_list(i), date_tag, vs_idx, ME.message);
    end
end

results = table(doy_list(:), date_tags(:), ok, ut_s, ut_h, max_n_he, max_alt_km, ...
    'VariableNames', {'doy', 'date_tag', 'ok', 'ut_s', 'ut_h', 'n_he_max_cm3', 'n_he_max_alt_km'});

fprintf('\nFinished: %d/%d succeeded.\n', nnz(ok), nday);
disp(results);

if ~isempty(z_km)
    results_csv = fullfile(save_root, sprintf('mshe_year_results_%d_vs%03d.csv', target_year, vs_idx));
    writetable(results, results_csv);
    fprintf('Saved summary: %s\n', results_csv);

    % Wide output: one row per day, one column per altitude.
    z_names = compose('z_%dkm', round(z_km));
    n_he_table = array2table(n_he_matrix, 'VariableNames', z_names);
    n_he_table = addvars(n_he_table, doy_list(:), date_tags(:), ok, ut_s, ut_h, 'Before', 1, ...
        'NewVariableNames', {'doy', 'date_tag', 'ok', 'ut_s', 'ut_h'});

    % Long output: (doy, date_tag, z_km, n_He_cm3)
    [doy_grid, z_grid] = ndgrid(doy_list(:), z_km(:));
    [ut_h_grid, ~] = ndgrid(ut_h, z_km(:));
    date_col = repmat(date_tags(:), numel(z_km), 1);
    n_he_profile = table(doy_grid(:), date_col, ut_h_grid(:), z_grid(:), n_he_matrix(:), ...
        'VariableNames', {'doy', 'date_tag', 'ut_h', 'z_km', 'n_He_cm3'});

    fprintf('\nCreated full-profile outputs:\n');
    fprintf('1) n_he_table   (wide: day x altitude)\n');
    fprintf('2) n_he_profile (long: doy, date_tag, ut_h, z_km, n_He_cm3)\n');

    data_mat = fullfile(save_root, sprintf('mshe_year_full_%d_vs%03d.mat', target_year, vs_idx));
    save(data_mat, 'target_year', 'fixed_lat', 'lon_tag', 'vs_idx', ...
        'doy_list', 'date_tags', 'ok', 'ut_s', 'ut_h', 'max_n_he', 'max_alt_km', ...
        'z_km', 'n_he_matrix', 'results', 'n_he_table', 'n_he_profile');
    fprintf('Saved data: %s\n', data_mat);

    valid_idx = find(ok);
    doy_valid = doy_list(valid_idx);
    n_he_valid = n_he_matrix(valid_idx, :);

    % Sort by DoY for smooth plotting.
    [doy_sorted, ord] = sort(doy_valid);
    n_he_sorted = n_he_valid(ord, :);

    % Plot 1: DoY-altitude heatmap of n_he.
    n_he_plot = n_he_sorted;
    n_he_plot(~isfinite(n_he_plot)) = nan;

    figure('Color', 'w', 'Position', [120, 120, 1100, 450]);
    [DOY, Z] = meshgrid(doy_sorted, z_km);
    contourf(DOY, Z, n_he_plot.', 30, 'LineStyle', 'none');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
    xlabel('Day of Year');
    ylabel('Altitude [km]');
    title(sprintf('He(2^3S) density at lat=%.2f, year=%d, vs=%03d', fixed_lat, target_year, vs_idx));
    cb = colorbar;
    cb.Label.String = 'n_{he} [cm^{-3}]';
    colormap(turbo);
    %caxis([0, 2]);
    grid on;
    heat_png = fullfile(save_root, sprintf('year_heatmap_nhe_%d_vs%03d.png', target_year, vs_idx));
    heat_fig = fullfile(save_root, sprintf('year_heatmap_nhe_%d_vs%03d.fig', target_year, vs_idx));
    saveas(gcf, heat_png);
    savefig(gcf, heat_fig);
    fprintf('Saved figure: %s\n', heat_png);
    fprintf('Saved figure: %s\n', heat_fig);

    % Plot 2: annual evolution of peak density and peak altitude.
    figure('Color', 'w', 'Position', [150, 150, 1000, 450]);
    yyaxis left;
    plot(doy_list(ok), max_n_he(ok), '-o', 'LineWidth', 1.2, 'MarkerSize', 3);
    ylabel('Peak n_{he} [cm^{-3}]');
    set(gca, 'YScale', 'log');

    yyaxis right;
    plot(doy_list(ok), max_alt_km(ok), '-s', 'LineWidth', 1.2, 'MarkerSize', 3);
    ylabel('Peak altitude [km]');

    xlabel('Day of Year');
    title(sprintf('Annual peak He(2^3S) metrics at lat=%.2f, vs=%03d', fixed_lat, vs_idx));
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
    grid on;
    box on;
    peak_png = fullfile(save_root, sprintf('year_peak_metrics_%d_vs%03d.png', target_year, vs_idx));
    peak_fig = fullfile(save_root, sprintf('year_peak_metrics_%d_vs%03d.fig', target_year, vs_idx));
    saveas(gcf, peak_png);
    savefig(gcf, peak_fig);
    fprintf('Saved figure: %s\n', peak_png);
    fprintf('Saved figure: %s\n', peak_fig);
end

