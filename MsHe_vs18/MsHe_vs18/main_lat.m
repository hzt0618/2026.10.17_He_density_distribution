% test_MsHe.m
% Batch smoke test for MsHe(glow_file, opts) from lat -60 to lat 60.

clear;
clc;

script_dir = fileparts(mfilename('fullpath'));
output_root = fullfile(script_dir, '..', '..', 'glow model', 'glow model', 'output');

doy_folder = '282';
date_tag = '20251009';
lon_tag = '109.00';
vs_idx = 5;
lat_list = -60:0.5:60;

opts = struct('make_plot', false);
nlat = numel(lat_list);
max_n_he = nan(nlat, 1);
max_alt_km = nan(nlat, 1);
ok = false(nlat, 1);
z_km = [];
n_he_matrix = [];

for i = 1:nlat
    lat = lat_list(i);
    lat_tag = sprintf('%.2f', lat);  % e.g., -60.00, 0.00, 60.00
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
        ok(i) = true;

        if isempty(z_km)
            z_km = r.grid.z_km(:);
            n_he_matrix = nan(nlat, numel(z_km));
        end
        n_he_matrix(i, :) = r.profile.n_He_cm3(:).';

        fprintf('lat=%6.2f: max_n_he=% .6e cm^-3, alt=%6.1f km\n', lat, max_n_he(i), max_alt_km(i));
    catch ME
        fprintf('lat=%6.2f failed: %s\n', lat, ME.message);
    end
end

results = table(lat_list(:), ok, max_n_he, max_alt_km, ...
    'VariableNames', {'lat', 'ok', 'n_he_max_cm3', 'n_he_max_alt_km'});

fprintf('\nFinished: %d/%d succeeded.\n', nnz(ok), nlat);
disp(results);

if ~isempty(z_km)
    % Wide output: one row per latitude, one column per altitude.
    z_names = compose('z_%dkm', round(z_km));
    n_he_table = array2table(n_he_matrix, 'VariableNames', z_names);
    n_he_table = addvars(n_he_table, lat_list(:), ok, 'Before', 1, ...
        'NewVariableNames', {'lat', 'ok'});

    % Long output: (lat, z_km, n_He_cm3), convenient for plotting/analysis.
    [lat_grid, z_grid] = ndgrid(lat_list(:), z_km(:));
    n_he_profile = table(lat_grid(:), z_grid(:), n_he_matrix(:), ...
        'VariableNames', {'lat', 'z_km', 'n_He_cm3'});

    fprintf('\nCreated full-profile outputs:\n');
    fprintf('1) n_he_table   (wide: latitude x altitude)\n');
    fprintf('2) n_he_profile (long: lat, z_km, n_He_cm3)\n');

    valid_idx = find(ok);
    lat_valid = lat_list(valid_idx);
    n_he_valid = n_he_matrix(valid_idx, :);

    % Plot 1: latitude-altitude heatmap of n_he
    n_he_plot = n_he_valid;
    n_he_plot(~isfinite(n_he_plot)) = nan;

    figure('Color', 'w', 'Position', [120, 120, 1100, 450]);
    [LAT, Z] = meshgrid(lat_valid, z_km);
    contourf(LAT, Z, n_he_plot.', 30, 'LineStyle', 'none');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
    xlabel('Latitude [deg]');
    ylabel('Altitude [km]');
    title(sprintf('He(2^3S) density, vs%03d, DoY %s', vs_idx, doy_folder));
    cb = colorbar;
    cb.Label.String = 'n_{he} [cm^{-3}]';
    colormap(turbo);
    caxis([0, 1]);
    grid on;

    % Plot 2: selected latitude profiles
    selected_lat = [-60, -30, 0, 30, 60];
    figure('Color', 'w', 'Position', [150, 150, 600, 700]);
    hold on;
    for j = 1:numel(selected_lat)
        [~, near_idx] = min(abs(lat_valid - selected_lat(j)));
        plot(n_he_valid(near_idx, :), z_km, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('lat=%.1f', lat_valid(near_idx)));
    end
    set(gca, 'XScale', 'log', 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
    xlabel('n_{he} [cm^{-3}]');
    ylabel('Altitude [km]');
    title(sprintf('He(2^3S) vertical profiles, vs%03d', vs_idx));
    legend('Location', 'best');
    grid on;
    box on;
    hold off;
end

