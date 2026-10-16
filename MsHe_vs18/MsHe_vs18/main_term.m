% main_term.m
% Unified MsHe post-processing script.
% Set plot_mode to:
%   'sphere'    -> spherical Cartesian 3D view
%   'lonlatalt' -> longitude-latitude-altitude 3D view

clear;
clc;

script_dir = fileparts(mfilename('fullpath'));

date_tag = '20250320';
run_suffix = 'mor_sth_05';
doy_folder = date_to_doy_folder(date_tag);
output_root = fullfile(script_dir, '..', '..', 'glow model', 'glow model', 'output', sprintf('%s_%s', date_tag, run_suffix));
conj_file = fullfile(script_dir, '..', '..', 'IGRF', 'output', date_tag, 'conjugate_morning_filled.txt');
vs_idx = 6;
file_tag = sprintf('vs%03d', vs_idx);
earth_radius_km = 1000;

plot_mode = 'sphere';   % 'sphere' or 'lonlatalt'
save_outputs = true;       % set false for quick look only

opts = struct('make_plot', false);

points = load_conjugate_points(conj_file);
npt = numel(points.idx);

max_n_he = nan(npt, 1);
max_alt_km = nan(npt, 1);
sza_selected_rad = nan(npt, 1);
sza_selected_deg = nan(npt, 1);
ut_selected_s = nan(npt, 1);
ut_selected_h = nan(npt, 1);
selected_vs_idx = nan(npt, 1);
selected_file = strings(npt, 1);
ok = false(npt, 1);
z_km = [];
n_he_matrix = [];
first_data_file = '';

qe_matrix = [];
qe_per_he_matrix = [];
qr_matrix = [];
q_matrix = [];
loss_matrix = [];
dif_matrix = [];
hn_matrix = [];
he_matrix = [];

for i = 1:npt
    lat = points.input_lat(i);
    lon = points.input_lon(i);

    lat_tag = sprintf('%.2f', lat);
    lon_tag = sprintf('%.2f', lon);
    site_tag = sprintf('Lat%s_Lon%s_%s', lat_tag, lon_tag, date_tag);
    file_name = sprintf('%s_vs%03d.txt', site_tag, vs_idx);
    glow_file = fullfile(output_root, site_tag, doy_folder, file_name);

    if ~isfile(glow_file)
        fprintf('Missing file (skip): idx=%d, %s\n', points.idx(i), glow_file);
        continue;
    end

    [sza_now, ut_now] = read_sza_ut_from_glow_file(glow_file);
    if isfinite(sza_now)
        sza_selected_rad(i) = sza_now;
        sza_selected_deg(i) = sza_now * 180 / pi;
    end
    if isfinite(ut_now)
        ut_selected_s(i) = ut_now;
        ut_selected_h(i) = ut_now / 3600;
    end
    selected_vs_idx(i) = vs_idx;
    selected_file(i) = string(get_file_name_only(glow_file));

    if isempty(first_data_file)
        first_data_file = glow_file;
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
            nz = numel(z_km);
            n_he_matrix = nan(npt, nz);
            qe_matrix = nan(npt, nz);
            qe_per_he_matrix = nan(npt, nz);
            qr_matrix = nan(npt, nz);
            q_matrix = nan(npt, nz);
            loss_matrix = nan(npt, nz);
            dif_matrix = nan(npt, nz);
            hn_matrix = nan(npt, nz);
            he_matrix = nan(npt, nz);
        end

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
        n_he_matrix(i, :) = r.profile.n_He_cm3(:).';

        fprintf('idx=%3d lat=%8.3f lon=%9.3f: max_n_he=% .6e cm^-3, alt=%6.1f km, file=%s\n', ...
            points.idx(i), lat, lon, max_n_he(i), max_alt_km(i), get_file_name_only(glow_file));
    catch ME
        fprintf('idx=%3d lat=%8.3f lon=%9.3f failed: %s\n', points.idx(i), lat, lon, ME.message);
    end
end

results = table(points.idx(:), points.input_lat(:), points.input_lon(:), ...
    points.input_x(:), points.input_y(:), points.input_z(:), ...
    ok, sza_selected_rad, sza_selected_deg, ut_selected_s, ut_selected_h, selected_vs_idx, ...
    selected_file, max_n_he, max_alt_km, ...
    'VariableNames', {'idx', 'input_lat', 'input_lon', 'input_x', 'input_y', 'input_z', ...
    'ok', 'sza_rad', 'sza_deg', 'ut_s', 'ut_h', 'vs_idx', 'selected_file', 'n_he_max_cm3', 'n_he_max_alt_km'});

fprintf('\nFinished: %d/%d succeeded.\n', nnz(ok), npt);
disp(results);

if isempty(z_km)
    return;
end

title_str = make_ut_title_from_first_file(first_data_file, date_tag);
valid_idx = find(ok);
lon_valid = points.input_lon(valid_idx);
lat_valid = points.input_lat(valid_idx);
[mlat_valid, mlon_valid] = geo_to_centered_dipole(points.input_lat(valid_idx), points.input_lon(valid_idx));
n_he_valid = n_he_matrix(valid_idx, :);

out_dir = '';
if save_outputs
    [~, input_data_name] = fileparts(output_root);
    out_dir = fullfile(script_dir, '..', 'output', input_data_name);
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    z_names = compose('z_%dkm', round(z_km));
    n_he_table = array2table(n_he_matrix, 'VariableNames', z_names);
    n_he_table = addvars(n_he_table, points.idx(:), points.input_lat(:), points.input_lon(:), ...
        points.input_x(:), points.input_y(:), points.input_z(:), ok, 'Before', 1, ...
        'NewVariableNames', {'idx', 'input_lat', 'input_lon', 'input_x', 'input_y', 'input_z', 'ok'});

    [idx_grid, z_grid] = ndgrid(points.idx(:), z_km(:));
    [lat_grid, ~] = ndgrid(points.input_lat(:), z_km(:));
    [lon_grid, ~] = ndgrid(points.input_lon(:), z_km(:));
    [x_grid, ~] = ndgrid(points.input_x(:), z_km(:));
    [y_grid, ~] = ndgrid(points.input_y(:), z_km(:));
    [zxyz_grid, ~] = ndgrid(points.input_z(:), z_km(:));
    n_he_profile = table(idx_grid(:), lat_grid(:), lon_grid(:), x_grid(:), y_grid(:), zxyz_grid(:), ...
        z_grid(:), n_he_matrix(:), ...
        'VariableNames', {'idx', 'input_lat', 'input_lon', 'input_x', 'input_y', 'input_z', 'z_km', 'n_He_cm3'});

    param_mat = fullfile(out_dir, sprintf('mshe_term_params_full_%s_%s.mat', date_tag, file_tag));
    save(param_mat, 'points', 'results', 'z_km', 'max_n_he', 'max_alt_km', ...
        'sza_selected_rad', 'sza_selected_deg', 'ut_selected_s', 'ut_selected_h', ...
        'selected_vs_idx', 'selected_file', 'n_he_matrix', 'qe_matrix', 'qe_per_he_matrix', ...
        'qr_matrix', 'q_matrix', 'loss_matrix', 'dif_matrix', 'hn_matrix', 'he_matrix', ...
        'n_he_table', 'n_he_profile', 'file_tag');
    fprintf('Saved data: %s\n', param_mat);

    results_csv = fullfile(out_dir, sprintf('mshe_term_results_%s_%s.csv', date_tag, file_tag));
    writetable(results, results_csv);
    fprintf('Saved summary: %s\n', results_csv);
end

switch lower(plot_mode)
    case 'sphere'
        plot_nhe_sphere(points, valid_idx, z_km, n_he_valid, earth_radius_km, title_str, date_tag, file_tag, out_dir);
    case 'lonlatalt'
        plot_nhe_lonlatalt(lon_valid, lat_valid, z_km, n_he_valid, title_str, date_tag, file_tag, out_dir);
        plot_nhe_maglonlatalt(mlon_valid, mlat_valid, z_km, n_he_valid, title_str, date_tag, file_tag, out_dir);
        if save_outputs
            plot_param_lonlatalt(lon_valid, lat_valid, z_km, qe_matrix(valid_idx, :), ...
                'qe', 'cm^{-3} s^{-1}', date_tag, file_tag, title_str, out_dir, false, 300);
            plot_param_lonlatalt(lon_valid, lat_valid, z_km, qe_per_he_matrix(valid_idx, :), ...
                'qe/He', 's^{-1}', date_tag, file_tag, title_str, out_dir, false, 300);
            plot_param_lonlatalt(lon_valid, lat_valid, z_km, qr_matrix(valid_idx, :), ...
                'qr', 'cm^{-3} s^{-1}', date_tag, file_tag, title_str, out_dir, false, 300);
            plot_param_lonlatalt(lon_valid, lat_valid, z_km, q_matrix(valid_idx, :), ...
                'q=qe+qr', 'cm^{-3} s^{-1}', date_tag, file_tag, title_str, out_dir, false, 300);
            plot_param_lonlatalt(lon_valid, lat_valid, z_km, he_matrix(valid_idx, :), ...
                'He', 'cm^{-3}', date_tag, file_tag, title_str, out_dir, true, 200);
            plot_param_lonlatalt(lon_valid, lat_valid, z_km, n_he_matrix(valid_idx, :), ...
                'n_he', 'cm^{-3}', date_tag, file_tag, title_str, out_dir, false, 300);
        end
    otherwise
        error('Unknown plot_mode: %s', plot_mode);
end

function plot_nhe_sphere(points, valid_idx, z_km, n_he_valid, earth_radius_km, title_str, date_tag, file_tag, out_dir)
n_he_plot = n_he_valid;
n_he_plot(~isfinite(n_he_plot)) = nan;
x_valid = points.input_x(valid_idx);
y_valid = points.input_y(valid_idx);
z_valid = points.input_z(valid_idx);
[x_valid, y_valid, z_valid, n_he_plot] = densify_spherical_track(x_valid, y_valid, z_valid, n_he_plot, 3);

[right_vec, up_vec, view_vec] = fit_track_plane_basis(x_valid, y_valid, z_valid);
track_u = x_valid * right_vec(1) + y_valid * right_vec(2) + z_valid * right_vec(3);
track_v = x_valid * up_vec(1) + y_valid * up_vec(2) + z_valid * up_vec(3);

[eq_u, eq_v] = projected_equator_uv(earth_radius_km, right_vec, up_vec, view_vec, false);
rot_angle = equator_alignment_angle(eq_u, eq_v);
[track_u, track_v] = rotate_projected_uv(track_u, track_v, rot_angle);

[north_u, north_v] = projected_north_pole_uv(earth_radius_km, right_vec, up_vec, view_vec);
[north_u, north_v] = rotate_projected_uv(north_u, north_v, rot_angle);
north_rot_angle = north_up_alignment_angle(north_u, north_v);
[track_u, track_v] = rotate_projected_uv(track_u, track_v, north_rot_angle);

theta_track = unwrap(atan2(track_v, track_u));
radial_km = earth_radius_km + z_km(:).';
[THETA, R3] = ndgrid(theta_track, radial_km);
u_plot = cos(THETA(:)) .* R3(:);
v_plot = sin(THETA(:)) .* R3(:);
c_nhe = n_he_plot(:);
keep = isfinite(u_plot) & isfinite(v_plot) & isfinite(c_nhe);

fig = figure('Color', 'w', 'Position', [120, 120, 980, 760], 'Renderer', 'painters');
hold on;

u_front = u_plot(keep);
v_front = v_plot(keep);
cv_front = c_nhe(keep);

theta = linspace(0, 2 * pi, 721);
disk_x = earth_radius_km * cos(theta);
disk_y = earth_radius_km * sin(theta);
patch(disk_x, disk_y, [0.85, 0.88, 0.92], 'EdgeColor', 'none', 'FaceAlpha', 1.0);
plot(disk_x, disk_y, 'Color', [0.18, 0.18, 0.18], 'LineWidth', 1.0);

[eq_u_full, eq_v_full] = projected_equator_uv(earth_radius_km, right_vec, up_vec, view_vec, false);
[eq_u_full, eq_v_full] = rotate_projected_uv(eq_u_full, eq_v_full, rot_angle);
[eq_u_full, eq_v_full] = rotate_projected_uv(eq_u_full, eq_v_full, north_rot_angle);
plot(eq_u_full, eq_v_full, '--', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.0);

plot_coastlines_on_projected_sphere(earth_radius_km, view_vec, right_vec, up_vec, rot_angle, north_rot_angle);
scatter(u_front, v_front, 12, cv_front, 'filled');

hainan_xyz = earth_radius_km * [-0.3068940506, 0.8912850400, 0.3338068592];
[uh, vh] = project_to_view_plane(hainan_xyz(1), hainan_xyz(2), hainan_xyz(3), right_vec, up_vec, view_vec);
[uh, vh] = rotate_projected_uv(uh, vh, rot_angle);
[uh, vh] = rotate_projected_uv(uh, vh, north_rot_angle);
plot(uh, vh, 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'w', 'LineWidth', 1.6);
text(uh + 40, vh + 20, 'hainan', ...
    'FontName', 'Times New Roman', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'r');

plot_extent = max(abs([u_front(:); v_front(:); earth_radius_km]));
plot_padding = max(80, 0.08 * plot_extent);

set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
xlabel('Projected horizontal distance [km]');
ylabel('Projected vertical distance [km]');
title(title_str);
cb = colorbar;
cb.Label.String = 'n_{he} [cm^{-3}]';
colormap(turbo);
auto_clim = auto_color_limits(cv_front);
if ~isempty(auto_clim)
    caxis(auto_clim);
end
axis equal;
xlim([-plot_extent - plot_padding, plot_extent + plot_padding]);
ylim([-plot_extent - plot_padding, plot_extent + plot_padding]);
box on;
grid off;
view(0, 90);

if ~isempty(out_dir)
    fig_png = fullfile(out_dir, sprintf('term_sphere_nhe_%s_%s.png', date_tag, file_tag));
    fig_fig = fullfile(out_dir, sprintf('term_sphere_nhe_%s_%s.fig', date_tag, file_tag));
    saveas(fig, fig_png);
    savefig(fig, fig_fig);
    fprintf('Saved figure: %s\n', fig_png);
    fprintf('Saved figure: %s\n', fig_fig);
end
hold off;
end

function clim = auto_color_limits(cvals)
cvals = cvals(isfinite(cvals));
if isempty(cvals)
    clim = [];
    return;
end
lo = prctile(cvals, 1);
hi = prctile(cvals, 99);
if ~isfinite(lo) || ~isfinite(hi) || hi <= lo
    lo = min(cvals);
    hi = max(cvals);
end
if hi <= lo
    pad = max(1e-6, 0.05 * max(abs(hi), 1));
    clim = [lo - pad, hi + pad];
else
    clim = [lo, hi];
end
end

function [ue, ve] = projected_equator_uv(radius_km, right_vec, up_vec, view_vec, front_only)
lon_eq = linspace(-180, 180, 721);
lat_eq = zeros(size(lon_eq));
[xe, ye, ze] = sph2cart(deg2rad(lon_eq), deg2rad(lat_eq), radius_km * ones(size(lat_eq)));
if front_only
    front_mask = (xe * view_vec(1) + ye * view_vec(2) + ze * view_vec(3)) >= 0;
    xe(~front_mask) = nan;
    ye(~front_mask) = nan;
    ze(~front_mask) = nan;
end
[ue, ve] = project_to_view_plane(xe, ye, ze, right_vec, up_vec, view_vec);
end

function [un, vn] = projected_north_pole_uv(radius_km, right_vec, up_vec, view_vec)
[un, vn] = project_to_view_plane(0, 0, radius_km, right_vec, up_vec, view_vec);
end

function ang = equator_alignment_angle(u, v)
mask = isfinite(u) & isfinite(v);
uv = [u(mask), v(mask)];
uv = uv - mean(uv, 1);
[~, ~, V] = svd(uv, 'econ');
major = V(:, 1);
ang = atan2(major(2), major(1));
end

function ang = north_up_alignment_angle(u, v)
ang = atan2(-u, v);
end

function [u_rot, v_rot] = rotate_projected_uv(u, v, ang)
ca = cos(ang);
sa = sin(ang);
u_rot = ca * u + sa * v;
v_rot = -sa * u + ca * v;
end

function [right_vec, up_vec, view_vec] = fit_track_plane_basis(x, y, z)
P = [x(:), y(:), z(:)];
P = P - mean(P, 1);
[~, ~, V] = svd(P, 'econ');
right_vec = V(:, 1).';
up_vec = V(:, 2).';
view_vec = V(:, 3).';

% Resolve the front/back ambiguity of the fitted plane by requiring the UK
% side of the Earth to be on the visible hemisphere.
uk_lat_deg = 54.0;
uk_lon_deg = -2.5;
uk_xyz = [cosd(uk_lat_deg) * cosd(uk_lon_deg), ...
          cosd(uk_lat_deg) * sind(uk_lon_deg), ...
          sind(uk_lat_deg)];
if dot(uk_xyz, view_vec) < 0
    view_vec = -view_vec;
    right_vec = -right_vec;
end

right_vec = right_vec / norm(right_vec);
up_vec = cross(view_vec, right_vec);
up_vec = up_vec / norm(up_vec);
end

function [x_dense, y_dense, z_dense, c_dense] = densify_spherical_track(x, y, z, cmat, n_insert)
if n_insert <= 0 || numel(x) < 2
    x_dense = x;
    y_dense = y;
    z_dense = z;
    c_dense = cmat;
    return;
end

n_base = numel(x);
nz = size(cmat, 2);
n_total = n_base * (n_insert + 1);
x_dense = nan(n_total, 1);
y_dense = nan(n_total, 1);
z_dense = nan(n_total, 1);
c_dense = nan(n_total, nz);
out_idx = 1;

for i = 1:n_base
    j = mod(i, n_base) + 1;
    p0 = [x(i), y(i), z(i)];
    p1 = [x(j), y(j), z(j)];
    c0 = cmat(i, :);
    c1 = cmat(j, :);
    for k = 0:n_insert
        t = k / (n_insert + 1);
        p = (1 - t) * p0 + t * p1;
        p = p / norm(p);
        x_dense(out_idx) = p(1);
        y_dense(out_idx) = p(2);
        z_dense(out_idx) = p(3);
        c_dense(out_idx, :) = (1 - t) * c0 + t * c1;
        out_idx = out_idx + 1;
    end
end
end

function [u, v, depth] = project_to_view_plane(x, y, z, right_vec, up_vec, view_vec)
u = x * right_vec(1) + y * right_vec(2) + z * right_vec(3);
v = x * up_vec(1) + y * up_vec(2) + z * up_vec(3);
depth = x * view_vec(1) + y * view_vec(2) + z * view_vec(3);
end

function plot_coastlines_on_projected_sphere(radius_km, view_vec, right_vec, up_vec, rot_angle, north_rot_angle)
coast = load_builtin_coastlines();
if isempty(coast)
    return;
end

lat = coast.lat(:);
lon = coast.lon(:);
is_break = isnan(lat) | isnan(lon);
start_idx = 1;

for i = 1:(numel(lat) + 1)
    at_break = (i > numel(lat)) || is_break(i);
    if at_break
        if i - start_idx >= 2
            seg_lat = lat(start_idx:i-1);
            seg_lon = lon(start_idx:i-1);
            [xc, yc, zc] = sph2cart(deg2rad(seg_lon), deg2rad(seg_lat), radius_km * ones(size(seg_lat)));
            front_mask = (xc * view_vec(1) + yc * view_vec(2) + zc * view_vec(3)) >= 0;
            xc(~front_mask) = nan;
            yc(~front_mask) = nan;
            zc(~front_mask) = nan;
            [uc, vc] = project_to_view_plane(xc, yc, zc, right_vec, up_vec, view_vec);
            [uc, vc] = rotate_projected_uv(uc, vc, rot_angle);
            [uc, vc] = rotate_projected_uv(uc, vc, north_rot_angle);
            plot(uc, vc, 'Color', [0.18, 0.18, 0.18], 'LineWidth', 0.6);
        end
        start_idx = i + 1;
    end
end
end

function coast = load_builtin_coastlines()
coast = [];

try
    script_dir = fileparts(mfilename('fullpath'));
    ne_path = fullfile(script_dir, 'data', 'coastlines_ne110m.mat');
    if isfile(ne_path)
        s = load(ne_path);
        if isfield(s, 'coastlat') && isfield(s, 'coastlon')
            coast.lat = s.coastlat;
            coast.lon = s.coastlon;
            return;
        end
    end
catch
end

try
    s = load('coastlines');
    if isfield(s, 'coastlat') && isfield(s, 'coastlon')
        coast.lat = s.coastlat;
        coast.lon = s.coastlon;
        return;
    end
catch
end

try
    s = load('coast');
    if isfield(s, 'lat') && isfield(s, 'long')
        coast.lat = s.lat;
        coast.lon = s.long;
        return;
    end
catch
end

coast = [];
end

function plot_nhe_lonlatalt(lon_valid, lat_valid, z_km, n_he_valid, title_str, date_tag, file_tag, out_dir)
n_he_plot = n_he_valid;
n_he_plot(~isfinite(n_he_plot)) = nan;
[LON, ALT] = ndgrid(lon_valid, z_km(:));
[LAT, ~] = ndgrid(lat_valid, z_km(:));
lon_plot = LON(:);
lat_plot = LAT(:);
alt_plot = ALT(:);
c_nhe = n_he_plot(:);
keep = isfinite(lon_plot) & isfinite(lat_plot) & isfinite(alt_plot) & isfinite(c_nhe);

figure('Color', 'w', 'Position', [120, 120, 1100, 560]);
scatter3(lon_plot(keep), lat_plot(keep), alt_plot(keep), 18, c_nhe(keep), 'filled');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
zlabel('Altitude [km]');
title(sprintf('%s (lon-lat-alt)', title_str));
cb = colorbar;
cb.Label.String = 'n_{he} [cm^{-3}]';
colormap(turbo);
caxis([0, 1]);
grid on;
box on;
view(45, 28);

if ~isempty(out_dir)
    fig3d_png = fullfile(out_dir, sprintf('term_lonlatalt_nhe_%s_%s.png', date_tag, file_tag));
    fig3d_fig = fullfile(out_dir, sprintf('term_lonlatalt_nhe_%s_%s.fig', date_tag, file_tag));
    saveas(gcf, fig3d_png);
    savefig(gcf, fig3d_fig);
    fprintf('Saved figure: %s\n', fig3d_png);
    fprintf('Saved figure: %s\n', fig3d_fig);
end
end

function plot_nhe_maglonlatalt(mlon_valid, mlat_valid, z_km, n_he_valid, title_str, date_tag, file_tag, out_dir)
n_he_plot = n_he_valid;
n_he_plot(~isfinite(n_he_plot)) = nan;
[MLON, ALT] = ndgrid(mlon_valid, z_km(:));
[MLAT, ~] = ndgrid(mlat_valid, z_km(:));
mlon_plot = MLON(:);
mlat_plot = MLAT(:);
alt_plot = ALT(:);
c_nhe = n_he_plot(:);
keep = isfinite(mlon_plot) & isfinite(mlat_plot) & isfinite(alt_plot) & isfinite(c_nhe);

figure('Color', 'w', 'Position', [120, 120, 1100, 560]);
scatter3(mlon_plot(keep), mlat_plot(keep), alt_plot(keep), 18, c_nhe(keep), 'filled');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
xlabel('Magnetic Longitude [deg]');
ylabel('Magnetic Latitude [deg]');
zlabel('Altitude [km]');
title(sprintf('%s (mag-lon-mag-lat-alt)', title_str));
cb = colorbar;
cb.Label.String = 'n_{he} [cm^{-3}]';
colormap(turbo);
caxis([0, 1]);
grid on;
box on;
view(45, 28);

if ~isempty(out_dir)
    fig3d_png = fullfile(out_dir, sprintf('term_maglonlatalt_nhe_%s_%s.png', date_tag, file_tag));
    fig3d_fig = fullfile(out_dir, sprintf('term_maglonlatalt_nhe_%s_%s.fig', date_tag, file_tag));
    saveas(gcf, fig3d_png);
    savefig(gcf, fig3d_fig);
    fprintf('Saved figure: %s\n', fig3d_png);
    fprintf('Saved figure: %s\n', fig3d_fig);
end
end

function points = load_conjugate_points(file_path)
fid = fopen(file_path, 'r');
if fid < 0
    error('Cannot open conjugate file: %s', file_path);
end
cleaner = onCleanup(@() fclose(fid));

idx = [];
input_lat = [];
input_lon = [];
input_x = [];
input_y = [];
input_z = [];

while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end

    line = strtrim(line);
    if isempty(line) || startsWith(line, '#') || startsWith(line, 'section=') || startsWith(line, 'columns:')
        continue;
    end

    vals = sscanf(line, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f');
    if numel(vals) < 7
        continue;
    end

    if numel(vals) >= 10
        input_lat_i = vals(2);
        input_lon_i = vals(3);
        input_x_i = vals(4);
        input_y_i = vals(5);
        input_z_i = vals(6);
        conj_lat_i = vals(7);
        conj_lon_i = vals(8);
    else
        input_lat_i = vals(2);
        input_lon_i = vals(3);
        input_x_i = nan;
        input_y_i = nan;
        input_z_i = nan;
        conj_lat_i = vals(4);
        conj_lon_i = vals(5);
    end

    if any(isnan([input_lat_i, input_lon_i, conj_lat_i, conj_lon_i]))
        continue;
    end

    idx(end+1, 1) = vals(1); %#ok<AGROW>
    input_lat(end+1, 1) = input_lat_i; %#ok<AGROW>
    input_lon(end+1, 1) = input_lon_i; %#ok<AGROW>
    input_x(end+1, 1) = input_x_i; %#ok<AGROW>
    input_y(end+1, 1) = input_y_i; %#ok<AGROW>
    input_z(end+1, 1) = input_z_i; %#ok<AGROW>
end

points = struct('idx', idx, 'input_lat', input_lat, 'input_lon', input_lon, ...
    'input_x', input_x, 'input_y', input_y, 'input_z', input_z);
end

function title_str = make_ut_title_from_first_file(file_path, fallback_date_tag)
title_str = sprintf('He(2^3S) UT,%s --:--:--', fallback_date_tag);
if isempty(file_path) || ~isfile(file_path)
    return;
end

fid = fopen(file_path, 'r');
if fid < 0
    return;
end
cleaner = onCleanup(@() fclose(fid));
line1 = fgetl(fid);
if ~ischar(line1)
    return;
end

nums = sscanf(strtrim(line1), '%f');
if numel(nums) < 5
    return;
end

idate_doy = round(nums(4));
ut_s = nums(5);
year = floor(idate_doy / 1000);
doy = mod(idate_doy, 1000);
if year < 1000 || doy < 1 || doy > 366
    date_tag = fallback_date_tag;
else
    date_dt = datetime(year, 1, 1) + days(doy - 1);
    date_tag = datestr(date_dt, 'yyyymmdd');
end

ut_s = mod(ut_s, 86400);
hh = floor(ut_s / 3600);
mm = floor(mod(ut_s, 3600) / 60);
ss = floor(mod(ut_s, 60));
title_str = sprintf('He(2^3S) UT,%s %02d:%02d:%02d', date_tag, hh, mm, ss);
end

function plot_param_lonlatalt(lon_values, lat_values, z_km, value_matrix, tag, unit_txt, date_tag, file_tag, title_str, out_dir, use_log_scale, clim_z_min_km)
fig = figure('Color', 'w', 'Position', [120, 120, 1100, 480]);
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

[LON, ALT] = ndgrid(lon_values(:), z_km(:));
[LAT, ~] = ndgrid(lat_values(:), z_km(:));
x_plot = LON(:);
y_plot = LAT(:);
z_plot = ALT(:);
c_plot = plot_matrix(:);
keep = isfinite(x_plot) & isfinite(y_plot) & isfinite(z_plot) & isfinite(c_plot);

scatter3(x_plot(keep), y_plot(keep), z_plot(keep), 18, c_plot(keep), 'filled');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1);
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
zlabel('Altitude [km]');
title(sprintf('%s | %s (lon-lat-alt)', title_str, title_suffix));
cb = colorbar;
cb.Label.String = color_label;
colormap(turbo);

if strcmp(tag, 'He')
    zlim([200, max(z_km)]);
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
box on;
view(45, 28);

safe_tag = regexprep(tag, '[^a-zA-Z0-9_]', '_');
png_path = fullfile(out_dir, sprintf('lonlatalt_%s_%s_%s.png', safe_tag, date_tag, file_tag));
fig_path = fullfile(out_dir, sprintf('lonlatalt_%s_%s_%s.fig', safe_tag, date_tag, file_tag));
saveas(fig, png_path);
savefig(fig, fig_path);
fprintf('Saved figure: %s\n', png_path);
fprintf('Saved figure: %s\n', fig_path);
end

function [sza, ut_s] = read_sza_ut_from_glow_file(file_path)
sza = nan;
ut_s = nan;
fid = fopen(file_path, 'r');
if fid < 0
    return;
end
cleaner = onCleanup(@() fclose(fid));
line1 = fgetl(fid);
if ~ischar(line1)
    return;
end

nums = sscanf(strtrim(line1), '%f');
if numel(nums) < 5
    return;
end

sza = nums(3);
ut_s = nums(5);
end

function name_only = get_file_name_only(file_path)
[~, name, ext] = fileparts(file_path);
name_only = [name, ext];
end

function [mlat_deg, mlon_deg] = geo_to_centered_dipole(glat_deg, glon_deg)
pole_lat_deg = 80.0;
pole_lon_deg = -72.2;

glat = deg2rad(glat_deg(:));
glon = deg2rad(glon_deg(:));
plat = deg2rad(pole_lat_deg);
plon = deg2rad(pole_lon_deg);
dlam = glon - plon;

sin_mlat = sin(glat) .* sin(plat) + cos(glat) .* cos(plat) .* cos(dlam);
sin_mlat = max(-1, min(1, sin_mlat));
mlat = asin(sin_mlat);

y = cos(glat) .* sin(dlam);
x = sin(glat) .* cos(plat) - cos(glat) .* sin(plat) .* cos(dlam);
mlon = atan2(y, x);

mlat_deg = rad2deg(mlat);
mlon_deg = rad2deg(mlon);
end

function doy_folder = date_to_doy_folder(date_tag)
date_dt = datetime(date_tag, 'InputFormat', 'yyyyMMdd');
doy_folder = sprintf('%03d', day(date_dt, 'dayofyear'));
end

