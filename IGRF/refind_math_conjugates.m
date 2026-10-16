function refind_math_conjugates()
% Re-find NaN endpoints as IGRF mathematical conjugates.
%
% This script is intentionally less conservative than refind_far_conjugates:
% it does not reject high B_ratio values, because high-latitude IGRF field
% lines may need a long mathematical trace even when the physical
% magnetospheric interpretation is questionable.
%
% Output files are written as *_math.txt and do not overwrite the physical
% conservative *_refind.txt files.

script_dir = fileparts(mfilename('fullpath'));
output_dir = fullfile(script_dir, 'output');

cfg.altitude_km = 0.1;
cfg.conj_alt_km = 0;
cfg.min_surface_distance_km = 500;
cfg.min_path_length_km = 1000;
cfg.max_path_length_km = 500000;
cfg.max_iter = 100000;
cfg.step.base = 100;
cfg.step.segment_remain_km = [5000, 2000, 1000, 500, 200, 50];
cfg.step.segment_multipliers = [4.0, 3.0, 2.0, 1.25, 1.0, 0.5, 0.25];
cfg.step.min_abs = 1;
cfg.step.max_abs = 500;
cfg.step_candidates = [1, -1, 2, -2, 0.5, -0.5, 4, -4, 0.25, -0.25];

repair_one(resolve_input(output_dir, 'conjugate_morning'), cfg);
repair_one(resolve_input(output_dir, 'conjugate_evening'), cfg);
end


function path_in = resolve_input(output_dir, base_name)
refind_path = fullfile(output_dir, [base_name, '_refind.txt']);
raw_path = fullfile(output_dir, [base_name, '.txt']);

if isfile(refind_path)
    path_in = refind_path;
elseif isfile(raw_path)
    path_in = raw_path;
else
    error('Input file not found for %s', base_name)
end
end


function repair_one(path_in, cfg)
[section_name, rows] = read_conjugate_txt(path_in);
repaired = 0;
failed = 0;
skipped = 0;

for i = 1:size(rows, 1)
    src_lat = rows(i, 2);
    src_lon = rows(i, 3);
    conj_lat = rows(i, 7);
    conj_lon = rows(i, 8);

    if any(~isfinite([src_lat, src_lon]))
        skipped = skipped + 1;
        continue
    end

    if all(isfinite([conj_lat, conj_lon]))
        skipped = skipped + 1;
        continue
    end

    [ok, best_lat, best_lon, best_ratio] = refind_math_endpoint(src_lat, src_lon, cfg);
    if ok
        rows(i, 7) = best_lat;
        rows(i, 8) = best_lon;
        rows(i, 9) = cfg.conj_alt_km;
        rows(i, 10) = best_ratio;
        repaired = repaired + 1;
    else
        rows(i, 7) = NaN;
        rows(i, 8) = NaN;
        rows(i, 9) = NaN;
        rows(i, 10) = NaN;
        failed = failed + 1;
    end
end

[folder, name, ext] = fileparts(path_in);
if endsWith(name, '_refind')
    name = extractBefore(name, strlength(name) - strlength('_refind') + 1);
end
path_out = fullfile(folder, [name, '_math', ext]);
write_conjugate_txt(path_out, section_name, rows);

fprintf('Math refind repaired %d rows, unresolved %d rows, skipped %d rows: %s -> %s\n', ...
    repaired, failed, skipped, path_in, path_out);
end


function [ok, best_lat, best_lon, best_ratio] = refind_math_endpoint(src_lat, src_lon, cfg)
ok = false;
best_lat = NaN;
best_lon = NaN;
best_ratio = NaN;
best_dist = -inf;

for k = 1:numel(cfg.step_candidates)
    trial = cfg.step;
    trial.base = cfg.step.base * cfg.step_candidates(k);
    try
        [cand_lat, cand_lon, cand_ratio, cand_dist] = trace_math_fieldline(src_lat, src_lon, cfg, trial);
        valid_candidate = isfinite(cand_dist) ...
            && cand_dist >= cfg.min_surface_distance_km ...
            && isfinite(cand_ratio);
        if valid_candidate && cand_dist > best_dist
            best_dist = cand_dist;
            best_lat = cand_lat;
            best_lon = cand_lon;
            best_ratio = cand_ratio;
            ok = true;
        end
    catch
        % Try next step candidate.
    end
end
end


function [conj_lat, conj_lon, B_ratio, dist_km] = trace_math_fieldline(src_lat, src_lon, cfg, step_cfg)
Re = 6371.2;
target_r = Re + cfg.conj_alt_km;
stop_tol_km = get_stop_tolerance(step_cfg);

theta1 = (90 - src_lat) * pi / 180;
lambda1 = src_lon * pi / 180;
r1 = Re + cfg.altitude_km;

[B_r, B_theta, B_lambda] = IGRF_vs3(r1, theta1, lambda1);
[Bx, By, Bz] = sph_to_cart(B_r, B_theta, B_lambda, theta1, lambda1);
Blocal = norm([Bx By Bz]);

Bm = [];
path_length_km = 0;
iter = 0;
prev_remain = r1 - target_r;

while true
    iter = iter + 1;
    if iter > cfg.max_iter
        error('math trace exceeded max_iter')
    end

    [B_r, B_theta, B_lambda] = IGRF_vs3(r1, theta1, lambda1);
    [B_x, B_y, B_z] = sph_to_cart(B_r, B_theta, B_lambda, theta1, lambda1);
    B = sqrt(B_x^2 + B_y^2 + B_z^2);
    if ~isfinite(B) || B <= 0
        error('invalid field magnitude')
    end
    Bm(end + 1) = B; %#ok<AGROW>

    [x1, y1, z1] = rtp_xyz(r1, theta1, lambda1);
    step_now = choose_step(step_cfg, r1 - target_r);

    x2 = x1 + step_now * B_x / B;
    y2 = y1 + step_now * B_y / B;
    z2 = z1 + step_now * B_z / B;

    seg_len = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);
    path_length_km = path_length_km + seg_len;
    if path_length_km > cfg.max_path_length_km
        error('math trace exceeded path length')
    end

    [r2, theta2, lambda2] = xyz_rtp(x2, y2, z2);
    remain = r2 - target_r;
    lat2 = 90 - theta2 * 180 / pi;
    lon2 = wrap_lon(lambda2 * 180 / pi);
    surf_dist = spherical_distance(src_lat, src_lon, lat2, lon2, Re);

    if path_length_km >= cfg.min_path_length_km && surf_dist >= cfg.min_surface_distance_km && remain <= stop_tol_km
        if prev_remain > 0 && remain <= 0
            alpha = prev_remain / (prev_remain - remain);
            alpha = min(1, max(0, alpha));
            xh = x1 + alpha * (x2 - x1);
            yh = y1 + alpha * (y2 - y1);
            zh = z1 + alpha * (z2 - z1);
            [~, theta2, lambda2] = xyz_rtp(xh, yh, zh);
            lat2 = 90 - theta2 * 180 / pi;
            lon2 = wrap_lon(lambda2 * 180 / pi);
            surf_dist = spherical_distance(src_lat, src_lon, lat2, lon2, Re);
        end
        conj_lat = lat2;
        conj_lon = lon2;
        dist_km = surf_dist;
        B_ratio = max(Bm) / Blocal;
        return
    end

    r1 = r2;
    theta1 = theta2;
    lambda1 = lambda2;
    prev_remain = remain;
end
end


function [section_name, rows] = read_conjugate_txt(path_in)
fid = fopen(path_in, 'r');
if fid < 0
    error('Cannot open file: %s', path_in)
end
cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

section_name = '';
rows = [];

while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end

    line = strtrim(line);
    if isempty(line)
        continue
    end
    if startsWith(line, 'section=')
        section_name = extractAfter(line, 'section=');
        continue
    end
    if startsWith(line, 'columns:')
        continue
    end

    vals = sscanf(line, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f');
    if numel(vals) == 10
        rows(end + 1, :) = vals.'; %#ok<AGROW>
    end
end
end


function write_conjugate_txt(path_out, section_name, rows)
fid = fopen(path_out, 'w');
if fid < 0
    error('Cannot write file: %s', path_out)
end
cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'section=%s\n', section_name);
fprintf(fid, 'columns: idx,input_lat,input_lon,input_x,input_y,input_z,conj_lat,conj_lon,conj_alt_km,B_ratio\n');
for i = 1:size(rows, 1)
    fprintf(fid, '%d,%.6f,%.6f,%.10f,%.10f,%.10f,%.6f,%.6f,%.6f,%.8f\n', ...
        round(rows(i, 1)), rows(i, 2), rows(i, 3), rows(i, 4), rows(i, 5), rows(i, 6), ...
        rows(i, 7), rows(i, 8), rows(i, 9), rows(i, 10));
end
end


function lon = wrap_lon(lon)
if lon > 180
    lon = lon - 360;
elseif lon < -180
    lon = lon + 360;
end
end


function dist = spherical_distance(lat1, lon1, lat2, lon2, Re)
lat1 = lat1 * pi / 180;
lon1 = lon1 * pi / 180;
lat2 = lat2 * pi / 180;
lon2 = lon2 * pi / 180;

dlat = lat2 - lat1;
dlon = lon2 - lon1;
a = sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));
dist = Re * c;
end


function step_now = choose_step(step_cfg, remain_km)
[thresholds, multipliers] = normalize_segment_arrays(step_cfg);
scale = multipliers(end);
for i = 1:numel(thresholds)
    if remain_km > thresholds(i)
        scale = multipliers(i);
        break
    end
end

step_abs = abs(step_cfg.base) * scale;
step_abs = min(step_cfg.max_abs, max(step_cfg.min_abs, step_abs));
step_now = sign(step_cfg.base) * step_abs;
end


function [thresholds, multipliers] = normalize_segment_arrays(step_cfg)
thresholds = step_cfg.segment_remain_km(:).';
multipliers = step_cfg.segment_multipliers(:).';
if numel(multipliers) ~= numel(thresholds) + 1
    error('segment_multipliers must have one more element than segment_remain_km')
end
[thresholds, order] = sort(thresholds, 'descend');
multipliers = [multipliers(order), multipliers(end)];
end


function stop_tol_km = get_stop_tolerance(step_cfg)
min_multiplier = min(abs(step_cfg.segment_multipliers));
stop_tol_km = abs(step_cfg.base) * min_multiplier;
stop_tol_km = min(step_cfg.max_abs, max(step_cfg.min_abs, stop_tol_km));
end


function [Bx, By, Bz] = sph_to_cart(Br, Btheta, Blambda, theta, lambda)
Bx = sin(theta) * cos(lambda) * Br ...
   + cos(theta) * cos(lambda) * Btheta ...
   - sin(lambda) * Blambda;

By = sin(theta) * sin(lambda) * Br ...
   + cos(theta) * sin(lambda) * Btheta ...
   + cos(lambda) * Blambda;

Bz = cos(theta) * Br ...
   - sin(theta) * Btheta;
end
