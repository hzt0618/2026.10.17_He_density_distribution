function debug_step_candidates()
lat = -0.5;
lon = 106;
altitude_km = 0.1;
rcon = 0;

step = struct( ...
    'base', 200, ...
    'segment_remain_km', [2000, 1000, 500, 200, 50], ...
    'segment_multipliers', [4.0, 2.5, 1.5, 1.0, 0.5, 0.25], ...
    'min_abs', 1, ...
    'max_abs', 1000, ...
    'retry_multipliers', [1, -1, 2, -2], ...
    'small_retry_multipliers', [0.5, -0.5, 0.25, -0.25]);

candidates = build_step_candidates_local(step);

for k = 1:numel(candidates)
    trial_step = candidates(k);
    trial_cfg = step;
    trial_cfg.base = trial_step;

    try
        [f_lat, f_lon, f_alt, f_ratio, f_dist] = trace_fieldline_local(lat, lon, altitude_km, rcon, trial_cfg);
        fprintf('STEP=%g FORWARD LAT=%.10f LON=%.10f ALT=%.10f BRATIO=%.10f DIST=%.10f\n', ...
            trial_step, f_lat, f_lon, f_alt, f_ratio, f_dist);
    catch ME
        fprintf('STEP=%g FORWARD FAIL MSG=%s\n', trial_step, ME.message);
        f_dist = NaN;
    end

    try
        rev_cfg = reverse_step_local(trial_cfg);
        [r_lat, r_lon, r_alt, r_ratio, r_dist] = trace_fieldline_local(lat, lon, altitude_km, rcon, rev_cfg);
        fprintf('STEP=%g REVERSE LAT=%.10f LON=%.10f ALT=%.10f BRATIO=%.10f DIST=%.10f\n', ...
            trial_step, r_lat, r_lon, r_alt, r_ratio, r_dist);
    catch ME
        fprintf('STEP=%g REVERSE FAIL MSG=%s\n', trial_step, ME.message);
        r_dist = NaN;
    end

    if isfinite(f_dist) && isfinite(r_dist)
        if r_dist > f_dist
            chosen = 'REVERSE';
        else
            chosen = 'FORWARD';
        end
    elseif isfinite(f_dist)
        chosen = 'FORWARD_ONLY';
    elseif isfinite(r_dist)
        chosen = 'REVERSE_ONLY';
    else
        chosen = 'NONE';
    end

    fprintf('STEP=%g CHOSEN=%s\n', trial_step, chosen);
end
end

function candidates = build_step_candidates_local(step)
base = step.base;
candidates = base .* step.retry_multipliers;
if isfield(step, 'small_retry_multipliers') && ~isempty(step.small_retry_multipliers)
    candidates = [candidates, base .* step.small_retry_multipliers];
end
candidates = unique(candidates, 'stable');
end

function [conj_lat, conj_lon, conj_alt, B_ratio, dist_km] = trace_fieldline_local(rlat, rlong, altitude_km, rcon, step)
Re = 6371.2;
step_cfg = normalize_step_config_local(step);
stop_tol_km = get_stop_tolerance_local(step_cfg);
target_r = Re + rcon;

theta1 = (90 - rlat) * pi / 180;
lambda1 = rlong * pi / 180;
r1 = Re + altitude_km;

[B_r, B_theta, B_lambda] = IGRF_vs3(r1, theta1, lambda1);
[Bx, By, Bz] = sph_to_cart_local(B_r, B_theta, B_lambda, theta1, lambda1);
Blocal = norm([Bx By Bz]);

Bm = [];
max_iter = 20000;
iter = 0;
prev_remain = inf;
stall_count = 0;
stall_limit = 2000;
remain = r1 - target_r;
path_length_km = 0;
max_path_length_km = 200000;

while true
    iter = iter + 1;
    if iter > max_iter
        error('Tracing failed: too many iterations')
    end

    [B_r, B_theta, B_lambda] = IGRF_vs3(r1, theta1, lambda1);

    B_x = sin(theta1) * cos(lambda1) * B_r + cos(theta1) * cos(lambda1) * B_theta - sin(lambda1) * B_lambda;
    B_y = sin(theta1) * sin(lambda1) * B_r + cos(theta1) * sin(lambda1) * B_theta + cos(lambda1) * B_lambda;
    B_z = cos(theta1) * B_r - sin(theta1) * B_theta;

    B = sqrt(B_x^2 + B_y^2 + B_z^2);
    if ~isfinite(B) || B <= 0
        error('Tracing failed: invalid magnetic field magnitude')
    end
    Bm(end + 1) = B; %#ok<AGROW>

    [x1, y1, z1] = rtp_xyz(r1, theta1, lambda1);
    step_now = choose_step_local(step_cfg, remain);

    x2 = x1 + step_now * B_x / B;
    y2 = y1 + step_now * B_y / B;
    z2 = z1 + step_now * B_z / B;

    path_length_km = path_length_km + sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);
    if path_length_km > max_path_length_km
        error('Tracing failed: path length limit exceeded')
    end

    prev_r = r1;
    prev_theta = theta1;
    prev_lambda = lambda1;
    prev_remain = remain;
    [r1, theta1, lambda1] = xyz_rtp(x2, y2, z2);
    remain = r1 - target_r;

    if iter > 1 && remain <= stop_tol_km
        if prev_remain > 0 && remain <= 0
            alpha = prev_remain / (prev_remain - remain);
            alpha = min(1, max(0, alpha));
            x_hit = x1 + alpha * (x2 - x1);
            y_hit = y1 + alpha * (y2 - y1);
            z_hit = z1 + alpha * (z2 - z1);
            [r1, theta1, lambda1] = xyz_rtp(x_hit, y_hit, z_hit);
        elseif abs(remain) > stop_tol_km && abs(prev_remain) < abs(remain)
            r1 = prev_r;
            theta1 = prev_theta;
            lambda1 = prev_lambda;
            remain = prev_remain;
        end
        break
    end

    if remain >= prev_remain - 1e-3
        stall_count = stall_count + 1;
    else
        stall_count = 0;
    end

    if stall_count > stall_limit
        error('Tracing failed: radial progress stalled')
    end
end

conj_lat = 90 - theta1 * 180 / pi;
conj_lon = lambda1 * 180 / pi;
conj_alt = rcon;

if conj_lon > 180
    conj_lon = conj_lon - 360;
elseif conj_lon < -180
    conj_lon = conj_lon + 360;
end

dist_km = spherical_distance_local(rlat, rlong, conj_lat, conj_lon, Re);

if isempty(Bm) || Blocal <= 0 || ~isfinite(Blocal)
    B_ratio = NaN;
else
    B_ratio = max(Bm) / Blocal;
end
end

function step_out = reverse_step_local(step_in)
step_out = step_in;
step_out.base = -step_in.base;
end

function [Bx, By, Bz] = sph_to_cart_local(B_r, B_theta, B_lambda, theta, lambda)
Bx = sin(theta) * cos(lambda) * B_r + cos(theta) * cos(lambda) * B_theta - sin(lambda) * B_lambda;
By = sin(theta) * sin(lambda) * B_r + cos(theta) * sin(lambda) * B_theta + cos(lambda) * B_lambda;
Bz = cos(theta) * B_r - sin(theta) * B_theta;
end

function step_cfg = normalize_step_config_local(step)
step_cfg = step;
if ~isfield(step_cfg, 'segment_remain_km')
    step_cfg.segment_remain_km = [];
end
if ~isfield(step_cfg, 'segment_multipliers')
    step_cfg.segment_multipliers = [];
end
if ~isfield(step_cfg, 'far_multiplier')
    step_cfg.far_multiplier = 2.5;
end
if ~isfield(step_cfg, 'mid_multiplier')
    step_cfg.mid_multiplier = 1.0;
end
if ~isfield(step_cfg, 'near_multiplier')
    step_cfg.near_multiplier = 0.25;
end
if ~isfield(step_cfg, 'far_remain_km')
    step_cfg.far_remain_km = 1000;
end
if ~isfield(step_cfg, 'near_remain_km')
    step_cfg.near_remain_km = 200;
end
if ~isfield(step_cfg, 'min_abs')
    step_cfg.min_abs = 1;
end
if ~isfield(step_cfg, 'max_abs')
    step_cfg.max_abs = 50;
end
end

function step_now = choose_step_local(step_cfg, remain_km)
if ~isempty(step_cfg.segment_remain_km) || ~isempty(step_cfg.segment_multipliers)
    [thresholds, multipliers] = normalize_segment_arrays_local(step_cfg);
    scale = multipliers(end);
    for i = 1:numel(thresholds)
        if remain_km > thresholds(i)
            scale = multipliers(i);
            break
        end
    end
else
    if remain_km > step_cfg.far_remain_km
        scale = step_cfg.far_multiplier;
    elseif remain_km > step_cfg.near_remain_km
        scale = step_cfg.mid_multiplier;
    else
        scale = step_cfg.near_multiplier;
    end
end

step_abs = abs(step_cfg.base) * scale;
step_abs = min(step_cfg.max_abs, max(step_cfg.min_abs, step_abs));
step_now = sign(step_cfg.base) * step_abs;
if step_now == 0
    step_now = step_abs;
end
end

function [thresholds, multipliers] = normalize_segment_arrays_local(step_cfg)
thresholds = step_cfg.segment_remain_km(:).';
multipliers = step_cfg.segment_multipliers(:).';

if isempty(thresholds) && isempty(multipliers)
    thresholds = [step_cfg.far_remain_km, step_cfg.near_remain_km];
    multipliers = [step_cfg.far_multiplier, step_cfg.mid_multiplier, step_cfg.near_multiplier];
    return
end

if numel(multipliers) ~= numel(thresholds) + 1
    error('segment_multipliers must have exactly one more element than segment_remain_km')
end

[thresholds, order] = sort(thresholds, 'descend');
multipliers = [multipliers(order), multipliers(end)];
end

function tol = get_stop_tolerance_local(step_cfg)
if ~isempty(step_cfg.segment_multipliers)
    min_multiplier = min(abs(step_cfg.segment_multipliers));
else
    min_multiplier = abs(step_cfg.near_multiplier);
end

tol = abs(step_cfg.base) * min_multiplier;
tol = min(step_cfg.max_abs, max(step_cfg.min_abs, tol));
end

function dist = spherical_distance_local(lat1, lon1, lat2, lon2, Re)
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
