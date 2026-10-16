function [conj_lat, conj_lon, conj_alt, B_ratio] = con_mag(rlat, rlong, altitude_km, rcon, step)

% Trace both directions and keep the farther valid endpoint.
[cand1.lat, cand1.lon, cand1.alt, cand1.ratio, cand1.dist] = ...
    trace_fieldline_with_fallback(rlat, rlong, altitude_km, rcon, step);

step_rev = reverse_step(step);
[cand2.lat, cand2.lon, cand2.alt, cand2.ratio, cand2.dist] = ...
    trace_fieldline_with_fallback(rlat, rlong, altitude_km, rcon, step_rev);

[best, reason] = pick_best_candidate(cand1, cand2);

if ~isempty(reason)
    fprintf('Warning: %s (forward=%.1f km, reverse=%.1f km)\n', ...
        reason, cand1.dist, cand2.dist);
end

conj_lat = best.lat;
conj_lon = best.lon;
conj_alt = best.alt;
B_ratio = best.ratio;

end


function [conj_lat, conj_lon, conj_alt, B_ratio, dist_km] = ...
    trace_fieldline_with_fallback(rlat, rlong, altitude_km, rcon, step)
try
    [conj_lat, conj_lon, conj_alt, B_ratio, dist_km] = ...
        trace_fieldline(rlat, rlong, altitude_km, rcon, step);
catch ME
    if ~should_use_polar_fallback(rlat, ME)
        rethrow(ME)
    end

    fprintf('Info: polar fallback enabled at lat=%.2f (%s)\n', rlat, ME.message);
    try
        [conj_lat, conj_lon, conj_alt, B_ratio, dist_km] = ...
            trace_fieldline_polar(rlat, rlong, altitude_km, rcon, step);
    catch ME2
        fprintf('Info: polar direct fallback enabled at lat=%.2f (%s)\n', rlat, ME2.message);
        [conj_lat, conj_lon, conj_alt, B_ratio, dist_km] = ...
            trace_fieldline_polar_direct(rlat, rlong, altitude_km, rcon, step);
    end
end
end


function [conj_lat, conj_lon, conj_alt, B_ratio, dist_km] = trace_fieldline(rlat, rlong, altitude_km, rcon, step)

Re = 6371.2;
step_cfg = normalize_step_config(step);
stop_tol_km = get_stop_tolerance(step_cfg);
target_r = Re + rcon;
use_local_footpoint_guard = abs(rlat) >= 65;

theta1 = (90 - rlat) * pi / 180;
lambda1 = mod(rlong * pi / 180, 2 * pi);
r1 = Re + altitude_km;

[B_r, B_theta, B_lambda] = IGRF_vs3(r1, theta1, lambda1);
[Bx, By, Bz] = sph_to_cart(B_r, B_theta, B_lambda, theta1, lambda1);
Blocal = norm([Bx By Bz]);

Bm = [];
max_iter = 20000;
iter = 0;
prev_remain = inf;
stall_count = 0;
stall_limit = 2000;
stall_shrink_trigger = 25;
remain = r1 - target_r;
path_length_km = 0;
max_path_length_km = 200000;
departure_threshold_km = max(5, min(100, abs(step_cfg.base)));
max_r_seen = r1;
has_departed_surface = false;

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
    step_now = choose_step(step_cfg, remain);
    step_now = apply_polar_step_limit(step_now, theta1, step_cfg);

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
    lambda1 = mod(lambda1, 2 * pi);
    max_r_seen = max(max_r_seen, r1);
    if max_r_seen - target_r >= departure_threshold_km
        has_departed_surface = true;
    end
    remain = r1 - target_r;

    % Always integrate at least one step for very low start altitudes.
    if iter > 1 && remain <= stop_tol_km && (~use_local_footpoint_guard || has_departed_surface)
        if prev_remain > 0 && remain <= 0
            alpha = prev_remain / (prev_remain - remain);
            alpha = min(1, max(0, alpha));
            x_hit = x1 + alpha * (x2 - x1);
            y_hit = y1 + alpha * (y2 - y1);
            z_hit = z1 + alpha * (z2 - z1);
            [r1, theta1, lambda1] = xyz_rtp(x_hit, y_hit, z_hit);
            lambda1 = mod(lambda1, 2 * pi);
        elseif abs(remain) > stop_tol_km && abs(prev_remain) < abs(remain)
            r1 = prev_r;
            theta1 = prev_theta;
            lambda1 = prev_lambda;
            remain = prev_remain;
        end
        break
    end

    if use_local_footpoint_guard && iter > 1 && remain <= 0 && ~has_departed_surface
        error('Tracing failed: local footpoint encountered before reaching conjugate branch')
    end

    if remain >= prev_remain - 1e-3
        stall_count = stall_count + 1;
    else
        stall_count = 0;
    end

    if stall_count >= stall_shrink_trigger
        step_cfg.base = sign(step_cfg.base) * max(step_cfg.min_abs, abs(step_cfg.base) / 2);
        stall_count = 0;
    end

    if stall_count > stall_limit
        error('Tracing failed: radial progress stalled')
    end

    prev_remain = remain;
end

conj_lat = 90 - theta1 * 180 / pi;
conj_lon = lambda1 * 180 / pi;
conj_alt = rcon;

if conj_lon > 180
    conj_lon = conj_lon - 360;
elseif conj_lon < -180
    conj_lon = conj_lon + 360;
end

dist_km = spherical_distance(rlat, rlong, conj_lat, conj_lon, Re);

if isempty(Bm) || Blocal <= 0 || ~isfinite(Blocal)
    B_ratio = NaN;
else
    B_ratio = max(Bm) / Blocal;
end

end


function tf = should_use_polar_fallback(rlat, ME)
tf = false;
if abs(rlat) < 65
    return
end

msg = ME.message;
tf = contains(msg, 'path length limit exceeded') || ...
    contains(msg, 'too many iterations') || ...
    contains(msg, 'radial progress stalled') || ...
    contains(msg, 'local footpoint encountered');
end


function [conj_lat, conj_lon, conj_alt, B_ratio, dist_km] = ...
    trace_fieldline_polar(rlat, rlong, altitude_km, rcon, step)
Re = 6371.2;
target_r = Re + rcon;
theta1 = (90 - rlat) * pi / 180;
lambda1 = mod(rlong * pi / 180, 2 * pi);
r1 = Re + altitude_km;
start_t = tic;
time_limit_s = get_step_time_limit(step, 'polar_time_limit_s', 12);

step_sign = sign(step.base);
if step_sign == 0
    step_sign = 1;
end
stop_tol_km = 0.5;
max_iter = 80000;
path_length_km = 0;
max_path_length_km = 400000;
Bm = [];
departure_threshold_km = max(5, min(50, abs(step.base)));
max_r_seen = r1;
has_departed_surface = false;
fast_departure_threshold_km = max(150, 3 * departure_threshold_km);
fast_step_abs = max(20, min(100, abs(step.base)));
fine_step_abs = max(5, min(20, abs(step.base) / 10));
has_entered_remote_branch = false;
apex_drop_tol_km = 0.05;
prev_r_raw = r1;
prev_theta_raw = theta1;
prev_lambda_raw = lambda1;

[B_r, B_theta, B_lambda] = IGRF_vs3(r1, theta1, lambda1);
[Bx, By, Bz] = sph_to_cart(B_r, B_theta, B_lambda, theta1, lambda1);
Blocal = norm([Bx By Bz]);

for iter = 1:max_iter
    if toc(start_t) > time_limit_s
        error('Tracing failed: polar fallback time limit exceeded')
    end
    [B_r, B_theta, B_lambda] = IGRF_vs3(r1, theta1, lambda1);

    B_x = sin(theta1) * cos(lambda1) * B_r + cos(theta1) * cos(lambda1) * B_theta - sin(lambda1) * B_lambda;
    B_y = sin(theta1) * sin(lambda1) * B_r + cos(theta1) * sin(lambda1) * B_theta + cos(lambda1) * B_lambda;
    B_z = cos(theta1) * B_r - sin(theta1) * B_theta;

    B = sqrt(B_x^2 + B_y^2 + B_z^2);
    if ~isfinite(B) || B <= 0
        error('Tracing failed: invalid magnetic field magnitude in polar fallback')
    end
    Bm(end + 1) = B; %#ok<AGROW>

    [x1, y1, z1] = rtp_xyz(r1, theta1, lambda1);
    remain = r1 - target_r;

    local_scale = max(0.3, min(1.0, abs(sin(theta1)) / 0.5));
    if ~has_entered_remote_branch
        trial_step_abs = fast_step_abs;
    else
        trial_step_abs = max(fine_step_abs, min(fast_step_abs / 3, max(abs(remain), fine_step_abs)));
    end
    trial_step = step_sign * max(0.1, trial_step_abs * local_scale);

    x2 = x1 + trial_step * B_x / B;
    y2 = y1 + trial_step * B_y / B;
    z2 = z1 + trial_step * B_z / B;

    path_length_km = path_length_km + sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);
    if path_length_km > max_path_length_km
        error('Tracing failed: polar fallback path length limit exceeded')
    end

    prev_r = r1;
    prev_theta = theta1;
    prev_lambda = lambda1;
    prev_remain = remain;

    [r1, theta1, lambda1] = xyz_rtp(x2, y2, z2);
    lambda1 = mod(lambda1, 2 * pi);
    max_r_seen = max(max_r_seen, r1);
    if max_r_seen - target_r >= departure_threshold_km
        has_departed_surface = true;
    end
    if has_departed_surface && (max_r_seen - target_r >= fast_departure_threshold_km) ...
            && (r1 < prev_r - apex_drop_tol_km || r1 < max_r_seen - apex_drop_tol_km)
        has_entered_remote_branch = true;
    end
    remain = r1 - target_r;

    if abs(remain) <= stop_tol_km && has_departed_surface && has_entered_remote_branch
        break
    end

    if prev_remain > 0 && remain <= 0
        if ~has_departed_surface || ~has_entered_remote_branch
            error('Tracing failed: local footpoint encountered before polar conjugate branch')
        end
        [r1, theta1, lambda1] = refine_surface_hit( ...
            x1, y1, z1, x2, y2, z2, target_r, stop_tol_km);
        break
    end

    if abs(remain) > abs(prev_remain)
        if ~has_entered_remote_branch
            fast_step_abs = max(1, fast_step_abs / 1.5);
        else
            fine_step_abs = max(2, fine_step_abs / 1.5);
        end
        r1 = prev_r;
        theta1 = prev_theta;
        lambda1 = prev_lambda;
    elseif has_entered_remote_branch && abs(remain) < 20
        fine_step_abs = max(2, fine_step_abs / 1.5);
    end

    prev_r_raw = prev_r;
    prev_theta_raw = prev_theta;
    prev_lambda_raw = prev_lambda;
end

if iter == max_iter
    error('Tracing failed: polar fallback exceeded maximum iterations')
end

conj_lat = 90 - theta1 * 180 / pi;
conj_lon = lambda1 * 180 / pi;
conj_alt = rcon;

if conj_lon > 180
    conj_lon = conj_lon - 360;
elseif conj_lon < -180
    conj_lon = conj_lon + 360;
end

dist_km = spherical_distance(rlat, rlong, conj_lat, conj_lon, Re);
if isempty(Bm) || Blocal <= 0 || ~isfinite(Blocal)
    B_ratio = NaN;
else
    B_ratio = max(Bm) / Blocal;
end
end


function [r_hit, theta_hit, lambda_hit] = refine_surface_hit( ...
    x_lo, y_lo, z_lo, x_hi, y_hi, z_hi, target_r, tol_km)
for i = 1:30
    x_mid = (x_lo + x_hi) / 2;
    y_mid = (y_lo + y_hi) / 2;
    z_mid = (z_lo + z_hi) / 2;
    r_mid = sqrt(x_mid^2 + y_mid^2 + z_mid^2);
    if abs(r_mid - target_r) <= tol_km
        [r_hit, theta_hit, lambda_hit] = xyz_rtp(x_mid, y_mid, z_mid);
        lambda_hit = mod(lambda_hit, 2 * pi);
        return
    end
    if r_mid > target_r
        x_lo = x_mid;
        y_lo = y_mid;
        z_lo = z_mid;
    else
        x_hi = x_mid;
        y_hi = y_mid;
        z_hi = z_mid;
    end
end
[r_hit, theta_hit, lambda_hit] = xyz_rtp(x_mid, y_mid, z_mid);
lambda_hit = mod(lambda_hit, 2 * pi);
end


function [conj_lat, conj_lon, conj_alt, B_ratio, dist_km] = ...
    trace_fieldline_polar_direct(rlat, rlong, altitude_km, rcon, step)
Re = 6371.2;
target_r = Re + rcon;
theta1 = (90 - rlat) * pi / 180;
lambda1 = mod(rlong * pi / 180, 2 * pi);
r1 = Re + altitude_km;
start_t = tic;
time_limit_s = get_step_time_limit(step, 'polar_direct_time_limit_s', 8);

[B_r, B_theta, B_lambda] = IGRF_vs3(r1, theta1, lambda1);
[Bx, By, Bz] = sph_to_cart(B_r, B_theta, B_lambda, theta1, lambda1);
Blocal = norm([Bx By Bz]);
[x1, y1, z1] = rtp_xyz(r1, theta1, lambda1);
B = norm([Bx By Bz]);

probe_step = max(10, min(40, abs(step.base)));
[rp, ~, ~] = xyz_rtp(x1 + probe_step * Bx / B, y1 + probe_step * By / B, z1 + probe_step * Bz / B);
[rm, ~, ~] = xyz_rtp(x1 - probe_step * Bx / B, y1 - probe_step * By / B, z1 - probe_step * Bz / B);
if rp >= rm
    step_sign = 1;
else
    step_sign = -1;
end

step_abs = probe_step;
departure_threshold_km = max(20, min(120, abs(step.base)));
apex_threshold_km = max(200, 4 * departure_threshold_km);
max_r_seen = r1;
has_departed_surface = false;
has_entered_remote_branch = false;
path_length_km = 0;
max_path_length_km = 500000;
stop_tol_km = 1;
Bm = [];

for iter = 1:30000
    if toc(start_t) > time_limit_s
        error('Tracing failed: polar direct fallback time limit exceeded')
    end
    [B_r, B_theta, B_lambda] = IGRF_vs3(r1, theta1, lambda1);
    [Bx, By, Bz] = sph_to_cart(B_r, B_theta, B_lambda, theta1, lambda1);
    B = norm([Bx By Bz]);
    if ~isfinite(B) || B <= 0
        error('Tracing failed: invalid magnetic field magnitude in polar direct fallback')
    end
    Bm(end + 1) = B; %#ok<AGROW>

    [x1, y1, z1] = rtp_xyz(r1, theta1, lambda1);
    prev_r = r1;
    prev_theta = theta1;
    prev_lambda = lambda1;
    prev_remain = r1 - target_r;

    x2 = x1 + step_sign * step_abs * Bx / B;
    y2 = y1 + step_sign * step_abs * By / B;
    z2 = z1 + step_sign * step_abs * Bz / B;
    path_length_km = path_length_km + step_abs;
    if path_length_km > max_path_length_km
        error('Tracing failed: polar direct fallback path length limit exceeded')
    end

    [r1, theta1, lambda1] = xyz_rtp(x2, y2, z2);
    lambda1 = mod(lambda1, 2 * pi);
    remain = r1 - target_r;

    max_r_seen = max(max_r_seen, r1);
    if max_r_seen - target_r >= departure_threshold_km
        has_departed_surface = true;
    end
    if has_departed_surface && max_r_seen - target_r >= apex_threshold_km && r1 < prev_r - 0.1
        has_entered_remote_branch = true;
    end

    if has_entered_remote_branch && prev_remain > 0 && remain <= 0
        [r1, theta1, lambda1] = refine_surface_hit( ...
            x1, y1, z1, x2, y2, z2, target_r, stop_tol_km);
        break
    end

    if has_entered_remote_branch && remain < 100
        step_abs = max(2, min(step_abs, 10));
    elseif has_entered_remote_branch
        step_abs = max(5, min(step_abs, 20));
    end
end

if ~(has_entered_remote_branch && abs(r1 - target_r) <= 2)
    error('Tracing failed: polar direct fallback did not reach remote surface crossing')
end

conj_lat = 90 - theta1 * 180 / pi;
conj_lon = lambda1 * 180 / pi;
conj_alt = rcon;
if conj_lon > 180
    conj_lon = conj_lon - 360;
elseif conj_lon < -180
    conj_lon = conj_lon + 360;
end
dist_km = spherical_distance(rlat, rlong, conj_lat, conj_lon, Re);
if isempty(Bm) || Blocal <= 0 || ~isfinite(Blocal)
    B_ratio = NaN;
else
    B_ratio = max(Bm) / Blocal;
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


function value = get_step_time_limit(step, field_name, default_value)
value = default_value;
if isstruct(step) && isfield(step, field_name)
    candidate = step.(field_name);
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate) && candidate > 0
        value = candidate;
    end
end
end


function [best, reason] = pick_best_candidate(cand1, cand2)
valid1 = isfinite(cand1.dist);
valid2 = isfinite(cand2.dist);

if valid1 && valid2
    if cand2.dist > cand1.dist
        best = cand2;
        reason = 'selected reverse trace as farther valid endpoint';
    else
        best = cand1;
        reason = '';
    end
    return
end

if valid1
    best = cand1;
    reason = 'selected forward trace because reverse trace is invalid';
    return
end

if valid2
    best = cand2;
    reason = 'selected reverse trace because forward trace is invalid';
    return
end

best = cand1;
reason = 'both trace directions are invalid';
end


function step_cfg = normalize_step_config(step)
if ~isstruct(step)
    step_cfg.base = step;
    step_cfg.segment_remain_km = [2000, 1000, 500, 200, 50];
    step_cfg.segment_multipliers = [4.0, 2.5, 1.5, 1.0, 0.5, 0.25];
    step_cfg.far_multiplier = 2.5;
    step_cfg.mid_multiplier = 1.0;
    step_cfg.near_multiplier = 0.25;
    step_cfg.far_remain_km = 1000;
    step_cfg.near_remain_km = 200;
    step_cfg.min_abs = 1;
    step_cfg.max_abs = 50;
    return
end

step_cfg = step;

if ~isfield(step_cfg, 'base')
    error('Step config must include "base"')
end
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


function step_now = choose_step(step_cfg, remain_km)
if ~isempty(step_cfg.segment_remain_km) || ~isempty(step_cfg.segment_multipliers)
    [thresholds, multipliers] = normalize_segment_arrays(step_cfg);
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
end


function step_now = apply_polar_step_limit(step_now, theta, step_cfg)
% High-latitude tracing is more sensitive because azimuth changes rapidly
% near the magnetic pole. Shrink steps progressively as theta approaches 0 or pi.
polar_scale = max(0.1, min(1.0, abs(sin(theta)) / 0.7));
step_abs = min(step_cfg.max_abs, max(step_cfg.min_abs, abs(step_now) * polar_scale));
step_now = sign(step_now) * step_abs;
end


function [thresholds, multipliers] = normalize_segment_arrays(step_cfg)
thresholds = step_cfg.segment_remain_km(:).';
multipliers = step_cfg.segment_multipliers(:).';

if isempty(thresholds) && isempty(multipliers)
    thresholds = [step_cfg.far_remain_km, step_cfg.near_remain_km];
    multipliers = [step_cfg.far_multiplier, step_cfg.mid_multiplier, step_cfg.near_multiplier];
    return
end

if numel(multipliers) ~= numel(thresholds) + 1
    error('Step config segment_multipliers must have exactly one more element than segment_remain_km')
end

[thresholds, order] = sort(thresholds, 'descend');
multipliers = [multipliers(order), multipliers(end)];
end


function step_out = reverse_step(step_in)
if isstruct(step_in)
    step_out = step_in;
    step_out.base = -step_out.base;
else
    step_out = -step_in;
end
end


function stop_tol_km = get_stop_tolerance(step_cfg)
% Use a tight fixed altitude tolerance so tracing reaches the conjugate
% endpoint instead of stopping tens of kilometers early when base steps are large.
stop_tol_km = 1;
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
