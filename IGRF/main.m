clear
clc

% Input terminator file exported by Python script
script_dir = fileparts(mfilename('fullpath'));
code_dir = fileparts(script_dir);
input_file = fullfile(code_dir, 'terminators', 'terminator_20260125_19.50_109.00.txt');

% IGRF parameters
alt = 0.1;      % start altitude (km)
conj_alt = 0;   % conjugate target altitude (km)
step.base = 200;               % base tracing step (km)
step.segment_remain_km = [2000, 1000, 500, 200, 50];             %N个分段
step.segment_multipliers = [4.0, 2.5, 1.5, 1.0, 0.5, 0.25];      %N+1个步长倍率
step.min_abs = 1;             % lower clamp for adaptive step magnitude (km)
step.max_abs = 1000;            % upper clamp for adaptive step magnitude (km)
step.retry_multipliers = [1, -1, 2, -2];
step.small_retry_multipliers = [0.5, -0.5, 0.25, -0.25];
step.low_lat_retry_multipliers = [1, -1, 2, -2, 0.5, -0.5, 0.25, -0.25];
step.mid_lat_retry_multipliers = [1, -1, 2, -2, 0.5, -0.5, 0.25, -0.25];
step.high_lat_retry_multipliers = [1, -1];
step.mid_lat_threshold = 55;
step.high_lat_threshold = 65;
step.fill_retry_lat_limit = 65;
sample_stride = 1;  % read one point every N points

[sections, section_names] = read_terminator_sections(input_file);

if numel(sections) < 2
    error('Need at least two sections (morning/evening) in: %s', input_file)
end

run_date = datestr(now, 'yyyymmdd');
out_dir = fullfile(script_dir, run_date);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fprintf('Parsed %d sections\n', numel(sections));
fprintf('Section 1 points: %d\n', size(sections{1}, 1));
fprintf('Section 2 points: %d\n', size(sections{2}, 1));

% Section #1 as morning, #2 as evening
morning_points = sections{1};
evening_points = sections{2};

% Interval sampling: keep one point every "sample_stride" points
morning_points = downsample_points(morning_points, sample_stride);
evening_points = downsample_points(evening_points, sample_stride);

fprintf('After sampling (stride=%d): morning=%d, evening=%d\n', ...
    sample_stride, size(morning_points, 1), size(evening_points, 1));
fprintf('Adaptive step config: base=%.2f km, thresholds=%s, multipliers=%s, retries=%s\n', ...
    step.base, mat2str(step.segment_remain_km), mat2str(step.segment_multipliers), mat2str(step.retry_multipliers));

fprintf('\nTracing morning section...\n');
morning_result = trace_points_with_conjugate(morning_points, alt, conj_alt, step, 'morning');
morning_result = fill_missing_conjugates(morning_result, morning_points, alt, conj_alt, step, 'morning');
save(fullfile(out_dir, 'conjugate_morning.mat'), 'morning_result');
write_result_txt(fullfile(out_dir, 'conjugate_morning.txt'), morning_result, section_names{1});

fprintf('\nTracing evening section...\n');
evening_result = trace_points_with_conjugate(evening_points, alt, conj_alt, step, 'evening');
evening_result = fill_missing_conjugates(evening_result, evening_points, alt, conj_alt, step, 'evening');

morning_mat = fullfile(out_dir, 'conjugate_morning.mat');
evening_mat = fullfile(out_dir, 'conjugate_evening.mat');
save(morning_mat, 'morning_result');
save(evening_mat, 'evening_result');

morning_txt = fullfile(out_dir, 'conjugate_morning.txt');
evening_txt = fullfile(out_dir, 'conjugate_evening.txt');
write_result_txt(morning_txt, morning_result, section_names{1});
write_result_txt(evening_txt, evening_result, section_names{2});

fprintf('\nDone. Output files:\n');
fprintf('  %s\n', morning_mat);
fprintf('  %s\n', evening_mat);
fprintf('  %s\n', morning_txt);
fprintf('  %s\n', evening_txt);


function [sections, section_names] = read_terminator_sections(path)
fid = fopen(path, 'r');
if fid < 0
    error('Cannot open input file: %s', path)
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
lines = raw{1};

sections = {};
section_names = {};
current = [];
current_name = '';
in_section = false;

for i = 1:numel(lines)
    line = strtrim(lines{i});
    if isempty(line)
        continue
    end

    % Section header examples:
    % [閺呫劎鍤嶿
    % [閺呫劎鍤嶿 points=720
    if startsWith(line, '[')
        if in_section && ~isempty(current)
            sections{end + 1} = current; %#ok<AGROW>
            section_names{end + 1} = current_name; %#ok<AGROW>
        end
        current = [];
        current_name = line;
        in_section = true;
        continue
    end

    if ~in_section
        continue
    end

    % Data line format:
    %   index,lon_deg,lat_deg
    %   index,lon_deg,lat_deg,x,y,z
    vals = sscanf(line, '%f,%f,%f,%f,%f,%f');
    if numel(vals) == 3 || numel(vals) == 6
        idx = vals(1);
        lon = vals(2);
        lat = vals(3);
        xyz = [NaN, NaN, NaN];
        if numel(vals) == 6
            xyz = vals(4:6).';
        end
        if isfinite(idx) && isfinite(lon) && isfinite(lat)
            current(end + 1, :) = [idx, lat, lon, xyz]; %#ok<AGROW>
        end
    end
end

if in_section && ~isempty(current)
    sections{end + 1} = current;
    section_names{end + 1} = current_name;
end

if isempty(sections)
    error('No curve points parsed from file: %s', path)
end
end


function result = trace_points_with_conjugate(points, alt, conj_alt, step, section_label)
n = size(points, 1);
result = zeros(n, 10);
ok_count = 0;
fail_count = 0;
update_interval = max(1, ceil(n / 20));  % about every 5%
t0 = tic;

fprintf('[%s] start: total points = %d\n', section_label, n);

for i = 1:n
    idx = points(i, 1);
    lat = points(i, 2);
    lon = points(i, 3);
    xyz = points(i, 4:6);

    [ok, clat, clon, calt, Bratio, used_step, err_msg] = solve_conjugate_safe(lat, lon, alt, conj_alt, step);

    if ok
        % [idx, input_lat, input_lon, input_x, input_y, input_z, conj_lat, conj_lon, conj_alt_km, B_ratio]
        result(i, :) = [idx, lat, lon, xyz, clat, clon, calt, Bratio];
        ok_count = ok_count + 1;
    else
        % Keep pipeline running even if one point fails.
        result(i, :) = [idx, lat, lon, xyz, NaN, NaN, NaN, NaN];
        fail_count = fail_count + 1;
        fprintf('[%s] warning: idx=%d lat=%.6f lon=%.6f failed (%s)\n', ...
            section_label, round(idx), lat, lon, err_msg);
    end

    if mod(i, update_interval) == 0 || i == n
        elapsed = toc(t0);
        avg = elapsed / i;
        remain = avg * (n - i);
        pct = 100.0 * i / n;
        fprintf('[%s] progress: %d/%d (%.1f%%), ok=%d, fail=%d, last_step=%.3f, elapsed=%.1fs, eta=%.1fs\n', ...
            section_label, i, n, pct, ok_count, fail_count, used_step, elapsed, remain);
    end
end

fprintf('[%s] done: ok=%d, fail=%d, total_elapsed=%.1fs\n', ...
    section_label, ok_count, fail_count, toc(t0));
end


function [ok, clat, clon, calt, Bratio, used_step, err_msg] = solve_conjugate_safe(lat, lon, alt, conj_alt, step)
ok = false;
clat = NaN;
clon = NaN;
calt = NaN;
Bratio = NaN;
used_step = get_step_base(step);
err_msg = '';

candidates = build_step_candidates_for_lat(step, lat);

for k = 1:numel(candidates)
    trial_step = candidates(k);
    if abs(trial_step) < 1e-6
        continue
    end

    try
        [clat, clon, calt, Bratio] = con_mag(lat, lon, alt, conj_alt, build_trial_step(step, trial_step));
        ok = true;
        used_step = trial_step;
        return
    catch ME
        err_msg = ME.message;
        used_step = trial_step;
    end
end
end


function base = get_step_base(step)
if isstruct(step)
    base = step.base;
else
    base = step;
end
end


function candidates = build_step_candidates(step)
base = get_step_base(step);

if isstruct(step) && isfield(step, 'retry_multipliers') && ~isempty(step.retry_multipliers)
    candidates = base .* step.retry_multipliers;
    if isfield(step, 'small_retry_multipliers') && ~isempty(step.small_retry_multipliers)
        candidates = [candidates, base .* step.small_retry_multipliers];
    end
else
    candidates = [base, -base, base / 2, -base / 2, base / 4, -base / 4, base * 2, -base * 2, base * 4, -base * 4];
end

candidates = unique(candidates, 'stable');
end


function candidates = build_step_candidates_for_lat(step, lat)
base = get_step_base(step);
abs_lat = abs(lat);

if ~isstruct(step)
    candidates = [base, -base];
    return
end

if isfield(step, 'high_lat_threshold') && abs_lat >= step.high_lat_threshold
    multipliers = get_step_multipliers(step, 'high_lat_retry_multipliers', [1, -1]);
elseif isfield(step, 'mid_lat_threshold') && abs_lat >= step.mid_lat_threshold
    multipliers = get_step_multipliers(step, 'mid_lat_retry_multipliers', [1, -1, 0.5, -0.5]);
else
    multipliers = get_step_multipliers(step, 'low_lat_retry_multipliers', [1, -1]);
end

candidates = unique(base .* multipliers, 'stable');
end


function multipliers = get_step_multipliers(step, field_name, default_value)
multipliers = default_value;
if isfield(step, field_name) && ~isempty(step.(field_name))
    multipliers = step.(field_name);
end
end


function step_out = build_trial_step(step_in, trial_base)
if isstruct(step_in)
    step_out = step_in;
    step_out.base = trial_base;
else
    step_out = trial_base;
end
end


function write_result_txt(path, result, section_name)
fid = fopen(path, 'w');
if fid < 0
    error('Cannot write file: %s', path)
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'section=%s\n', section_name);
fprintf(fid, 'columns: idx,input_lat,input_lon,input_x,input_y,input_z,conj_lat,conj_lon,conj_alt_km,B_ratio\n');

for i = 1:size(result, 1)
    fprintf(fid, '%d,%.6f,%.6f,%.10f,%.10f,%.10f,%.6f,%.6f,%.6f,%.8f\n', ...
        round(result(i, 1)), result(i, 2), result(i, 3), ...
        result(i, 4), result(i, 5), result(i, 6), ...
        result(i, 7), result(i, 8), result(i, 9), result(i, 10));
end
end


function points_out = downsample_points(points_in, stride)
if stride <= 1
    points_out = points_in;
    return
end

idx = 1:stride:size(points_in, 1);

% Ensure the last point is preserved.
if idx(end) ~= size(points_in, 1)
    idx = [idx, size(points_in, 1)];
end

points_out = points_in(idx, :);
end


function result = fill_missing_conjugates(result, points, alt, conj_alt, step, section_label)
missing = isnan(result(:, 7)) | isnan(result(:, 8)) | isnan(result(:, 9)) | isnan(result(:, 10));
if any(missing)
    retry_step = step;
    retry_step.base = max(get_step_base(step) / 4, step.min_abs);
    if isfield(step, 'small_retry_multipliers') && ~isempty(step.small_retry_multipliers)
        retry_step.retry_multipliers = step.small_retry_multipliers;
    else
        retry_step.retry_multipliers = [1, -1, 0.5, -0.5];
    end

    retried = 0;
    missing_idx = find(missing);
    for j = 1:numel(missing_idx)
        i = missing_idx(j);
        lat = points(i, 2);
        if should_skip_fill_retry(lat, step)
            continue
        end
        lon = points(i, 3);
        xyz = points(i, 4:6);
        [ok, clat, clon, calt, Bratio] = solve_conjugate_safe(lat, lon, alt, conj_alt, retry_step);
        if ok
            result(i, :) = [points(i, 1), lat, lon, xyz, clat, clon, calt, Bratio];
            retried = retried + 1;
        end
    end

    if retried > 0
        fprintf('[%s] recovered %d missing rows with small-step retry\n', section_label, retried);
    end
end

missing = isnan(result(:, 7)) | isnan(result(:, 8)) | isnan(result(:, 9)) | isnan(result(:, 10));
if ~any(missing)
    return
end
fprintf('[%s] warning: %d rows remain unresolved after small-step retry\n', section_label, nnz(missing));
end


function tf = should_skip_fill_retry(lat, step)
tf = false;
if isstruct(step) && isfield(step, 'fill_retry_lat_limit')
    tf = abs(lat) >= step.fill_retry_lat_limit;
end
end
