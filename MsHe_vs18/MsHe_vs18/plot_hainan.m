% Plot Hainan 2025-01-19 NDMH_L2A heatmap (time vs altitude).
clc; clear;

script_dir = fileparts(mfilename('fullpath'));
data_dir = fullfile(script_dir, '..', '..', 'hainan_data', '20250323', 'ODAZH_ALAL01_NDMH_L2A');

if ~isfolder(data_dir)
    error('Data folder not found: %s', data_dir);
end

files = dir(fullfile(data_dir, '*.TXT'));
if isempty(files)
    files = dir(fullfile(data_dir, '*.txt'));
end
if isempty(files)
    error('No TXT files found in: %s', data_dir);
end

name0 = files(1).name;
tok = regexp(name0, '^(.*)_\d{14}_V\d+\.\d+\.txt$', 'tokens', 'once', 'ignorecase');
if isempty(tok)
    [plot_name, ~] = fileparts(name0);
else
    plot_name = tok{1};
end
date_label = char(regexp(plot_name, '(\d{8})$', 'match', 'once'));

nf = numel(files);
time_vec = NaT(nf, 1);
z_all = [];
z_cell = cell(nf, 1);
den_cell = cell(nf, 1);

for i = 1:nf
    fpath = fullfile(files(i).folder, files(i).name);
    [time_i, z_i, den_i] = i_read_one_ndmh_txt(fpath);
    time_vec(i) = time_i;
    z_cell{i} = z_i;
    den_cell{i} = den_i;
    z_all = union(z_all, z_i);
end

% sort by observation time
[~, ord] = sort(time_vec);
time_vec = time_vec(ord);
z_cell = z_cell(ord);
den_cell = den_cell(ord);
if any(~isnat(time_vec))
    date_label = datestr(time_vec(find(~isnat(time_vec), 1, 'first')), 'yyyy-mm-dd');
end

z_km = z_all(:);
nt = numel(time_vec);
nz = numel(z_km);
den_2d = nan(nt, nz);

for i = 1:nt
    [tf, idx] = ismember(z_cell{i}, z_km);
    den_2d(i, idx(tf)) = den_cell{i}(tf);
end

% Keep original matrix, then build a filled matrix for plotting.
den_2d_raw = den_2d;
den_2d_fill = den_2d_raw;

% Fill NaN by local mean (time direction first, then altitude direction).
den_2d_fill = fillmissing(den_2d_fill, 'movmean', 3, 1, 'EndValues', 'nearest');
den_2d_fill = fillmissing(den_2d_fill, 'movmean', 3, 2, 'EndValues', 'nearest');

% Fallback for any remaining NaN: use global mean of valid values.
if any(isnan(den_2d_fill), 'all')
    m = mean(den_2d_fill(~isnan(den_2d_fill)), 'omitnan');
    den_2d_fill(isnan(den_2d_fill)) = m;
end

figure('Color', 'w', 'Position', [120, 100, 980, 520]);
t_hour = hour(time_vec) + minute(time_vec) / 60 + second(time_vec) / 3600;
z_num = z_km(:);
den_grid = den_2d_fill.'; % [z x t]
[UT, Z] = meshgrid(t_hour, z_num);
contourf(UT, Z, den_grid,10, 'LineStyle', 'none');
set(gca, 'YDir', 'normal');
xlim([0, 24]);
xticks(0:2:24);
xlabel('UT Time');
ylabel('Altitude (km)');
title(date_label, 'Interpreter', 'none');
colormap(turbo);
cb = colorbar;
ylabel(cb, 'Density (cm^{-3})');

grid on;
box on;
set(gca, 'FontSize', 11, 'LineWidth', 1);

function [obs_time, z_km, den_cm3] = i_read_one_ndmh_txt(fpath)
obs_time = NaT;
fid = fopen(fpath, 'r');
if fid < 0
    error('Cannot open file: %s', fpath);
end

cleaner = onCleanup(@() fclose(fid));
z_km = [];
den_cm3 = [];
in_table = false;

while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end

    if startsWith(line, '#DataTime:')
        tstr = strtrim(strrep(line, '#DataTime:', ''));
        try
            obs_time = datetime(tstr, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
        catch
            obs_time = NaT;
        end
    end

    if contains(line, 'Elevation') && contains(line, 'Density') && contains(line, 'DenEr')
        in_table = true;
        continue;
    end

    if in_table
        vals = sscanf(line, '%f %f %f');
        if numel(vals) == 3
            z_km(end+1, 1) = vals(1); %#ok<AGROW>
            den_cm3(end+1, 1) = vals(2); %#ok<AGROW>
        end
    end
end

% replace file-defined missing values with NaN
z_km(z_km <= -9.9e7) = NaN;
den_cm3(den_cm3 <= -9.9e11) = NaN;
den_cm3(den_cm3 < 0) = NaN;

% keep only 200-1000 km
keep = isfinite(z_km) & (z_km >= 200) & (z_km <= 1000);
z_km = z_km(keep);
den_cm3 = den_cm3(keep);
end
