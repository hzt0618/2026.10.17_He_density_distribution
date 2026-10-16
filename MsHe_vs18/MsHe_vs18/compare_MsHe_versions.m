function comparison = compare_MsHe_versions(glow_file)
%COMPARE_MSHE_VERSIONS Compare modern and legacy MsHe results on one GLOW file.
%   comparison = compare_MsHe_versions(glow_file)
%
% This script does not modify either source MsHe.m file. It:
%   1) runs the modern function-based MsHe.m
%   2) replays the legacy calculation chain with the legacy folder's functions
%   3) compares qe, qr, q, Loss, Dif, Hn, and n_he on the same input file
%
% If glow_file is omitted, the script searches a few likely locations.

clearvars -except glow_file
clc;

script_dir = fileparts(mfilename('fullpath'));
modern_dir = script_dir;
legacy_dir = 'D:\1CH2\1study\DXK\zhanshi\鏄熻浇\15_绋嬪簭闆嗗悎\MsHe_vs18\MsHe_vs18';

if nargin < 1 || isempty(glow_file)
    glow_file = find_default_glow_file(script_dir);
end

if ~isfile(glow_file)
    error('Glow file not found: %s', glow_file);
end

opts = struct( ...
    'usehu', 1000, ...
    'usehd', 100, ...
    'interval', 1, ...
    'whpp', 0, ...
    'whre', 1, ...
    'make_plot', false);

fprintf('Comparing MsHe versions with input file:\n%s\n\n', glow_file);

modern_result = run_modern_wrapper(modern_dir, legacy_dir, glow_file, opts);
modern_core = run_core_chain_in_folder(modern_dir, glow_file, opts);
legacy_core = run_core_chain_in_folder(legacy_dir, glow_file, opts);

comparison = struct();
comparison.glow_file = glow_file;
comparison.options = opts;
comparison.modern_result = modern_result;
comparison.modern_core = modern_core;
comparison.legacy_core = legacy_core;
comparison.wrapper_vs_modern_core = compare_results(modern_result.profile, modern_core.profile);
comparison.modern_vs_legacy = compare_results(modern_core.profile, legacy_core.profile);

fprintf('Modern wrapper vs modern core replay:\n');
print_metrics(comparison.wrapper_vs_modern_core);

fprintf('\nModern core vs legacy core:\n');
print_metrics(comparison.modern_vs_legacy);

fprintf('\nPeak n_he comparison:\n');
fprintf('  modern: %.12e cm^-3 at %.1f km\n', ...
    modern_core.profile.n_he_max_cm3, modern_core.profile.n_he_max_alt_km);
fprintf('  legacy: %.12e cm^-3 at %.1f km\n', ...
    legacy_core.profile.n_he_max_cm3, legacy_core.profile.n_he_max_alt_km);

assignin('base', 'comparison', comparison);
fprintf('\nComparison struct has been assigned to workspace variable: comparison\n');
end

function result = run_modern_wrapper(modern_dir, legacy_dir, glow_file, opts)
    cleanup = prioritize_modern_mshe(modern_dir, legacy_dir);
    result = MsHe(glow_file, opts);
    clear MsHe
    clear cleanup
end

function glow_file = find_default_glow_file(script_dir)
    candidates = {
        fullfile(script_dir, '..', '..', 'glow model', 'glow model', 'output', ...
            'HaiNan20240106', '6Module', 'output_vs50.txt')
        fullfile(script_dir, '..', '..', '..', '15_绋嬪簭闆嗗悎', 'glow model', 'glow model', ...
            'HaiNan2021_all', '1', 'output_vs50.txt')
        fullfile(script_dir, '..', '..', '..', '15_绋嬪簭闆嗗悎', 'glow model', 'glow model', ...
            'HaiNan2021_all', '1', 'output_vs1.txt')
        };

    glow_file = '';
    for i = 1:numel(candidates)
        if isfile(candidates{i})
            glow_file = candidates{i};
            return;
        end
    end

    error('No default GLOW file found. Please pass glow_file explicitly.');
end

function result = run_core_chain_in_folder(folder, glow_file, opts)
    cleanup = push_folder_to_path(folder);

    z = opts.usehd:opts.interval:opts.usehu;
    nz = numel(z);

    [heion, elecden, electem, jmax, sza, idate, ut, glat, glong, ...
        f107a, f107, f107p, ap, ef, ec, ener, del, eflx] = read_glowdata(z, glow_file);

    [DoY, Yr] = date2doy(idate);
    lt_hour = UtToLt(glong, ut / 3600);
    [Tem, He, O1, N2, O2, den] = gain_atmdentem(z, glat, glong, Yr, DoY, ut, f107a, f107, ap);

    qe = he_produc1_vs2(eflx, ener, He);
    qr = he_produc2(heion, elecden, electem);
    q = qe + qr;

    Loss = he_loss(O1, O2, N2, Tem);
    Dif = he_dif(He, O1, Tem);
    [Hn, a, gstar, gam] = he_Hn1vs3(z, nz, Tem, sza, glat, opts.whpp);
    n_he = he_denvs6(z, Dif, Hn, Loss, q);

    [n_he_max_cm3, max_idx] = max(n_he);

    result = struct();
    result.meta = struct( ...
        'jmax', jmax, ...
        'sza_rad', sza, ...
        'sza_deg', sza * 180 / pi, ...
        'idate', idate, ...
        'ut_s', ut, ...
        'lt_hour', lt_hour, ...
        'glat', glat, ...
        'glong', glong, ...
        'f107a', f107a, ...
        'f107', f107, ...
        'f107p', f107p, ...
        'ap', ap, ...
        'ef', ef, ...
        'ec', ec, ...
        'doy', DoY, ...
        'year', Yr);
    result.grid = struct('z_km', z, 'nz', nz);
    result.raw = struct( ...
        'heion', heion, ...
        'elecden', elecden, ...
        'electem', electem, ...
        'ener_ev', ener, ...
        'del_ev', del, ...
        'eflx', eflx, ...
        'Tem', Tem, ...
        'He', He, ...
        'O1', O1, ...
        'N2', N2, ...
        'O2', O2, ...
        'den', den);
    result.profile = struct( ...
        'qe', qe, ...
        'qr', qr, ...
        'q', q, ...
        'Loss', Loss, ...
        'Dif', Dif, ...
        'Hn_cm', Hn, ...
        'a', a, ...
        'gstar', gstar, ...
        'gam', gam, ...
        'n_He_cm3', n_he, ...
        'n_he_max_cm3', n_he_max_cm3, ...
        'n_he_max_alt_km', z(max_idx));

    clear_functions_for_switch();
    clear cleanup
end

function cleanup = push_folder_to_path(folder)
    clear_functions_for_switch();
    addpath(folder, '-begin');
    rehash;
    cleanup = onCleanup(@() pop_folder_from_path(folder));
end

function cleanup = prioritize_modern_mshe(modern_dir, legacy_dir)
    clear MsHe

    had_legacy = contains(path, [legacy_dir pathsep]);
    legacy_mshe = fullfile(legacy_dir, 'MsHe.m');
    modern_mshe = fullfile(modern_dir, 'MsHe.m');
    if strcmp(which('MsHe'), legacy_mshe)
        had_legacy = true;
    end

    if had_legacy
        rmpath(legacy_dir);
    end

    addpath(modern_dir, '-begin');
    rehash;

    cleanup = onCleanup(@() restore_mshe_path(modern_dir, legacy_dir, had_legacy));
end

function restore_mshe_path(modern_dir, legacy_dir, had_legacy)
    clear MsHe
    modern_mshe = fullfile(modern_dir, 'MsHe.m');
    if contains(path, [modern_dir pathsep]) || strcmp(which('MsHe'), modern_mshe)
        rmpath(modern_dir);
    end
    if had_legacy
        addpath(legacy_dir, '-begin');
    end
    rehash;
end

function pop_folder_from_path(folder)
    rmpath(folder);
    clear_functions_for_switch();
    rehash;
end

function clear_functions_for_switch()
    clear read_glowdata date2doy UtToLt gain_atmdentem
    clear he_produc1_vs2 he_produc2 he_loss he_dif he_Hn1vs3 he_denvs6
end

function metrics = compare_results(lhs_profile, rhs_profile)
    field_map = {
        'qe', 'qe'
        'qr', 'qr'
        'q', 'q'
        'Loss', 'Loss'
        'Dif', 'Dif'
        'Hn_cm', 'Hn_cm'
        'n_He_cm3', 'n_He_cm3'
        };

    metrics = struct();
    for i = 1:size(field_map, 1)
        label = field_map{i, 1};
        lhs = lhs_profile.(field_map{i, 2});
        rhs = rhs_profile.(field_map{i, 2});
        delta = lhs - rhs;
        denom = max(abs(lhs), abs(rhs));
        denom(denom == 0) = 1;

        metrics.(label) = struct( ...
            'max_abs_diff', max(abs(delta(:))), ...
            'max_rel_diff', max(abs(delta(:)) ./ denom(:)), ...
            'rms_diff', sqrt(mean(delta(:).^2)));
    end
end

function print_metrics(metrics)
    labels = fieldnames(metrics);
    for i = 1:numel(labels)
        item = metrics.(labels{i});
        fprintf('  %-8s max_abs=% .6e   max_rel=% .6e   rms=% .6e\n', ...
            labels{i}, item.max_abs_diff, item.max_rel_diff, item.rms_diff);
    end
end

