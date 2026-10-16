function result = MsHe(glow_file, opts)
%MSHE Run metastable helium forward model from one GLOW output file.
%   result = MsHe(glow_file, opts)
%   - glow_file: full/relative path of one output_vsXXX.txt file
%   - opts: optional struct
%       opts.usehu (default 1000 km)
%       opts.usehd (default 100 km)
%       opts.interval (default 1 km)
%       opts.whpp (default 0)
%       opts.whre (default 1)
%       opts.make_plot (default true)
%
%   Backward-compatible usage:
%       MsHe
%   runs the historical default case under ../glow_data/HaiNan20240106/006.
%
%   物理流程：
%       GLOW输入 -> 背景大气 -> 产生/损失/扩散
%       -> 亚稳态氦尺度高度 -> He(2^3S)稳态密度廓线

    if nargin < 2 || isempty(opts)
        opts = struct();
    end

    opts = apply_defaults(opts, struct(...
        'usehu', 1000, ...
        'usehd', 100, ...
        'interval', 1, ...
        'whpp', 0, ...
        'whre', 1, ...
        'make_plot', true));

    if nargin < 1 || isempty(glow_file)
        script_dir = fileparts(mfilename('fullpath'));
        glow_file = fullfile(script_dir, '..', 'glow_data', 'HaiNan20240106', '006', 'output_vs050.txt');
    end

    if ~isfile(glow_file)
        error('Glow data file not found: %s', glow_file);
    end

    % Height grid (km)
    nz = (opts.usehu - opts.usehd) / opts.interval + 1;
    z = linspace(opts.usehd, opts.usehu, nz);

    [heion, elecden, electem, jmax, sza, idate, ut, glat, glong, ...
        f107a, f107, f107p, ap, ef, ec, ener, del, eflx] = read_glowdata(z, glow_file);

    % 将时间信息转换为年积日和地方时，供大气参数计算与结果诊断使用。
    [DoY, Yr] = date2doy(idate);
    lt_hour = UtToLt(glong, ut / 3600);

    % 由经验大气模型获取背景中性大气成分和温度。
    [Tem, He, O1, N2, O2, den] = gain_atmdentem(z, glat, glong, Yr, DoY, ut, f107a, f107, ap);

    % 亚稳态氦产生项：
    % qe 为光电子碰撞基态氦产生的速率，
    % qr 为氦离子与电子复合产生的速率。
    qe = he_produc1_vs2(eflx, ener, He);
    qr = he_produc2(heion, elecden, electem);
    q = qe + qr;

    % 损失率、扩散系数和有效尺度高度共同决定
    % He(2^3S) 稳态平衡方程中的输运与耗散项。
    Loss = he_loss(O1, O2, N2, Tem);
    Dif = he_dif(He, O1, Tem);
    [Hn, a, gstar, gam] = he_Hn1vs3(z, nz, Tem, sza, glat, opts.whpp);

    % 求解一维稳态扩散-损失方程，得到 He(2^3S) 密度廓线。
    n_he = he_denvs6(z, Dif, Hn, Loss, q);

    [max_n_he, max_idx] = max(n_he);
    max_altitude = z(max_idx);

    result = struct();
    result.input = struct('glow_file', glow_file, 'options', opts);
    % 保存复现实验所需的运行条件和输入元信息。
    result.meta = struct(...
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
    % 保存求解过程中用到的原始输入和中间大气场。
    result.raw = struct(...
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

    % 保存正演得到的各类廓线结果和统计量。
    result.profile = struct(...
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
        'n_he_max_cm3', max_n_he, ...
        'n_he_max_alt_km', max_altitude);

    if opts.make_plot
        result.figure_handle = plot_diagnostics(z, eflx, ener, q, Dif, n_he, result.meta, max_n_he, max_altitude);
    else
        result.figure_handle = [];
    end
end

function fig = plot_diagnostics(z, eflx, ener, q, Dif, n_he, meta, max_n_he, max_altitude)
    % 诊断图：电子通量谱、总产生率、扩散系数、
    % 以及最终的亚稳态氦密度廓线。
    fig = figure('Position', [100, 100, 1000, 800], 'Color', 'w');

    subplot(2,2,1);
    hold on;
    colors = jet(46);
    for i = 35:80
        plot(eflx(:, i), z, 'Color', colors(i - 34, :));
    end
    xlabel('eflx [cm^{-2}s^{-1}ev^{-1}]');
    ylabel('z [km]');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'FontWeight', 'Bold', 'linewidth', 1);
    set(gca, 'XScale', 'log');
    colormap(colors);
    tick_indices = round(linspace(35, 80, 10));
    tick_values = ener(tick_indices);
    c = colorbar;
    c.Ticks = linspace(0, 1, 10);
    c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), tick_values, 'UniformOutput', false);
    c.Label.String = 'eV';
    c.Label.FontSize = 10;
    c.Label.FontWeight = 'bold';
    c.Label.Rotation = 0;
    c.Label.VerticalAlignment = 'bottom';
    c.Label.Position = [0.5 1.02 0];
    text(-0.1, 1.1, '(a)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    hold off;

    subplot(2,2,2);
    plot(q, z, 'linewidth', 1.2);
    xlabel('metstable excitation rate:qe [cm^{-3}s^{-1}]');
    ylabel('z [km]');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'FontWeight', 'Bold', 'linewidth', 1);
    text(-0.1, 1.1, '(b)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

    subplot(2,2,3);
    plot(Dif, z, 'linewidth', 1.2);
    xlabel('The diffusion coefficient:D [cm^{2}s^{-1}]');
    ylabel('z [km]');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'FontWeight', 'Bold', 'linewidth', 1);
    text(-0.1, 1.1, '(c)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

    subplot(2,2,4);
    plot(n_he, z, 'b-', 'linewidth', 1.2);
    hold on;
    plot(max_n_he, max_altitude, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    xline(0, '--k', 'LineWidth', 0.8);
    xlim([1e-20, 1.0]);
    ylim([100, 1000]);
    xlabel('He(2^3S) numebr density [cm^{-3}]');
    ylabel('Altitude [km]');
    text(0.95 * max(xlim), 0.8 * max(ylim), sprintf(['sza=%.2fdeg\n' ...
        'Ap:%.2f\n' ...
        'F10.7:%.2f\n' ...
        'F10.7a:%.2f\n' ...
        'F10.7p:%.2f'], ...
        meta.sza_deg, meta.ap, meta.f107, meta.f107a, meta.f107p), ...
        'HorizontalAlignment', 'right', ...
        'FontSize', 10, ...
        'FontName', 'Times New Roman', ...
        'EdgeColor', 'black', ...
        'LineWidth', 1, ...
        'BackgroundColor', 'white');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'FontWeight', 'Bold', 'linewidth', 1);
    text(-0.1, 1.1, '(d)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
    hold off;
end

function s = apply_defaults(s, defaults)
    % 仅为缺失或为空的参数填入默认值，保留调用者显式传入的设置。
    keys = fieldnames(defaults);
    for k = 1:numel(keys)
        key = keys{k};
        if ~isfield(s, key) || isempty(s.(key))
            s.(key) = defaults.(key);
        end
    end
end

