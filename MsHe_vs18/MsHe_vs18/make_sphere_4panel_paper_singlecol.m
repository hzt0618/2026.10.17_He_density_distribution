% make_sphere_4panel_paper_singlecol.m
% Single-column version of the four seasonal projected-sphere figure.

clear;
clc;

script_dir = fileparts(mfilename('fullpath'));
out_root = fullfile(script_dir, '..', 'output');

panels = struct( ...
    'fig_path', { ...
        fullfile(out_root, '20250320_mor_sth_05', 'term_sphere_nhe_20250320_vs006.fig'), ...
        fullfile(out_root, '20250620_mor_sth_05', 'term_sphere_nhe_20250620_vs006.fig'), ...
        fullfile(out_root, '20250920_mor_sth_05', 'term_sphere_nhe_20250920_vs006.fig'), ...
        fullfile(out_root, '20251221_mor_sth_05', 'term_sphere_nhe_20251221_vs006.fig')}, ...
    'title1', {'Spring equinox', 'Summer solstice', 'Autumn equinox', 'Winter solstice'}, ...
    'title2', {'20 Mar 2025', '20 Jun 2025', '20 Sep 2025', '21 Dec 2025'}, ...
    'label', {'(a)', '(b)', '(c)', '(d)'} ...
);

for k = 1:numel(panels)
    if ~isfile(panels(k).fig_path)
        error('Missing source figure: %s', panels(k).fig_path);
    end
end

fig_font = 'Times New Roman';
label_font_pt = 9;
title_font_pt = 9;
panel_font_pt = 9;
tick_font_pt = 8;

src = struct([]);
global_clim = [inf, -inf];
for k = 1:numel(panels)
    f = openfig(panels(k).fig_path, 'invisible');
    ax = find_main_axes(f);
    delete(findall(ax, 'Type', 'text', 'String', 'hainan'));

    src(k).fig = f; %#ok<SAGROW>
    src(k).ax = ax; %#ok<SAGROW>
    src(k).xlim = xlim(ax); %#ok<SAGROW>
    src(k).ylim = ylim(ax); %#ok<SAGROW>
    src(k).clim = caxis(ax); %#ok<SAGROW>
    global_clim(1) = min(global_clim(1), src(k).clim(1));
    global_clim(2) = max(global_clim(2), src(k).clim(2));
end
if ~all(isfinite(global_clim)) || global_clim(2) <= global_clim(1)
    global_clim = [0, 1];
end
global_clim(1) = min(0, global_clim(1));

fig_w = 3.46;
fig_h = 4.65;
paper_fig = figure( ...
    'Color', 'w', ...
    'Units', 'inches', ...
    'Position', [0.5, 0.5, fig_w, fig_h], ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0, 0, fig_w, fig_h], ...
    'PaperSize', [fig_w, fig_h], ...
    'Renderer', 'painters', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none');

t = tiledlayout(paper_fig, 2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for k = 1:numel(panels)
    ax_new = nexttile(t, k);
    hold(ax_new, 'on');

    copyobj(allchild(src(k).ax), ax_new);
    delete(findall(ax_new, 'Type', 'text', 'String', 'hainan'));
    restore_projected_sphere_stack(ax_new);
    replace_density_scatter_with_surface(ax_new);
    add_inner_altitude_gap(ax_new, 1000, 1100);
    style_hainan_marker(ax_new);

    set(ax_new, ...
        'FontName', fig_font, ...
        'FontSize', tick_font_pt, ...
        'LineWidth', 0.75, ...
        'TickDir', 'out', ...
        'Box', 'on', ...
        'Layer', 'top', ...
        'CLim', global_clim);
    colormap(ax_new, turbo);
    axis(ax_new, 'equal');
    xlim(ax_new, src(k).xlim);
    ylim(ax_new, src(k).ylim);
    grid(ax_new, 'off');

    if k > 2
        xlabel(ax_new, 'Projected x distance [km]', 'FontName', fig_font, 'FontSize', label_font_pt);
    else
        xlabel(ax_new, '');
        ax_new.XTickLabel = [];
    end
    if mod(k, 2) == 1
        ylabel(ax_new, 'Projected y distance [km]', 'FontName', fig_font, 'FontSize', label_font_pt);
    else
        ylabel(ax_new, '');
        ax_new.YTickLabel = [];
    end
    title(ax_new, {panels(k).title1, panels(k).title2}, ...
        'FontName', fig_font, 'FontSize', title_font_pt, 'FontWeight', 'normal');

    panel_label = text(ax_new, 0.015, 0.985, panels(k).label, ...
        'Units', 'normalized', ...
        'FontName', fig_font, ...
        'FontSize', panel_font_pt, ...
        'FontWeight', 'normal', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'Color', 'k', ...
        'BackgroundColor', 'none', ...
        'Clipping', 'off', ...
        'Margin', 0.5);
    uistack(panel_label, 'bottom');
    hold(ax_new, 'off');
end

cb = colorbar;
cb.Layout.Tile = 'south';
cb.Label.String = '{\it n}_{\rm He} [cm^{-3}]';
cb.Label.Interpreter = 'tex';
cb.Label.FontName = fig_font;
cb.Label.FontSize = label_font_pt;
cb.FontName = fig_font;
cb.FontSize = tick_font_pt;
cb.LineWidth = 0.75;

for k = 1:numel(src)
    close(src(k).fig);
end

out_png = fullfile(out_root, 'sphere_4panel_20250320_20250620_20250920_20251221_paper_singlecol.png');
print(paper_fig, out_png, '-dpng', '-r600');
fprintf('Saved figure: %s\n', out_png);

out_pdf = fullfile(out_root, 'sphere_4panel_20250320_20250620_20250920_20251221_paper_singlecol.pdf');
exportgraphics(paper_fig, out_pdf, 'ContentType', 'vector');
fprintf('Saved figure: %s\n', out_pdf);

function ax = find_main_axes(fig)
axes_list = findall(fig, 'Type', 'axes');
keep = true(size(axes_list));
for i = 1:numel(axes_list)
    tag = get(axes_list(i), 'Tag');
    if strcmpi(tag, 'Colorbar') || strcmpi(tag, 'legend')
        keep(i) = false;
    end
end
axes_list = axes_list(keep);
if isempty(axes_list)
    error('No main axes found.');
end
areas = zeros(size(axes_list));
for i = 1:numel(axes_list)
    pos = get(axes_list(i), 'Position');
    areas(i) = pos(3) * pos(4);
end
[~, idx] = max(areas);
ax = axes_list(idx);
end

function restore_projected_sphere_stack(ax)
patches = findall(ax, 'Type', 'patch');
for i = 1:numel(patches)
    fc = get(patches(i), 'FaceColor');
    if isnumeric(fc) && numel(fc) == 3 && norm(fc - [0.85, 0.88, 0.92]) < 0.08
        uistack(patches(i), 'bottom');
    end
end

red_marks = findall(ax, 'Type', 'line', 'MarkerFaceColor', 'r');
if ~isempty(red_marks)
    uistack(red_marks, 'top');
end
end

function replace_density_scatter_with_surface(ax)
scatters = findall(ax, 'Type', 'scatter');
for i = 1:numel(scatters)
    x = get(scatters(i), 'XData');
    y = get(scatters(i), 'YData');
    c = get(scatters(i), 'CData');
    x = x(:);
    y = y(:);
    c = c(:);
    keep = isfinite(x) & isfinite(y) & isfinite(c);
    x = x(keep);
    y = y(keep);
    c = c(keep);
    if numel(c) < 10
        continue;
    end

    theta = atan2(y, x);
    theta(theta < 0) = theta(theta < 0) + 2*pi;
    r = sqrt(x.^2 + y.^2);
    r_key = round(r, 6);
    theta_key = round(theta, 6);
    r_vals = unique(r_key);
    theta_vals = unique(theta_key);
    if numel(r_vals) < 2 || numel(theta_vals) < 2
        continue;
    end

    c_grid = nan(numel(theta_vals), numel(r_vals));
    [~, r_idx] = ismember(r_key, r_vals);
    [~, theta_idx] = ismember(theta_key, theta_vals);
    lin_idx = sub2ind(size(c_grid), theta_idx, r_idx);
    c_grid(lin_idx) = c;
    [theta_vals, c_grid] = close_periodic_theta(theta_vals, c_grid);

    [theta_grid, r_grid] = ndgrid(theta_vals, r_vals);
    x_grid = r_grid .* cos(theta_grid);
    y_grid = r_grid .* sin(theta_grid);
    surf(ax, x_grid, y_grid, zeros(size(x_grid)), c_grid, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'interp', ...
        'HandleVisibility', 'off');
    delete(scatters(i));
end
end

function [theta_vals_closed, c_grid_closed] = close_periodic_theta(theta_vals, c_grid)
theta_vals_closed = theta_vals;
c_grid_closed = c_grid;
if isempty(theta_vals) || isempty(c_grid)
    return;
end

theta_step = median(diff(theta_vals));
if ~isfinite(theta_step) || theta_step <= 0
    return;
end

if theta_vals(1) > 0.5 * theta_step
    theta_vals_closed = [0; theta_vals_closed];
    c_grid_closed = [c_grid_closed(end, :); c_grid_closed];
end

if (2*pi - theta_vals_closed(end)) > 0.5 * theta_step
    theta_vals_closed = [theta_vals_closed; 2*pi];
    c_grid_closed = [c_grid_closed; c_grid_closed(1, :)];
else
    theta_vals_closed(end) = 2*pi;
    c_grid_closed(end, :) = c_grid_closed(1, :);
end
end

function style_hainan_marker(ax)
red_marks = findall(ax, 'Type', 'line', 'MarkerFaceColor', 'r');
for i = 1:numel(red_marks)
    set(red_marks(i), ...
        'MarkerSize', 3.0, ...
        'MarkerEdgeColor', 'r', ...
        'LineWidth', 0.4);
    uistack(red_marks(i), 'top');
end
end

function add_inner_altitude_gap(ax, earth_radius, data_inner_radius)
theta = linspace(0, 2*pi, 360);
x_outer = data_inner_radius * cos(theta);
y_outer = data_inner_radius * sin(theta);
x_inner = earth_radius * cos(fliplr(theta));
y_inner = earth_radius * sin(fliplr(theta));
patch(ax, [x_outer, x_inner], [y_outer, y_inner], 'w', ...
    'EdgeColor', 'none', ...
    'FaceColor', 'w', ...
    'HandleVisibility', 'off');
end
