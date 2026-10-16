% make_sphere_4panel_paper.m
% Recompose the four seasonal projected-sphere .fig files as a paper-ready
% 2-by-2 figure with one shared colorbar and common color limits.

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
    'title', {'Spring equinox, 20 Mar 2026', ...
        'Summer solstice, 20 Jun 2025', ...
        'Autumn equinox, 20 Sep 2025', ...
        'Winter solstice, 21 Dec 2025'}, ...
    'label', {'(a)', '(b)', '(c)', '(d)'} ...
);

for k = 1:numel(panels)
    if ~isfile(panels(k).fig_path)
        error('Missing source figure: %s', panels(k).fig_path);
    end
end

fig_font = 'Times New Roman';
axis_font_pt = 9;
label_font_pt = 9;
title_font_pt = 9;
panel_font_pt = 9;

% First pass: collect source axes, children, and common color range.
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

paper_fig = figure( ...
    'Color', 'w', ...
    'Units', 'inches', ...
    'Position', [0.5, 0.5, 7.25, 6.25], ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0, 0, 7.25, 6.25], ...
    'PaperSize', [7.25, 6.25], ...
    'Renderer', 'painters', ...
    'MenuBar', 'none', ...
    'ToolBar', 'none');

t = tiledlayout(paper_fig, 2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'loose');

for k = 1:numel(panels)
    ax_new = nexttile(t, k);
    hold(ax_new, 'on');

    src_children = allchild(src(k).ax);
    copyobj(src_children, ax_new);
    delete(findall(ax_new, 'Type', 'text', 'String', 'hainan'));
    restore_projected_sphere_stack(ax_new);

    set(ax_new, ...
        'FontName', fig_font, ...
        'FontSize', axis_font_pt, ...
        'LineWidth', 0.7, ...
        'TickDir', 'out', ...
        'Box', 'on', ...
        'Layer', 'top', ...
        'CLim', global_clim);
    colormap(ax_new, turbo);
    axis(ax_new, 'equal');
    xlim(ax_new, src(k).xlim);
    ylim(ax_new, src(k).ylim);
    grid(ax_new, 'off');

    xlabel(ax_new, 'Projected x distance [km]', 'FontName', fig_font, 'FontSize', label_font_pt);
    ylabel(ax_new, 'Projected y distance [km]', 'FontName', fig_font, 'FontSize', label_font_pt);
    title(ax_new, panels(k).title, 'FontName', fig_font, 'FontSize', title_font_pt, 'FontWeight', 'normal');

    text(ax_new, 0.03, 0.96, panels(k).label, ...
        'Units', 'normalized', ...
        'FontName', fig_font, ...
        'FontSize', panel_font_pt, ...
        'FontWeight', 'normal', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'Color', 'k', ...
        'BackgroundColor', 'w', ...
        'Margin', 1);
    hold(ax_new, 'off');
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'n(He(2^3S)) (cm^{-3})';
cb.Label.Interpreter = 'tex';
cb.Label.FontName = fig_font;
cb.Label.FontSize = label_font_pt;
cb.FontName = fig_font;
cb.FontSize = axis_font_pt;
cb.LineWidth = 0.7;

for k = 1:numel(src)
    close(src(k).fig);
end

out_png = fullfile(out_root, 'sphere_4panel_20260320_20250620_20250920_20251221_paper.png');
print(paper_fig, out_png, '-dpng', '-r600');
fprintf('Saved figure: %s\n', out_png);

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
