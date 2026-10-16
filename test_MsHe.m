
%% 1. 设置项目根目录
project_root = pwd;

%% 2. 添加 MsHe 模块路径
mshe_dir = fullfile(project_root, 'MsHe_vs18', 'MsHe_vs18');
addpath(mshe_dir);

%% 3. 检查 MATLAB 是否能找到主函数
which MsHe -all

%% 4. 指定一个现有的 GLOW 结果文件
glow_file = fullfile( ...
    project_root, ...
    'glow model', 'glow model', ...
    'output', ...
    'hainan_2025_sth_05', ...
    'Lat19.50_Lon109.00_20250101', ...
    '001', ...
    'Lat19.50_Lon109.00_20250101_vs001.txt');

%% 5. 检查输入文件是否存在
assert(isfile(glow_file), ...
    '找不到 GLOW 文件，请检查项目根目录和文件路径。');

%% 6. 运行亚稳态氦密度计算
r = MsHe(glow_file, struct('make_plot', false));

%% 7. 提取高度和 He* 密度
z_km = r.grid.z_km;
n_he_cm3 = r.profile.n_He_cm3;

%% 8. 绘制密度垂直廓线
figure;
plot(n_he_cm3, z_km, 'b-', 'LineWidth', 1.8);
xlabel('He(2^3S) density (cm^{-3})');
ylabel('Altitude (km)');
title('Metastable Helium Density Profile');
grid on;

%% 1. 查看结果结构体包含哪些字段
fieldnames(r)

% 2. 查看高度数组和密度数组的大小
size(r.grid.z_km)
size(r.profile.n_He_cm3)

% 3. 找出密度峰值
[n_max, idx] = max(r.profile.n_He_cm3);

% 4. 查看峰值及其对应高度
fprintf('最大密度 = %.4f cm^-3\n', n_max);
fprintf('对应高度 = %.1f km\n', r.grid.z_km(idx));

% 5. 抽查几个高度的密度
z_check = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000];

for z = z_check
    [~, idx] = min(abs(r.grid.z_km - z));
    fprintf('%4.0f km: %.5f cm^-3\n', ...
        r.grid.z_km(idx), r.profile.n_He_cm3(idx));
end