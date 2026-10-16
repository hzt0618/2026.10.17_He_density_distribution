%% 磁力线可视化，在三维图中绘制磁力线和晨昏圈
% 新增 rot_angle 参数：控制晨昏线绕地球自转轴（Z轴）的旋转角度（度）
% 调整该参数可让晨昏线覆盖磁力线区域

figure(1);

% 绘制磁力线
plot3(line(:,1), line(:,2), line(:,3), 'LineWidth', 1.5);
hold on;

%% 绘制地球
[X, Y, Z] = sphere(100);
r = 6371.2;  % 地球半径（km）
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;

s = surf(X2, Y2, Z2);
s.EdgeColor = 'none';  % 去除球面网格线
axis equal
hold on

%% 绘制晨昏圈（可旋转+覆盖磁力线区域）
% 1. 核心可调参数：晨昏线绕Z轴旋转角度（度），根据磁力线位置调整
rot_angle = 13;  % 旋转参数,匹配磁力线经度（rlong=109），可±调整

% 黄赤交角：晨昏面与地球自转轴（Z轴）的夹角（23.6°），修正之前的计算错误
alpha = (90-23.6) * pi / 180;  % 交角（弧度）
rot = rot_angle * pi / 180;  % 旋转角度转弧度

% 3. 晨昏面法向量（绕Z轴旋转后，确保过地心+可转动）
% 原始法向量（沿X轴）→ 绕Z轴旋转rot_angle → 适配磁力线区域
nx_origin = sin(alpha);  % 原始X分量（黄赤交角决定）
ny_origin = 0;           % 原始Y分量
nz_origin = cos(alpha);  % 原始Z分量（黄赤交角决定）

% 法向量绕Z轴旋转公式（保持Z分量不变，X/Y分量旋转）
nx = nx_origin * cos(rot) - ny_origin * sin(rot);
ny = nx_origin * sin(rot) + ny_origin * cos(rot);
nz = nz_origin;

% 4. 生成晨昏圈上的点（满足过地心平面方程：nx*x + ny*y + nz*z = 0）
phi = linspace(0, 2*pi, 360);  % 方位角范围（0到360度）

% 从平面方程推导极角theta（确保点严格在晨昏面上）
% 平面方程：nx*sinθ*cosφ + ny*sinθ*sinφ + nz*cosθ = 0 → tanθ = -nz/(nx*cosφ + ny*sinφ)
theta = atan2(-nz, nx*cos(phi) + ny*sin(phi));

% 转换为直角坐标
x_circle = r * sin(theta) .* cos(phi);
y_circle = r * sin(theta) .* sin(phi);
z_circle = r * cos(theta);

% 5. 筛选严格位于晨昏面上的点（过滤计算误差）
mask = abs(nx*x_circle + ny*y_circle + nz*z_circle) < 1e-3;
x_circle = x_circle(mask);
y_circle = y_circle(mask);
z_circle = z_circle(mask);

% 6. 绘制晨昏圈（红色虚线，加粗便于观察）
plot3(x_circle, y_circle, z_circle, 'r--', 'LineWidth', 1.5);

% 坐标轴标签
xlabel('(km)');
ylabel('(km)');
zlabel('(km)');