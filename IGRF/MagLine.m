% 输入经纬度及高度，计算得到地磁共轭点的经纬度

close all
clear

rlat=19.5;
rlong=109;

altitude_km = 700;        %本地高度(km)
rcon = 0;               %地磁共轭处输出高度(km)

%% 常量
Re = 6371.2; % 地球半径(km)
i = 1;
%% 处理输入坐标
% 确定纬度方向
if rlat >= 0
    lat_dir = 'N';
    lat_deg = rlat;
else
    lat_dir = 'S';
    lat_deg = -rlat;
end

% 确定经度方向
if rlong >= 0
    lon_dir = 'E';
    lon_deg = rlong;
else
    lon_dir = 'W';
    lon_deg = -rlong;
end

% 转换为弧度
theta0=(90-rlat)*pi/180;
lambda0=rlong*pi/180;
r0 = Re + altitude_km;                  % 初始半径(km)

%% 磁力线追踪
dB = 0.1; % 步长
r1 = r0;
theta1 = theta0;
lambda1 = lambda0;

while r1 <= Re + altitude_km && r1 >= Re + rcon
    % 获取当前点的磁场分量
    [B_r, B_theta, B_lambda] = IGRF_vs3(r1, theta1, lambda1);
    
    % 转换为直角坐标系分量
    B_x = sin(theta1)*cos(lambda1)*B_r + cos(theta1)*cos(lambda1)*B_theta - sin(lambda1)*B_lambda;
    B_y = sin(theta1)*sin(lambda1)*B_r + cos(theta1)*sin(lambda1)*B_theta + cos(lambda1)*B_lambda;
    B_z = cos(theta1)*B_r - sin(theta1)*B_theta;
    
    B = sqrt(B_x^2 + B_y^2 + B_z^2); % 磁场强度
    
    
    % 计算下一步位置
    [x1, y1, z1] = rtp_xyz(r1, theta1, lambda1);
    x2 = x1 - dB*B_x/B;
    y2 = y1 - dB*B_y/B;
    z2 = z1 - dB*B_z/B;

    line(i,1:3)=[x1,y1,z1];
    % 转换为球坐标
    [r1, theta1, lambda1] = xyz_rtp(x2, y2, z2);
    
    i = i+1;
end

% %% 精确计算交点
% a = B^2;
% b = 2*x1*B_x + 2*y1*B_y + 2*z1*B_z;
% c = x1^2 + y1^2 + z1^2 - Re^2;
% t = (-b + sqrt(b^2 - 4*a*c))/(2*a);
% 
% x2 = x1 + t*B_x;
% y2 = y1 + t*B_y;
% z2 = z1 + t*B_z;

% 转换为球坐标
[r1, theta1, lambda1] = xyz_rtp(x2, y2, z2);

%% 转换为地理坐标
conj_lat = 90 - theta1*180/pi; % 共轭点纬度(度)
conj_lon = lambda1*180/pi;     % 共轭点经度(度)

% 将经度调整到[-180,180]范围
if conj_lon > 180
    conj_lon = conj_lon - 360;
elseif conj_lon < -180
    conj_lon = conj_lon + 360;
end

%% 处理输出坐标
% 确定共轭点纬度方向
if conj_lat >= 0
    conj_lat_dir = 'N';
    conj_lat_deg = conj_lat;
else
    conj_lat_dir = 'S';
    conj_lat_deg = -conj_lat;
end

% 确定共轭点经度方向
if conj_lon >= 0
    conj_lon_dir = 'E';
    conj_lon_deg = conj_lon;
else
    conj_lon_dir = 'W';
    conj_lon_deg = -conj_lon;
end

%% 显示输入输出结果
% disp('=====  当前位置  =====');
% disp(['纬度: ', num2str(rlat), '° ', lat_dir]);
% disp(['经度: ', num2str(rlong), '° ', lon_dir]);
% disp(['高度: ', num2str(altitude_km), ' km']);

disp('===== 地磁共轭点 =====');
disp(['纬度: ', num2str(conj_lat_deg), '° ', conj_lat_dir]);
disp(['经度: ', num2str(conj_lon_deg), '° ', conj_lon_dir]);
disp(['高度: ', num2str(rcon), ' km']);

% conj = [conj_lat; conj_lon; rcon];
% writematrix(conj, 'C:\Users\31217\Desktop\氦激光雷达\IGRF\conj.txt');

