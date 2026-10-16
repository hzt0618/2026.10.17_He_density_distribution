% 输入经纬度及高度，计算得到地磁共轭点的经纬度和高度

close all
clear

rlat=-75.657694;  %纬度19.5
rlong=82.830639;  %经度109

altitude_km = 0.1;        %本地高度(km) 900
rcon = 0;               %地磁共轭处输出高度(km) 500
dB = -200;                   %磁力线追踪步长+-
Bm=[];                      %磁场最值
%% 常量
Re = 6371.2; % 地球平均半径(km)
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
r1 = r0;
theta1 = theta0;
lambda1 = lambda0;

% 计算初始磁场强度
[B_r_init, B_theta_init, B_lambda_init] = IGRF_vs3(r0, theta0, lambda0);
% 转换为直角坐标系分量（与循环内的转换逻辑一致）
B_x_init = sin(theta0)*cos(lambda0)*B_r_init + cos(theta0)*cos(lambda0)*B_theta_init - sin(lambda0)*B_lambda_init;
B_y_init = sin(theta0)*sin(lambda0)*B_r_init + cos(theta0)*sin(lambda0)*B_theta_init + cos(lambda0)*B_lambda_init;
B_z_init = cos(theta0)*B_r_init - sin(theta0)*B_theta_init;
% 计算初始磁场强度（赋值给local_B）
local_B = sqrt(B_x_init^2 + B_y_init^2 + B_z_init^2);

while r1 >= Re + rcon
    % 获取当前点的磁场分量
    [B_r, B_theta, B_lambda] = IGRF_vs3(r1, theta1, lambda1);
    
    % 转换为直角坐标系分量
    B_x = sin(theta1)*cos(lambda1)*B_r + cos(theta1)*cos(lambda1)*B_theta - sin(lambda1)*B_lambda;
    B_y = sin(theta1)*sin(lambda1)*B_r + cos(theta1)*sin(lambda1)*B_theta + cos(lambda1)*B_lambda;
    B_z = cos(theta1)*B_r - sin(theta1)*B_theta;
    
    B = sqrt(B_x^2 + B_y^2 + B_z^2); % 磁场强度
    Bm(end+1) = B;
    
    % 计算下一步位置
    [x1, y1, z1] = rtp_xyz(r1, theta1, lambda1);
    x2 = x1 + dB*B_x/B;
    y2 = y1 + dB*B_y/B;
    z2 = z1 + dB*B_z/B;

    line(i,1:3)=[x1,y1,z1];
    % 转换为球坐标
    [r1, theta1, lambda1] = xyz_rtp(x2, y2, z2);
    
    i = i+1;
end



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

% 显示输入输出结果
disp('=====  当前位置  =====');
disp(['纬度: ', num2str(rlat), '° ', lat_dir]);
disp(['经度: ', num2str(rlong), '° ', lon_dir]);
disp(['高度: ', num2str(altitude_km), ' km']);

disp('===== 地磁共轭点 =====');
disp(['纬度: ', num2str(conj_lat_deg), '° ', conj_lat_dir]);
disp(['经度: ', num2str(conj_lon_deg), '° ', conj_lon_dir]);
disp(['高度: ', num2str(rcon), ' km']);

% conjN = [conj_lat; conj_lon];
% writematrix(conjN, 'C:\Users\31217\Desktop\氦激光雷达\IGRF\conjN.txt');
% conjS = [conj_lat; conj_lon];
% writematrix(conjS, 'C:\Users\31217\Desktop\氦激光雷达\IGRF\conjS.txt');

max_B = max(Bm);
min_B = min(Bm);
disp(['磁场最大值与最小值的比值: ', num2str(max_B/min_B)]);
% disp([num2str(min_B), 'nT']);
disp([num2str(max_B), 'nT']);
disp([num2str(local_B), 'nT']);
disp(['磁场当前值与最大值的比值: ', num2str(local_B/max_B)]);
% Save key results to a txt file in the current folder
output_file = fullfile(pwd, 'B_conjugate_result.txt');
fid = fopen(output_file, 'w');
if fid == -1
    warning('Failed to create result file: %s', output_file);
else
    fprintf(fid, 'Input_Latitude_deg: %.10f\n', rlat);
    fprintf(fid, 'Input_Longitude_deg: %.10f\n', rlong);
    fprintf(fid, 'Input_Altitude_km: %.10f\n', altitude_km);
    fprintf(fid, 'Conjugate_Latitude_deg: %.10f\n', conj_lat);
    fprintf(fid, 'Conjugate_Longitude_deg: %.10f\n', conj_lon);
    fprintf(fid, 'Conjugate_Altitude_km: %.10f\n', rcon);
    fprintf(fid, 'Bmax_nT: %.10f\n', max_B);
    fprintf(fid, 'Bmin_nT: %.10f\n', min_B);
    fprintf(fid, 'LocalB_nT: %.10f\n', local_B);
    fprintf(fid, 'Bmax_to_Bmin: %.10f\n', max_B/min_B);
    fprintf(fid, 'LocalB_to_Bmax: %.10f\n', local_B/max_B);
    fclose(fid);
    disp(['Result saved to: ', output_file]);
end
