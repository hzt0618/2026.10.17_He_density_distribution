%% 给出给定点附近点的磁力线数据，并把磁力线数据记录下来
% 并给出该定点的磁共轭点的地心坐标

close all
clear
Re=6371.2;
altitude_km = 800;        %本地高度(km)
rcon = 0;               %地磁共轭处输出高度(km)


% SanYa
% theta=(90-(18+9/60+34/3600))*pi/180;
% lambda=(108+37/60+27/3600)*pi/180;

%海南儋州
rlat=19.5*pi/180;
rlong=109*pi/180;

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

theta=pi/2-rlat;
lambda=rlong;


% Millstone Hill
% theta=(90-42.6)*pi/180;
% lambda=(360-71.5)*pi/180;

% Arecibo
% theta=(90-18.35)*pi/180;
% lambda=(360-66.75)*pi/180;

% Kakioka
% theta=(90-36.2)*pi/180;
% lambda=140.2*pi/180;

% num=20;
% theta=linspace(30*pi/180,pi-30*pi/180,num);
% lambda=[pi/6,pi/6+pi];

% num=200;
% theta=(90-30)*pi/180;
% lambda=linspace(0,2*pi,num);


% MeCheng
% theta=(90-33)*pi/180;
% lambda=116.5*pi/180; 


% HeFei
% theta=(90-31-52/60)*pi/180;
% lambda=(117+17/60)*pi/180; 


% Huancayo
% theta=(90+12)*pi/180;
% lambda=286.7*pi/180;


% maglat corresponding to geolat30
% load maglat_corresponding_to_geolat30.mat
% theta=pi/2-geolat*pi/180;
% lambda=0:pi/6:2*pi-pi/6;

jk=1;j=1;k=1;
% for j=1:length(theta)
%     for k=1:length(lambda)
        
%         if k~=j
%             continue;
%         end
        
        %% 给出初始地点的地心坐标
        r0=Re + altitude_km;
        theta0=theta(j);
        lambda0=lambda(k);  
        
        %% 给出每一小线段延伸的长度
        dB=0.2;
        i=1;
        
        %% 把初始地点的地心坐标赋值给地点1的地心坐标
        r1=r0;
        theta1=theta0;
        lambda1=lambda0;
        
        %% 防止开始时磁场方向朝向里面
        flag=1;
        
        %% 不断循环给出整条磁力线，直到磁力线再次达到地球表面
        while r1>=Re + rcon
            %% 根据地点1的地心坐标求出球心坐标系中的磁场分量
            [B_r,B_theta,B_lambda]=IGRF_vs2(r1,theta1,lambda1);
            
            %% 把球心坐标系中的磁场分量转换为直角坐标系中的磁场分量
            B_x=sin(theta1)*cos(lambda1)*B_r+cos(theta1)*cos(lambda1)*B_theta-...
                sin(lambda1)*B_lambda;
            B_y=sin(theta1)*sin(lambda1)*B_r+cos(theta1)*sin(lambda1)*B_theta+...
                cos(lambda1)*B_lambda;
            B_z=cos(theta1)*B_r-sin(theta1)*B_theta;
            
            %% 求出地点1的总磁场强度
            B=sqrt(B_x^2+B_y^2+B_z^2);
            
            %% 把地点1的地心坐标转换为直角坐标
            [x1,y1,z1]=rtp_xyz(r1,theta1,lambda1);
            
            %% 记录地点1的直角坐标
            line{jk}(i,1:3)=[x1,y1,z1];
            
            %% 如果从给定地点开始磁场的方向朝向地球里面，进行反向
            if i==1 && x1*B_x+y1*B_y+z1*B_z<0
                flag=-1;
            end
            
            %% 根据地点1的直角坐标求出地点2的直角坐标
            x2=x1+dB*B_x/B*flag;
            y2=y1+dB*B_y/B*flag;
            z2=z1+dB*B_z/B*flag;
            
            %% 把地点2的直角坐标转换成地心坐标并赋值给地点1的地心坐标
            [r1,theta1,lambda1]=xyz_rtp(x2,y2,z2);
            
            i=i+1;
            
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

        disp('===== 地磁共轭点 =====');
        disp(['纬度: ', num2str(conj_lat_deg), '° ', conj_lat_dir]);
        disp(['经度: ', num2str(conj_lon_deg), '° ', conj_lon_dir]);
        disp(['高度: ', num2str(rcon), ' km']);
        
        conj = [conj_lat; conj_lon; rcon];
        writematrix(conj, 'C:\Users\31217\Desktop\氦激光雷达\IGRF\conj.txt');
    
        %% 求出地点1的磁共轭点的直角坐标
        t=(-2*(x1*B_x+y1*B_y+z1*B_z)-flag*2*sqrt((x1*B_x+y1*B_y+z1*B_z)^2-...
            B^2*(x1^2+y1^2+z1^2-Re^2)))/(2*B^2);
        x3=x1+B_x*t;
        y3=y1+B_y*t;
        z3=z1+B_z*t;
        
        %% 把共轭点的直角坐标转换为地心坐标并记录下来
        [r3,theta3,lambda3]=xyz_rtp(x3,y3,z3);
        line{jk}(i,1:3)=[x3,y3,z3];
        conjugate{jk}(i,1:3)=[r3,theta3,lambda3];
        
        %% 输出共轭点的经纬度
        disp([90-theta3*180/pi,lambda3*180/pi]);
        
        jk=jk+1;
%     end
% end


save line_glat_30_vs2.mat line

