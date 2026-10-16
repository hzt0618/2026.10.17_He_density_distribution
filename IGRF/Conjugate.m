%% 根据某点经纬度给出共轭位置的经纬度
function [rlat,rlong]=Conjugate(rlat,rlong,year,height)
% height   km
% ralt     rad
% rlong    rad
% year

Re=6371.2+height;

%% 给出初始地点的地心坐标
theta0=pi/2-rlat;
lambda0=rlong+(1-sign(rlong))*pi;

%% 给出每一小线段延伸的长度
dB=1;
i=1;

%% 把初始地点的地心坐标赋值给地点1的地心坐标
r0=Re;
r1=r0;
theta1=theta0;
lambda1=lambda0;

%% 防止开始时磁场方向朝向里面
flag=1;

%% 不断循环给出整条磁力线，直到磁力线再次达到地球表面
while r1>=Re
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
    
%     %% 记录地点1的直角坐标
%     line{jk}(i,1:3)=[x1,y1,z1];
    
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
%% 求出地点1的磁共轭点的直角坐标
t=(-2*(x1*B_x+y1*B_y+z1*B_z)-flag*2*sqrt((x1*B_x+y1*B_y+z1*B_z)^2-...
    B^2*(x1^2+y1^2+z1^2-Re^2)))/(2*B^2);
x3=x1+B_x*t;
y3=y1+B_y*t;
z3=z1+B_z*t;

% %% 把共轭点的直角坐标转换为地心坐标并记录下来
[r3,theta3,lambda3]=xyz_rtp(x3,y3,z3);
% line{jk}(i,1:3)=[x3,y3,z3];
% conjugate{jk}(i,1:3)=[r3,theta3,lambda3];

%% 输出共轭点的经纬度
rlat=90-theta3*180/pi;
rlong=lambda3*180/pi;
if rlong > 180
    rlong = rlong - 360;
elseif rlong < -180
    rlong = rlong + 360;
end
end



