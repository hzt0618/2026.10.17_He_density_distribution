%计算地磁场磁偏角随高度的值

%% 给出给定点的磁场参数

close all
clear
Re=6371.2;


nr=1000;

% HeFei
r1=linspace(Re,Re+1000,nr);
lambda1=(117+13/60+18.25/3600)*pi/180;
theta1=(90-31-49/60-21.30/3600)*pi/180;


% num=20;
% theta=linspace(30*pi/180,pi-30*pi/180,num);
% lambda=[pi/6,pi/6+pi];


De=zeros(1,nr);
for j=1:nr

    %% 根据地点1的地心坐标求出球心坐标系中的磁场分量
    [B_r,B_theta,B_lambda]=IGRF_vs2(r1(j),theta1,lambda1);
    
    %% 计算地磁场磁偏角随高度的值
    De(j)=atan(B_lambda/B_theta);
    
%     %% 把球心坐标系中的磁场分量转换为直角坐标系中的磁场分量
%     B_x=sin(theta1)*cos(lambda1)*B_r+cos(theta1)*cos(lambda1)*B_theta-...
%         sin(lambda1)*B_lambda;
%     B_y=sin(theta1)*sin(lambda1)*B_r+cos(theta1)*sin(lambda1)*B_theta+...
%         cos(lambda1)*B_lambda;
%     B_z=cos(theta1)*B_r-sin(theta1)*B_theta;
%     
%     %% 求出地点1的总磁场强度
%     B=sqrt(B_x^2+B_y^2+B_z^2);
%     
%     %% 把地点1的地心坐标转换为直角坐标
%     [x1,y1,z1]=rtp_xyz(r1,theta1,lambda1);

end
figure(1)
plot(r1,De*180/pi)

% save line_ziwu_vs2.mat line
