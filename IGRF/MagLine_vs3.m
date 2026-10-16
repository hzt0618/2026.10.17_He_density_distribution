%% 本程序画出400km处，磁场与当地垂直方向的夹角随纬度的变化


close all
clear
Re=6371.2;



dB=1;


theta=linspace(0,pi,100);
lambda=linspace(0,2*pi,200);
za=zeros(length(theta),length(lambda));

for i=1:length(theta)
    for j=1:length(lambda)
        %% 
        r0=Re;
        theta0=theta(i);
        lambda0=lambda(j);
        
        
        r1=r0+400;
        theta1=theta0;
        lambda1=lambda0;
        
        [B_r,B_theta,B_lambda]=IGRF_vs2(r1,theta1,lambda1);
        
        
        B_x=sin(theta1)*cos(lambda1)*B_r+cos(theta1)*cos(lambda1)*B_theta-...
            sin(lambda1)*B_lambda;
        B_y=sin(theta1)*sin(lambda1)*B_r+cos(theta1)*sin(lambda1)*B_theta+...
            cos(lambda1)*B_lambda;
        B_z=cos(theta1)*B_r-sin(theta1)*B_theta;
        
        [x0,y0,z0]=rtp_xyz(r0,theta0,lambda0);
        
        za(i,j)=acos((x0*B_x+y0*B_y+z0*B_z)/...
            (sqrt(x0^2+y0^2+z0^2)*sqrt(B_x^2+B_y^2+B_z^2)));
    end
end



figure(2)
s=contourf(lambda*180/pi,...
    90-theta*180/pi,abs(za*180/pi-(sign(za*180/pi-90)+1)/2*180),20,'Linestyle','none');
colorbar
xlabel('经度')
ylim([-90,90])
ylabel('纬度')
title('400km处磁场与当地垂直方向的夹角')














