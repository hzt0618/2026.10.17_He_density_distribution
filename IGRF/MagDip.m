%% 画出地面磁偏角等高图

close all
clear

Re=6371.2;

theta=linspace(pi/6,pi*5/6,20);
lambda=linspace(0,2*pi-pi/1000,40);
dip=zeros(length(theta),length(lambda));

for i=1:length(theta)
    for j=1:length(lambda)
        %% 
        r0=Re;
        theta0=theta(i);
        lambda0=lambda(j);

        [B_r,B_theta,B_lambda]=IGRF_vs2(r0,theta0,lambda0);
        if B_theta<0 
            disp(sign(B_theta));
        end
        dip(i,j)=atan(B_lambda./B_theta);
    end
end

glong=lambda*180/pi-(sign(lambda*180/pi-180)+1)*180;
glat=90-theta*180/pi;
[glong,I]=sort(glong);
for i=1:length(glat)
    dip(i,:)=dip(i,I);
end


levels=-20:5:20;
figure(2)
s=contour(glong,glat,-dip*180/pi,levels,'Linestyle','-','ShowText','on');
colorbar
xlabel('经度')
ylim([-30,30])
ylabel('纬度')
title('400km处磁场与当地垂直方向的夹角')
axis equal









