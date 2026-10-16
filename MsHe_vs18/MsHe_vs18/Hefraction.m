% He数密度
clear
close all

z=100:1000;
glat=19.5;
glong=109;
Yr=2020;
DoY=180;
ut=24*3600;
f107a=150;
f107=150;
ap=8;

[T,DEN]=atmosnrlmsise00(z*1000,glat,glong,Yr,DoY,ut,f107a,f107,ap);
He=DEN(:,1);
TDen=sum(DEN(:,[1,2,3,4,5,7,8,9]),2);


figure(1)
plot(He./TDen*100,z,'LineWidth',2)
xlabel('He Number Density Fraction [%]')
ylabel('Height [m]')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);


figure(2)
for i=[1,2,3,4,5,7,8]
    semilogx(DEN(:,i),z,'LineWidth',2)
    hold on
end
legend('He','O','N_2','O_2','Ar','H','N')
xlim([10^(-1),10^(20)])
xlabel('Number Density [m^{-3}]')
ylabel('Height [km]')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);