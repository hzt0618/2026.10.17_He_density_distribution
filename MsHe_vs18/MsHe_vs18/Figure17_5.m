%% 17_5是作为把MsHe_vs17_5的运算结果，画成图
% HaiNanNoConEF2021_2.mat代表30eV电子450km处的电子通量
% HaiNanNoConEF2021_3.mat代表30eV电子所有高度处电子通量的和
% HaiNanNoConEF2021_4.mat代表20-50eV电子所有高度处电子通量的和

close all
clear

ii=2;

load HaiNanNoConEF2021_4.mat
nl=size(Doy,1);
ndoy=size(Doy,2);

figure(10)
mesh(Doy(ii:nl/2,:),LT(ii:nl/2,:),EF(ii:nl/2,:))
xlabel('Doy')
ylabel('local time')
zlabel('30eV electran flux [cm^{-2} s^{-1} eV^{-1}]')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);
figure(12)
mesh(Doy(nl/2+ii:nl,:),LT(nl/2+ii:nl,:),EF(nl/2+ii:nl,:))
xlabel('Doy')
ylabel('local time')
zlabel('30eV electran flux [cm^{-2} s^{-1} eV^{-1}]')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);


figure(15)
mesh(Doy(ii:nl/2,:),sza0(ii:nl/2,:)*180/pi,EF(ii:nl/2,:))
xlabel('Doy')
xlim([1,365])
ylabel('SZA')
ylim([90,140])
zlabel('30eV electran flux [cm^{-2} s^{-1} eV^{-1}]')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);
figure(16)
mesh(Doy(nl/2+ii:nl,:),sza0(nl/2+ii:nl,:)*180/pi,EF(nl/2+ii:nl,:))
xlabel('Doy')
ylabel('SZA')
zlabel('30eV electran flux [cm^{-2} s^{-1} eV^{-1}]')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);


szat=[105,110,115];
ItDawn=zeros(ndoy,3);
ItDusk=zeros(ndoy,3);
for i=1:length(szat)
    for j=1:ndoy
        ItDawn(j,i)=interp1(sza0(1:nl/2,j),EF(1:nl/2,j),szat(i)*pi/180,'pchip');
        ItDusk(j,i)=interp1(sza0(nl/2+1:nl,j),EF(nl/2+1:nl,j),szat(i)*pi/180,'pchip');
    end
end
figure(20)
plot(Doy(1,:),ItDawn(:,1),'color','r','Linestyle','-','Linewidth',1.5)
hold on
plot(Doy(1,:),ItDawn(:,2),'color','b','Linestyle','-','Linewidth',1.5)
hold on
plot(Doy(1,:),ItDawn(:,3),'color','g','Linestyle','-','Linewidth',1.5)
hold on
plot(Doy(1,:),ItDusk(:,1),'color','r','Linestyle','--','Linewidth',1.5)
hold on
plot(Doy(1,:),ItDusk(:,2),'color','b','Linestyle','--','Linewidth',1.5)
hold on
plot(Doy(1,:),ItDusk(:,3),'color','g','Linestyle','--','Linewidth',1.5)
hold on
xlabel('doy')
xlim([1,365])
ylabel('30eV electran flux [cm^{-2} s^{-1} eV^{-1}]')
legend('Dawn, SZA=105^{\circ}','Dawn, SZA=110^{\circ}','Dawn, SZA=115^{\circ}',...
    'Dusk, SZA=105^{\circ}','Dusk, SZA=110^{\circ}','Dusk, SZA=115^{\circ}',...
    'Location','best')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);
% title('HaiNanDanZhou(19.5^{\circ}N, 109^{\circ}E)2021')
% title('Huancayo(12^{\circ}S, 73.3^{\circ}W)1973')
% title('HeFei(31.87^{\circ}N, 117.28^{\circ}E)2021')
title('HaiNanDanZhou(19.5^{\circ}N, 109^{\circ}E)2021WithoutConjugate')




