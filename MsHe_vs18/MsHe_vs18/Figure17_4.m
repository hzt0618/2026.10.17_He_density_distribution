%% 17_4是作为把MsHe_vs17_4的运算结果，画成图
close all
clear

ii=0;

load HuangHeStation2021_CD.mat
nl=size(Doy,1);
ndoy=size(Doy,2);

It=It1;

figure(11)
contourf(Doy(2:nl/2-ii,:),LT(2:nl/2-ii,:),It(2:nl/2-ii,:),'Levels',max(max(It))/100,...
    'LineColor','none')
hold on
contourf(Doy(nl/2+1+ii:nl,:),LT(nl/2+1+ii:nl,:),It(nl/2+1+ii:nl,:),'Levels',max(max(It))/100,...
    'LineColor','none')
xlabel('Day of year')
xlim([1,365])
ylabel('Local time, hour')
ylim([0,24])
hold on
c=colorbar;
c.Label.String = '1083 nm emisson intensity (R)';
title('The north the Yellow River station China');
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);


It=It2;

figure(12)
contourf(Doy(2:nl/2-ii,:),LT(2:nl/2-ii,:),It(2:nl/2-ii,:),'Levels',max(max(It))/100,...
    'LineColor','none')
hold on
contourf(Doy(nl/2+1+ii:nl,:),LT(nl/2+1+ii:nl,:),It(nl/2+1+ii:nl,:),'Levels',max(max(It))/100,...
    'LineColor','none')
xlabel('Day of year')
xlim([1,365])
ylabel('Local time, hour')
ylim([0,24])
hold on
c=colorbar;
c.Label.String = '1083 nm emisson intensity (R)';
title('The north the Yellow River station China');
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);


% figure(10)
% contour(Doy,LT,It)
% xlabel('Day of year')
% ylabel('Local time, hour')
% zlabel('1083 nm emisson intensity (R)')
% set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);
% 
% 
% 
% figure(12)
% mesh(Doy(2:nl/2-ii,:),LT(2:nl/2-ii,:),It(2:nl/2-ii,:))
% hold on
% mesh(Doy(nl/2+1+ii:nl,:),LT(nl/2+1+ii:nl,:),It(nl/2+1+ii:nl,:))
% xlabel('Doy')
% xlim([1,365])
% ylabel('local time')
% ylim([0,24])
% zlabel('4\piI(R)')
% set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);














%% 特定天顶角时的值

% szat=[105,110,115];
% ItDawn=zeros(ndoy,3);
% ItDusk=zeros(ndoy,3);
% for i=1:length(szat)
%     for j=1:ndoy
%         ItDawn(j,i)=interp1(sza0(1:nl/2,j),It(1:nl/2,j),szat(i)*pi/180,'pchip');
%         ItDusk(j,i)=interp1(sza0(nl/2+1:nl,j),It(nl/2+1:nl,j),szat(i)*pi/180,'pchip');
%     end
% end
% figure(20)
% plot(Doy(1,:),ItDawn(:,1),'color','r','Linestyle','-','Linewidth',1.5)
% hold on
% plot(Doy(1,:),ItDawn(:,2),'color','b','Linestyle','-','Linewidth',1.5)
% hold on
% plot(Doy(1,:),ItDawn(:,3),'color','g','Linestyle','-','Linewidth',1.5)
% hold on
% plot(Doy(1,:),ItDusk(:,1),'color','r','Linestyle','--','Linewidth',1.5)
% hold on
% plot(Doy(1,:),ItDusk(:,2),'color','b','Linestyle','--','Linewidth',1.5)
% hold on
% plot(Doy(1,:),ItDusk(:,3),'color','g','Linestyle','--','Linewidth',1.5)
% hold on
% xlabel('doy')
% xlim([1,365])
% ylabel('I [R]')
% legend('Dawn, SZA=105^{\circ}','Dawn, SZA=110^{\circ}','Dawn, SZA=115^{\circ}',...
%     'Dusk, SZA=105^{\circ}','Dusk, SZA=110^{\circ}','Dusk, SZA=115^{\circ}',...
%     'Location','best')
% set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);
% % title('HaiNanDanZhou(19.5^{\circ}N, 109^{\circ}E)2021')
% % title('Huancayo(12^{\circ}S, 73.3^{\circ}W)1973')
% % title('HeFei(31.87^{\circ}N, 117.28^{\circ}E)2021')
% title('HaiNanDanZhou(19.5^{\circ}N, 109^{\circ}E)2021WithoutConjugate')

