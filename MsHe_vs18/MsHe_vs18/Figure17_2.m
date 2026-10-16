% 17_2是作为把MsHe_vs17_2的运算结果，画成图
close all
clear

ii=2;

Local='HaiNan';


matname=[Local,'2021_CD_iriLatest.mat'];
load(matname)
nl=size(Doy,1);
ndoy=size(Doy,2);

It=It2;

figure(10)
mesh(Doy(ii:nl/2,:),LT(ii:nl/2,:),It(ii:nl/2,:))
hold on
mesh(Doy(nl/2+ii:nl,:),LT(nl/2+ii:nl,:),It(nl/2+ii:nl,:))
xlabel('Doy')
xlim([1,365])
ylabel('local time')
ylim([0,24])
zlabel('4\piI(R)')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);

figure(11)
contourf(Doy(ii:nl/2,:),LT(ii:nl/2,:),It(ii:nl/2,:),'Levels',max(max(It))/100,...
    'LineColor','none')
hold on
contourf(Doy(nl/2+ii:nl,:),LT(nl/2+ii:nl,:),It(nl/2+ii:nl,:),'Levels',max(max(It))/100,...
    'LineColor','none')

c=colorbar;
c.Label.String = '1083 nm emisson intensity (R)';
title(Local);
xlabel('Day of year')
xlim([1,365])
ylabel('Local time, hours')
ylim([0,24])
zlabel('4\piI(R)')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);


It=It2-It1;
figure(13)
contourf(Doy(ii:nl/2,:),LT(ii:nl/2,:),It(ii:nl/2,:),'Levels',max(max(It))/100,...
    'LineColor','none')
hold on
contourf(Doy(nl/2+ii:nl,:),LT(nl/2+ii:nl,:),It(nl/2+ii:nl,:),'Levels',max(max(It))/100,...
    'LineColor','none')
c=colorbar;
c.Label.String = '1083 nm emisson intensity (R)';
title(Local);
xlabel('Day of year')
xlim([1,365])
ylabel('Local time, hours')
ylim([0,24])
zlabel('4\piI(R)')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);


It=(It2-It1)./(It1+1e-9);
figure(12)
contourf(Doy(ii:nl/2,:),LT(ii:nl/2,:),It(ii:nl/2,:),'Levels',max(max(It))/100,...
    'LineColor','none')
hold on
contourf(Doy(nl/2+ii:nl,:),LT(nl/2+ii:nl,:),It(nl/2+ii:nl,:),'Levels',max(max(It))/100,...
    'LineColor','none')
c=colorbar;
c.Label.String = '1083 nm emisson intensity (R)';
title(Local);
xlabel('Day of year')
xlim([1,365])
ylabel('Local time, hours')
ylim([0,24])
zlabel('4\piI(R)')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);


figure(15)
mesh(Doy(ii:nl/2,:),sza0(ii:nl/2,:)*180/pi,It(ii:nl/2,:))
xlabel('Doy')
xlim([1,365])
ylabel('SZA')
ylim([90,140])
zlabel('4\piI(R)')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);
figure(16)
mesh(Doy(nl/2+ii:nl,:),sza0(nl/2+ii:nl,:)*180/pi,It(nl/2+ii:nl,:))
xlabel('Doy')
ylabel('SZA')
zlabel('4\piI(R)')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);


szat=[100,105,110];
ItDawn=zeros(ndoy,3);
ItDusk=zeros(ndoy,3);
for i=1:length(szat)
    for j=1:ndoy
        ItDawn(j,i)=interp1(sza0(1:nl/2,j),It(1:nl/2,j),szat(i)*pi/180,'pchip');
        ItDusk(j,i)=interp1(sza0(nl/2+1:nl,j),It(nl/2+1:nl,j),szat(i)*pi/180,'pchip');
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
ylabel('I [R]')
legend(['Dawn, SZA=',num2str(szat(1)),'^{\circ}'],...
    ['Dawn, SZA=',num2str(szat(2)),'^{\circ}'],...
    ['Dawn, SZA=',num2str(szat(3)),'^{\circ}'],...
    ['Dawn, SZA=',num2str(szat(1)),'^{\circ}'],...
    ['Dawn, SZA=',num2str(szat(2)),'^{\circ}'],...
    ['Dawn, SZA=',num2str(szat(3)),'^{\circ}'],...
    'Location','best')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);
title('HaiNanDanZhou(19.5^{\circ}N, 109^{\circ}E)2021')
title('Huancayo(12^{\circ}S, 73.3^{\circ}W)1973')
title('HeFei(31.87^{\circ}N, 117.28^{\circ}E)2021')
title('HaiNanDanZhou(19.5^{\circ}N, 109^{\circ}E)2021WithoutConjugate')
title('HuangHeStation(78.92^{\circ}N, 11.93^{\circ}E)2021')



