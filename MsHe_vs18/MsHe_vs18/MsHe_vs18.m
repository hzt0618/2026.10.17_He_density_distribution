%vs18
% 相对于vs15，vs17修正了对碰撞截面的计算。修正过后，更方便对横截面进行修改。
% 相对于vs17，vs17_2对某地全年的亚稳态氦密度进行模拟，
% 20220508,vs18把改变氦离子的来源为最新版的IRI模型
% 并得出被动观测强度随全年的变化，2022.02.25
% HaiNan2021_2相比于HaiNan2021相比，没有考虑qr
% HaiNan2021_vs3相比于HaiNan2021_vs2，把heion乘以1/5
% HaiNan2021_vs4相比于HaiNan2021_vs2，把heion乘以1/10

clear;
close all

%% 设置一些标记参数
whpp=0;  %  如果whpp=0，则需要考虑光压；如果whpp=1，则无需考虑光压；photon press
whre=0;  %  如果whre=0，则考虑共振发光，是按照被动观测的公式计算的；
         %  如果whre=1，则无需考虑共振，是为激光雷达观测准备的；resonance
%% variable indefine
usehu=1000;
usehd=100;
interval=1;
nz=(usehu-usehd)/interval+1;
z=linspace(usehd,usehu,nz);  %km
exz=1500;      % extension height, km


%% 设置一些随时间变化的参数
list=1:70;
nl=2*length(list);
doy=1:5:365;
ndoy=length(doy);

sza0=zeros(nl,ndoy);
It1=zeros(nl,ndoy);
It2=zeros(nl,ndoy);
LT=zeros(nl,ndoy);
Doy=zeros(nl,ndoy);
% Qe=zeros(nl,365);
%% 定义颜色
co=parula(256);
DawnDusk=["Dawn","Dusk"];

%% heion
filename='HaiNanHeionTwiFact.txt';
fileID=fopen(filename,'r');
data=textscan(fileID,repmat('%f',[1,8]));
fclose(fileID);
%% 正演
for doyi=doy
    disp(doyi)
    for DD=1:2
for inst=list
    color=interp1(1:256,co,inst*255/144+1,'pchip','extrap');
    %% read data
    filename=['./HaiNanNoCon2021_2/',num2str(doyi),DawnDusk{DD},'/output_vs',num2str(inst),'.txt'];
    [heion,elecden,electem,jmax,sza,idate,ut,glat,glong, ...
        f107a,f107,f107p,ap,ef,ec,ener,del,eflx]=read_glowdata2(z,filename);
    
    
    [heion,Tn,Te,Ti]=IRILatest(z,2020,doyi,UtToLt(glong,ut/3600),data);
    
    
    
    if ef>30
        break;
    end
    sza0(find(list==inst)+(DD-1)*nl/2,doy==doyi)=sza;
    LT(find(list==inst)+(DD-1)*nl/2,doy==doyi)=UtToLt(glong,ut/3600);
    %% get DoY and Yr
    [DoY,Yr]=date2doy(idate);
    Doy(find(list==inst)+(DD-1)*nl/2,doy==doyi)=DoY;
    %% get T and Density following altitude from atmosbrlmsise00
    % 获得大气参数
    [Tem,He,O1,N2,O2,den]=gain_atmdentem(z,glat,glong,Yr,DoY,ut,f107a,f107,ap);
    
    % 把大气温度，氧原子、氧气、氮气、电子温度外推到exz的高度，并生成ze
    [ze,Tem,O1,N2,O2,electem]=ExtensionHeight1(z,exz,Tem,O1,N2,O2,electem,glat);
    nze=length(ze);

    %% (1) the meta He production rate due to photonelectrons,
    % 计算电子碰撞基态氦产生亚稳态氦的速率轮廓      cm^(-3)s^(-1), that is, alpha * He
    qe=he_produc1_vs2(eflx,ener,He);
    
    % 把速率轮廓外推到ze高度
    qe=ExtensionHeight2(z,ze,qe,Tem,glat); 
    
    %% (2) the meta He production rate due to recombination, cm^(-3)s^(-1)
    % 计算氦离子与电子重组产生亚稳态氦的速率轮廓
    qr=he_produc2(heion,elecden,electem);
    
    % 将速率轮廓外推到ze高度
    qr=ExtensionHeight3(z,ze,qr,Tem,glat);
    
%     qe=qe+qr;
    %% ;; (3) Loss rate,  s^(-1)
    Loss=he_loss(O1,O2,N2,Tem);
    
    %% (4)  diffusion coefficiant and scale height
    % 把基态氦密度外推到ze的高度
    He=ExtensionHeight4(z,ze,He,Tem,glat);
    
    % 计算亚稳态氦的扩散系数
    Dif=he_dif(He,O1,Tem);
    
    % calculate metastable He scale rhight
    [Hn,a,gstar,gam]=he_Hn1vs3(ze,nz,Tem,sza,glat,whpp);   %cm
    %% (5) 解扩散方程
    
    % 考虑光压的影响，改变方程中的高度
    zstar=Zstar(ze,sza,gam,whpp);
    
    % 根据上下边界条件计算亚稳态氦密度轮廓
    result=he_denvs5(zstar,Dif,Hn,Loss,qe);
    
    % 计算解的相对误差
    per=SolveError(ze,qe,Loss,result);

    %% 求1083nm处不同高度的太阳光强度
    [den_e,ze2]=ExtensionHeight5(ze,interval,sza,Tem,result);
    
%     It(find(list==inst)+(DD-1)*nl/2,doy==doyi)=GlowIntensity(ze2,den_e,sza,whre);  % Rayleigh
%     
%     if DD==2
%     figure(11)
%     semilogx(heion,z)
%     hold on
%     end
    [It1(find(list==inst)+(DD-1)*nl/2,doy==doyi),...
        It2(find(list==inst)+(DD-1)*nl/2,doy==doyi)]...
        =GlowIntensity2(ze2,den_e,sza);  % Rayleigh
end
    end
end



% load Huancayo1973.mat
nl=size(Doy,1);
ndoy=size(Doy,2);

figure(10)
mesh(Doy(1:nl/2,:),LT(1:nl/2,:),It1(1:nl/2,:))
xlabel('Doy')
ylabel('local time')
zlabel('4\piI(R)')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);
figure(12)
mesh(Doy(nl/2+1:nl,:),LT(nl/2+1:nl,:),It1(nl/2+1:nl,:))
xlabel('Doy')
ylabel('local time')
zlabel('4\piI(R)')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);


figure(15)
mesh(Doy(1:nl/2,:),sza0(1:nl/2,:)*180/pi,It1(1:nl/2,:))
xlabel('Doy')
xlim([1,365])
ylabel('SZA')
ylim([90,140])
zlabel('4\piI(R)')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);
figure(16)
mesh(Doy(nl/2+1:nl,:),sza0(nl/2+1:nl,:)*180/pi,It1(nl/2+1:nl,:))
xlabel('Doy')
ylabel('SZA')
zlabel('4\piI(R)')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);


szat=[105,110,115];
ItDawn=zeros(ndoy,3);
ItDusk=zeros(ndoy,3);
for i=1:length(szat)
    for j=1:ndoy
        ItDawn(j,i)=interp1(sza0(1:nl/2,j),It1(1:nl/2,j),szat(i)*pi/180,'pchip');
        ItDusk(j,i)=interp1(sza0(nl/2+1:nl,j),It1(nl/2+1:nl,j),szat(i)*pi/180,'pchip');
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
legend('Dawn, SZA=105^{\circ}','Dawn, SZA=110^{\circ}','Dawn, SZA=115^{\circ}',...
    'Dusk, SZA=105^{\circ}','Dusk, SZA=110^{\circ}','Dusk, SZA=115^{\circ}',...
    'Location','best')
set(gca,'FontName','Times New Roman','FontSize',12,'FontWeight','Bold','linewidth',1);
title('HaiNanDanZhou(19.5^{\circ}N, 109^{\circ}E)2021WithoutConjugate')
% title('Huancayo(12^{\circ}S, 73.3^{\circ}W)1973')

save HaiNanNoCon2021_CD_iriLatest.mat Doy LT sza0 It1 It2
