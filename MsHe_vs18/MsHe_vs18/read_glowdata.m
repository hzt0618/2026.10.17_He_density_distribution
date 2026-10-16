% clear;
% close all;
% inst=1;
function [heion,elecden,electem,jmax,sza,idate,ut,glat,glong, ...
    f107a,f107,f107p,ap,ef,ec,ener,del,eflx]=read_glowdata(z0,filename)
fileID=fopen(filename,'r');
data4=textscan(fileID,repmat('%f',[1,14]),1);
jmax=data4{1,1};
nbins=data4{1,2};
data5=textscan(fileID,['%f',repmat('%f',[1,nbins])],jmax);
data6=textscan(fileID,['%f',repmat('%f',[1,nbins])],jmax);
data7=textscan(fileID,repmat('%f',[1,nbins]),1);
data8=textscan(fileID,repmat('%f',[1,nbins]),1,'ReturnOnError',1);

data1=textscan(fileID,repmat('%f',[1,jmax]),1);
data2=textscan(fileID,repmat('%f',[1,jmax]),1);
data3=textscan(fileID,repmat('%f',[1,jmax]),1);
fclose(fileID);

heion=cell2mat(data1);
elecden=cell2mat(data2);
electem=cell2mat(data3);

sza=data4{1,3};
idate=data4{1,4};
ut=data4{1,5};
glat=data4{1,6};
glong=data4{1,7};
f107a=data4{1,8};
f107=data4{1,9};
f107p=data4{1,10};
ap=data4{1,11};
ef=data4{1,12};
ec=data4{1,13};

data5mat=cell2mat(data5);
uflx=data5mat(:,2:191);        %as cm^-2 s-1 eV-1
data6mat=cell2mat(data6);
dflx=data6mat(:,2:191);        %as cm^-2 s-1 eV-1
z=data6mat(:,1);              %km
ener=cell2mat(data7);
del=cell2mat(data8);

eflx=uflx+dflx;

heion=interp1(z,heion,z0,'pchip','extrap');
elecden=interp1(z,elecden,z0,'pchip','extrap');
electem=interp1(z,electem,z0,'pchip','extrap');
eflx=interp1(z,eflx,z0,'pchip','extrap');

% energy=linspace(19.8,100,101);
% eflx1=interp1(ener,eflx',energy,'pchip','extrap');
% eflx=eflx1';
% % ii=find(ener>19.8 & ener<100);
% % ener=ener(ii);
% % eflx=eflx(:,ii);

end
