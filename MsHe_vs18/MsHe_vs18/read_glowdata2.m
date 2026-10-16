% clear;
% close all;
% inst=1;
function [jmax,sza,idate,ut,glat,glong, ...
    f107a,f107,f107p,ap,ef,ec,ener,del,eflx]=read_glowdata2(z0,filename)
fileID=fopen(filename,'r');
data4=textscan(fileID,repmat('%f',[1,14]),1);
jmax=data4{1,1};    %138
nbins=data4{1,2};   %190
data5=textscan(fileID,['%f',repmat('%f',[1,nbins])],jmax);
data6=textscan(fileID,['%f',repmat('%f',[1,nbins])],jmax);
data7=textscan(fileID,repmat('%f',[1,nbins]),1);
data8=textscan(fileID,repmat('%f',[1,nbins]),1,'ReturnOnError',1);
fclose(fileID);


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

eflx=interp1(z,eflx,z0,'pchip','extrap');

% energy=linspace(19.8,100,101);
% eflx1=interp1(ener,eflx',energy,'pchip','extrap');
% eflx=eflx1';
% % ii=find(ener>19.8 & ener<100);
% % ener=ener(ii);
% % eflx=eflx(:,ii);

end
