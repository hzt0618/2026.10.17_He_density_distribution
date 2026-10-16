function [heion,Tn,Te,Ti]=IRILatest(z0,iy,imd,lt,data)
% filename='HaiNanHeionTwiFact.txt';
% 
% fileID=fopen(filename,'r');
% data=textscan(fileID,repmat('%f',[1,8]));
% fclose(fileID);


% year=data{1}(data{2}==imd & data{1}==iy);
% doy=data{2}(data{2}==imd & data{1}==iy);
% lt=data{3}(data{2}==imd & data{1}==iy);
z=data{4}(data{2}==imd & data{1}==iy & abs(data{3}-lt)<0.01);
Tn=data{5}(data{2}==imd & data{1}==iy & abs(data{3}-lt)<0.01);
Te=data{6}(data{2}==imd & data{1}==iy & abs(data{3}-lt)<0.01);
Ti=data{7}(data{2}==imd & data{1}==iy & abs(data{3}-lt)<0.01);
heion=data{8}(data{2}==imd & data{1}==iy & abs(data{3}-lt)<0.01);
% disp(Tn)

Tn=interp1(z,Tn,z0,'pchip','extrap');
Te=interp1(z,Te,z0,'pchip','extrap');
Ti=interp1(z,Ti,z0,'pchip','extrap');
heion=interp1(z,heion,z0,'pchip','extrap');

% energy=linspace(19.8,100,101);
% eflx1=interp1(ener,eflx',energy,'pchip','extrap');
% eflx=eflx1';
% % ii=find(ener>19.8 & ener<100);
% % ener=ener(ii);
% % eflx=eflx(:,ii);

end