function [qe]=he_produc1_vs2(eflx,ener,He)
% cros      cm^2
% He        cm^-3
% eflx      cm^-2 s^-1 ev^-1
% qe        cm^-3 s^-1
[cros]=cross2();
fn=length(cros);
del_e=zeros(1,fn);

del_e(2:fn-1)=(100-19.82)/(fn-1);
del_e([1,fn])=(100-19.82)/(fn-1)/2;

energy=linspace(19.8,100,fn);
eflx1=interp1(ener,eflx',energy,'pchip','extrap');
eflx=eflx1';

for i=1:length(He)
    eflx(i,:)=eflx(i,:).*del_e;
end
qe=(eflx*(cros'))'.*He;

end
