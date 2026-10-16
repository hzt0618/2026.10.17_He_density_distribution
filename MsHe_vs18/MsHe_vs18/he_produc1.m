function [qe]=he_produc1(eflx,He)
% cros      cm^2
% He        cm^-3
% eflx      cm^-2 s^-1 ev^-1
% qe        cm^-3 s^-1
[cros]=cross();
fn=length(cross);
del_e=zeros(1,fn);

del_e(2:100)=(100-19.82)/(fn-1);
del_e([1,101])=(100-19.82)/(fn-1)/2;

for i=1:length(He)
    eflx(i,:)=eflx(i,:).*del_e;
end

qe=(eflx*(cros'))'.*He;
end