function qe=ExtensionHeight2(z,ze,qe,Tem,glat)
% some contants
kb = 1.38064852*10.^(-23);    %m2 kg s-2 K-1
me=9.10938215*10^(-31);   %kg电子质量
mp=1.672621637*10^(-27);   %kg质子(Proton)质量
mn=1.674927211*10^(-27);   %kg中子(neutron)质量
Re=6.371393*10^(6);    %m地球半径


nze=length(ze);
nz=length(z);

% calculate g
rlat=glat*pi/180;
g  = 978.0318*(1+0.0053024*sin(2*rlat)-0.0000059*sin(2*rlat)^2)...
    *10^(-2)*(Re/1000)^2./(Re/1000+ze).^2;    %m/s^2

% calclate extension heigh density
for i=nz+1:nze
    qe(i)=qe(i-1)*exp(-(ze(i)-ze(i-1)).*1000.*g(i)*2*...
        (mp+mn+me)./(kb*Tem(i)));
end

end