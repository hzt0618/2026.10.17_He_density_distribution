function [Hn,a,gstar,gam]=he_Hn1vs2(z,nz,Tem,sza,glat)
% 在vs1的基础上，vs2计算g^*的时候考虑了阴影高度的影响
%     d he_meta/dz = - He_meta/ H  - phi/ Dif
%     d phi/dz     = Qe - He_meta*(k1 + J2 + k3[O1] + k4[O2] + k5[N2]

% Scale Height  Bishop 1999 JGR
k_const = 1.38064852*10.^(-23);    %m2 kg s-2 K-1
me=9.10938215*10^(-31);   %kg电子质量
mp=1.672621637*10^(-27);   %kg质子(Proton)质量
mn=1.674927211*10^(-27);   %kg中子(neutron)质量
Re=6.371393*10^(6);    %m地球半径

rlat=glat*pi/180;
g  = 978.0318*(1+0.0053024*sin(2*rlat)-0.0000059*sin(2*rlat)^2)...
    *10^(-2)*(Re/1000)^2./(Re/1000+z).^2;    %m/s^2
a=1.68;        %m s^-2

h_shadow=(Re/1000+30)/sin(pi-sza)-Re/1000;
gstar=sqrt(g.^2+a^2+2*g*a*cos(sza));
gstar(z<h_shadow)=g(z<h_shadow);
gam=acos((g+a*cos(sza))./gstar);
gam(z<h_shadow)=0;

m_he=2*(mn+mp+me);

Hp = k_const .* Tem ./(m_he .* gstar);   %m

dz=(z(2)-z(1))*1000;   %m
dTdz=zeros(1,nz-1);
for i=1:nz-1
  dTdz(i)=(Tem(i+1)-Tem(i))/dz;    %K/m
end
dTdz=interp1(z(1:nz-1)*1000+dz/2,dTdz,z*1000,'pchip','extrap');
Hn=(dTdz./Tem+1./Hp).^(-1);  %m
Hn=Hn*100;            %cm
end