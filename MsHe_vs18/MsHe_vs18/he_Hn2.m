function [Hn,g,a,gstar]=he_Hn2(z,Tem,sza)
%     d he_meta/dz = - He_meta/ H  - phi/ Dif
%     d phi/dz     = Qe - He_meta*(k1 + J2 + k3[O1] + k4[O2] + k5[N2]

% Scale Height  Bishop 1999 JGR

%% some constant
k_const = 1.38064852*10.^(-23);    %m2 kg s-2 K-1
m_he    = 2* 1.6726231*10^(-27)+2*1.6749286*10^(-27)+2*9.10956*10^(-31);  % kg
G=6.67259*10^(-11);    %Nm^2/kg^2牛顿引力常量
M_E=5.965*10^(24);     %kg
Re=6.371393*10^(6);    %m地球半径

%% calculate g
g  = G*M_E./(Re+z*1000).^2;    %m/s^2

%% calculate g^*
a=1.68;        %m s^-2
gstar=sqrt(g.^2+a^2+2*g*a*cos(sza));

%% calculate Hp
Hp = k_const .* Tem ./(m_he .* gstar);   %m

%% calculate Hn
Hn=Hp;  %m
Hn=Hn*100;            %cm
end