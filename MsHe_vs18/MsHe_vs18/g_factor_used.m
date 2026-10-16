clear;
close all;

%% 设置一些可能用到的物理常数
h=6.62607015*10^(-34);   %Js普朗克常量
c0=2.99792458*10^(8);    %m/s真空光速
kb=1.380649*10^(-23);     %J/K玻尔兹曼常数（Boltzmann constant）
mp=1.6726231*10^(-27);   %kg质子(Proton)质量
mn=1.6749286*10^(-27);   %kg中子(neutron)质量
Re=6.371393*10^(6);    %m地球半径

%% 求1083nm处不同高度的太阳光强度

% 先做一个简单的假设
I_sun1082=0.6182;    % W m-2 nm-1
I_sun1083=0.6168;    % W m-2 nm-1
I_sun1084=0.6154;    % W m-2 nm-1


%% 求高度z处的共振散射光强
% 首先，对太阳光谱进行一个简单的线性插值
lam_num=1001;
lambda=linspace(1082.8,1083.2,lam_num);    % nm
Isun=interp1([1082,1083,1084],[I_sun1082,I_sun1083,I_sun1084],...
    lambda,'linear');    % W m-2 nm-1
dlambda=(max(lambda)-min(lambda))/(lam_num-1);  % nm
lambda=lambda*10^(-9);    % m


% 假设一个温度
T=1000;    % K

% 爱因斯坦A系数
Aki=1.0216*10^7;   % s-1


% 对J=0的共振，计算原子共振散射截面
sigma0=0;
lambda0=1082.909114*10^(-9);     % m
sigmaD=c0/lambda0*sqrt(kb*T/(2*mp+2*mn)/c0^2);   %s-1
for i=2:lam_num-1
    sigma0=sigma0+lambda0^2*...
        exp(-(c0/lambda(i)-c0/lambda0)^2/2/sigmaD^2)*...
        Isun(i)/h/c0*lambda(i)*dlambda;
end
sigma0=sigma0*Aki/8/pi/sqrt(2*pi)/sigmaD;



% 对J=1的共振，计算原子共振散射截面
sigma1=0;
lambda1=1083.025010*10^(-9);     % m
sigmaD=c0/lambda1*sqrt(kb*T/(2*mp+2*mn)/c0^2);   %s-1
for i=2:lam_num-1
    sigma1=sigma1+lambda1^2*...
        exp(-(c0/lambda(i)-c0/lambda1)^2/2/sigmaD^2)*...
        Isun(i)/h/c0*lambda(i)*dlambda;
end
sigma1=sigma1*Aki/8/pi/sqrt(2*pi)/sigmaD;


% 对J=2的共振，计算原子共振散射截面
sigma2=0;
lambda2=1083.034977*10^(-9);     % m
sigmaD=c0/lambda2*sqrt(kb*T/(2*mp+2*mn)/c0^2);   % s-1
for i=2:lam_num-1
    sigma2=sigma2+lambda2^2*...
        exp(-(c0/lambda(i)-c0/lambda2)^2/2/sigmaD^2)*...
        Isun(i)/h/c0*lambda(i)*dlambda;
end
sigma2=sigma2*Aki/8/pi/sqrt(2*pi)/sigmaD;



% 求总散射截面
sigmat=sigma0*3/33+sigma1*10/33+sigma2*20/33;




Isun0=interp1([1082,1083,1084],[I_sun1082,I_sun1083,I_sun1084],...
    lambda0*10^9,'linear');   % W m-2 nm-1
Isun1=interp1([1082,1083,1084],[I_sun1082,I_sun1083,I_sun1084],...
    lambda1*10^9,'linear');    % W m-2 nm-1
Isun2=interp1([1082,1083,1084],[I_sun1082,I_sun1083,I_sun1084],...
    lambda2*10^9,'linear');    % W m-2 nm-1
g0=1.33*10^(-21)*Isun0*10^(-4)/h/c0*lambda0*(lambda0*10^6)^4*Aki*0.95;
g1=1.33*10^(-21)*Isun1*10^(-4)/h/c0*lambda1*(lambda1*10^6)^4*Aki*0.95;
g2=1.33*10^(-21)*Isun2*10^(-4)/h/c0*lambda2*(lambda2*10^6)^4*Aki*0.95;
g=g1*10/33+g2*20/33;
g=3*g;



