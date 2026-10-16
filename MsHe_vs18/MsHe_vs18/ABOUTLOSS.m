clear
close all

T=linspace(500,3000,1000);
kN2=1.9*10^(-9)*exp(-(1.04+685./T));
kO2=1.9*10^(-9)*exp(-(0.63+418./T));
kO=1.0*10^(-9)*(T/1500).^(1/6);

figure(1)
plot(T,kN2,T,kO2,T,kO)
legend('kN2','kO2','kO')


