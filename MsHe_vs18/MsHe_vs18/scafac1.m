%% The program calculate scattering factor of metastable helium



function p=scafac1(theta)
r=0.1667;
a=(1-r)/(1-2*r/3);
b=r/(1-2*r/3);
p=a+b.*cos(theta).^2;
end








