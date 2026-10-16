function qr=he_produc2(heion,elecden,electem)
% from Verner and Ferland [1996]
% heion       cm^-3 from IRI
% elecden     cm^-3 from IRI
% electem     K, electron temperature from IRI
% alpha_r     cm^3 s^-1 from Verner and Ferland [1996]
% Te          K
% qr          cm^-3 s^-1, production rate of recombination 
Te=[100,316,1000,3160,10000];
alpha_r=[17,9.6,5.6,3.2,1.6]*10^(-13);
alpha_r=interp1(Te,alpha_r,electem,'pchip','extrap');
qr=(alpha_r.*heion).*elecden;
% qr(750:850) = qr(750:850) * 1e1;

end