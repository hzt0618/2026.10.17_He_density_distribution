function Dif=he_dif(He,O1,Tem)
% T0  = 750 ; %750K for test b1, b2
% b1_p0 = 1.4 * 10.^18 * sqrt(T0) ; %/(cm s) Parterson
% b2_p0 = 1.1 * 10.^18 * sqrt(T0)  ; %/(cm s) Parterson
% 
% b1_p1 = 1.4 * 10.^18 * sqrt(Tem)  ; %/(cm s)  Parterson
% b2_p1 = 1.1 * 10.^18 * sqrt(Tem) ; %/(cm s)  Parterson 

b1_w  = 9. * 10.^17 * sqrt(Tem); %/(cm s)  Waldrop
b2_w  = 3. * 10.^18 * sqrt(Tem);

b_he = b1_w;
b_O  = b2_w;

Dif = 1./(He./b_he + O1./b_O);    % cm^2*s^-1
% Dif(750:850) = Dif(750:850) / 15;

end