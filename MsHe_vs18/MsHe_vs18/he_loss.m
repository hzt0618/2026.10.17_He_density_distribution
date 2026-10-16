function Loss=he_loss(O1,O2,N2,Tem)
%/s  from Waldrop for He(23S) -> He(11S), not consider in Patterson 1967
k1 = 1.1 * 10.^(-4);    % s^-1

%/s  from Waldrop, 2005
J2 = 1.9 *  10.^(-3);    % s^-1

%J2 = 1.61 * 10.^(-3) ; /s  from Patterson 1967

k3 = 1.03 * 10.^(-9) * (Tem/750).^(1./6.);  %cm3/s for O1 from Waldrop, 2005
%k3 = 1.25 * 10.^(-10) ;cm3/s from Patterson 1967

k4 = 1.90 * 10.^(-9) * exp(-0.63-418./Tem);    %cm3/s for O2 from Waldrop, 2005
%k4 = 3.50 * 10.^(-10) ; cm3/s from Patterson 1967

k5 = 1.90 * 10.^(-9) *exp(-1.04-685./Tem);   %cm3/s for N2 from Waldrop, 2005
%k5 = 1.75 * 10.^(-10) ; cm3/s from Patterson 1967

Loss = k1+ J2 + k3.*O1 + k4.*O2 + k5.*N2;
end