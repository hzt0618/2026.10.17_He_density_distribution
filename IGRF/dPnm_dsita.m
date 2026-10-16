function res=dPnm_dsita(n,m,sita)

res=(n*cos(sita).*Pnm(n,m,sita)-sqrt((n-m)*(n+m))*Pnm(n-1,m,sita)) ...
    ./sin(sita);

end
