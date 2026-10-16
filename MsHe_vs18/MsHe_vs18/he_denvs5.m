function He_den=he_denvs5(z,Dif,Hn,Loss,qe)
% 在第一、二和最后一、二个值之间插值，作为边界条件
% 在vs3的基础上，进一步的更改上边界条件
% 并研究解对上边界值的敏感性
% 结果有影响，上边界条件对结果的影响很大
% 与前人工作相比，本次改进明显更符合

nz=length(z);

f = zeros(nz, nz);
g = zeros(nz,1);

for n =3:nz-2
    dz2=(z(n)-z(n-2))*10^5;
    dz3=(z(n+1)-z(n-1))*10^5;
    dz1=(z(n+2)-z(n))*10^5;
    f(n, n-2) = -Dif(n-1)/dz2;
    f(n, n-1) = Dif(n-1)/Hn(n-1);
    f(n, n  ) = Loss(n)*dz3+Dif(n+1)/dz1+Dif(n-1)/dz2;
    f(n, n+1) = -Dif(n+1)/Hn(n+1);
    f(n, n+2) = -Dif(n+1)/dz1;
    g(n) = qe(n)*dz3;
end



f(1,1)=Loss(1);
g(1)=qe(1);


f(2,2)=Loss(2);
g(2)=qe(2);





% f(nz,nz-2)=-Dif(nz-1);
% dz1=(z(nz)-z(nz-2))*10^5;   % cm
% f(nz,nz-1)=Dif(nz-1)/Hn(nz-1)*dz1;
% dz2=(z(nz)-z(nz-1))*2*10^5;   % cm
% f(nz,nz)=Dif(nz-1)+Loss(nz)*dz1*dz2;
% g(nz)=qe(nz)*dz1*dz2;
% 
% 
% 
% f(nz-1,nz-3)=-Dif(nz-2);
% dz1=(z(nz-1)-z(nz-3))*10^5;   % cm
% f(nz-1,nz-2)=Dif(nz-2)/Hn(nz-2)*dz1;
% dz2=(z(nz)-z(nz-2))*10^5;   % cm
% f(nz-1,nz-1)=Dif(nz-2)+Loss(nz-1)*dz1*dz2;
% g(nz-1)=qe(nz-1)*dz1*dz2;


dz1=(z(nz)-z(nz-2))/2*10^5;   % cm
dz2=(z(nz)-z(nz-1))*10^5;   % cm
f(nz,nz-2)=-Dif(nz-1)/2/dz1;
f(nz,nz-1)=Dif(nz-1)/Hn(nz-1)+(Loss(nz)+Loss(nz-1))/2/2*dz2;
f(nz,nz)=Dif(nz-1)/2/dz1+(Loss(nz)+Loss(nz-1))/2/2*dz2;
g(nz)=(qe(nz)+qe(nz-1))/2*dz2;


dz1=(z(nz-1)-z(nz-3))*10^5;   % cm
dz2=(z(nz)-z(nz-2))*10^5;   % cm
f(nz-1,nz-3)=-Dif(nz-2)/dz1;
f(nz-1,nz-2)=Dif(nz-2)/Hn(nz-2);
f(nz-1,nz-1)=Dif(nz-2)/dz1+Loss(nz-1)*dz2;
g(nz-1)=qe(nz-1)*dz2;



He_den=f\g;%((f'*f)^-1)*f'*g';
end