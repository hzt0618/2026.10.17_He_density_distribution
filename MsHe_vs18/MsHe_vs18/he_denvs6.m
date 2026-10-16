function He_den=he_denvs6(z,Dif,Hn,Loss,qe)
% 在第一、二和最后一、二个值之间插值，作为边界条件
% 在vs5的基础上，进一步的更改上边界条件
% 并研究解对上边界值的敏感性
% 结果有影响，上边界条件对结果的影响很大
% 与前人工作相比，本次改进明显更符合


nz=length(z);



f = zeros(nz, nz);
g = zeros(nz,1);

for i =2:nz-1
    dz=(z(i+1)-z(i-1))/2*10^5;
    f(i, i-1) = (  (Dif(i+1)-Dif(i-1))/2/dz+Dif(i)/Hn(i)  )/2/dz-Dif(i)/dz^2;
    f(i, i  ) = -(  1/Hn(i)*(Dif(i+1)-Dif(i-1))/2/dz - Dif(i)/Hn(i)^2*(Hn(i+1)-Hn(i-1))/2/dz  )+2*Dif(i)/dz^2+Loss(i);
    f(i, i+1) = -(  (Dif(i+1)-Dif(i-1))/2/dz+Dif(i)/Hn(i)  )/2/dz-Dif(i)/dz^2;
    g(i) = qe(i);
end


% 下边界条件：局部化学平衡
f(1,1)=Loss(1);
g(1)=qe(1);


% 上边界条件：零净通量边界条件
dz=(z(nz)-z(nz-1))*10^5;   % cm
f(nz,nz-1)=-2*Dif(nz)/dz^2;
f(nz,nz)=(  Dif(nz)/Hn(nz)+Dif(nz)/Hn(nz)*(Hn(nz)-Hn(nz-1))/dz  )/Hn(nz) + 2*Dif(nz)/dz/Hn(nz) + 2*Dif(nz)/dz^2+Loss(nz);
g(nz)=qe(nz);




He_den=f\g;%((f'*f)^-1)*f'*g';
end