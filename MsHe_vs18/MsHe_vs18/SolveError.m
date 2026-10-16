function per=SolveError(ze,qe,Loss,result)
%% 本程序用来计算亚稳态氦分布方程解的相对误差的
dze=ze(2)-ze(1);
nze=length(ze);
LossTotal=Loss(1)*result(1)*dze/2;
for i=2:nze-1
    LossTotal=LossTotal+Loss(i)*result(i)*dze;
end
LossTotal=LossTotal+Loss(end)*result(end)*dze/2;
PTotal=sum(qe)*dze-(qe(1)+qe(end))*dze/2;
per=abs(LossTotal-PTotal)*2/(LossTotal+PTotal);
end