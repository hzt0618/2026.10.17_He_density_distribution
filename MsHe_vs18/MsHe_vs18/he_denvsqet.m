function n_he_all = he_denvsqet(z, Dif, Hn, Loss, qe, dt, nt, n_he_init, z1, z2, dqet, dqe)
% 输入：z(1×nz, km)、Dif(1×nz, cm²/s)、Hn(1×nz, cm)、Loss(1×nz, s⁻¹)、qe(1×nz, cm⁻³s⁻¹)
%       dt(时间步长, s)、nt(总时间步数)、n_he_init(初始密度, 1×nz)
%       z1/z2(源项增强区间, km)、dqet(增强持续时间, s)、dqe(增强系数)
% 输出：n_he_all(nz×(nt+1), cm⁻³)：各时间层密度矩阵

% 变量预处理
nz = length(z);
Dif = Dif(:);       % nz×1
Loss = Loss(:);     % nz×1
Hn = Hn(:);         % nz×1
qe = qe(:);         % nz×1
n_he_init = n_he_init(:);  % nz×1
n_he_all = zeros(nz, nt+1);
n_he_all(:, 1) = n_he_init;  % 初始条件

dz_km = z(2) - z(1);
dz = dz_km * 1e5;  % 转换为cm

% 预计算导数项
dDdz = zeros(nz, 1);
dHdz = zeros(nz, 1);
for i = 2:nz-1
    dDdz(i) = (Dif(i+1) - Dif(i-1)) / (2*dz);
    dHdz(i) = (Hn(i+1) - Hn(i-1)) / (2*dz);
end
dDdz(1) = (Dif(2) - Dif(1)) / dz;
dHdz(1) = (Hn(2) - Hn(1)) / dz;
dDdz(nz) = (Dif(nz) - Dif(nz-1)) / dz;
dHdz(nz) = (Hn(nz) - Hn(nz-1)) / dz;

M = eye(nz);
B = eye(nz);

% 1. 内部点（2≤i≤nz-1）
for i = 2:nz-1
    C0 = (dDdz(i)/Hn(i)) - (Dif(i)/(Hn(i)^2))*dHdz(i) - Loss(i);
    C1 = dDdz(i) + Dif(i)/Hn(i);
    C2 = Dif(i);
    
    % 左端矩阵M（n^(m+1)项系数）
    M(i, i-1) = (dt*C1)/(4*dz) - (dt*C2)/(2*dz^2);
    M(i, i) = 1 - (dt*C0)/2 + (dt*C2)/dz^2;
    M(i, i+1) = - (dt*C1)/(4*dz) - (dt*C2)/(2*dz^2);
    
    % 右端矩阵B（n^m项系数）
    B(i, i-1) = - (dt*C1)/(4*dz) + (dt*C2)/(2*dz^2);
    B(i, i) = 1 + (dt*C0)/2 - (dt*C2)/dz^2;
    B(i, i+1) = (dt*C1)/(4*dz) + (dt*C2)/(2*dz^2);
end

% 下边界：局部化学平衡
    M(1, 1) = 1 + dt*Loss(1)/2;
    B(1, 1) = 1 - dt*Loss(1)/2;

% 上边界：零净通量
    % term1 = Dif(nz)/Hn(nz) + (Dif(nz)/Hn(nz))*(Hn(nz)-Hn(nz-1))/dz;
    % f_nz_nz = term1/Hn(nz) + 2*Dif(nz)/(dz*Hn(nz)) + 2*Dif(nz)/dz^2 + Loss(nz);
    % f_nz_nz1 = -2*Dif(nz)/dz^2;
    % M(nz, nz-1) = dt*f_nz_nz1 / 2;
    % M(nz, nz) = 1 + dt*f_nz_nz / 2;
    % B(nz, nz-1) = -dt*f_nz_nz1 / 2;
    % B(nz, nz) = 1 - dt*f_nz_nz / 2;
    C0(nz) = (dDdz(nz)/Hn(nz)) - (Dif(nz)/(Hn(nz)^2))*dHdz(nz) - Loss(nz);
    C1(nz) = dDdz(nz) + Dif(nz)/Hn(nz);
    C2(nz) = Dif(nz);
    M(nz, nz-2) = -dt*C2(nz)/dz^2/2;
    M(nz, nz-1) = dt*C2(nz)/dz^2;
    M(nz, nz)   = 1 + dt*(-C0(nz)+C1(nz)/Hn(nz)-C2(nz)/dz^2)/2;
    B(nz, nz-2) = dt*C2(nz)/dz^2/2;
    B(nz, nz-1) = -dt*C2(nz)/dz^2;
    B(nz, nz)   = 1 + dt*(C0(nz)-C1(nz)/Hn(nz)+C2(nz)/dz^2)/2;
    
% 4. 隐式时间步进
for j = 1:nt
    n_he_prev = n_he_all(:, j);  % 上一时刻密度
    
    % 源项P（局部增强）
    qqee = qe;
    if dt*j >= 0 && dt*j <= dqet
        qqee(z1-99 : z2-99) = qe(z1-99 : z2-99)*dqe;
    end
    
    % 右端项：B*n_prev + dt*P
    term_rhs = B*n_he_prev + dt*qqee;
    
    % 求解当前时间密度
    n_he_curr = M \ term_rhs;
    n_he_all(:, j+1) = n_he_curr;
end

% 保存总密度数据
n_total = [dt*(1:nt+1)', sum(n_he_all(:, 1:nt+1), 1)'];
writematrix(n_total, 'n_total.txt');
end