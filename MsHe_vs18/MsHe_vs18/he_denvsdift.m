function n_he_all = he_denvsdift(z, Dif, Hn, Loss, qe, dt, nt, n_he_init, z1, z2, ddift, ddif)
% 输入：z(1×nz, km)、Dif(1×nz, cm²/s)、Hn(1×nz, cm)、Loss(1×nz, s⁻¹)、qe(1×nz, cm⁻³s⁻¹)
%       dt(时间步长, s)、nt(总时间步数)、n_he_init(初始密度, 1×nz)
%       z1/z2(区间, km)、ddift(持续时间, s)、ddif(系数)
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
dHdz = zeros(nz, 1);
for i = 2:nz-1
    dHdz(i) = (Hn(i+1) - Hn(i-1)) / (2*dz);
end
dHdz(1) = (Hn(2) - Hn(1)) / dz;
dHdz(nz) = (Hn(nz) - Hn(nz-1)) / dz;

for j = 1:nt
    n_he_prev = n_he_all(:, j);  % 上一时刻密度

    % 1. 当前时刻扩散系数
    Dif_curr = Dif;
    if dt*j >= 0 && dt*j <= ddift
        Dif_curr(z1-99:z2-99) = Dif(z1-99:z2-99)/ddif;
    end

    % 2. 计算扩散系数导数dDdz
    dDdz = zeros(nz, 1);
    for i = 2:nz-1
        dDdz(i) = (Dif_curr(i+1) - Dif_curr(i-1)) / (2*dz);  % 中心差分
    end
    dDdz(1) = (Dif_curr(2) - Dif_curr(1)) / dz;              % 下边界前向差分
    dDdz(nz) = (Dif_curr(nz) - Dif_curr(nz-1)) / dz;         % 上边界后向差分

    % 3. 实时构建M矩阵（左端：n^(m+1)项系数）和B矩阵（右端：n^m项系数）
    M = eye(nz);
    B = eye(nz);
    C0 = zeros(nz, 1);
    C1 = zeros(nz, 1);
    C2 = zeros(nz, 1);

    % 内部点
    for i = 2:nz-1
        C0(i) = (dDdz(i)/Hn(i)) - (Dif_curr(i)/(Hn(i)^2))*dHdz(i) - Loss(i);
        C1(i) = dDdz(i) + Dif_curr(i)/Hn(i);
        C2(i) = Dif_curr(i);

        % M矩阵（左端：n^(m+1)系数）
        M(i, i-1) = (dt*C1(i))/(4*dz) - (dt*C2(i))/(2*dz^2);
        M(i, i)   = 1 - (dt*C0(i))/2 + (dt*C2(i))/dz^2;
        M(i, i+1) = - (dt*C1(i))/(4*dz) - (dt*C2(i))/(2*dz^2);

        % B矩阵（右端：n^m系数）
        B(i, i-1) = - (dt*C1(i))/(4*dz) + (dt*C2(i))/(2*dz^2);
        B(i, i)   = 1 + (dt*C0(i))/2 - (dt*C2(i))/dz^2;
        B(i, i+1) = (dt*C1(i))/(4*dz) + (dt*C2(i))/(2*dz^2);
    end

    % 下边界：局部化学平衡
    M(1, 1) = 1 + dt*Loss(1)/2;
    B(1, 1) = 1 - dt*Loss(1)/2;

    % 上边界：零净通量
    C0(nz) = (dDdz(nz)/Hn(nz)) - (Dif_curr(nz)/(Hn(nz)^2))*dHdz(nz) - Loss(nz);
    C1(nz) = dDdz(nz) + Dif_curr(nz)/Hn(nz);
    C2(nz) = Dif_curr(nz);
    M(nz, nz-2) = -dt*C2(nz)/dz^2/2;
    M(nz, nz-1) = dt*C2(nz)/dz^2;
    M(nz, nz)   = 1 + dt*(-C0(nz)+C1(nz)/Hn(nz)-C2(nz)/dz^2)/2;
    B(nz, nz-2) = dt*C2(nz)/dz^2/2;
    B(nz, nz-1) = -dt*C2(nz)/dz^2;
    B(nz, nz)   = 1 + dt*(C0(nz)-C1(nz)/Hn(nz)+C2(nz)/dz^2)/2;

    % 4. 计算右端项
    term_rhs = B * n_he_prev + qe*dt;

    % 5. 隐式求解线性方程组
    n_he_curr = M \ term_rhs;

    % 6. 当前时刻密度
    n_he_all(:, j+1) = n_he_curr;
end

% 保存总密度数据
n_total = [dt*(1:nt+1)', sum(n_he_all(:, 1:nt+1), 1)'];
writematrix(n_total, 'n_total.txt');
end