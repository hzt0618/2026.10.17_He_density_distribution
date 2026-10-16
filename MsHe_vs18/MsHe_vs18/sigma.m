function cross_section = sigma(ener)
    % 输入参数：ener - 1×190的能量数组（单位：eV）
    % 输出参数：cross_section - 1×190的截面数组（单位：cm²）
    
    % 原始数据 (Bishop1999)
    EeV = [19.82,19.9,20.1,20.3,20.5,20.7,20.9,21.1,21.3,21.5, ...
           21.7,21.9,22.1,22.3,22.5,22.7,22.9,30,50,100];    % 能量点(eV)
    cross_raw = [0.00,2.15,5.09,6.10,5.40,3.85,4.88,5.01,4.73,4.80, ...
                 4.92,5.09,5.27,5.23,6.32,5.46,5.53,8.80,4.20,0.51]*1e-18;   % 截面(cm²)
    
    % 初始化截面数组
    cross_section = zeros(size(ener));
    
    % 1. 处理能量低于19.82eV的情况（已初始化为0，无需额外操作）
    
    % 2. 处理能量高于1000eV的情况：先计算1000eV对应的截面
    E_ref = 100;                  % 参考能量点
    sigma_ref = cross_raw(end);   % 100eV对应的截面
    sigma_1000 = sigma_ref * (E_ref / 1000)^2 * exp(-(1000 - E_ref)/500);
    idx_above = ener > 1000;
    cross_section(idx_above) = sigma_1000;
    
    % 3. 处理19.82eV至1000eV之间的能量
    % 3.1 低能段(19.82-22.7eV)：样条插值
    idx_low = (ener >= 19.82) & (ener <= 22.7);
    if any(idx_low)
        cross_section(idx_low) = interp1(EeV(1:16), cross_raw(1:16), ener(idx_low), 'spline');
    end
    
    % 3.2 中能段(22.5-100eV)：PCHIP插值（保证单调性）
    idx_mid = (ener >= 22.5) & (ener <= 100);
    if any(idx_mid)
        cross_section(idx_mid) = interp1(EeV(15:20), cross_raw(15:20), ener(idx_mid), 'pchip');
    end
    
    % 3.3 高能段(100-1000eV)：物理模型计算
    idx_high = (ener > 100) & (ener <= 1000);
    if any(idx_high)
        high_ener = ener(idx_high);
        cross_section(idx_high) = sigma_ref .* (E_ref ./ high_ener).^2 .* exp(-(high_ener - E_ref)/500);
    end
end