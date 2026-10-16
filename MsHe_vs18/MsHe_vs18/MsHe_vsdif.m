%final
% 相对于vs15，vs17修正了对碰撞截面的计算。修正过后，更方便对横截面进行修改。
% 相对于vs17，vs17_2对某地全年的亚稳态氦密度进行模拟，
% 并得出被动观测强度随全年的变化，2022.02.25
% HaiNan2021_2相比于HaiNan2021相比，没有考虑qr
% HaiNan2021_vs3相比于HaiNan2021_vs2，把heion乘以1/5
% HaiNan2021_vs4相比于HaiNan2021_vs2，把heion乘以1/10
% 更改扩散系数dif

clear;
close all

%% 设置一些标记参数
whpp=0;  %  如果whpp=0，则需要考虑光压；如果whpp=1，则无需考虑光压；photon press
whre=1;  %  如果whre=0，则考虑共振发光，是按照被动观测的公式计算的；
         %  如果whre=1，则无需考虑共振，是为激光雷达观测准备的；resonance
% k_D = 1.5e10;  %考虑电子碰撞使扩散系数减小的系数/s
usehu=1000;
usehd=100;
interval=1;  %高度间隔(km)
nz=(usehu-usehd)/interval+1;
z=linspace(usehd,usehu,nz);  %km
exz=1500;      % extension height, km

%% 时间参数
% 时间
dt = 1e0;    % 时间步长(s) 不能太小
nt = 1e2;   % 步数(总时间dt*nt)
% 初始条件
z1 = 950; z2 = 990;  %初始变化高度(km)
ddift = 1e5;  %初始dif降低持续时间(s)
ddif = 10;  %初始dif降低倍数
n_he_init = zeros(1, nz);
t_key = [0, dt*nt/10, dt*nt/10*2, dt*nt/10*3, dt*nt/10*4, dt*nt/10*5,...
    dt*nt/10*6, dt*nt/10*7, dt*nt/10*8, dt*nt/10*9, dt*nt];  % 5个关键时刻（s）
t_idx = round(t_key / dt) + 1;  % 对应n_he_all的列索引

%% 设置一些随日期变化的参数
list=1:50;
% inst = 50;  %按需调整
nl=length(list);
doy = 6 ;              %1:5:31;
doyi=doy;
disp(doy)
ndoy=length(doy);
sza0=zeros(nl,ndoy);
It1=zeros(nl,ndoy);
It2=zeros(nl,ndoy);
LT=zeros(nl,ndoy);
Doy=zeros(nl,ndoy);
% Qe=zeros(nl,365);
% 定义颜色
co=parula(256);
DawnModule=["Dawn","Module"];
%% 正演
for inst=list(end)  %list
    %for DD=1:2
    for DM = 1     %1,2对应不含/含He+重组
    color=interp1(1:256,co,inst*255/144+1,'pchip','extrap');
    %% read data
    filename=['C:\Users\31217\Desktop\氦激光雷达\地磁场光电子模型\glow model\glow model\HaiNan20240106/',...
                num2str(doyi),DawnModule{DM},'/output_vs',num2str(inst),'.txt'];
    [jmax,sza,idate,ut,glat,glong, ...
        f107a,f107,f107p,ap,ef,ec,ener,del,eflx]=read_glowdata2(z,filename);

    % filename=['C:\Users\admin\Desktop\氦激光雷达\地磁场光电子模型\glow model\glow model\HaiNan2024/',...
    %             num2str(doyi),DawnModule{1},'/output_vs',num2str(inst),'.txt'];
    % [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, eflxinit] = read_glowdata2(z, filename);

    if ef>30
        break;
    end
    sza0(find(list==inst)+(DM-1)*nl/2,doy==doyi)=sza;
    LT(find(list==inst)+(DM-1)*nl/2,doy==doyi)=UtToLt(glong,ut/3600);
    %% get DoY and Yr
    [DoY,Yr]=date2doy(idate);
    Doy(find(list==inst)+(DM-1)*nl/2,doy==doyi)=DoY;
    %% get T and Density following altitude from atmosbrlmsise00
    % 获得大气参数
    [Tem,He,O1,N2,O2,den]=gain_atmdentem(z,glat,glong,Yr,DoY,ut,f107a,f107,ap);

    %% (1) the meta He production rate due to photonelectrons,
    % 计算电子碰撞基态氦产生亚稳态氦的速率轮廓      cm^(-3)s^(-1), that is, alpha * He
    qe=he_produc1_vs2(eflx,ener,He);

    %% ;; (3) Loss rate,  s^(-1)
    Loss=he_loss(O1,O2,N2,Tem);
    
    %% (4)  diffusion coefficiant and scale height
    
    % 计算亚稳态氦的扩散系数
    Dif=he_dif(He,O1,Tem);
    
    % calculate metastable He scale rhight
    [Hn,a,gstar,gam]=he_Hn1vs3(z,nz,Tem,sza,glat,whpp);   %cm

    %% (5) 非稳态扩散方程  
    n_he_init = he_denvs6(z, Dif, Hn, Loss, qe);
    n_he = he_denvsdift(z, Dif, Hn, Loss, qe, dt, nt, n_he_init, z1, z2, ddift, ddif);
    % n_he_steady_trans = n_he(:, end);  % 最终时刻
    
    % 原稳态解
    % n_he_original = he_denvs6(z, Dif, Hn, Loss, qe);  % 原程序输出
    % 最大误差
    % err = max(abs(n_he_steady_trans - n_he_original));
    % fprintf('稳态还原误差：%.2e cm⁻³\n', err);

    %% (5) 扩散方程
    % % 根据上下边界条件计算亚稳态氦密度轮廓
    % result=he_denvs6(z,Dif,Hn,Loss,qe);
    % % 计算解的相对误差
    % per=SolveError(z,qe,Loss,result);

    end
end

%% He(23S)随时图
figure('Position', [150, 150, 1000, 700]);  % 整体图尺寸（宽×高）
colors = [
    0.85, 0.11, 0.12;  % 1. 红色
    0.12, 0.56, 0.12;  % 2. 绿色
    0.13, 0.29, 0.80;  % 3. 蓝色
    0.63, 0.27, 0.08;  % 4. 橙色
    0.51, 0.00, 0.50;  % 5. 紫色
    0.00, 0.00, 0.00;  % 6. 黑色/白色
    0.92, 0.75, 0.00;  % 7. 黄色
    0.00, 0.70, 0.73;  % 8. 青色
    0.80, 0.40, 0.80;  % 9. 粉色
    0.60, 0.60, 0.60;  % 10. 深灰色
    0.20, 0.80, 0.20]; % 11. 亮绿色
subplot(1, 2, 1);
hold on; grid on;
% 绘制关键时刻
for k = 1:length(t_key)
    idx = t_idx(k);
    if idx > size(n_he, 2), idx = size(n_he, 2); end  % 防止索引溢出
    plot(n_he(:, idx), z, 'Color', colors(k,:), ...
         'LineWidth', 1.5, 'DisplayName', sprintf('t=%.3fs', t_key(k)));
end
xlabel('$He(2^3S)$ Density (cm$^{-3}$)', 'Interpreter', 'Latex', 'FontSize', 14);
ylabel('Altitude (km)', 'FontSize', 14);
xlim([0, 1]);
ylim([100, 1000]);
title('亚稳态氦密度随时间分布', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 12);
box on;

subplot(1, 2, 2);
hold on; grid on;
z_zoom_idx = z >= 800 & z <= 1000;  % 800-1000km的索引
z_zoom = z(z_zoom_idx);             % 放大后的高度
% 循环绘制对应时刻的局部密度
for k = 1:length(t_key)
    idx = t_idx(k);
    if idx > size(n_he, 2), idx = size(n_he, 2); end
    % 提取800-1000km的密度数据
    n_zoom = n_he(z_zoom_idx, idx);
    plot(n_zoom, z_zoom, 'Color', colors(k,:), ...
         'LineWidth', 1.5, 'DisplayName', sprintf('t=%.3fs', t_key(k)));
end
xlabel('$He(2^3S)$ Density (cm$^{-3}$)', 'Interpreter', 'Latex', 'FontSize', 14);
ylabel('Altitude (km)', 'FontSize', 14);
ylim([800, 1000]);
xlim auto;
title('800-1000km局部放大图', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 12);
box on;
sgtitle(sprintf('%.0f-%.0fkm前%.0fs扩散系数Dif减小%.1f倍时亚稳态氦密度随时间分布', z1, z2, ddift, ddif), 'FontSize', 18, 'FontWeight', 'bold');
set(gcf, 'Position', [150, 150, 1200, 500]);
set(findall(gcf, 'Type', 'Axes'), 'FontSize', 12);  % 统一坐标轴字体大小
hold off;
%% 视频
% figure('Position', [150, 150, 1000, 800]);
% video = VideoWriter('C:\Users\admin\Desktop\氦激光雷达\地磁场光电子模型\模拟视频\1.mp4', 'MPEG-4');
% video.FrameRate = 30;
% open(video);
% for j = 1:nt
%     clf;
%     hold on; grid on; box on;
%     plot(n_he(:,j), z, 'Color', [0.13, 0.29, 0.80], ...
%          'LineWidth', 2.5, 'DisplayName', '亚稳态氦密度');
%     xlabel('$He(2^3S)$ Density (cm$^{-3}$)', 'Interpreter', 'Latex', 'FontSize', 14);
%     ylabel('Altitude (km)', 'FontSize', 14);
%     xlim([0, 1]);
%     ylim([100, 1000]);
%     title(sprintf('%.1f-%.1fkm前%.1fs扩散系数dif减小%.1f倍时亚稳态氦密度随时间分布,t=%.1fs', z1, z2, dqet, dqe, j*dt), 'FontSize', 16, 'FontWeight', 'bold');
%     legend('Location', 'southeast', 'FontSize', 12);
%     frame = getframe(gcf);
%     writeVideo(video, frame);
% end
% close(video);
%% dif随时图
cmap = parula(nt);  %parula viridis
diff = he_dif(He,O1,Tem);
figure;
hold on;grid on;
xlabel('Dif(cm$^{2}$s$^{-1}$)', 'Interpreter', 'Latex');
ylabel('z(km)', 'Interpreter', 'Latex');
for j = 1:nt
    diff_d = diff;
    if dt*j >= 0 && dt*j <= ddift  %(s) 局部减小扩散系数
        diff_d(z1-99:z2-99) = diff(z1-99:z2-99)/ddif;
    end
    color = cmap(j, :);
    plot(diff_d, z, 'Color', color, 'LineWidth', 1.2);
end
c = colorbar;
c.Label.String = 't(s)';  % 颜色条标签
c.Label.FontSize = 10;
c.FontSize = 9;
clim([0, dt*nt]);
c.Ticks = linspace(0, dt*nt, 6);
c.TickLabels = sprintf('%.1f\n', linspace(0, dt*nt, 6));
hold off;
%% n_total随时图
n_t = readmatrix('n_total.txt');
tt = n_t(:,1);
n_total = n_t(:,2);
figure;
plot(tt, n_total, 'b-', 'LineWidth', 1.2);
% ylim([115,120]);
xlabel('t(s)', 'Interpreter', 'Latex');  
ylabel('n$_{total}$(cm$^{-3}$)', 'Interpreter', 'Latex');
%% 对比图
% figure;
% plot(n_he_original, z, 'b-', 'LineWidth', 1.5);
% hold on;
% plot(n_he_steady_trans, z, 'r-', 'LineWidth', 1.5);
% xlabel('密度(cm⁻³)'); ylabel('高度 (km)');
% xlim([0, 1]);
% legend('原稳态解', sprintf('t=%.1fs时的非稳态解', dt*nt));
% title(sprintf('原稳态解与非稳态%.1fs后的对比', dt*nt));

%% 各参数图
% figure('Position', [100, 100, 1000, 800], 'Color', 'w');
% 
% %电子通量
% subplot(2,2,1);
% pos = get(gca, 'Position');
% set(gca, 'Position', [pos(1), pos(2)-0.02, pos(3), pos(4)]);  %调整子图间距
% hold on;
% %获取颜色映射（80-35+1=46）
% colors = jet(46);
% %循环绘制第35到80列的曲线
% for i = 35:80
%     %使用i-34作为颜色索引（因为i从35开始）
%     plot(eflx(:, i), z, 'Color', colors(i-34, :)); 
% end
% %添加标签和标题
% xlabel('eflx [cm^{-2}s^{-1}ev^{-1}]');
% ylabel('z [km]');
% set(gca,'FontName','Times New Roman','FontSize',10,'FontWeight','Bold','linewidth',1);
% set(gca, 'XScale', 'log');
% %添加颜色条
% colormap(colors);
% %调整刻度位置和标签（对应35到80列的能量值）
% tick_positions = linspace(1, 46, 10);  % 46个数据点
% tick_indices = round(linspace(35, 80, 10));  % 对应的列索引
% tick_values = ener(tick_indices);  % 获取对应的能量值
% c = colorbar;
% c.Ticks = linspace(0, 1, 10);
% c.TickLabels = arrayfun(@(x) sprintf('%.2f', x), tick_values, 'UniformOutput', false);
% c.Label.String = 'eV';
% c.Label.FontSize = 10;
% c.Label.FontWeight = 'bold';
% c.Label.Rotation = 0;
% c.Label.VerticalAlignment = 'bottom';
% c.Label.Position = [0.5 1.02 0];
% text(-0.1, 1.1, '(a)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% hold off;
% 
% %产生率
% subplot(2,2,2);
% pos = get(gca, 'Position');
% set(gca, 'Position', [pos(1)-0.01, pos(2)-0.02, pos(3), pos(4)]);
% plot(qe,z, 'linewidth',1.2)
% xlabel('metstable excitation rate:qe [cm^{-3}s^{-1}]')
% ylabel('z [km]')
% set(gca,'FontName','Times New Roman','FontSize',10,'FontWeight','Bold','linewidth',1);
% text(-0.1, 1.1, '(b)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% 
% % %损失率
% % subplot(2,2,3);
% % plot(Loss,z, 'linewidth',1.2)
% % xlabel('metastable He loss rates:Loss [s^{-1}]')
% % ylabel('z [km]')
% % set(gca,'FontName','Times New Roman','FontSize',10,'FontWeight','Bold','linewidth',1);
% % text(-0.1, 1.1, '(c)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% 
% %扩散系数
% subplot(2,2,3);
% plot(Dif,z, 'linewidth',1.2)
% xlabel('The diffusion coefficient:D [cm^{2}s^{-1}]')
% ylabel('z [km]')
% set(gca,'FontName','Times New Roman','FontSize',10,'FontWeight','Bold','linewidth',1);
% text(-0.1, 1.1, '(c)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% 
% %亚稳态氦密度
% subplot(2,2,4);
% pos = get(gca, 'Position');
% set(gca, 'Position', [pos(1)-0.01, pos(2), pos(3), pos(4)]);
% % figure;
% fileID=fopen('C:\Users\admin\Desktop\氦激光雷达\地磁场光电子模型\sza.txt','r');
% szaa=fscanf(fileID, '%f');
% fclose(fileID);
% [max_n_he, max_idx] = max(n_he_init);
% max_altitude = z(max_idx);
% plot(n_he_init,z, 'b-', 'linewidth',1.2);
% hold on;
% plot(max_n_he, max_altitude, 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
% xline(0, '--k', 'LineWidth', 0.8); 
% xlim([10^(-20),1.0])
% ylim([100,1000])
% xlabel('He(2^3S) numebr density [cm^{-3}]')
% ylabel('Altitude [km]')
% text(0.95 * max(xlim), 0.8 * max(ylim), sprintf(['sza=%.2f°\n' ...
%     'Ap:3\n' ...
%     'F10.7:159\n' ...
%     'F10.7a:149.3\n' ...
%     'F10.7p:153'], szaa), ...
%     'HorizontalAlignment', 'right', ...
%     'FontSize', 10, ...
%     'FontName', 'Times New Roman', ...
%     'EdgeColor', 'black', ...
%     'LineWidth', 1, ...
%     'BackgroundColor', 'white');
% % set(gcf, 'Position', [500, 400, 600, 500]);
% set(gca,'FontName','Times New Roman','FontSize',10,'FontWeight','Bold','linewidth',1);
% text(-0.1, 1.1, '(d)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% 
% %调整整体布局
% set(gcf, 'Color', 'w');
% ha = findobj(gcf,'type','axes');
% for i = 1:length(ha)
%     ha(i).Title.FontSize = 12;
%     ha(i).Title.FontWeight = 'bold';
% end

% %添加总标题
% annotation('textbox', [0.1, 0.96, 0.8, 0.01], ...  % 位置：[左, 下, 宽, 高]
%            'String', '', ...  % 总标题文本
%            'FontName', 'Times New Roman', ...
%            'FontSize', 14, ...
%            'FontWeight', 'bold', ...
%            'HorizontalAlignment', 'center', ...  % 水平居中
%            'VerticalAlignment', 'middle', ...    % 垂直居中
%            'EdgeColor', 'none', ...              % 无边框
%            'BackgroundColor', 'w');              % 背景色与图形一致

% % sig = sigma(ener);
% % figure(2)
% % plot(ener,sig)
% % xlim([19, 100])

%%
% extra_metaHe = [result, z'];
% writematrix(extra_metaHe, 'C:\Users\admin\Desktop\氦激光雷达\地磁场光电子模型\result.txt', ...
%                             'Delimiter', 'tab'); %result,conjugate
% saveas(gcf, 'C:\Users\admin\Desktop\氦激光雷达\本科毕业论文\MsHe.png'); 
