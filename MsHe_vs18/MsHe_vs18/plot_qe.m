%% 产生率、氦离子绘图
close all
% figure;
% subplot(2,1,1);
% plot(qe, z, 'b-', 'linewidth', 1.2, 'DisplayName', 'qe');
% xlabel('qe [cm^{-3}s^{-1}]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
% ylabel('z [km]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
% set(gca, 'FontName','Times New Roman','FontSize',10,'FontWeight','Bold','linewidth',1);
% text(-0.1, 1.1, '(a)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% legend('Location','best','FontName','Times New Roman','FontSize',10);
% grid on; box on;
% subplot(2,1,2);
% plot(qr, z, 'r-', 'linewidth', 1.2, 'DisplayName', 'qr');
% xlabel('qr [cm^{-3}s^{-1}]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
% ylabel('z [km]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
% set(gca, 'FontName','Times New Roman','FontSize',10,'FontWeight','Bold','linewidth',1);
% text(-0.1, 1.1, '(b)', 'Units', 'normalized', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
% legend('Location','best','FontName','Times New Roman','FontSize',10);
% grid on; box on;

figure;
plot(qe, z, 'b-', 'linewidth', 1.2, 'DisplayName', 'qe');
hold on;
plot(qr, z, 'r-', 'linewidth', 1.2, 'DisplayName', 'qr');
xlabel('metstable helium production rate:q [cm^{-3}s^{-1}]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
ylabel('z [km]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
set(gca, 'FontName','Times New Roman','FontSize',10,'FontWeight','Bold','linewidth',1);
legend('Location','best','FontName','Times New Roman','FontSize',10);
grid on; box on;
hold off;

% figure
% plot(heion, z, 'b-', 'linewidth', 1.2, 'DisplayName', 'He^{+}');
% xlabel('He^{+} number density [cm^{-3}]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
% ylabel('z [km]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
% set(gca, 'FontName','Times New Roman','FontSize',10,'FontWeight','Bold','linewidth',1);
% 
% figure
% plot(elecden, z, 'b-', 'linewidth', 1.2, 'DisplayName', 'elecden');
% xlabel('electron density: n_{e} [cm^{-3}]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
% ylabel('z [km]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
% set(gca, 'FontName','Times New Roman','FontSize',10,'FontWeight','Bold','linewidth',1);
% 
% figure
% plot(electem, z, 'b-', 'linewidth', 1.2, 'DisplayName', 'electem');
% xlabel('electron temperature: T_{e} [K]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
% ylabel('z [km]','FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
% set(gca, 'FontName','Times New Roman','FontSize',10,'FontWeight','Bold','linewidth',1);
