%% see /home/sergio/MATLABCODE/~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/contrast_05_10_15_20.m

load fig3.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf
plot(data_fig3a.vchan,data_fig3a.q_mean_BT); xlabel('Wavenumber cm-1'); ylabel('BT(K)');
xlim([640 1640]); plotaxis2; hl = legend('Q50','Q80','Q90','Q95','Q97','location','best','fontsize',8);
set(gca,'fontsize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf
ta = tiledlayout(2,1);
tafov(1) = nexttile; plot(data_fig3b.vchan,data_fig3b.q01.trend,data_fig3b.vchan,data_fig3b.q02.trend,data_fig3b.vchan,data_fig3b.q03.trend,data_fig3b.vchan,data_fig3b.q04.trend,data_fig3b.vchan,data_fig3b.q05.trend); 
  plotaxis2;
  xlim([640 1620]); ylabel('d(BT)/dt [K yr^{-1}]');
  %legend('Q50','Q80','Q90','Q95','Q97','fontsize',10,'location','south');
  %legend('Q50','Q80','Q90','Q95','Q97','fontsize',10,'Position',[0.45 0.45 0.25 0.25])
  set(gca,'fontsize',14)

tafov(2) = nexttile; plot(data_fig3b.vchan,data_fig3b.q01.unc,data_fig3b.vchan,data_fig3b.q02.unc,data_fig3b.vchan,data_fig3b.q03.unc,data_fig3b.vchan,data_fig3b.q04.unc,data_fig3b.vchan,data_fig3b.q05.unc); 
  plotaxis2;
  xlim([640 1620]);
  %legend('Q50','Q80','Q90','Q95','Q97','fontsize',8,'location','southeast');
  legend('Q50','Q80','Q90','Q95','Q97','fontsize',10,'Position',[0.45 0.4625 0.2 0.2])
  xlabel('Wavenumber [cm^{-1}]'); ylabel('d(BT)/dt [K yr^{-1}]');
  ylim([0 0.031])
  set(gca,'fontsize',14)

ta.Padding = 'none';
ta.TileSpacing = 'compact';
ta.Padding = 'compact';
ta.TileSpacing = 'tight';
tafov(1).XTickLabel = '';
