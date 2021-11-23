addpath /home/sergio/MATLABCODE/PLOTMISC   %% for shadedErrorBar
addpath /home/sergio/MATLABCODE/           %% for plotaxis2

figure(21);
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
plot(f,airsobs,'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'),'r',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.umbc_spectral_rates'),'m',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.umbc_spectral_rates'),'m')
hold on
shadedErrorBar(f,airsobs,airsobs_unc,'b',0.3);
plot(f,airsobs,'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'),'r')
hold off
plotaxis2;

tafov(2) = nexttile;
plot(f,airsobs,'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates'),'r',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.airsL3_spectral_rates'),'m',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.airsL3_spectral_rates'),'m')
hold on
shadedErrorBar(f,airsobs,airsobs_unc,'b',0.3);
plot(f,airsobs,'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates'),'r')
hold off
plotaxis2;

tafov(3) = nexttile;
plot(f,airsobs,'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'),'r',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.era5_spectral_rates'),'m',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.era5_spectral_rates'),'m')
hold on
shadedErrorBar(f,airsobs,airsobs_unc,'b',0.3);
plot(f,airsobs,'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'),'r')
hold off
plotaxis2;

tafov(4) = nexttile;
plot(f,airsobs,'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates'),'r',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.cmip6_spectral_rates'),'m',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.cmip6_spectral_rates'),'m')
hold on
shadedErrorBar(f,airsobs,airsobs_unc,'b',0.3);
plot(f,airsobs,'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates'),'r')
hold off
plotaxis2;

%xlabel(tafov,'Wavenumber cm-1');
%ylabel(tafov,'dBT(K)/dt');
tafov(1).YLabel.String = 'd(BT)/dt K/yr'; tafov(1).YLabel.FontSize = 10;
tafov(3).YLabel.String = 'd(BT)/dt K/yr'; tafov(3).YLabel.FontSize = 10;
tafov(3).XLabel.String = 'Wavenumber cm^{-1}'; tafov(3).XLabel.FontSize = 10;
tafov(4).XLabel.String = 'Wavenumber cm^{-1}'; tafov(4).XLabel.FontSize = 10;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;
tafov(4).FontSize = 10;

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
for ii = [2 4]
   tafov(ii).YTickLabel = '';
end
% Remove all xtick labels except for 2nd row
for ii = [1 2]
   tafov(ii).XTickLabel = '';
end

%% put titles
title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(3),'ERA5', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);
title(tafov(4),'CMIP6', 'Units', 'normalized', 'Position', [0.5, +1.05, 0]);

set(tafov,'Xlim',[640 1640]);
set(tafov,'Ylim',[-0.1 +0.1]);

figure(22); clf
plot(f,airsobs,f,airsobs_unc,...
     f,abs((nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.umbc_spectral_rates')-nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'))),...
     f,abs((nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.airsL3_spectral_rates')-nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates'))),...
     f,abs((nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.era5_spectral_rates')-nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'))),...
     f,abs((nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.cmip6_spectral_rates')-nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates'))),...
     'linewidth',0.5)
plotaxis2; xlim([640 1640]);
hl = legend('AIRS obs dBT/dt','AIRS unc dBT/dt','UMBC','AIRS L3','ERA5','CMIP6','location','south','fontsize',8);
ylabel('dBT/dt (K/yr)'); xlabel('Wavenumber cm^{-1}');

figure(22); clf
plot(f,airsobs_unc,...
     f,abs((nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.airsL3_spectral_rates')-nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates'))),...
     f,abs((nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.era5_spectral_rates')-nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'))),...
     'linewidth',0.5)
plotaxis2; xlim([1240 1640]);
hl = legend('AIRS obs','AIRS L3','ERA5','location','best','fontsize',8);
ylabel('dBT/dt unc (K/yr)'); xlabel('Wavenumber cm^{-1}');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
addpath /asl/matlib/plotutils
figure(21); figname = 'spectral_plot_with_unc.pdf';  figure(21); aslprint(figname)
figure(22); figname = 'spectral_plot_with_unc2.pdf'; figure(22); aslprint(figname)


save strow_rates_with_unc.mat f *nwp_spectral_trends_cmip6_era5_airsL3_umbc* airsobs airsobs_unc
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(23); clf;

plot(f,airsobs,'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'),'r',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.umbc_spectral_rates'),'m',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.umbc_spectral_rates'),'g')
hold on
shadedErrorBar(f,airsobs,airsobs_unc,'b',0.3);
plot(f,airsobs,'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'),'r')
hold off
plotaxis2;
title('UMBC')
