addpath /home/sergio/MATLABCODE/PLOTMISC
figure(21);
ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
plot(f,nanmean(rates'),'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'),'r',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.umbc_spectral_rates'),'m',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.umbc_spectral_rates'),'m')
hold on
shadedErrorBar(f,nanmean(rates'),mean_desc_trend_unc,'b',0.3);
plot(f,nanmean(rates'),'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'),'r')
hold off

tafov(2) = nexttile;
plot(f,nanmean(rates'),'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates'),'r',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.airsL3_spectral_rates'),'m',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.airsL3_spectral_rates'),'m')
hold on
shadedErrorBar(f,nanmean(rates'),mean_desc_trend_unc,'b',0.3);
plot(f,nanmean(rates'),'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates'),'r')
hold off

tafov(3) = nexttile;
plot(f,nanmean(rates'),'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'),'r',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.era5_spectral_rates'),'m',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.era5_spectral_rates'),'m')
hold on
shadedErrorBar(f,nanmean(rates'),mean_desc_trend_unc,'b',0.3);
plot(f,nanmean(rates'),'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'),'r')
hold off

tafov(4) = nexttile;
plot(f,nanmean(rates'),'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates'),'r',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.cmip6_spectral_rates'),'m',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.cmip6_spectral_rates'),'m')
hold on
shadedErrorBar(f,nanmean(rates'),mean_desc_trend_unc,'b',0.3);
plot(f,nanmean(rates'),'b',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates'),'r')
hold off

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
