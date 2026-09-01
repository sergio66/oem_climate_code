clear all64x72_trends
if ~exist('all64x72_trends')
  all64x72_trends = load('iType_2_convert_sergio_clearskygrid_obsonly_Q16.mat');
end

all64x72_trends.b_desc = reshape(all64x72_trends.b_desc,64*72,2645)';
all64x72_trends.b_err_desc = reshape(all64x72_trends.b_err_desc,64*72,2645)';

mean_desc_trend     = nanmean(all64x72_trends.b_desc');
mean_desc_trend_unc = nanmean(all64x72_trends.b_err_desc');

num   = zeros(size(mean_desc_trend'));
denom = zeros(size(mean_desc_trend'));
sig   = zeros(size(mean_desc_trend'));
for ii = 1 : 4608
  num = num + all64x72_trends.b_desc(:,ii)./(all64x72_trends.b_err_desc(:,ii).^2);
  denom = denom + 1./(all64x72_trends.b_err_desc(:,ii).^2);
  sig = sig + (all64x72_trends.b_err_desc(:,ii).^2);;
end
mean_wgt_desc_trend = num./denom;
mean_wgt_desc_trend_unc = 1/sqrt(4608)*sqrt(sig);
mean_wgt_desc_trend_unc = 1/4608*sqrt(sig);

figure(12)
plot(f,mean_desc_trend,'b',f,mean_wgt_desc_trend,'r',f,mean_desc_trend_unc,'c',f,mean_wgt_desc_trend_unc,'m');
  hl = legend('mean','std','mean unbiased','correct sigma','location','best','fontsize',10);

airsobs = mean_desc_trend;
airsobs_unc = mean_desc_trend_unc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

make_UMBC_ERA5_AIRSL3_CMIP6_trends_64x72_withunc

resultsAIRS_unc.resultsTG = resultsunc;
resultsAIRS_unc.resultsWV_err = resultsWVunc;
resultsAIRS_unc.resultsT_err  = resultsTunc;
resultsAIRS_unc.resultsO3_err = resultsO3unc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nwp_spectral_trends_cmip6_era5_airsL3_umbc0 = make_profile_spectral_trends_64x72(cmip6_64x72,era5_64x72,airsL3_64x72,results,resultsWV,resultsT,resultsO3,resultsAIRS_unc,fits,rates,pavg,plays,f,2,0,1);
figure(13); plot(f,nanmean(rates'),f,nanmean(fits'),f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc0.umbc_spectral_rates'))
  plotaxis2; xlim([640 1640]); hl = legend('AIRS OBS','OEM FIT','UMBC','location','best'); grid;

figure(14); plot(f,nanmean(fits'),'x-',...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc0.airsL3_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc0.cmip6_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc0.era5_spectral_rates'),...
                 f,nanmean(rates'));
figure(14); plot(f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc0.umbc_spectral_rates'),'x-',...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc0.airsL3_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc0.cmip6_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc0.era5_spectral_rates'),...
                 f,nanmean(rates'));
  plotaxis2; xlim([640 1640]); hl = legend('UMBC','AIRSL3','CMIP6','ERA5','AIRS Obs','location','best','fontsize',10); grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror = make_profile_spectral_trends_64x72(cmip6_64x72,era5_64x72,airsL3_64x72,results,resultsWV,resultsT,resultsO3,resultsAIRS_unc,fits,rates,pavg,plays,f,2,1,1);

figure(15); plot(f,nanmean(rates'),f,nanmean(fits'),f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.umbc_spectral_rates'))
  plotaxis2; xlim([640 1640]); hl = legend('AIRS OBS','OEM FIT','UMBC','location','best'); grid;

figure(16); plot(f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.umbc_spectral_rates'),'x-',...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.airsL3_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.cmip6_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.era5_spectral_rates'),...
                 f,nanmean(rates'));
  plotaxis2; xlim([640 1640]); hl = legend('UMBC','AIRSL3','CMIP6','ERA5','AIRS Obs','location','best','fontsize',10); grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror = make_profile_spectral_trends_64x72(cmip6_64x72,era5_64x72,airsL3_64x72,results,resultsWV,resultsT,resultsO3,resultsAIRS_unc,fits,rates,pavg,plays,f,2,-1,1);

figure(17); plot(f,nanmean(rates'),f,nanmean(fits'),f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.umbc_spectral_rates'))
  plotaxis2; xlim([640 1640]); hl = legend('AIRS OBS','OEM FIT','UMBC','location','best'); grid;

figure(18); plot(f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.umbc_spectral_rates'),'x-',...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.airsL3_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.cmip6_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.era5_spectral_rates'),...
                 f,nanmean(rates'));
  plotaxis2; xlim([640 1640]); hl = legend('UMBC','AIRSL3','CMIP6','ERA5','AIRS Obs','location','best','fontsize',10); grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_fig21_spectral_closure
