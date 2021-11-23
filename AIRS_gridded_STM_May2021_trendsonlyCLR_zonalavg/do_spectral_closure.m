all64_trends = load('all64_airs_dbt_trends.mat');
mean_desc_trend     = nanmean(all64_trends.airs_dbt_desc');
mean_desc_trend_unc = nanmean(all64_trends.airs_dbt_desc_unc');

num   = zeros(size(mean_desc_trend'));
denom = zeros(size(mean_desc_trend'));
sig   = zeros(size(mean_desc_trend'));
for ii = 1 : 64
  num = num + all64_trends.airs_dbt_desc(:,ii)./(all64_trends.airs_dbt_desc_unc(:,ii).^2);
  denom = denom + 1./(all64_trends.airs_dbt_desc_unc(:,ii).^2);
  sig = sig + (all64_trends.airs_dbt_desc_unc(:,ii).^2);;
end
mean_wgt_desc_trend = num./denom;
mean_wgt_desc_trend_unc = 1/8*sqrt(sig);

figure(12)
plot(f,mean_desc_trend,'b',f,mean_wgt_desc_trend,'r',f,mean_desc_trend_unc,'c',f,mean_wgt_desc_trend_unc,'m');
  hl = legend('mean','std','mean unbiased','correct sigma','location','best','fontsize',10);

airsobs = mean_desc_trend;
airsobs_unc = mean_desc_trend_unc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultsAIRS_unc.resultsTG = resultsunc;
resultsAIRS_unc.resultsWV_err = resultsWVunc;
resultsAIRS_unc.resultsT_err  = resultsTunc;
resultsAIRS_unc.resultsO3_err = resultsO3unc;

nwp_spectral_trends_cmip6_era5_airsL3_umbc = make_profile_spectral_trends_64(cmip6_64,era5_64,airsL3_64,results,resultsWV,resultsT,resultsO3,resultsAIRS_unc,fits,rates,pavg,plays,f,2,0,1);

figure(13); plot(f,nanmean(rates'),f,nanmean(fits'),f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'))
  plotaxis2; xlim([640 1640]); hl = legend('AIRS OBS','OEM FIT','UMBC','location','best'); grid; 

figure(14); plot(f,nanmean(fits'),'x-',...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'),...
                 f,nanmean(rates'));
figure(14); plot(f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'),'x-',...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'),...
                 f,nanmean(rates'));
  plotaxis2; xlim([640 1640]); hl = legend('UMBC','AIRSL3','CMIP6','ERA5','AIRS Obs','location','best','fontsize',10); grid; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror = make_profile_spectral_trends_64(cmip6_64,era5_64,airsL3_64,results,resultsWV,resultsT,resultsO3,resultsAIRS_unc,fits,rates,pavg,plays,f,2,1,1);

figure(15); plot(f,nanmean(rates'),f,nanmean(fits'),f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.umbc_spectral_rates'))
  plotaxis2; xlim([640 1640]); hl = legend('AIRS OBS','OEM FIT','UMBC','location','best'); grid; 

figure(16); plot(f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.umbc_spectral_rates'),'x-',...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.airsL3_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.cmip6_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Perror.era5_spectral_rates'),...
                 f,nanmean(rates'));
  plotaxis2; xlim([640 1640]); hl = legend('UMBC','AIRSL3','CMIP6','ERA5','AIRS Obs','location','best','fontsize',10); grid; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror = make_profile_spectral_trends_64(cmip6_64,era5_64,airsL3_64,results,resultsWV,resultsT,resultsO3,resultsAIRS_unc,fits,rates,pavg,plays,f,2,-1,1);

figure(17); plot(f,nanmean(rates'),f,nanmean(fits'),f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.umbc_spectral_rates'))
  plotaxis2; xlim([640 1640]); hl = legend('AIRS OBS','OEM FIT','UMBC','location','best'); grid; 

figure(18); plot(f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.umbc_spectral_rates'),'x-',...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.airsL3_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.cmip6_spectral_rates'),...
                 f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc_Merror.era5_spectral_rates'),...
                 f,nanmean(rates'));
  plotaxis2; xlim([640 1640]); hl = legend('UMBC','AIRSL3','CMIP6','ERA5','AIRS Obs','location','best','fontsize',10); grid; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noTG_nwp_spectral_trends_cmip6_era5_airsL3_umbc = make_profile_spectral_trends_64(cmip6_64,era5_64,airsL3_64,results,resultsWV,resultsT,resultsO3,resultsAIRS_unc,fits,rates,pavg,plays,f,2,0,0);

ratesOut = noTG_nwp_spectral_trends_cmip6_era5_airsL3_umbc.ratesOut;
figure(19); plot(f,nanmean(ratesOut'),f,nanmean(fits'),f,nanmean(noTG_nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'))
  plotaxis2; xlim([640 2780]); hl = legend('AIRS OBS','OEM FIT','UMBC','location','best'); grid; 

figure(19); plot(f,nanmean(fits'),'x-',...
                 f,nanmean(noTG_nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates'),...
                 f,nanmean(noTG_nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates'),...
                 f,nanmean(noTG_nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'),...
                 f,nanmean(ratesOut'));
figure(20); plot(f,nanmean(noTG_nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'),'x-',...
                 f,nanmean(noTG_nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates'),...
                 f,nanmean(noTG_nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates'),...
                 f,nanmean(noTG_nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'),...
                 f,nanmean(ratesOut'));
  plotaxis2; xlim([640 2780]); hl = legend('UMBC','AIRSL3','CMIP6','ERA5','AIRS Obs','location','best','fontsize',10); grid; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_fig21_spectral_closure
