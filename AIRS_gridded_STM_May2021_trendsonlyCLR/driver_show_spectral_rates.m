driver_get_the_model_trends

obs = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX100_50fatlayers_CLIMCAPS_MERRA2_AMIP6_feedback.mat','rates');
obs = obs.rates;

load h2645structure.mat

subplot(211); 
plot(airsL3.fchanx,nanmean(obs'),airsL3.fchanx,nanmean(airsL3.airsL3_spectral_rates'),...
     airsL3.fchanx,nanmean(era5.era5_spectral_rates'),airsL3.fchanx,nanmean(cmip6.cmip6_spectral_rates'),'linewidth',2);
plotaxis2;
hl = legend('AIRS OBS','AIRS L3','ERA5','CMIP6','location','best','fontsize',10); ylabel('dBT/dt (K/yr)'); 
xlim([645 1620])

subplot(212);
plot(airsL3.fchanx,nanmean(obs'-obs'),airsL3.fchanx,nanmean(obs'-airsL3.airsL3_spectral_rates'),...
     airsL3.fchanx,nanmean(obs'-era5.era5_spectral_rates'),airsL3.fchanx,nanmean(obs'-cmip6.cmip6_spectral_rates'),'linewidth',2);
plotaxis2;
hl = legend('AIRS OBS','AIRS L3','ERA5','CMIP6','location','best','fontsize',10); ylabel('\delta dBT/dt (K/yr)'); xlabel('Wavenumber cm-1')
xlim([645 1620])
ylim([-1 +1]*0.05)
