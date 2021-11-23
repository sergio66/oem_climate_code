umbc_64.trend_RH    = deltaRHlat';
umbc_64.trend_gas_1 = fracWVlat';
umbc_64.trend_gas_3 = fracO3lat';
umbc_64.trend_ptemp = deltaTlat';
umbc_64.trend_stemp = results(:,6)';

cmip6_64.trend_RH = squeeze(nanmean(reshape(cmip6.trend_RH,100,72,64),2));
cmip6_64.trend_gas_1 = squeeze(nanmean(reshape(cmip6.trend_gas_1,100,72,64),2));
cmip6_64.trend_gas_3 = squeeze(nanmean(reshape(cmip6.trend_gas_3,100,72,64),2));
cmip6_64.trend_ptemp = squeeze(nanmean(reshape(cmip6.trend_ptemp,100,72,64),2));
cmip6_64.trend_stemp = squeeze(nanmean(reshape(cmip6.trend_stemp,72,64),1));

era5_64.trend_RH = squeeze(nanmean(reshape(era5.trend_RH,100,72,64),2));
era5_64.trend_gas_1 = squeeze(nanmean(reshape(era5.trend_gas_1,100,72,64),2));
era5_64.trend_gas_3 = squeeze(nanmean(reshape(era5.trend_gas_3,100,72,64),2));
era5_64.trend_ptemp = squeeze(nanmean(reshape(era5.trend_ptemp,100,72,64),2));
era5_64.trend_stemp = squeeze(nanmean(reshape(era5.trend_stemp,72,64),1));

airsL3_64.trend_RH = squeeze(nanmean(airsL3.thestats64x72.RHrate,1))';
airsL3_64.trend_gas_1 = squeeze(nanmean(airsL3.thestats64x72.waterrate,1))';
airsL3_64.trend_gas_3 = squeeze(nanmean(airsL3.thestats64x72.ozonerate,1))';
airsL3_64.trend_ptemp = squeeze(nanmean(airsL3.thestats64x72.ptemprate,1))';
airsL3_64.trend_stemp = squeeze(nanmean(airsL3.thestats64x72.stemprate,1));
airsL3_64.Qlevs = airsL3.Qlevs;
airsL3_64.Tlevs = airsL3.Tlevs;
