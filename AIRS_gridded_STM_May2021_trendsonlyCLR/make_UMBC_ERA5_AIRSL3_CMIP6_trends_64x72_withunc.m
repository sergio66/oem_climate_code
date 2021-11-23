umbc_64x72.trend_RH    = deltaRH;
umbc_64x72.trend_gas_1 = fracWV;
umbc_64x72.trend_gas_3 = fracO3;
umbc_64x72.trend_ptemp = deltaT;
umbc_64x72.trend_stemp = results(:,6)';
%
umbc_64x72.trend_RH_unc    = deltaRHunc;
umbc_64x72.trend_gas_1_unc = fracWVunc;
umbc_64x72.trend_gas_3_unc = fracO3unc;
umbc_64x72.trend_ptemp_unc = deltaTunc;
umbc_64x72.trend_stemp_unc = resultsunc(:,6)';

%%%%%%%%%%%%%%%%%%%%%%%%%

cmip6_64x72.trend_RH    = cmip6.trend_RH;
cmip6_64x72.trend_gas_1 = cmip6.trend_gas_1;;
cmip6_64x72.trend_gas_3 = cmip6.trend_gas_3;;
cmip6_64x72.trend_ptemp = cmip6.trend_ptemp;;
cmip6_64x72.trend_stemp = cmip6.trend_stemp;
%
cmip6_64x72.trend_RH_unc    = cmip6.trend_RH_err;;
cmip6_64x72.trend_gas_1_unc = cmip6.trend_gas_1_err;;
cmip6_64x72.trend_gas_3_unc = cmip6.trend_gas_3_err;;
cmip6_64x72.trend_ptemp_unc = cmip6.trend_ptemp_err;;
cmip6_64x72.trend_stemp_unc = cmip6.trend_stemp_err;

%%%%%%%%%%%%%%%%%%%%%%%%%

era5_64x72.trend_RH    = era5.trend_RH;;
era5_64x72.trend_gas_1 = era5.trend_gas_1;;
era5_64x72.trend_gas_3 = era5.trend_gas_3;;
era5_64x72.trend_ptemp = era5.trend_ptemp;;
era5_64x72.trend_stemp = era5.trend_stemp;
%
era5_64x72.trend_RH_unc    = era5.trend_RH_err;;
era5_64x72.trend_gas_1_unc = era5.trend_gas_1_err;;
era5_64x72.trend_gas_3_unc = era5.trend_gas_3_err;;
era5_64x72.trend_ptemp_unc = era5.trend_ptemp_err;;
era5_64x72.trend_stemp_unc = era5.trend_stemp_err;

%%%%%%%%%%%%%%%%%%%%%%%%%

airsL3_64x72.trend_RH    = reshape(airsL3.thestats64x72.RHrate,64*72,12)';
airsL3_64x72.trend_gas_1 = reshape(airsL3.thestats64x72.waterrate,64*72,12)';
airsL3_64x72.trend_gas_3 = reshape(airsL3.thestats64x72.ozonerate,64*72,24)';
airsL3_64x72.trend_ptemp = reshape(airsL3.thestats64x72.ptemprate,64*72,24)';
airsL3_64x72.trend_stemp = reshape(airsL3.thestats64x72.stemprate,64*72,1)';
%
airsL3_64x72.trend_RH_unc    = reshape(airsL3.thestats64x72.RHratestd,64*72,12)';
airsL3_64x72.trend_gas_1_unc = reshape(airsL3.thestats64x72.waterratestd,64*72,12)';
airsL3_64x72.trend_gas_3_unc = reshape(airsL3.thestats64x72.ozoneratestd,64*72,24)';
airsL3_64x72.trend_ptemp_unc = reshape(airsL3.thestats64x72.ptempratestd,64*72,24)';
airsL3_64x72.trend_stemp_unc = reshape(airsL3.thestats64x72.stempratestd,64*72,1)';
%
airsL3_64x72.Qlevs = airsL3.Qlevs;
airsL3_64x72.Tlevs = airsL3.Tlevs;

