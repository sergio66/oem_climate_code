umbc_spectral_olr = struct;    %% so it has no fields
umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,-1,rlat65,rlon73,'UMBC');

%%%%%%%%%%%%%%%%%%%%%%%%%

era5_spectral_olr = struct;    %% so it has no fields
era5_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,era5.trend_ptemp,era5.trend_gas_1,era5.trend_gas_3,era5_spectral_olr,-1,'ERA5');

aL3trend.stemp = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
aL3trend.ptemp = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp;
aL3trend.gas_1 = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1;
aL3trend.gas_3 = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_3;
airsL3_spectral_olr = struct;    %% so it has no fields
airsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,aL3trend.stemp,aL3trend.ptemp,aL3trend.gas_1,aL3trend.gas_3,airsL3_spectral_olr,-1,'AIRS L3');

x6trend.stemp = cmip6.trend_stemp;
x6trend.ptemp = cmip6.trend_ptemp;
x6trend.gas_1 = cmip6.trend_gas_1;
x6trend.gas_3 = cmip6.trend_gas_3;
cmip6_spectral_olr = struct;    %% so it has no fields
cmip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,x6trend.stemp,x6trend.ptemp,x6trend.gas_1,x6trend.gas_3,cmip6_spectral_olr,-1,'CMIP6');

%%%%%%%%%%%%%%%%%%%%%%%%%

merra2_spectral_olr = struct;    %% so it has no fields
merra2_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,merra2.trend_stemp,merra2.trend_ptemp,merra2.trend_gas_1,merra2.trend_gas_3,merra2_spectral_olr,-1,'MERRA2');

aL3trend.stemp = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.climcapsL3_100_layertrends.stemp;
aL3trend.ptemp = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.climcapsL3_100_layertrends.ptemp;
aL3trend.gas_1 = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.climcapsL3_100_layertrends.gas_1;
aL3trend.gas_3 = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.climcapsL3_100_layertrends.gas_3;
climcapsL3_spectral_olr = struct;    %% so it has no fields
climcapsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,aL3trend.stemp,aL3trend.ptemp,aL3trend.gas_1,aL3trend.gas_3,climcapsL3_spectral_olr,-1,'CLIMCAPS L3');

x6trend.stemp = amip6.trend_stemp;
x6trend.ptemp = amip6.trend_ptemp;
x6trend.gas_1 = amip6.trend_gas_1;
x6trend.gas_3 = amip6.trend_gas_3;
amip6_spectral_olr = struct;    %% so it has no fields
amip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,x6trend.stemp,x6trend.ptemp,x6trend.gas_1,x6trend.gas_3,amip6_spectral_olr,-1,'AMIP6');
