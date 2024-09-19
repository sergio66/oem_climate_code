 addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/CERES_L3_TRENDS
addpath /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR

disp('look at ../AIRS_gridded_STM_May2021_trendsonlyCLR/show_fluxband_anomaly_timeseries.m')
disp('look at ../AIRS_gridded_STM_May2021_trendsonlyCLR/modis_ceres_era5_anoms.m')
disp('look at ../AIRS_gridded_STM_May2021_trendsonlyCLR/do_ceres_vs_umbc_anomalies_zonalavg.m')
disp('look at ../FIND_NWP_MODEL_TRENDS/driver_compare_read64x72flux72_save64zonalflux_ERA5_dRH_anom.m')

disp('to ccompare ERA5,AIRSL3 monthly OLR against that from SARTA, run ../AIRS_gridded_STM_May2021_trendsonlyCLR/compare_OLR_from_monthly_2_OLR_from_sarta.m')
disp('to ccompare ERA5,AIRSL3 monthly OLR against that from SARTA, run ../AIRS_gridded_STM_May2021_trendsonlyCLR/compare_OLR_from_monthly_2_OLR_from_sarta.m')
disp('to ccompare ERA5,AIRSL3 monthly OLR against that from SARTA, run ../AIRS_gridded_STM_May2021_trendsonlyCLR/compare_OLR_from_monthly_2_OLR_from_sarta.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('radanomD_fatbins')
  radanomD_fatbins = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/anomalyD_zonalavg_globalavg_and_28_averages_timeseries_Q03_numyears_22_iNumAnomTimeSteps_500.mat');
  radanomD_fatbins.latbins = [-90 -75 -60 [-55:5:+55] +60 +75 +90];
  radanomD_fatbins.globalavg = radanomD_fatbins.btavgAnomFinal(:,1:500);
  if ~exist('f_ind')
    f_ind - instr_chans2645;
  end
  i0667 = find(f_ind >= 667.3,1);
  i0723 = find(f_ind >= 723,1);
  i1231 = find(f_ind >= 1231,1);
  i1419 = find(f_ind >= 1419,1);

  %junk = btavg_cos;
  junk = radanomD_fatbins.globalavg;
  daysSince2002A = change2days(yyD,mmD,ddD,2002);

  %% see stand_alone_make_globalavg_and_N_average_anomalies_zonalavg.m
  iboo = 0;
  iboo = iboo + 1; ijunk = find(f_ind >= 0667.019,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(f_ind >= 0681.457,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(f_ind >= 0704.436,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(f_ind >= 0723.029,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(f_ind >= 0801.099,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(f_ind >= 0900.000,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(f_ind >= 0960.000,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(f_ind >= 1040.083,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(f_ind >= 1231.300,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(f_ind >= 1419.100,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(f_ind >= 1519.070,1); iaChan(iboo) = ijunk;
  iboo = iboo + 1; ijunk = find(f_ind >= 1598.490,1); iaChan(iboo) = ijunk;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iYS = 2002; iYE = 2024;
yy = []; mm = []; dd = [];
for ii = iYS : iYE
  clear yyx mmx ddx
  if ii == iYS
    inum = 4;
    yyx(1:inum) = ii;
    mmx = 9:12;
    ddx = ones(size(mmx)) * 15;
%  elseif ii == 2021
%    %inum = 7;
%    %yyx(1:inum) = ii;
%    %mmx = 1 : 7;
%    inum = 8;
%    yyx(1:inum) = ii;
%    mmx = 1 : 8;
%    ddx = ones(size(mmx)) * 15;
  elseif ii == iYE
    inum = 8;
    yyx(1:inum) = ii;
    mmx = 1 : 8;
    ddx = ones(size(mmx)) * 15;
  else
    inum = 12;
    yyx(1:inum) = ii;
    mmx = 1:12;
    ddx = ones(size(mmx)) * 15;
  end
  fprintf(1,'%4i %2i \n',[ii inum])
  yy = [yy yyx];
  mm = [mm mmx];
  dd = [dd ddx];
end

yymm = yy + (mm-1)/12;
daysSince2002 = change2days(yy,mm,dd,2002);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('ceres_trend')
  %% see modis_ceres_era5_anoms.m
  ceres_trend = quick_ceres_flux_cloud_anomaly_look();
  modis_cloud = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/MODIS_L3_AEROSOL_TRENDS/modis_L3_cloud_trends_2002_0009_2024_0006.mat');
    aod_anom = reshape(modis_cloud.anom_cldfrac,4608,262);
    aod_anom_cosavg = sum(aod_anom.*(cos(reshape(YY,4608,1)*pi/180)*ones(1,262)))./sum(cos(reshape(YY,4608,1)*pi/180)*ones(1,262));
    deepblue_anom = reshape(modis_cloud.anom_cldtop,4608,262);
    deepblue_anom_cosavg = sum(deepblue_anom.*(cos(reshape(YY,4608,1)*pi/180)*ones(1,262)))./sum(cos(reshape(YY,4608,1)*pi/180)*ones(1,262));

  ohc = ocean_heat_content();
  solar = annual_solar_irradiance_data();
  wind_era5 = wind_speed_changes();
  wind_merra2 = wind_speed_changes();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ht = 2022 + 15/30/12;  %% hunga tonga
plot_clouds_vs_AIRSchans_anomalies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('airsL3_direct_flux')
  airsL3_direct_flux = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2024_22yr_desc_oopOLRanom.mat');

  wah = squeeze(airsL3_direct_flux.thestats64x72_other.olranom);
  wah = squeeze(nanmean(wah,1));
  coswah = cos(rlat*pi/180)*ones(1,262);
  wahmean_airsL3_direct = nansum(wah.*coswah,1)./nansum(coswah,1);
  plot(yymm(1:262),smooth(wahmean_airsL3_direct,7)); plotaxis2;

  wah = squeeze(airsL3_direct_flux.thestats64x72_other.clr_olranom);
  wah = squeeze(nanmean(wah,1));
  coswah = cos(rlat*pi/180)*ones(1,262);
  wahmean_airsL3_direct_clr = nansum(wah.*coswah,1)./nansum(coswah,1);
  plot(yymm(1:262),smooth(wahmean_airsL3_direct,7)); plotaxis2;

  %%%%%

  airsL3_sarta_flux = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_all64_anomflux_14RRTMbands.mat');
  wah = squeeze(airsL3_sarta_flux.anomflux(:,1,:));
  coswah = cos(rlat*pi/180)*ones(1,262);
  wahmean_airsL3 = nansum(wah.*coswah,1)./nansum(coswah,1);
  plot(smooth(squeeze(nanmean(airsL3_sarta_flux.anomflux(:,1,:),1)),7)*3.5788); plotaxis2;
  plot(1:262,smooth(wahmean_airsL3,7)*3.5788); plotaxis2;

 %%%%%%%%%%%%%%%%%%%%%%%%%

  era5_direct_flux = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_2024_08_trends_desc_surf_anomaly.mat');

  wah = era5_direct_flux.anom_olr;
  wah = reshape(wah,72,64,264);
  wah = squeeze(nanmean(wah,1));
  coswah = cos(rlat*pi/180)*ones(1,264);
  wahmean_era5_direct = nansum(wah.*coswah,1)./nansum(coswah,1);

  wah = era5_direct_flux.anom_olr_clr;
  wah = reshape(wah,72,64,264);
  wah = squeeze(nanmean(wah,1));
  coswah = cos(rlat*pi/180)*ones(1,264);
  wahmean_era5_direct_clr = nansum(wah.*coswah,1)./nansum(coswah,1);

  wah = era5_direct_flux.anom_stemp;
  wah = reshape(wah,72,64,264);
  wah = squeeze(nanmean(wah,1));
  coswah = cos(rlat*pi/180)*ones(1,264);
  wahmean_era5_direct_stemp = nansum(wah.*coswah,1)./nansum(coswah,1);

  plot(yymm(1:264),smooth(wahmean_era5_direct,7),yymm(1:264),smooth(wahmean_era5_direct_stemp,7)); plotaxis2;

  %%%%%%%%%%%%%%%%%%%%%%%%%

  figure(5)
  plot(yymm(1:262),smooth(wahmean_airsL3_direct_clr,7),'b',yymm(1:262),smooth(wahmean_airsL3_direct,7),'c',...
       yymm(1:264),smooth(wahmean_era5_direct_clr,7),'r',yymm(1:264),smooth(wahmean_era5_direct,7),'m',yymm(1:264),smooth(wahmean_era5_direct_stemp,7),'linewidth',2); plotaxis2;
    hl = legend('AIRS L3 CLR','AIRS L3 ALL','ERA5 CLR','ERA5 ALL','ERA5 STEMP','location','best','fontsize',8);
  set(gca,'fontsize',10);
  ylim([-1 +1]*1.5)
  xlim([2002 2025])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('fluxanomA')
  load /asl/s1/sergio/JUNK/anomaly_iQAX_3_zonalavg_ALL_Q03_numyears_22.00_iNumAnomTimeSteps_500_A.mat  % asc, Q90
  load /asl/s1/sergio/JUNK/anomaly_iQAX_3_zonalavg_ALL_Q03_numyears_22.00_iNumAnomTimeSteps_500_D.mat  % desc, Q90
  load /asl/s1/sergio/JUNK/anomaly_iQAX_4_zonalavg_ALL_Q01_numyears_22.00_iNumAnomTimeSteps_500_D.mat  % desc, allsky
  yymmL1C = yyA + (mmA-1)/12 + (ddA-1)/30/12;
  do_XX_YY_from_X_Y
end

disp('ret to continue to OLR anomalies'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_XX_YY_from_X_Y
coslat    = cos(rlat*pi/180)* ones(1,12*22);
coslatL1C = cos(rlat*pi/180)* ones(1,500);

smn = 3;
smn = 7;

a0 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/era5_64x72_Sept2002_Aug2024_22yr_desc_btanom_all64_anomflux_14RRTMbands.mat');
a1 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/era5_64x72_Sept2002_Aug2024_22yr_desc_btanom_all64_idRH_1_anomflux_14RRTMbands.mat');
a2 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/era5_64x72_Sept2002_Aug2024_22yr_desc_btanom_all64_idRH_2_anomflux_14RRTMbands.mat');
a3 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/era5_64x72_Sept2002_Aug2024_22yr_desc_btanom_all64_idRH_3_anomflux_14RRTMbands.mat');
a4 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/era5_64x72_Sept2002_Aug2024_22yr_desc_btanom_all64_idRH_4_anomflux_14RRTMbands.mat');
a5 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/era5_64x72_Sept2002_Aug2024_22yr_desc_btanom_all64_idRH_5_anomflux_14RRTMbands.mat');
a6 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/era5_64x72_Sept2002_Aug2024_22yr_desc_btanom_all64_idRH_6_anomflux_14RRTMbands.mat');

era5_direct_flux = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_2024_08_trends_desc_surf_anomalyat');
era5.olr     = reshape(era5_direct_flux.anom_olr,72,64,12*22);     era5.olr64 = squeeze(nanmean(era5.olr,1));         era5.olr = squeeze(nanmean(era5.olr,1));          era5.olr = sum(coslat.*era5.olr)./sum(coslat);
era5.olr_clr = reshape(era5_direct_flux.anom_olr_clr,72,64,12*22); era5.olr_clr64 = squeeze(nanmean(era5.olr_clr,1)); era5.olr_clr = squeeze(nanmean(era5.olr_clr,1));  era5.olr_clr = sum(coslat.*era5.olr_clr)./sum(coslat);

figure(1)
plot(yymm,smooth(era5.olr_clr,smn),'b',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw_clr,smn),'c',yymm,smooth(era5.olr,smn),'r',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw,smn),'m','linewidth',2);
plotaxis2; title('CosAvg LW/MW IR Flux Anomaly W/m2'); ylabel('Flux Anomaly [W/m2]')
legend('ERA5 CLR','CERES LW CLR','ERA5 CLD', 'CERES LW CLD','location','best','fontsize',8);
ylim([-1 +1]*1.5)
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])
set(gca,'fontsize',10)

figure(2)
fluxA = sum(coslatL1C.*fluxanomA(:,:,1),1)./sum(coslatL1C,1) * 3.5585; 
fluxD = sum(coslatL1C.*fluxanomD(:,:,1),1)./sum(coslatL1C,1) * 3.5585; 
plot(yymm,smooth(era5.olr_clr,smn),'b',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw_clr,smn),'c',yymm,smooth(era5.olr,smn),'r',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw,smn),'m',...
     yymmL1C,smooth(squeeze(nanmean(coslatL1C.*fluxanomD(:,:,1),1)),smn),'k',yymmL1C,smooth(squeeze(nanmean(coslatL1C.*fluxanomD(:,:,1),1)),smn),'g','linewidth',2);
plot(yymm,smooth(era5.olr_clr,smn),'b',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw_clr,smn),'c',yymm,smooth(era5.olr,smn),'r',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw,smn),'m',...
     yymmL1C,smooth(fluxA,smn*2+1),'k',yymmL1C,smooth(fluxD,smn*2+1),'g','linewidth',2);
plotaxis2; title('CosAvg LW/MW IR Flux Anomaly W/m2'); ylabel('Flux Anomaly [W/m2]')
legend('ERA5 CLR','CERES LW CLR','ERA5 CLD', 'CERES LW CLD','AIRS L1C Clrsky','AIRS L1C AllSky','location','best','fontsize',8);
ylim([-1 +1]*1.5)
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])
set(gca,'fontsize',10)

figure(3)
plot(yymm,smooth(era5.olr_clr,smn),'b',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw_clr,smn),'k',yymm(1:262),smooth(wahmean_airsL3_direct_clr,smn),'g',yymmL1C,smooth(fluxA,smn*2+1),'r.-','linewidth',2);
plotaxis2; ylabel('Flux Anomaly [W/m2]'); %title('CLRSKY CosAvg LW/MW IR Flux Anomaly W/m2'); 
legend('ERA5 CLR','CERES LW CLR','AIRS L3 CLR','AIRS L1C ClrSky','location','best','fontsize',8);
ylim([-1 +1]*1.5)
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])
set(gca,'fontsize',10)

figure(4)
plot(yymm,smooth(era5.olr,smn),'b',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw,smn),'k',yymm(1:262),smooth(wahmean_airsL3_direct,smn),'g',yymmL1C,smooth(fluxD,smn*2+1),'r.-','linewidth',2);
plotaxis2;  ylabel('Flux Anomaly [W/m2]'); %title('ALLSKY CosAvg LW/MW IR Flux Anomaly W/m2');
legend('ERA5 CLD', 'CERES LW CLD','AIRS L3 CLD','AIRS L1C Allsky','location','best','fontsize',8);
ylim([-1 +1]*1.5)
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])
set(gca,'fontsize',10)

%{
figure(3); sergioprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/clrskyOLR_era5_airsL3_ceres_umbc');
figure(4); sergioprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/allskyOLR_era5_airsL3_ceres_umbc');
%}

disp('ret to continue to OLR trends'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now do fluxes
for rr = 1 : 64
  [B, stats, err] = Math_tsfit_lin_robust(daysSince2002,era5.olr64(rr,:),1);     era5.trendolr64(rr) = B(2);     era5.trendolr64_unc(rr) = stats.se(2);
  [B, stats, err] = Math_tsfit_lin_robust(daysSince2002,era5.olr_clr64(rr,:),1); era5.trendolr_clr64(rr) = B(2); era5.trendolr_clr64_unc(rr) = stats.se(2);
    
  [B, stats, err] = Math_tsfit_lin_robust((yymmL1C-2002)*365.25,squeeze(fluxanomD(rr,:,1))*3.5578,1); umbc.trendolr64(rr) = B(2);     umbc.trendolr64_unc(rr) = stats.se(2);
  [B, stats, err] = Math_tsfit_lin_robust((yymmL1C-2002)*365.25,squeeze(fluxanomA(rr,:,1))*3.5578,1); umbc.trendolr_clr64(rr) = B(2); umbc.trendolr_clr64_unc(rr) = stats.se(2);
end

figure(5); 
  plot(rlat,era5.trendolr_clr64,'b',rlat,nanmean(reshape(ceres_trend.trend_toa_lw_clr_t_4608,72,64),1),'k',...
       rlat,nanmean(airsL3_direct_flux.thestats64x72_other.clrolrrate,1),'g',rlat,umbc.trendolr_clr64,'r','linewidth',2)
  plotaxis2;  ylabel('Flux Trend [W/m2/Yr]'); %title('ALLSKY CosAvg LW/MW IR Flux Anomaly W/m2');
  legend('ERA5 CLR', 'CERES LW CLR','AIRS L3 CLR','AIRS L1C Clrsky','location','best','fontsize',8);
  xlim([-1 +1]*90); xlabel('Latitude')
  set(gca,'fontsize',10)

figure(6); 
  plot(rlat,era5.trendolr64,'b',rlat,nanmean(reshape(ceres_trend.trend_toa_lw_all_4608,72,64),1),'k',...
       rlat,nanmean(airsL3_direct_flux.thestats64x72_other.olrrate,1),'g',rlat,umbc.trendolr64,'r','linewidth',2)
  plotaxis2;  ylabel('Flux Trend [W/m2/Yr]'); %title('ALLSKY CosAvg LW/MW IR Flux Anomaly W/m2');
  legend('ERA5 CLD', 'CERES LW CLD','AIRS L3 CLD','AIRS L1C Allsky','location','best','fontsize',8);
  xlim([-1 +1]*90); xlabel('Latitude')
  set(gca,'fontsize',10)


%{
figure(5); sergioprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/clrskyOLR_era5_airsL3_ceres_umbc_trends');
figure(6); sergqgrep -in '/airsL3_direct_vs_sarta' *.mioprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/allskyOLR_era5_airsL3_ceres_umbc_trends');
%}

disp('ret to continue to OLR anomalies, and zoom into 2020-2025'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tropics = find(abs(rlat) < 30);
RRTM_bands0 = [10 250 500 630 700 820 980 1080 1180 1390 1480 1800 2080 2250 2380 2600 3000];
RRTM_bands = RRTM_bands0(4:end);

compare_read64x72flux72_save64zonalflux_ERA5_dRH_anom_basic
disp('ret to continue to RRTM band OLR anomalies, and zoom into 2020-2025'); pause

compare_read64x72flux72_save64zonalflux_ERA5_dRH_anom_loopbasic
disp('ret to continue to RRTM band OLR anomalies, and zoom into 2020-2025, smarter'); pause

compare_read64x72flux72_save64zonalflux_ERA5_dRH_anom_loopadv
disp('ret to continue to AIRS L3 vs OBS')

compare_read64x72flux72_save64zonalflux_ERA5_AIRL3_anom_loopadv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
