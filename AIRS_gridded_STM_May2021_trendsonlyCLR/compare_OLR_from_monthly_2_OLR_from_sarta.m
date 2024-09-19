figure(21); clf

if ~exist('rlat')
  do_XX_YY_from_X_Y
end

ixjunk = 0;
for yyjunk = 2002 : 2024
  mmSjunk = 1; mmEjunk = 12;
  if yyjunk == 2002
    mmSjunk = 9;
  elseif yyjunk == 2002
    mmEjunk = 8;
  end
  for mmjunk = mmSjunk : mmEjunk
    ixjunk = ixjunk + 1;
    yy_month(ixjunk) = yyjunk;
    mm_month(ixjunk) = mmjunk;
  end
end
yymm_monthly = yy_month + (mm_month-1)/12;
%% factor of 3.5788 from stand_alone_make_globalavg_and_N_average_anomalies_zonalavg.m

%%%%%%%%%%%%%%%%%%%%%%%%%
%% see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_computeERA5_monthly_trends_desc_or_asc.m
era5_direct_flux = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_2024_08_trends_desc_surf_anomaly.mat');
  wah = era5_direct_flux.anom_olr_clr;
  wah = reshape(wah,72,64,264);
  wah = squeeze(nanmean(wah,1));
  coswah = cos(rlat*pi/180)*ones(1,264);
  wahmean_era5_direct = nansum(wah.*coswah,1)./nansum(coswah,1);

  wah = era5_direct_flux.anom_stemp;
  wah = reshape(wah,72,64,264);
  wah = squeeze(nanmean(wah,1));
  coswah = cos(rlat*pi/180)*ones(1,264);
  wahmean_era5_direct_stemp = nansum(wah.*coswah,1)./nansum(coswah,1);

  plot(1:264,smooth(wahmean_era5_direct,7),1:264,smooth(wahmean_era5_direct_stemp,7)); plotaxis2;

era5_sarta_flux = load('SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/era5_64x72_Sept2002_Aug2024_22yr_desc_btanom_all64_anomflux_14RRTMbands.mat');
  wah = squeeze(era5_sarta_flux.anomflux(:,1,:));
  coswah = cos(rlat*pi/180)*ones(1,264);
  wahmean_era5 = nansum(wah.*coswah,1)./nansum(coswah,1);
  plot(smooth(squeeze(nanmean(era5_sarta_flux.anomflux(:,1,:),1)),7)*3.5788); plotaxis2;
  plot(1:264,smooth(wahmean_era5,7)*3.5788); plotaxis2;

%%%%%%%%%%%%%%%%%%%%%%%%%
%% this should really be in /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2024_22yr_desc.mat
%% see driver_compute_AIRSL3_trends_desc_or_ascNOQuestioN.m
%% but I forgot to do it, so look at driver_compute_AIRSL3_trends_desc_or_ascNOQuestioN_oopsOLRanom.m

wah = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2024_22yr_desc.mat');
wah.Qlevs = [1000         925         850         700         600         500         400         300         250         200         150         100];
%pcolor(meanvaluebin(rlat),wah.Qlevs,squeeze(nanmean(wah.thestats64x72.waterrate,1))'); shading interp; colormap(usa2); caxis([-1 +1]*0.015); set(gca,'ydir','reverse');
pcolor(rlat',wah.Qlevs,squeeze(nanmean(wah.thestats64x72.waterrate,1))'); shading interp; colormap(usa2); caxis([-1 +1]*0.015); set(gca,'ydir','reverse'); colorbar
title('AIRS L3 fracWV trends /yr'); xlabel('Latitude'); ylabel('Pressure')

airsL3_direct_flux = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2024_22yr_desc_oopOLRanom.mat');
  wah = squeeze(airsL3_direct_flux.thestats64x72_other.clr_olranom);
  wah = squeeze(nanmean(wah,1));
  coswah = cos(rlat*pi/180)*ones(1,262);
  wahmean_airsL3_direct = nansum(wah.*coswah
,1)./nansum(coswah,1);
  plot(1:262,smooth(wahmean_airsL3_direct,7)); plotaxis2;

%% see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_read64x72flux72_save64zonalflux_ERA5_AIRSL3.m
airsL3_sarta_flux = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_all64_anomflux_14RRTMbands.mat');
  wah = squeeze(airsL3_sarta_flux.anomflux(:,1,:));
  coswah = cos(rlat*pi/180)*ones(1,262);
  wahmean_airsL3 = nansum(wah.*coswah,1)./nansum(coswah,1);
  plot(smooth(squeeze(nanmean(airsL3_sarta_flux.anomflux(:,1,:),1)),7)*3.5788); plotaxis2;
  plot(1:262,smooth(wahmean_airsL3,7)*3.5788); plotaxis2;

%%%%%%%%%%%%%%%%%%%%%%%%%

plot(yymm_monthly(1:262),smooth(squeeze(nanmean(airsL3_sarta_flux.anomflux(:,1,:),1)),7)*3.5788,...
       yymm_monthly(1:264),smooth(squeeze(nanmean(era5_sarta_flux.anomflux(:,1,:),1)),7)*3.5788,'linewidth',2); 
plot(yymm_monthly(1:262),smooth(wahmean_airsL3,7)*3.5788,'c',...
       yymm_monthly(1:264),smooth(wahmean_era5,7)*3.5788,'m',...
       yymm_monthly(1:262),smooth(wahmean_airsL3_direct,7),'b',...
       yymm_monthly(1:264),smooth(wahmean_era5_direct,7),'r','linewidth',2); 
plotaxis2;
hl = legend('AIRS L3 --> SARTA','ERA5 --> SARTA','AIRS L3 direct','ERA5 direct','location','best');
title('CosAvg MONTHLY fields --> SARTA or direct Flux Anomaly W/m2')
set(gca,'fontsize',10)
xlabel('Time'); ylabel('Anomaly Flux [W/m2]'); axis([2002.75 2024.5 -1.5 +1.5]);

%{
plot(yymm_monthly(1:262),smooth(wahmean_airsL3,7)*3.5788,'b',...
       yymm_monthly(1:262),smooth(wahmean_airsL3_direct,7),'r','linewidth',2)
plotaxis2;
hl = legend('AIRS L3 --> SARTA','AIRS L3 direct','location','best');
%title('CosAvg MONTHLY fields --> SARTA or direct Flux Anomaly W/m2')
set(gca,'fontsize',10)
xlabel('Time'); ylabel('Anomaly Flux [W/m2]'); axis([2002.75 2024.5 -1.5 +1.5]);
sergioprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/airsL3_direct_vs_sarta');

plot(yymm_monthly(1:264),smooth(wahmean_era5,7)*3.5788,'b',...
     yymm_monthly(1:264),smooth(wahmean_era5_direct,7),'r','linewidth',2); 
plotaxis2; axis([2002.75 2024.5 -1.5 +1.5]);
hl = legend('ERA5 --> SARTA','ERA5 direct','location','best');
%title('CosAvg MONTHLY fields --> SARTA or direct Flux Anomaly W/m2')
set(gca,'fontsize',10)
xlabel('Time'); ylabel('Anomaly Flux [W/m2]'); 
sergioprint('/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/era5_direct_vs_sarta');
%}
