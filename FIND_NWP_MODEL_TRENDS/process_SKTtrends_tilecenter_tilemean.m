if ~exist('trend')
  disp('loading in SKT center-of-tile and mean-over-tile trended time series from skt_trends_tilecenter_tilemean.mat')
  load skt_trends_tilecenter_tilemean.mat
end
%% new to test pushing smoothn

clear newz*
%% see compare_SKT_trends_Day_vs_Night.m
  disp('loading in earlier SKT trends from center-of-tile')
newz6 = load('L3_SKT_TIMESERIES_2002_09_to_2022_08/compare_SKT_trends_Day_vs_Night_fig22.mat');

if ~exist('landfrac')
  addpath /asl/matlib/h4tools
  addpath /asl/matlib/science/
  [h,ha,p,pa] = rtpread('summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.op.rtp');
  [salti,landfrac] =  usgs_deg10_dem(p.rlat,p.rlon);
  p.landfrac = landfrac;
end

load llsmap5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using ORIG from Jan 2025 submission
clc; disp('using ORIG/CNTR data')
newz11 = newz6.newz11;  %% umbc
newz12 = newz6.newz12;  %% v7
newz13 = newz6.newz13;  %% climcaps
newz21 = newz6.newz21;  %% giss
newz22 = newz6.newz22;  %% era5
newz23 = newz6.newz23;  %% merra2
make_table3('ORIG',newz11,newz12,newz13,newz21,newz22,newz23);
disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

%% using CNTR
clc; disp('using NEW/CNTR data')
%newz11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  
newz11 = newz6.newz11;  
  newz12 = (trend.v7.cntrA + trend.v7.cntrD)*0.5;
  newz13 = (trend.climcaps.cntrA + trend.climcaps.cntrD)*0.5;
newz21 = trend.giss.cntr;                                   newz22 = (trend.ERA5.cntrD + trend.ERA5.cntrA)*0.5;                             newz23 = trend.merra2.cntr;
make_table3('CNTR',newz11,newz12,newz13,newz21,newz22,newz23);
disp('ret to continue'); pause

%% using MEAN
clc; disp('using NEW/MEAN data')
%newz11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  
newz11 = newz6.newz11;  
  newz12 = (trend.v7.meanA + trend.v7.meanD)*0.5;
  newz13 = (trend.climcaps.meanA + trend.climcaps.meanD)*0.5;
newz21 = trend.giss.mean;                                   newz22 = (trend.ERA5.meanD + trend.ERA5.meanA)*0.5;                             newz23 = trend.merra2.mean;
make_table3('MEAN',newz11,newz12,newz13,newz21,newz22,newz23);
disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

%newz11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  
newz11 = newz6.newz11;  
  newz12 = (trend.v7.meanA + trend.v7.meanD)*0.5;
  newz13 = (trend.climcaps.meanA + trend.climcaps.meanD)*0.5;
newz21 = trend.giss.mean;                                   newz22 = (trend.ERA5.meanD + trend.ERA5.meanA)*0.5;                             newz23 = trend.merra2.mean;

iFig = 22;
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
  plotoptions.smooth = 1; 
  figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);

iFig = 23;
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
  plotoptions.smooth = -1; 
  figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);
disp('these were using MEAN SMOOTH and MEAN NOSMOOTH ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%newz11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  
newz11 = newz6.newz11;  
  newz12 = (trend.v7.cntrA + trend.v7.cntrD)*0.5;
  newz13 = (trend.climcaps.cntrA + trend.climcaps.cntrD)*0.5;
newz21 = trend.giss.cntr;                                   newz22 = (trend.ERA5.cntrD + trend.ERA5.cntrA)*0.5;                             newz23 = trend.merra2.cntr;

iFig = 24;
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
  plotoptions.smooth = 1; 
  figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);

iFig = 25;
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
  plotoptions.smooth = -1; 
  figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);
disp('thse were using TILE CENTER SMOOTH and TILE CENTER NOSMOOTH ret to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iFig = 20; 
figure(iFig); clf; 

xnewz11 = newz6.newz11 - newz11;
xnewz12 = newz6.newz12 - newz12;
xnewz13 = newz6.newz13 - newz13;
xnewz21 = newz6.newz21 - newz21;
xnewz22 = newz6.newz22 - newz22;
xnewz23 = newz6.newz23 - newz23;
xplotoptions = plotoptions;
xplotoptions.cx = [-1 +1]*0.151/10;
xplotoptions.barstr = 'DELTA OLD-NEW dSKT/dt [K/yr]';
figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(xnewz11,xnewz12,xnewz13,xnewz21,xnewz22,xnewz23,iFig,xplotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_XX_YY_from_X_Y
figure(12); sizefig; ; clf; aslmap(12,rlat65,rlon73,smoothn(reshape(newz6.newz11,72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT NIGHT SMOOTHN');
  caxis([-1 +1]*0.151); colormap(llsmap5);
window = 7;
window = 4;
window = 5;
figure(13); sizefig; ; clf; aslmap(13,rlat65,rlon73,smoothdata(reshape(newz6.newz11,72,64)',"movmean",window),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT NIGHT SMOOTH MEAN');
  caxis([-1 +1]*0.151); colormap(llsmap5);
figure(14); sizefig; ; clf; aslmap(14,rlat65,rlon73,reshape(newz6.newz11,72,64)',[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT NIGHT NO SMOOTH');
  caxis([-1 +1]*0.151); colormap(llsmap5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

newz11 = newz6.newz11;  
  newz12 = (trend.v7.meanA + trend.v7.meanD)*0.5;
  newz13 = (trend.climcaps.meanA + trend.climcaps.meanD)*0.5;
newz21 = trend.giss.mean;                                   newz22 = (trend.ERA5.meanD + trend.ERA5.meanA)*0.5;                             newz23 = trend.merra2.mean;

data7.umbc = newz11;    data7.airsv7 = newz12;    data7.climcaps = newz13;
data7.giss = newz21;    data7.era5   = newz22;    data7.merra2   = newz23;  
data7.comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_SKTtrends_tilecenter_tilemean.m';

  iFig = 123;
    figure(iFig); sizefig; ; clf;
    clear plotoptions
    plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
    plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
    plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
    plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
    plotoptions.barstr = 'dSKT/dt [K/yr]';
    plotoptions.smooth = -1; 
    figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);
    figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(data7.umbc,data7.airsv7,data7.climcaps,data7.giss,data7.era5,data7.merra2,iFig,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data8 = data7;
  data8.comment2 = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/compare_RHsurf_trends_Day_vs_Night.m';
  data8.rlat = rlat;

  data8.landfrac = reshape(landfrac,72,64);
  data8.ocean_umbc = newz11;    data8.ocean_airsv7 = newz12;    data8.ocean_climcaps = newz13;
  data8.ocean_giss = newz21;    data8.ocean_era5   = newz22;    data8.ocean_merra2   = newz23;  
  data8.ocean_umbc(data8.landfrac > eps) = nan;   data8.ocean_airsv7(data8.landfrac > eps) = nan;   data8.ocean_climcaps(data8.landfrac > eps) = nan;
  data8.ocean_giss(data8.landfrac > eps) = nan;   data8.ocean_era5(data8.landfrac > eps)   = nan;   data8.ocean_merra2(data8.landfrac > eps) = nan;

  figure(181); clf
  plot(data8.rlat,nanmean(data8.umbc,1),'k',data8.rlat,nanmean(data8.airsv7,1),'b',data8.rlat,nanmean(data8.climcaps,1),'g',data8.rlat,nanmean(data8.era5,1),'r',data8.rlat,nanmean(data8.merra2),'m',data8.rlat,nanmean(data8.giss),'c','linewidth',2);
  xlim([-1 +1]*90); ylim([-0.06 +0.11]); 
  plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','GISS','location','best','fontsize',8);

  %plot(data8.rlat,nanmean(data8.umbc,1),'k',data8.rlat,nanmean(data8.airsv7,1),'b',data8.rlat,nanmean(data8.climcaps,1),'g',data8.rlat,nanmean(data8.era5,1),'r',data8.rlat,nanmean(data8.merra2),'m','linewidth',2);
  %plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);

  xlabel('Latitude [deg]'); ylabel('dST/dt [K yr^{-1}]');

  %%%%%%%%%%%%%%%%%%%%%%%%%

  figure(182); clf
  plot(data8.rlat,nanmean(data8.ocean_umbc,1),'k',data8.rlat,nanmean(data8.ocean_airsv7,1),'b',data8.rlat,nanmean(data8.ocean_climcaps,1),'g',...
      data8.rlat,nanmean(data8.ocean_era5,1),'r',data8.rlat,nanmean(data8.ocean_merra2),'m',data8.rlat,nanmean(data8.ocean_giss),'c','linewidth',2);
  xlim([-1 +1]*90); ylim([-0.06 +0.11]); 
  plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','GISS','location','best','fontsize',8);

  %plot(data8.rlat,nanmean(data8.ocean_umbc,1),'k',data8.rlat,nanmean(data8.ocean_airsv7,1),'b',data8.rlat,nanmean(data8.ocean_climcaps,1),'g',data8.rlat,nanmean(data8.ocean_era5,1),'r',data8.rlat,nanmean(data8.ocean_merra2),'m','linewidth',2);
  %plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);

  xlabel('Latitude [deg]'); ylabel('dST_{ocean}/dt [K yr^{-1}]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iSaveJGR = -1;
if iSaveJGR > 0
  disp('want to save off <<fig 23>> : fig8new_mean_nosmooth  and also the zonal dST/dt')

  save /home/sergio/MATLABCODE/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig7.mat data7
  save /home/sergio/MATLABCODE/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig8.mat data8
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(22); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_mean_smoothnn');
figure(23); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_mean_nosmooth');
figure(24); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_cntr_smoothnn');
figure(25); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_cntr_nosmooth');

figure(181); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/nosmooth_dST_dt_land_ocean_zonal_allmodels_and_GISS');
figure(182); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/nosmooth_dST_dt_ocean_zonal_allmodels_and_GISS');

figure(1); clf; generic_subtiles_fig_into_surface_matr('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_mean_smoothnn.fig'); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8xnew_mean_smoothnn');
figure(1); clf; generic_subtiles_fig_into_surface_matr('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_mean_nosmooth.fig'); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8xnew_mean_nosmooth');
figure(1); clf; generic_subtiles_fig_into_surface_matr('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_cntr_smoothnn.fig'); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8xnew_cntr_smoothnn');
figure(1); clf; generic_subtiles_fig_into_surface_matr('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8new_cntr_nosmooth.fig'); sergioprintfig('L3_SKT_TIMESERIES_2002_09_to_2022_08/fig8xnew_cntr_nosmooth');
%}
