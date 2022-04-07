addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS

%% see driver_gather_spectralrates_AIRSL3_NWP_XMIP6.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airsL3_geo_rates         = getdata_AIRSL3vsCLIMCAPSL3(1);
airsCLIMCAPSL3_geo_rates = getdata_AIRSL3vsCLIMCAPSL3(-1);
merra2_geo_rates         = getdata_NWP(2);
era5_geo_rates           = getdata_NWP(5);
cmip6_geo_rates          = getdata_XMIP6(-1);
amip6_geo_rates          = getdata_XMIP6(+1);

airsL3_geo_rates         = airsL3_geo_rates.thestats64x72.stemprate;         
  airsL3_geo_rates = airsL3_geo_rates(:)';
airsCLIMCAPSL3_geo_rates = airsCLIMCAPSL3_geo_rates.thestats64x72.stemprate;
  airsCLIMCAPSL3_geo_rates = airsCLIMCAPSL3_geo_rates(:)';
merra2_geo_rates         = merra2_geo_rates.trend_stemp;
era5_geo_rates           = era5_geo_rates.trend_stemp;
amip6_geo_rates          = amip6_geo_rates.trend_stemp;
cmip6_geo_rates          = cmip6_geo_rates.trend_stemp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load('llsmap5');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear plotoptions;
plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'ERA5';    plotoptions.str12 = 'MERRA2';    
plotoptions.str21 = 'AIRS L3'; plotoptions.str22 = 'CLIMCAPS2'; 
plotoptions.str31 = 'UMBC';    plotoptions.str32 = 'CMIP6';     
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
figure(5); clf; aslmap_6tiledlayout(era5_geo_rates,merra2_geo_rates,airsL3_geo_rates,airsCLIMCAPSL3_geo_rates,umbcst,cmip6_geo_rates,5,plotoptions);
%save FIGS/Figs_JPL_Apr2022/strow_jpl_Apr2022_sstrates era5_geo_rates merra2_geo_rates airsL3_geo_rates airsCLIMCAPSL3_geo_rates cmip6_geo_rates umbcst

