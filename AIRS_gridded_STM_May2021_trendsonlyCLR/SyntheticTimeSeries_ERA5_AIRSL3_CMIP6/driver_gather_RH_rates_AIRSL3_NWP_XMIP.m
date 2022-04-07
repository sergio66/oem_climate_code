addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS

%% see driver_gather_spectralrates_AIRSL3_NWP_XMIP6.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airsL3         = getdata_AIRSL3vsCLIMCAPSL3(1);
airsCLIMCAPSL3 = getdata_AIRSL3vsCLIMCAPSL3(-1);
merra2_geo_rates         = getdata_NWP(2);
era5_geo_rates           = getdata_NWP(5);
cmip6_geo_rates          = getdata_XMIP6(-1);
amip6_geo_rates          = getdata_XMIP6(+1);

merra2_geo_rates         = merra2_geo_rates.trend_RH;
era5_geo_rates           = era5_geo_rates.trend_RH;
amip6_geo_rates          = amip6_geo_rates.trend_RH;
cmip6_geo_rates          = cmip6_geo_rates.trend_RH;

savename = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithERA5_uncX3.mat';
savename = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX3_50fatlayers_AIRSL3_ERA5_CMIP6_feedback.mat';

plays = load(savename,'plays'); plays = plays.plays;
pavg  = load(savename,'pavg');  pavg = pavg.pavg;
junk    = load(savename,'resultsWV');
  resultsWV = junk.resultsWV;
junk    = load(savename,'resultsT');
  resultsT = junk.resultsT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end

  ind = (ii-1)*72 + (1:72);

  for jjj = 1 : 72

    %%%%%%%%%%%%%%%%%%%%%%%%%
    Qlevs = airsL3.Qlevs;
    Tlevs = airsL3.Tlevs;
    %junkrate = airsL3.thestats64x72.stemprate(:,ii)';   
    %  airsL3_100_layertrends.stemp(ind) = junkrate;
    %xjunkrate = airsL3.thestats64x72.ptemprate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
    %  for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
    %  airsL3_100_layertrends.ptemp(:,ind) = junkrate;
    xjunkrate = airsL3.thestats64x72.RHrate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      airsL3_100_layertrends.RH(:,ind) = junkrate;
    %xjunkrate = airsL3.thestats64x72.waterrate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
    %  for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Qlevs) | plays >= max(Qlevs)); junkrate(:,jjj) = 0;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
    %  airsL3_100_layertrends.gas_1(:,ind) = junkrate;
    %xjunkrate = airsL3.thestats64x72.ozonerate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
    %  for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0;       
    %  airsL3_100_layertrends.gas_3(:,ind) = junkrate;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    Qlevs = airsCLIMCAPSL3.Qlevs;
    Tlevs = airsCLIMCAPSL3.Tlevs;
    %junkrate = airsCLIMCAPSL3.thestats64x72.stemprate(:,ii)';   
    %  airsCLIMCAPSL3_100_layertrends.stemp(ind) = junkrate;
    %xjunkrate = airsCLIMCAPSL3.thestats64x72.ptemprate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
    %  for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
    %  airsCLIMCAPSL3_100_layertrends.ptemp(:,ind) = junkrate;
    xjunkrate = airsCLIMCAPSL3.thestats64x72.RHrate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      airsCLIMCAPSL3_100_layertrends.RH(:,ind) = junkrate;
    %xjunkrate = airsCLIMCAPSL3.thestats64x72.waterrate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
    %  for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Qlevs) | plays >= max(Qlevs)); junkrate(:,jjj) = 0;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
    %  airsCLIMCAPSL3_100_layertrends.gas_1(:,ind) = junkrate;
    %xjunkrate = airsCLIMCAPSL3.thestats64x72.ozonerate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
    %  for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0;       
    %  airsCLIMCAPSL3_100_layertrends.gas_3(:,ind) = junkrate;
  
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    Qlevs = pavg;
    Tlevs = pavg;
    %  umbc_20_layertrends.stemp(ind) = results(ind,6);
    xjunkrate = resultsT(ind,:); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      umbc_20_layertrends.ptemp(:,ind) = junkrate;
    %xjunkrate = resultsRH(ind,:); clear junkrate
    %  for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
    %  umbc_20_layertrends.RH(:,ind) = junkrate;
    xjunkrate = resultsWV(ind,:); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Qlevs) | plays >= max(Qlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      umbc_20_layertrends.gas_1(:,ind) = junkrate;
    %xjunkrate = resultsO3(ind,:); clear junkrate
    %  junkrate = xjunkrate;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
    %  umbc_20_layertrends.gas_3(:,ind) = junkrate;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now need to get in mean profiles
%% see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/plot_driver_gather_gridded_retrieval_results.m
[hMean17years,ha,pMean17years,pa]     = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
iLoad = 1;
  iDorA = 1;
  if iDorA > 0
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/DESC/era_tile_center_timestep_' num2str(iLoad,'%03d') '.mat'];
  else
    fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA/Tile_Center/ASC/era_tile_center_timestep_' num2str(iLoad,'%03d') '.mat'];
  end
  era_prof = load(fin);
  hTimeStep1 = era_prof.hnew_op;
  pTimeStep1 = era_prof.pnew_op;
%% now see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/find_T_RH_trends.m
h = hMean17years; p = pMean17years;
h = hTimeStep1; p = pTimeStep1;

pjunkN = p.plevs(1:100,:)-p.plevs(2:101,:);
pjunkD = log(p.plevs(1:100,:)./p.plevs(2:101,:));
pavgLAY = pjunkN./pjunkD;

pert = p;
pert.ptemp = pert.ptemp(1:100,:) + umbc_20_layertrends.ptemp;
pert.gas_1 = pert.gas_1(1:100,:).*(1 + umbc_20_layertrends.gas_1);

[xRH0,xRH1km0,xcolwater0] = layeramt2RH(h,p);
[xRHpert,xRH1kmpert,xcolwaterpert] = layeramt2RH(h,pert);

umbc_20_layertrends.RH = xRHpert - xRH0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
umbc_geo_rates     = umbc_20_layertrends.RH;
airsL3_geo_rates   = airsL3_100_layertrends.RH;
climcaps_geo_rates = airsCLIMCAPSL3_100_layertrends.RH;

junk = load('../../FIND_NWP_MODEL_TRENDS/MLS_atm_data_2004_09_to_2020_08_trends.mat');
mls_geo_rates = junk.trend_RH;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load('llsmap5');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);
end

rlat = junk.trend_rlat64;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear plotoptions;
plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dRH/dt'; plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'ERA5';    plotoptions.str12 = 'MERRA2';    
plotoptions.str21 = 'AIRS L3'; plotoptions.str22 = 'CLIMCAPS2'; 
plotoptions.str31 = 'CMIP6';   plotoptions.str32 = 'AMIP6';     
plotoptions.str41 = 'UMBC';    plotoptions.str42 = 'MLS L3';     
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = -1;
plotoptions.yReverseDir = +1;
plotoptions.yLimits = [100 1000];

airsL3_geo_ratesXY   = squeeze(nanmean(reshape(airsL3_geo_rates,100,72,64),2));
climcaps_geo_ratesXY = squeeze(nanmean(reshape(climcaps_geo_rates,100,72,64),2));
amip6_geo_ratesXY    = squeeze(nanmean(reshape(amip6_geo_rates,100,72,64),2));
cmip6_geo_ratesXY    = squeeze(nanmean(reshape(cmip6_geo_rates,100,72,64),2));
merra2_geo_ratesXY   = squeeze(nanmean(reshape(merra2_geo_rates,100,72,64),2));
era5_geo_ratesXY     = squeeze(nanmean(reshape(era5_geo_rates,100,72,64),2));
umbc_geo_ratesXY     = squeeze(nanmean(reshape(umbc_geo_rates,100,72,64),2));
mls_geo_ratesXY      = squeeze(nanmean(reshape(mls_geo_rates,100,72,64),2));

iFig = 3; figure(iFig); clf; profile_plots_8tiledlayout(rlat,plays,era5_geo_ratesXY,merra2_geo_ratesXY,airsL3_geo_ratesXY,climcaps_geo_ratesXY,cmip6_geo_ratesXY,amip6_geo_ratesXY,umbc_geo_ratesXY,mls_geo_ratesXY,iFig,plotoptions);
%save FIGS/Figs_JPL_Apr2022/strow_jpl_Apr2022_RHrates.mat era5_geo_ratesXY merra2_geo_ratesXY airsL3_geo_ratesXY climcaps_geo_ratesXY cmip6_geo_ratesXY amip6_geo_ratesXY umbc_geo_ratesXY mls_geo_ratesXY

