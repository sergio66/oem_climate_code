addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS

%% see driver_gather_spectralrates_AIRSL3_NWP_XMIP6.m

load('llsmap5');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);
end

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
junk = load('../../FIND_NWP_MODEL_TRENDS/MLS_atm_data_2004_09_to_2020_08_trends.mat');
mls_geo_rates = junk.trend_RH;
rlat = junk.trend_rlat64;

set_gather_savename_rates_AIRSL3_NWP_XMIP

plays = load(savename,'plays'); plays = plays.plays;
pavg  = load(savename,'pavg');  pavg = pavg.pavg;
junk    = load(savename,'results');
  resultsSST = junk.results(:,6);
junk  = load(savename,'resultsWV');
  resultsWV = junk.resultsWV;
junk  = load(savename,'resultsWVunc');
  resultsWVunc = junk.resultsWVunc;
junk  = load(savename,'resultsT');
  resultsT = junk.resultsT;
junk  = load(savename,'resultsTunc');
  resultsTunc = junk.resultsTunc;

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
    xjunkrate = resultsTunc(ind,:); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      umbc_20_layertrends.ptempunc(:,ind) = junkrate;
    %xjunkrate = resultsRH(ind,:); clear junkrate
    %  for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
    %  umbc_20_layertrends.RH(:,ind) = junkrate;
    xjunkrate = resultsWV(ind,:); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Qlevs) | plays >= max(Qlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      umbc_20_layertrends.gas_1(:,ind) = junkrate;
    xjunkrate = resultsWVunc(ind,:); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Qlevs) | plays >= max(Qlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      umbc_20_layertrends.gas_1unc(:,ind) = junkrate;
    %xjunkrate = resultsO3(ind,:); clear junkrate
    %  junkrate = xjunkrate;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
    %  umbc_20_layertrends.gas_3(:,ind) = junkrate;
  end
end
fprintf(1,'\n');

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
h = hTimeStep1; p = pTimeStep1;      %% I been using this in eg /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_gather_gridded_retrieval_results
h = hMean17years; p = pMean17years;  %% I think I should use this

pjunkN = p.plevs(1:100,:)-p.plevs(2:101,:);
pjunkD = log(p.plevs(1:100,:)./p.plevs(2:101,:));
pavgLAY = pjunkN./pjunkD;

pert = p;
pert.ptemp(1:100,:) = pert.ptemp(1:100,:) + umbc_20_layertrends.ptemp(1:100,:);
pert.gas_1(1:100,:) = pert.gas_1(1:100,:).*(1 + umbc_20_layertrends.gas_1(1:100,:));
[xRH0,xRH1km0,xcolwater0] = layeramt2RH(h,p);
[xRHpert,xRH1kmpert,xcolwaterpert] = layeramt2RH(h,pert);
umbc_20_layertrends.RH = xRHpert - xRH0;

pertLLS = p;
pertLLS.ptemp(1:95,:) = pertLLS.ptemp(1:95,:) + umbc_20_layertrends.ptemp(1:95,:);
pertLLS.ptemp(96:100,:) = pertLLS.ptemp(96:100,:) + ones(5,1)*resultsSST';
pertLLS.gas_1(1:100,:) = pertLLS.gas_1(1:100,:).*(1 + umbc_20_layertrends.gas_1(1:100,:));
[xRH0,xRH1km0,xcolwater0] = layeramt2RH(h,p);
[xRHpertLLS,xRH1kmpertLLS,xcolwaterpertLLS] = layeramt2RH(h,pertLLS);
umbc_20_layertrends.RH = xRHpertLLS - xRH0;

%%%%%%%%%%%%%%%%%%%%%%%%%

umbc_geo_ratesT_unc  = umbc_20_layertrends.ptempunc;
umbc_geo_ratesWV_unc = umbc_20_layertrends.gas_1unc;

perty = p;
perty.ptemp = pert.ptemp(1:100,:) + umbc_20_layertrends.ptemp + umbc_geo_ratesT_unc;
perty.gas_1 = pert.gas_1(1:100,:).*(1 + umbc_20_layertrends.gas_1 + umbc_geo_ratesWV_unc);
[yRHpert,yRH1kmpert,ycolwaterpert] = layeramt2RH(h,perty);
RHunc1 = abs(yRHpert - xRH0);

perty = p;
perty.ptemp = pert.ptemp(1:100,:) + umbc_20_layertrends.ptemp - umbc_geo_ratesT_unc;
perty.gas_1 = pert.gas_1(1:100,:).*(1 + umbc_20_layertrends.gas_1 - umbc_geo_ratesWV_unc);
[yRHpert,yRH1kmpert,ycolwaterpert] = layeramt2RH(h,perty);
RHunc2 = yRHpert - xRH0;

perty = p;
perty.ptemp = pert.ptemp(1:100,:) + umbc_20_layertrends.ptemp + umbc_geo_ratesT_unc;
perty.gas_1 = pert.gas_1(1:100,:).*(1 + umbc_20_layertrends.gas_1 - umbc_geo_ratesWV_unc);
[yRHpert,yRH1kmpert,ycolwaterpert] = layeramt2RH(h,perty);
RHunc3 = yRHpert - xRH0;

perty = p;
perty.ptemp = pert.ptemp(1:100,:) + umbc_20_layertrends.ptemp - umbc_geo_ratesT_unc;
perty.gas_1 = pert.gas_1(1:100,:).*(1 + umbc_20_layertrends.gas_1 + umbc_geo_ratesWV_unc);
[yRHpert,yRH1kmpert,ycolwaterpert] = layeramt2RH(h,perty);
RHunc4 = yRHpert - xRH0;

RHunc = sqrt(RHunc1.^2 + RHunc2.^2 + RHunc3.^2 + RHunc4.^2)/4;
RHunc = nanmax(RHunc1,nanmax(RHunc2,nanmax(RHunc3,RHunc4)));

umbc_20_layertrends.RHunc = RHunc;

%{
umbc_20_layertrends.RHunc = RHunc;
umbc_geo_rates_unc = umbc_20_layertrends.RHunc;
plotoptions.yLimits = [100 1000];
  umbc_geo_ratesXY_unc = abs(umbc_geo_ratesXY_unc);  
iFig = 26; figure(iFig); subplot(223); pcolor(rlat,plays,umbc_geo_ratesXY_unc); colormap jet; colorbar; title('RH unc')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim(plotoptions.yLimits); shading interp
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
umbc_geo_rates_unc = umbc_20_layertrends.RHunc;
umbc_geo_rates     = umbc_20_layertrends.RH;
airsL3_geo_rates   = airsL3_100_layertrends.RH;
climcaps_geo_rates = airsCLIMCAPSL3_100_layertrends.RH;

boo = min(airsCLIMCAPSL3.Qlevs); boo = find(plays <= boo); climcaps_geo_rates(boo,:) = 0;
boo = min(airsL3.Qlevs);         boo = find(plays <= boo); airsL3_geo_rates(boo,:) = 0;
fprintf(1,'min AIRS WV Qlevs = %8.6f  min CLIMCAPS WV Qlev = %8.6f \n',min(airsL3.Qlevs),min(airsCLIMCAPSL3.Qlevs))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear plotoptions;
plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dRH/dt'; plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'ERA5';    plotoptions.str12 = 'MERRA2';    
plotoptions.str21 = 'AIRS L3'; plotoptions.str22 = 'CLIMCAPS2'; 
plotoptions.str31 = 'CMIP6';   plotoptions.str32 = 'AMIP6';     
plotoptions.str41 = 'UMBC';    plotoptions.str42 = 'MLS L3';     
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = -1;
%plotoptions.yLinearOrLog = +1;
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
umbc_geo_ratesXY_unc = squeeze(nanmean(reshape(umbc_geo_rates_unc,100,72,64),2))/sqrt(72);
  umbc_geo_ratesXY_unc = abs(umbc_geo_ratesXY_unc);  

iFig = 24; figure(iFig); clf; profile_plots_8tiledlayout(rlat,plays,era5_geo_ratesXY,merra2_geo_ratesXY,airsL3_geo_ratesXY,climcaps_geo_ratesXY,cmip6_geo_ratesXY,amip6_geo_ratesXY,umbc_geo_ratesXY,mls_geo_ratesXY,iFig,plotoptions);

iFig = 26; figure(iFig); subplot(223); pcolor(rlat,plays,umbc_geo_ratesXY_unc); colormap jet; colorbar; title('RH unc')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim(plotoptions.yLimits); shading interp

iFig = 27; figure(iFig); clf
subplot(221); zT = p.ptemp(1:100,:);                                                                     zT = squeeze(nanmean(reshape(zT,100,72,64),2));  pcolor(rlat,plays,zT); colormap jet; colorbar; title('T (K)'); caxis([200 300])
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); shading interp
subplot(222); zWV = layers2ppmv(h,p,1:length(p.stemp),1); 
  [mm,nn] = size(zWV);
  for jj = mm+1:100
    zWV(jj,:) = zWV(mm,:);
    zWV(jj,:) = NaN;
  end
  zWV = squeeze(nanmean(reshape(zWV,100,72,64),2));  pcolor(rlat,plays,zWV); colormap jet; colorbar; title('WV ppmv'); caxis([0 3e4])
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); shading interp
subplot(223); zRH = xRH0;                                                                                zRH = squeeze(nanmean(reshape(zRH,100,72,64),2));  pcolor(rlat,plays,zRH); colormap jet; colorbar; title('RH'); caxis([0 100])
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); shading interp
subplot(224); zO3 = layers2ppmv(h,p,1:length(p.stemp),3); 
  [mm,nn] = size(zO3);
  for jj = mm+1:100
    zO3(jj,:) = zO3(mm,:);
    zO3(jj,:) = NaN;
  end  
  zO3 = squeeze(nanmean(reshape(zO3,100,72,64),2));  pcolor(rlat,plays,zO3); colormap jet; colorbar; title('O3 ppmv'); caxis([0 10])
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([0.1 200]); shading interp

ind_layer_rates

if iSave > 0
  saver = ['save FIGS/Figs_JPL_Apr2022/strow_jpl_Apr2022_RHrates' savestr '.mat rlat plays'];
  saver = [saver ' era5_geo_ratesXY merra2_geo_ratesXY airsL3_geo_ratesXY climcaps_geo_ratesXY cmip6_geo_ratesXY amip6_geo_ratesXY umbc_geo_ratesXY mls_geo_ratesXY umbc_geo_ratesXY_unc zT zWV zRH zO3'];
  eval(saver)
end

