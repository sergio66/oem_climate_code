addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
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

merra2_geo_rates         = merra2_geo_rates.trend_ptemp;
era5_geo_rates           = era5_geo_rates.trend_ptemp;
amip6_geo_rates          = amip6_geo_rates.trend_ptemp;
cmip6_geo_rates          = cmip6_geo_rates.trend_ptemp;
junk = load('../../FIND_NWP_MODEL_TRENDS/MLS_atm_data_2004_09_to_2020_08_trends.mat');
mls_geo_rates = junk.trend_ptemp;
rlat = junk.trend_rlat64;

set_gather_savename_rates_AIRSL3_NWP_XMIP

plays = load(savename,'plays'); plays = plays.plays;
pavg  = load(savename,'pavg');  pavg = pavg.pavg;
junk  = load(savename,'resultsT');
  resultsT    = junk.resultsT;
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
    xjunkrate = airsL3.thestats64x72.ptemprate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      airsL3_100_layertrends.ptemp(:,ind) = junkrate;
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
    xjunkrate = airsCLIMCAPSL3.thestats64x72.ptemprate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      airsCLIMCAPSL3_100_layertrends.ptemp(:,ind) = junkrate;
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
    %xjunkrate = resultsWV(ind,:); clear junkrate
    %  junkrate = xjunkrate; 
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
    %  umbc_20_layertrends.gas_1(:,ind) = junkrate;
    %xjunkrate = resultsO3(ind,:); clear junkrate
    %  junkrate = xjunkrate;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
    %  umbc_20_layertrends.gas_3(:,ind) = junkrate;
  end
end
fprintf(1,'\n');

umbc_geo_rates     = umbc_20_layertrends.ptemp;
umbc_geo_rates_unc = umbc_20_layertrends.ptempunc;
airsL3_geo_rates   = airsL3_100_layertrends.ptemp;
climcaps_geo_rates = airsCLIMCAPSL3_100_layertrends.ptemp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear plotoptions;
plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dT/dt'; plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'ERA5';    plotoptions.str12 = 'MERRA2';    
plotoptions.str21 = 'AIRS L3'; plotoptions.str22 = 'CLIMCAPS2'; 
plotoptions.str31 = 'CMIP6';   plotoptions.str32 = 'AMIP6';     
plotoptions.str41 = 'UMBC';    plotoptions.str42 = 'MLS L3';     
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = -1;
plotoptions.yReverseDir = +1;
plotoptions.yLimits = [10 1000];

airsL3_geo_ratesXY   = squeeze(nanmean(reshape(airsL3_geo_rates,100,72,64),2));
climcaps_geo_ratesXY = squeeze(nanmean(reshape(climcaps_geo_rates,100,72,64),2));
amip6_geo_ratesXY    = squeeze(nanmean(reshape(amip6_geo_rates,100,72,64),2));
cmip6_geo_ratesXY    = squeeze(nanmean(reshape(cmip6_geo_rates,100,72,64),2));
merra2_geo_ratesXY   = squeeze(nanmean(reshape(merra2_geo_rates,100,72,64),2));
era5_geo_ratesXY     = squeeze(nanmean(reshape(era5_geo_rates,100,72,64),2));
umbc_geo_ratesXY     = squeeze(nanmean(reshape(umbc_geo_rates,100,72,64),2));
mls_geo_ratesXY      = squeeze(nanmean(reshape(mls_geo_rates,100,72,64),2));
umbc_geo_ratesXY_unc = squeeze(nanmean(reshape(umbc_geo_rates_unc,100,72,64),2))/sqrt(72);

iFig = 22; figure(iFig); clf; profile_plots_8tiledlayout(rlat,plays,era5_geo_ratesXY,merra2_geo_ratesXY,airsL3_geo_ratesXY,climcaps_geo_ratesXY,cmip6_geo_ratesXY,amip6_geo_ratesXY,umbc_geo_ratesXY,mls_geo_ratesXY,iFig,plotoptions);

iFig = 26; figure(iFig); subplot(221); pcolor(rlat,plays,umbc_geo_ratesXY_unc); colormap jet; colorbar; title('T unc')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim(plotoptions.yLimits); shading interp

if iSave > 0
  saver = ['save FIGS/Figs_JPL_Apr2022/strow_jpl_Apr2022_Trates'  savestr '..mat rlat plays era5_geo_ratesXY merra2_geo_ratesXY airsL3_geo_ratesXY climcaps_geo_ratesXY cmip6_geo_ratesXY amip6_geo_ratesXY umbc_geo_ratesXY mls_geo_ratesXY umbc_geo_ratesXY_unc'];
  eval(saver);
end
