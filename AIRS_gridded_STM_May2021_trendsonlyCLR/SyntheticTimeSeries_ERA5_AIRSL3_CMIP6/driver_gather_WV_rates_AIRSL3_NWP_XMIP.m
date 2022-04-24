addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
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

merra2_geo_rates         = merra2_geo_rates.trend_gas_1;
era5_geo_rates           = era5_geo_rates.trend_gas_1;
amip6_geo_rates          = amip6_geo_rates.trend_gas_1;
cmip6_geo_rates          = cmip6_geo_rates.trend_gas_1;
junk = load('../../FIND_NWP_MODEL_TRENDS/MLS_atm_data_2004_09_to_2020_08_trends.mat');
mls_geo_rates = junk.trend_gas_1;
rlat = junk.trend_rlat64;

set_gather_savename_rates_AIRSL3_NWP_XMIP

plays = load(savename,'plays'); plays = plays.plays;
pavg  = load(savename,'pavg');  pavg = pavg.pavg;
junk  = load(savename,'resultsWV');
  resultsWV = junk.resultsWV;
junk  = load(savename,'resultsWVunc');
  resultsWVunc = junk.resultsWVunc;

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
    xjunkrate = airsL3.thestats64x72.waterrate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Qlevs) | plays >= max(Qlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      airsL3_100_layertrends.gas_1(:,ind) = junkrate;
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
    xjunkrate = airsCLIMCAPSL3.thestats64x72.waterrate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Qlevs) | plays >= max(Qlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      airsCLIMCAPSL3_100_layertrends.gas_1(:,ind) = junkrate;
    %xjunkrate = airsCLIMCAPSL3.thestats64x72.ozonerate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
    %  for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0;       
    %  airsCLIMCAPSL3_100_layertrends.gas_3(:,ind) = junkrate;
  
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    Qlevs = pavg;
    Tlevs = pavg;
    %  umbc_20_layertrends.stemp(ind) = results(ind,6);
    %xjunkrate = resultsT(ind,:); clear junkrate
    %  for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
    %  junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
    %  umbc_20_layertrends.ptemp(:,ind) = junkrate;
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

umbc_geo_rates     = umbc_20_layertrends.gas_1;
umbc_geo_rates_unc = umbc_20_layertrends.gas_1unc;
airsL3_geo_rates   = airsL3_100_layertrends.gas_1;
climcaps_geo_rates = airsCLIMCAPSL3_100_layertrends.gas_1;

boo = min(airsCLIMCAPSL3.Qlevs); boo = find(plays <= boo); climcaps_geo_rates(boo,:) = 0;
boo = min(airsL3.Qlevs);         boo = find(plays <= boo); airsL3_geo_rates(boo,:) = 0;
fprintf(1,'min AIRS WV Qlevs = %8.6f  min CLIMCAPS WV Qlev = %8.6f \n',min(airsL3.Qlevs),min(airsCLIMCAPSL3.Qlevs))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear plotoptions;
plotoptions.cx = [-1 +1]*0.15/10; plotoptions.maintitle = 'dWVfrac/dt'; plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'ERA5';    plotoptions.str12 = 'MERRA2';    
plotoptions.str21 = 'AIRS L3'; plotoptions.str22 = 'CLIMCAPS2'; 
plotoptions.str31 = 'CMIP6';   plotoptions.str32 = 'AMIP6';     
plotoptions.str41 = 'UMBC';    plotoptions.str42 = 'MLS L3';     
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = -1;
plotoptions.yLinearOrLog = +1;
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

iFig = 23; figure(iFig); clf; profile_plots_8tiledlayout(rlat,plays,era5_geo_ratesXY,merra2_geo_ratesXY,airsL3_geo_ratesXY,climcaps_geo_ratesXY,cmip6_geo_ratesXY,amip6_geo_ratesXY,umbc_geo_ratesXY,mls_geo_ratesXY,iFig,plotoptions);

iFig = 26; figure(iFig); subplot(222); pcolor(rlat,plays,umbc_geo_ratesXY_unc); colormap jet; colorbar; title('WV frac unc')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim(plotoptions.yLimits); shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/
get_the_mean_profiles  %% this is if you want to change from fracWV or fracO3 to eg ppmv

zWV = layers2ppmv(h,p,1:length(p.stemp),1); 
  [mm,nn] = size(zWV);
  for jj = mm+1:100
    zWV(jj,:) = zWV(mm,:);
    zWV(jj,:) = NaN;
  end
figure(30); clf
  zWV0 = zWV; zWV = squeeze(nanmean(reshape(zWV,100,72,64),2));  pcolor(rlat,plays,zWV); colormap jet; colorbar; title('WV ppmv'); caxis([0 3e4])
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); shading interp

iFig = 30; figure(iFig); clf; 
plotoptions.cx = [-1 +1]*100; plotoptions.maintitle = 'dWV(ppmv)/dt'; plotoptions.plotcolors = llsmap5;
profile_plots_8tiledlayout(rlat,plays,era5_geo_ratesXY.*zWV,merra2_geo_ratesXY.*zWV,airsL3_geo_ratesXY.*zWV,...
                           climcaps_geo_ratesXY.*zWV,cmip6_geo_ratesXY.*zWV,amip6_geo_ratesXY.*zWV,umbc_geo_ratesXY.*zWV,mls_geo_ratesXY.*zWV,iFig,plotoptions);

iFig = 31; figure(iFig); clf; 
plotoptions.cx = [-1 +1]*4; plotoptions.maintitle = 'dWV(ppmv)/dt'; plotoptions.plotcolors = llsmap5;
profile_plots_8tiledlayout(rlat,plays,abslog(era5_geo_ratesXY.*zWV),abslog(merra2_geo_ratesXY.*zWV),abslog(airsL3_geo_ratesXY.*zWV),...
                                      abslog(climcaps_geo_ratesXY.*zWV),abslog(cmip6_geo_ratesXY.*zWV),abslog(amip6_geo_ratesXY.*zWV),abslog(umbc_geo_ratesXY.*zWV),abslog(mls_geo_ratesXY.*zWV),iFig,plotoptions);

ind_layer_rates

if iSave > 0
  saver = ['save FIGS/Figs_JPL_Apr2022/strow_jpl_Apr2022_WVrates'  savestr '.mat rlat plays era5_geo_ratesXY merra2_geo_ratesXY airsL3_geo_ratesXY climcaps_geo_ratesXY cmip6_geo_ratesXY amip6_geo_ratesXY umbc_geo_ratesXY mls_geo_ratesXY umbc_geo_ratesXY_unc'];
  eval(saver)
end


