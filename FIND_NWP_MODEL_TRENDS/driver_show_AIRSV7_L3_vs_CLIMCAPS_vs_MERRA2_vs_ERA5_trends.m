function [era5,merra2,airsL3,climcapsL3,umbc,thecorr,amp] = driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends(strUMBC,iNumYears,iPentagonPlot);

%{
  strUMBC = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithMLSL3_uncX100_50fatlayers_AIRSL3_ERA5_CMIP6_globalSSTfeedback.mat';
  strUMBC = '/asl/s1/sergio/JUNK/test9_guessstartWV_Vers1_march22_2023.mat';
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q05_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_feedback.mat';              %% not too bad at lower atm/polar!!!!
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q05_newERA5_2021jacs_startwith_MLSL3_TOA_guessWV_dRH_zero_bot_50fatlayers.mat';  %% too much oomph at gnd : use MLS L3 TOA and dRH/dt = 0 at bottom

  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q16_newERA5_2021jacs_startwith0_50fatlayers_ERA5calcs_spectraltrends.mat';       iNumYears = 20; %% compare geophysical trends from ERA5     spectral trends
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q16_newERA5_2021jacs_startwith0_50fatlayers_MERRA2calcs_spectraltrends.mat';     iNumYears = 20; %% compare geophysical trends from MERRA2   spectral trends
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q16_newERA5_2021jacs_startwith0_50fatlayers_AIRSL3calcs_spectraltrends.mat';     iNumYears = 20; %% compare geophysical trends from AIRS L3  spectral trends
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q16_newERA5_2021jacs_startwith0_50fatlayers_CLIMCAPSL3calcs_spectraltrends.mat'; iNumYears = 20; %% compare geophysical trends from CLIMCAPS spectral trends

  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_MLS.mat'; iNumYears = 20;     %% use CarbonTracker CO2 trends, MLS a priori
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 20;         %% use CarbonTracker CO2 trends

  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset10_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 05;         %% use CarbonTracker CO2 trends
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset11_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 10;         %% use CarbonTracker CO2 trends
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset12_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 15;         %% use CarbonTracker CO2 trends

  driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends(strUMBC,iNumYears);  
%}

set(0,'DefaultaxesLineWidth',1);
set(0,'DefaultaxesFontSize',16);

amp = [];

if nargin == 0

  strUMBC = [];
  iUMBC = -1;
  umbc = [];
  iNumYears = -1;
  iPentagonPlot = +1;

  thecorr = [];

elseif nargin == 1
  iNumYears = 20;
  iPentagonPlot = +1;
  iUMBC = +1;
elseif nargin == 2
  iPentagonPlot = +1;
  iUMBC = +1;
else
  iUMBC = +1;
end

%% see How closely do changes in surface and column water vapor follow Clausius-Clapeyron scaling in climate-change simulations?
%% P A O’Gorman, C J Muller, https://core.ac.uk/download/pdf/4426849.pdf, or saved in PDF/change_of_RH_with_stemp_GCM_PGorman.pdf
%% P A O'Gorman and C J Muller 2010 Environ. Res. Lett. 5 025207 DOI 10.1088/1748-9326/5/2/025207

addpath /home/sergio/MATLABCODE/matlib/science/
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/NANROUTINES
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools
addpath /asl/matlib/plotutils
addpath /home/sergio/MATLABCODE/matlib/science/
addpath ../AIRS_gridded_STM_May2021_trendsonlyCLR
addpath /home/sergio/MATLABCODE/PLOTMISC

load llsmap5

plays100 = load('ERA5_atm_data_2002_09_to_2022_08_trends_desc.mat','trend_plays');
plays100 = plays100.trend_plays;

%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2; 
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

addpath /home/sergio/MATLABCODE/matlib/science/            %% for usgs_deg10_dem.m that has correct paths
[salti, landfrac] = usgs_deg10_dem(Y(:),X(:));
Ylat = Y(:);
Xlon = X(:);
mu = cos(Ylat*pi/180); mu = mu';

%%%%%%%%%%%%%%%%%%%%%%%%%

[h,ha,p,pa] = rtpread('summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.op.rtp');
mmw0 = mmwater_rtp(h,p);
[salti, landfrac] = usgs_deg10_dem(p.rlat, p.rlon);
lfmaskA = ones(1,4608);
lfmaskL = (landfrac > 0);
lfmaskO = (landfrac == 0);

for ii = 1 : 12;
  figure(ii); clf
end

%iX = input('fastgrib (+1/default) or slow calc (-1) : ');
%if length(iX) == 0
%  iX = 1;
%end
iX = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iX > 0
  %load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat
  fin = ['/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug' num2str(2002+iNumYears) '_' num2str(iNumYears) 'yr_desc.mat'];
  iNumYearsPowWow = iNumYears;
  if ~exist(fin)
    iaJunk = 20;
    moo = abs(iaJunk-iNumYears);
    moo = find(moo == min(moo),1);    
    fin = ['/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug' num2str(2002+iaJunk(moo)) '_' num2str(iaJunk(moo)) 'yr_desc.mat'];    
    iNumYearsPowWow = 2002+iaJunk(moo);
  end
  loader = ['load ' fin];
  eval(loader)
  fprintf(1,'loaded %2i years AIRS     geophysical trends from %s \n',iNumYearsPowWow,fin)
else
  fin = ['/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug' num2str(2002+iNumYears) '_' num2str(iNumYears) 'yr_desc.mat'];
  iNumYearsPowWow = iNumYears;
  if ~exist(fin)
    iaJunk = [5 10 15 20 12 18 19];
    moo = abs(iaJunk-iNumYears);
    moo = find(moo == min(moo),1);    
    iNumYearsPowWow = 2002+iaJunk(moo);
    fin = ['/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug' num2str(2002+iaJunk(moo)) '_' num2str(iaJunk(moo)) 'yr_desc.mat'];    
  end
  loader = ['load ' fin];
  eval(loader)
  fprintf(1,'loaded %2i years AIRS     geophysical trends from %s \n',iNumYearsPowWow,fin)
end

z11 = thestats64x72.stemprate(:);
z11 = z11';

for ii = 1 : 72
  for jj = 1 : 64
    junk = squeeze(thestats64x72.RHrate(ii,jj,:));
    airsL3.RHrate(ii,jj,:) = interp1(log10(Qlevs),junk,log10(plays100),[],'extrap');
    junk = squeeze(thestats64x72.waterrate(ii,jj,:));
    airsL3.waterrate(ii,jj,:) = interp1(log10(Qlevs),junk,log10(plays100),[],'extrap');
    junk = squeeze(thestats64x72.ptemprate(ii,jj,:));
    airsL3.ptemprate(ii,jj,:) = interp1(log10(Tlevs),junk,log10(plays100),[],'extrap');
  end
end
airsL3.stemprate = thestats64x72.stemprate(:)';
bad = find(isnan(airsL3.stemprate) | isinf(airsL3.stemprate)); airsL3.stemprate(bad) = 0;
bad = find(isnan(airsL3.ptemprate) | isinf(airsL3.ptemprate)); airsL3.ptemprate(bad) = 0;
bad = find(isnan(airsL3.waterrate) | isinf(airsL3.waterrate)); airsL3.waterrate(bad) = 0;

pert = p;
junk = reshape(permute(airsL3.ptemprate,[3 1 2]),100,4608); pert.ptemp(1:100,:) = pert.ptemp(1:100,:) + junk;
junk = reshape(permute(airsL3.waterrate,[3 1 2]),100,4608); pert.gas_1(1:100,:) = pert.gas_1(1:100,:) .* (1 + junk);
pert.stemp = pert.stemp + airsL3.stemprate;
mmw_airsL3 = mmwater_rtp(h,pert);
airsL3.trend_mmw = mmw_airsL3 - mmw0;
figure(13); dmmw_dST_airsL3 = dmmw_dsst_VS_lat(lfmaskA,lfmaskL,lfmaskO,mmw0,airsL3.trend_mmw,airsL3.stemprate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fin = ['/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug' num2str(2002+iNumYears) '_' num2str(iNumYears) 'yr_desc.mat'];
iNumYearsPowWow = iNumYears;
if ~exist(fin)
  iaJunk = [5 10 15 20];
  moo = abs(iaJunk-iNumYears);
  moo = find(moo == min(moo),1);    
  fin = ['/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug' num2str(2002+iaJunk(moo)) '_' num2str(iaJunk(moo)) 'yr_desc.mat'];    
  iNumYearsPowWow = 2002+iaJunk(moo);
end
loader = ['load ' fin];
eval(loader)
fprintf(1,'loaded %2i years CLIMCAPS geophysical trends from %s \n',iNumYearsPowWow,fin)

z12 = thestats64x72.stemprate(:);
z12 = z12';

for ii = 1 : 72
  for jj = 1 : 64
    junk = squeeze(thestats64x72.RHrate(ii,jj,:));
    climcapsL3.RHrate(ii,jj,:) = interp1(log10(Qlevs/100),100*junk,log10(plays100),[],'extrap');
    junk = squeeze(thestats64x72.waterrate(ii,jj,:));
    climcapsL3.waterrate(ii,jj,:) = interp1(log10(Qlevs/100),junk,log10(plays100),[],'extrap');
    junk = squeeze(thestats64x72.ptemprate(ii,jj,:));
    climcapsL3.ptemprate(ii,jj,:) = interp1(log10(Tlevs/100),junk,log10(plays100),[],'extrap');
  end
end
climcapsL3.stemprate = thestats64x72.stemprate(:)';
bad = find(isnan(climcapsL3.stemprate) | isinf(climcapsL3.stemprate)); climcapsL3.stemprate(bad) = 0;
bad = find(isnan(climcapsL3.ptemprate) | isinf(climcapsL3.ptemprate)); climcapsL3.ptemprate(bad) = 0;
bad = find(isnan(climcapsL3.waterrate) | isinf(climcapsL3.waterrate)); climcapsL3.waterrate(bad) = 0;

pert = p;
junk = reshape(permute(climcapsL3.ptemprate,[3 1 2]),100,4608); pert.ptemp(1:100,:) = pert.ptemp(1:100,:) + junk;
junk = reshape(permute(climcapsL3.waterrate,[3 1 2]),100,4608); pert.gas_1(1:100,:) = pert.gas_1(1:100,:) .* (1 + junk);
pert.stemp = pert.stemp + climcapsL3.stemprate;
mmw_climcapsL3 = mmwater_rtp(h,pert);
climcapsL3.trend_mmw = mmw_climcapsL3 - mmw0;
figure(13); clf; dmmw_dST_climcapsL3 = dmmw_dsst_VS_lat(lfmaskA,lfmaskL,lfmaskO,mmw0,climcapsL3.trend_mmw,climcapsL3.stemprate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fin = ['MERRA2_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc.mat'];
iNumYearsPowWow = iNumYears;
if ~exist(fin)
  iaJunk = [5 10 15 20];
  moo = abs(iaJunk-iNumYears);
  moo = find(moo == min(moo),1);    
  fin = ['MERRA2_atm_data_2002_09_to_' num2str(2002+iaJunk(moo)) '_08_trends_desc.mat']; 
  iNumYearsPowWow = 2002+iaJunk(moo);
end
loader = ['load ' fin];
eval(loader)
fprintf(1,'loaded %2i years MERRA2   geophysical trends from %s \n',iNumYearsPowWow,fin)

z21 = trend_stemp;
merra2.stemprate = trend_stemp;
merra2.RHrate   = trend_RH;
merra2.waterrate = trend_gas_1;
merra2.ptemprate = trend_ptemp;
merra2.trend_mmw = trend_mmw;

bad = find(isnan(merra2.stemprate) | isinf(merra2.stemprate)); merra2.stemprate(bad) = 0;
bad = find(isnan(merra2.ptemprate) | isinf(merra2.ptemprate)); merra2.ptemprate(bad) = 0;
bad = find(isnan(merra2.waterrate) | isinf(merra2.waterrate)); merra2.waterrate(bad) = 0;

iFig = 14;
iFig = 1;
figure(iFig); clf
[h,ha,p,pa] = rtpread('summary_19years_all_lat_all_lon_2002_2021_monthlyMERRA2.op.rtp');
figure(13); clf; dmmw_dST_MERRA2 = dmmw_dsst_VS_lat(lfmaskA,lfmaskL,lfmaskO,mmw0,merra2.trend_mmw,merra2.stemprate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iNumYears > 0
  eraX = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc.mat'];
  eraX = ['ERA5_atm_N_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc.mat'];
else
  eraX = ['ERA5_atm_data_2002_09_to_' num2str(2022) '_08_trends_desc.mat'];
end
loader = ['load ' eraX];
fprintf(1,'loaded %2i years ERA5     geophysical trends from %s \n',iNumYears,eraX)
eval(loader);

z22 = trend_stemp;
era5.stemprate = trend_stemp;
era5.RHrate    = trend_RH;
era5.waterrate = trend_gas_1;
era5.ptemprate = trend_ptemp;
era5.trend_mmw = trend_mmw;

bad = find(isnan(era5.stemprate) | isinf(era5.stemprate)); era5.stemprate(bad) = 0;
bad = find(isnan(era5.ptemprate) | isinf(era5.ptemprate)); era5.ptemprate(bad) = 0;
bad = find(isnan(era5.waterrate) | isinf(era5.waterrate)); era5.waterrate(bad) = 0;

figure(iFig); clf
[h,ha,p,pa] = rtpread('summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.op.rtp');
mmw0 = mmwater_rtp(h,p);
figure(13); clf; dmmw_dST_ERA5 = dmmw_dsst_VS_lat(lfmaskA,lfmaskL,lfmaskO,mmw0,era5.trend_mmw,era5.stemprate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(13); clf
plot(trend_rlat64,dmmw_dST_MERRA2.all,'b',trend_rlat64,smooth(dmmw_dST_MERRA2.all,10),'c',trend_rlat64,dmmw_dST_ERA5.all,'r',trend_rlat64,smooth(dmmw_dST_ERA5.all,10),'m','linewidth',2); 
plot(trend_rlat64,smooth(dmmw_dST_MERRA2.all,10),'b',trend_rlat64,smooth(dmmw_dST_ERA5.all,10),'r','linewidth',2); 
plotaxis2; hl = legend('MERRA2','ERA5','location','best'); ylabel('d%mmw/dST'); xlabel('Latitude');

figure(13); clf
plot(trend_rlat64,smooth(dmmw_dST_MERRA2.all,10),'g',trend_rlat64,smooth(dmmw_dST_ERA5.all,10),'r',trend_rlat64,smooth(dmmw_dST_airsL3.all,10),'b',trend_rlat64,smooth(dmmw_dST_climcapsL3.all,10),'c','linewidth',2); 
plotaxis2; hl = legend('MERRA2','ERA5','AIRS L3','CLIMCAPSL3','location','best','fontsize',8); ylabel('d%mmw/dST'); xlabel('Latitude'); ylim([-10 +20]); title('Land + Ocean')

figure(14); clf
plot(trend_rlat64,smooth(dmmw_dST_MERRA2.land,10),'g',trend_rlat64,smooth(dmmw_dST_ERA5.land,10),'r',trend_rlat64,smooth(dmmw_dST_airsL3.land,10),'b',trend_rlat64,smooth(dmmw_dST_climcapsL3.land,10),'c','linewidth',2); 
plotaxis2; hl = legend('MERRA2','ERA5','AIRS L3','CLIMCAPSL3','location','best','fontsize',8); ylabel('d%mmw/dST'); xlabel('Latitude'); ylim([-10 +20]); title('Land')

figure(15); clf
plot(trend_rlat64,smooth(dmmw_dST_MERRA2.ocean,10),'g',trend_rlat64,smooth(dmmw_dST_ERA5.ocean,10),'r',trend_rlat64,smooth(dmmw_dST_airsL3.ocean,10),'b',trend_rlat64,smooth(dmmw_dST_climcapsL3.ocean,10),'c','linewidth',2); 
plotaxis2; hl = legend('MERRA2','ERA5','AIRS L3','CLIMCAPSL3','location','best','fontsize',8); ylabel('d%mmw/dST'); xlabel('Latitude'); ylim([-10 +20]); title('Ocean')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = iFig + 1;
figure(iFig); clf; clear plotoptions
maskLF = ones(size(z22));
plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'AIRS L3';
plotoptions.str12 = 'CLIMCAPS L3';
plotoptions.str21 = 'MERRA2';
plotoptions.str22 = 'ERA5';
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = +1;
aslmap_2x2tiledlayout(z11,z12,z21,z22,iFig,plotoptions);

disp('SKT trend K/yr  | AIRS   CLIMCAPS   MERRA2      ERA5    ');
disp('----------------|---------------------------------------');
junk = [nanmean(z11)  nanmean(z12) nanmean(z21)  nanmean(z22)];
fprintf(1,'(ALL  )         |  %8.6f %8.6f %8.6f %8.6f \n',junk);
boo = find(landfrac == 0);  
junk = [nanmean(z11(boo))  nanmean(z12(boo)) nanmean(z21(boo))  nanmean(z22(boo))];
fprintf(1,'(OCEAN)         |  %8.6f %8.6f %8.6f %8.6f \n',junk);
boo = find(landfrac > 0.99);  
junk = [nanmean(z11(boo))  nanmean(z12(boo)) nanmean(z21(boo))  nanmean(z22(boo))];
fprintf(1,'(LAND )         |  %8.6f %8.6f %8.6f %8.6f \n',junk);
disp('----------------|---------------------------------------');

disp('SKT trend K/yr  | AIRS   CLIMCAPS   MERRA2      ERA5    ');
disp('area weighted with cos(lat)')
disp('----------------|---------------------------------------');
junk = [nansum(z11.*mu)/nansum(mu)  nansum(z12.*mu)/nansum(mu) nansum(z21.*mu)/nansum(mu)  nansum(z22.*mu)/nansum(mu)];
fprintf(1,'(ALL  )         |  %8.6f %8.6f %8.6f %8.6f \n',junk);
boo = find(landfrac == 0);  
junk = [nansum(z11(boo).*mu(boo))/nansum(mu(boo))  nansum(z12(boo).*mu(boo))/nansum(mu(boo)) nansum(z21(boo).*mu(boo))/nansum(mu(boo))  nansum(z22(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'(OCEAN)         |  %8.6f %8.6f %8.6f %8.6f \n',junk);
boo = find(landfrac > 0.99);  
junk = [nansum(z11(boo).*mu(boo))/nansum(mu(boo))  nansum(z12(boo).*mu(boo))/nansum(mu(boo)) nansum(z21(boo).*mu(boo))/nansum(mu(boo))  nansum(z22(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'(LAND )         |  %8.6f %8.6f %8.6f %8.6f \n',junk);
disp('----------------|---------------------------------------');

clear plotoptions
plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'AIRS L3';
plotoptions.str12 = 'CLIMCAPS L3';
plotoptions.str21 = 'ERA5';
plotoptions.str22 = 'MERRA2';
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = +1;
plotoptions.yReverseDir  = +1;
plotoptions.yLimits = [100 1000];
plotoptions.xaxis_metric = 'sine';

iFig = iFig + 1;
figure(iFig); clf; 
plotoptions.maintitle = 'RH trends [%/yr]'; 
plotoptions.cx = [-1 +1]*0.5; 
miaow11 = squeeze(nanmean(permute(airsL3.RHrate,[3 1 2]),2));
miaow12 = squeeze(nanmean(permute(climcapsL3.RHrate,[3 1 2]),2));
miaow21 = squeeze(nanmean(reshape(era5.RHrate,100,72,64),2));
miaow22 = squeeze(nanmean(reshape(merra2.RHrate,100,72,64),2));
profile_plots_2x2tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

iFig = iFig + 1;
figure(iFig); clf; 
plotoptions.maintitle = 'd(log(WV))/dt [1/yr]'; 
plotoptions.cx = [-1 +1]*0.015; 
miaow11 = squeeze(nanmean(permute(airsL3.waterrate,[3 1 2]),2));
miaow12 = squeeze(nanmean(permute(climcapsL3.waterrate,[3 1 2]),2));
miaow21 = squeeze(nanmean(reshape(era5.waterrate,100,72,64),2));
miaow22 = squeeze(nanmean(reshape(merra2.waterrate,100,72,64),2));
profile_plots_2x2tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

iFig = iFig + 1;
figure(iFig); clf; 
plotoptions.maintitle = 'Temperature Trends [K/yr]'; 
plotoptions.cx = [-1 +1]*0.151; 
plotoptions.yLinearOrLog = -1;
plotoptions.yReverseDir  = +1;
plotoptions.yLimits = [10 1000];
miaow11 = squeeze(nanmean(permute(airsL3.ptemprate,[3 1 2]),2));
miaow12 = squeeze(nanmean(permute(climcapsL3.ptemprate,[3 1 2]),2));
miaow21 = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
miaow22 = squeeze(nanmean(reshape(merra2.ptemprate,100,72,64),2));
profile_plots_2x2tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%iUMBC = input('Load UMBC (-1 default no /  +1 yes) : ');
%if length(iUMBC) == 0
%  iUMBC = -1;
%end

if iUMBC  > 0
  %strUMBC = input('Enter UMBC results fname eg ''/asl/s1/sergio/JUNK/test5_march9_2023.mat'' : ');
  fprintf(1,'<<< loading in %s >>> \n',strUMBC);
  umbcX        = load(strUMBC,'deltaRH','deltaT','fracWV','results');

  strGISS = ['ChrisHTrends/giss_trends_2002_' num2str(2002+iNumYears) '.mat'];
  iNumYearsPowWow = iNumYears;
  if ~exist(strGISS)
    iaJunk = [05 10 15 818 19 20];
    moo = abs(iaJunk-iNumYears);
    moo = find(moo == min(moo),1);    
    strGISS = ['ChrisHTrends/giss_trends_2002_' num2str(2002+iaJunk(moo)) '.mat'];
    iNumYearsPowWow = 2002+iaJunk(moo);
  end

  fprintf(1,'loaded %2i years %s \n',iNumYearsPowWow,strGISS);
  loader = ['load ' strGISS];
  eval(loader);

  if isfield(umbcX,'fracWV')
    umbcX.waterrate = umbcX.fracWV(1:100,:);
    umbcX.deltaT = umbcX.deltaT(1:100,:);
  else
    umbcX = load(strUMBC,'resultsT','resultsWV','pjunk20','results','deltaRH');
    for ii = 1 : 4608
      junk = umbcX.resultsWV(ii,:);
      umbcX.fracWV(ii,:) = interp1(log10(umbcX.pjunk20),junk,log10(plays100),[],'extrap')';
      junk = umbcX.resultsT(ii,:);
      umbcX.deltaT(ii,:) = interp1(log10(umbcX.pjunk20),junk,log10(plays100),[],'extrap')';
    end
    umbcX.fracWV = umbcX.fracWV';
    umbcX.deltaT = umbcX.deltaT';
  end
  umbcX.fracWV = umbcX.fracWV(1:100,:);
  umbcX.stemprate = umbcX.results(:,6)';

  umbc.stemprate  = umbcX.stemprate;
  umbc.ptemprate  = umbcX.deltaT;
  umbc.RHrate     = umbcX.deltaRH;
  umbc.waterrate  = umbcX.fracWV;

  clear umbcX

  %%%%%%%%%%%%%%%%%%%%%%%%%
  figure(13); clf
  pert = p;
  junk = umbc.ptemprate; pert.ptemp(1:100,:) = pert.ptemp(1:100,:) + junk;
  junk = umbc.waterrate; pert.gas_1(1:100,:) = pert.gas_1(1:100,:) .* (1 + junk);
  pert.stemp = pert.stemp + umbc.stemprate;
  mmw_umbc = mmwater_rtp(h,pert);
  umbc.trend_mmw = mmw_umbc - mmw0;
  dmmw_dST_umbc = dmmw_dsst_VS_lat(lfmaskA,lfmaskL,lfmaskO,mmw0,umbc.trend_mmw,umbc.stemprate);

  figure(13); clf
  plot(trend_rlat64,smooth(dmmw_dST_MERRA2.ocean,10),'g',trend_rlat64,smooth(dmmw_dST_ERA5.ocean,10),'r',trend_rlat64,smooth(dmmw_dST_airsL3.ocean,10),'b',trend_rlat64,smooth(dmmw_dST_climcapsL3.ocean,10),'c',...
       trend_rlat64,smooth(dmmw_dST_umbc.ocean,10),'k','linewidth',2); 
  plotaxis2; hl = legend('MERRA2','ERA5','AIRS L3','CLIMCAPSL3','UMBC','location','best','fontsize',8); ylabel('d%mmw/dST'); xlabel('Latitude'); ylim([-10 +20]); title('Ocean')

  figure(14); clf
  plot(trend_rlat64,smooth(dmmw_dST_MERRA2.land,10),'g',trend_rlat64,smooth(dmmw_dST_ERA5.land,10),'r',trend_rlat64,smooth(dmmw_dST_airsL3.land,10),'b',trend_rlat64,smooth(dmmw_dST_climcapsL3.land,10),'c',...
       trend_rlat64,smooth(dmmw_dST_umbc.land,10),'k','linewidth',2); 
  plotaxis2; hl = legend('MERRA2','ERA5','AIRS L3','CLIMCAPSL3','UMBC','location','best','fontsize',8); ylabel('d%mmw/dST'); xlabel('Latitude'); ylim([-10 +20]); title('Land')

  figure(15); clf
  plot(trend_rlat64,smooth(dmmw_dST_MERRA2.all,10),'g',trend_rlat64,smooth(dmmw_dST_ERA5.all,10),'r',trend_rlat64,smooth(dmmw_dST_airsL3.all,10),'b',trend_rlat64,smooth(dmmw_dST_climcapsL3.all,10),'c',...
       trend_rlat64,smooth(dmmw_dST_umbc.all,10),'k','linewidth',2); 
  plotaxis2; hl = legend('MERRA2','ERA5','AIRS L3','CLIMCAPSL3','UMBC','location','best','fontsize',8); ylabel('d%mmw/dST'); xlabel('Latitude'); ylim([-10 +20]); title('All')

  %%%%%%%%%%%%%%%%%%%%%%%%%

  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  maskLF = ones(size(z22));
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'UMBC';
  plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';
  plotoptions.str22 = 'ERA5';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  z11x = umbc.stemprate;
  aslmap_2x2tiledlayout(z11x,z12,z21,z22,iFig,plotoptions);
  
  clear plotoptions
  plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'UMBC';
  plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'ERA5';
  plotoptions.str22 = 'MERRA2';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  plotoptions.yReverseDir  = +1;
  plotoptions.xaxis_metric = 'sine';
  
  iFig = iFig + 1;
  figure(iFig); clf; 
  plotoptions.maintitle = 'RH trends [%/yr]'; 
  plotoptions.cx = [-1 +1]*0.5; 
  plotoptions.yLimits = [100 1000];
  miaow11 = squeeze(nanmean(reshape(umbc.RHrate,100,72,64),2));
  miaow12 = squeeze(nanmean(permute(climcapsL3.RHrate,[3 1 2]),2));
  miaow21 = squeeze(nanmean(reshape(era5.RHrate,100,72,64),2));
  miaow22 = squeeze(nanmean(reshape(merra2.RHrate,100,72,64),2));
  profile_plots_2x2tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

  iFig = iFig + 1;
  figure(iFig); clf; 
  plotoptions.maintitle = 'd(log(WV))/dt [1/yr]'; 
  plotoptions.cx = [-1 +1]*0.015; 
  plotoptions.yLimits = [100 1000];
  miaow11 = squeeze(nanmean(reshape(umbc.waterrate,100,72,64),2));
  miaow12 = squeeze(nanmean(permute(climcapsL3.waterrate,[3 1 2]),2));
  miaow21 = squeeze(nanmean(reshape(era5.waterrate,100,72,64),2));
  miaow22 = squeeze(nanmean(reshape(merra2.waterrate,100,72,64),2));
  profile_plots_2x2tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

  iFig = iFig + 1;
  figure(iFig); clf; 
  plotoptions.maintitle = 'Temperature Trends [K/yr]'; 
  plotoptions.cx = [-1 +1]*0.151; 
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits = [10 1000];
  miaow11 = squeeze(nanmean(reshape(umbc.ptemprate,100,72,64),2));
  miaow12 = squeeze(nanmean(permute(climcapsL3.ptemprate,[3 1 2]),2));
  miaow21 = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
  miaow22 = squeeze(nanmean(reshape(merra2.ptemprate,100,72,64),2));
  profile_plots_2x2tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

  iFig = iFig + 1;
  figure(iFig); clf; 
  plotoptions.maintitle = 'dT/dt UA'; 
  plotoptions.cx = [-1 +1]*0.151; 
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits = [0.1 100];
  miaow11 = squeeze(nanmean(reshape(umbc.ptemprate,100,72,64),2));
  miaow12 = squeeze(nanmean(permute(climcapsL3.ptemprate,[3 1 2]),2));
  miaow21 = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
  miaow22 = squeeze(nanmean(reshape(merra2.ptemprate,100,72,64),2));
  profile_plots_2x2tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

  %{
  % strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q05_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_feedback.mat';
  % driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends(strUMBC);
  %  quick_print_figs_compare_trends.m
  figure(9); aslprint('QuickFigs/Strow_May2023/Trates_UMBC_4panels.pdf')
  figure(8); aslprint('QuickFigs/Strow_May2023/WVfracrates_UMBC_4panels.pdf')  
  figure(7); aslprint('QuickFigs/Strow_May2023/RHrates_UMBC_4panels.pdf')
  figure(16); aslprint('QuickFigs/Strow_May2023/stemprates_UMBC_6panels.pdf')  
  %}

  %%%%%%%%%%%%%%%%%%%%%%%%%
  iFig = iFig + 1;
  figure(iFig); clf; 

  plotoptions.str11 = 'RH 10(UMBC-ERA5)';
  plotoptions.str12 = 'WVfrac (UMBC/ERA5-1)/100';
  plotoptions.str21 = 'T(z) UMBC-ERA5';
  plotoptions.str22 = 'Stemp UMBC-ERA5';

  plotoptions.maintitle = 'Compare UA UMBC/ERA5 or UMBC-ERA5'; 
  plotoptions.cx = [-1 +1]*0.151; 
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits = [0.1 100];
  plotoptions.yLimits = [0.1 100];

  miaowA = squeeze(nanmean(reshape(umbc.RHrate,100,72,64),2));
  miaowB = squeeze(nanmean(reshape(era5.RHrate,100,72,64),2));
  miaow11 = 10*(miaowA - miaowB);

  miaowA = squeeze(nanmean(reshape(umbc.waterrate,100,72,64),2));
  miaowB = squeeze(nanmean(reshape(era5.waterrate,100,72,64),2));
  miaow12 = (miaowA ./ miaowB - 1)/100;

  miaowA = squeeze(nanmean(reshape(umbc.ptemprate,100,72,64),2));
  miaowB = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
  miaow21 = miaowA - miaowB;

  miaow22 = umbc.stemprate - era5.stemprate;
  miaow22 = ones(100,1) * miaow22;
  miaow22 = squeeze(nanmean(reshape(miaow22,100,72,64),2));

  profile_plots_2x2tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

  %%%%%%%%%%%%%%%%%%%%%%%%%
  iFig = iFig + 1;
  figure(iFig); clf; 

  plotoptions.str11 = 'RH (UMBC-ERA5)';
  plotoptions.str12 = 'WVfrac (UMBC/ERA5-1)/100';
  plotoptions.str21 = 'T(z) UMBC-ERA5';
  plotoptions.str22 = 'Stemp UMBC-ERA5';

  plotoptions.maintitle = 'Compare LA UMBC/ERA5 or UMBC-ERA5'; 
  plotoptions.cx = [-1 +1]*0.151; 
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits = [0.1 100];
  plotoptions.yLimits = [100 1000];

  miaowA = squeeze(nanmean(reshape(umbc.RHrate,100,72,64),2));
  miaowB = squeeze(nanmean(reshape(era5.RHrate,100,72,64),2));
  miaow11 = 1*(miaowA - miaowB);

  miaowA = squeeze(nanmean(reshape(umbc.waterrate,100,72,64),2));
  miaowB = squeeze(nanmean(reshape(era5.waterrate,100,72,64),2));
  miaow12 = (miaowA ./ miaowB - 1)/100;

  miaowA = squeeze(nanmean(reshape(umbc.ptemprate,100,72,64),2));
  miaowB = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
  miaow21 = miaowA - miaowB;

  miaow22 = umbc.stemprate - era5.stemprate;
  miaow22 = ones(100,1) * miaow22;
  miaow22 = squeeze(nanmean(reshape(miaow22,100,72,64),2));

  profile_plots_2x2tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

  %%%%%
  iFig = iFig + 1;
  figure(iFig)
  miaowA = squeeze(nanmean(reshape(umbc.ptemprate,100,72,64),2));
  miaowB = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
  miaow21 = miaowA - miaowB;
  pcolor(trend_rlat64,plays100,miaow12); colormap(llsmap5); ylim([10 1000])
  title('dTz/dt : UMBC - ERA5'); shading interp; colorbar
  caxis([-1 +1]*0.05)
  set(gca,'ydir','reverse')

  %%%%%
  iFig = iFig + 1;
  figure(iFig)
  miaowA = squeeze(nanmean(reshape(umbc.waterrate,100,72,64),2));
  miaowB = squeeze(nanmean(reshape(era5.waterrate,100,72,64),2));
  miaow12 = (miaowA ./ miaowB - 1)/100;
  miaow12 = (miaowA ./ miaowB) - 1;
  pcolor(trend_rlat64,plays100,miaow12); caxis([0 2]); colormap(llsmap5); ylim([100 1000])
  miaow12 = (miaowA ./ miaowB) - 1;
  pcolor(trend_rlat64,plays100,miaow12); caxis([-1 +1]); colormap(llsmap5); ylim([100 1000])
  title('dWVfrac/dt : UMBC ./ ERA5'); shading interp; colorbar
  set(gca,'ydir','reverse')

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% do correlations

  [r,chisqr,P] = nanlinearcorrelation(umbc.stemprate,era5.stemprate);        thecorr.ST(1) = r;
  [r,chisqr,P] = nanlinearcorrelation(umbc.stemprate,merra2.stemprate);      thecorr.ST(2) = r;
  [r,chisqr,P] = nanlinearcorrelation(umbc.stemprate,airsL3.stemprate);      thecorr.ST(3) = r;
  [r,chisqr,P] = nanlinearcorrelation(umbc.stemprate,climcapsL3.stemprate);  thecorr.ST(4) = r;
  [r,chisqr,P] = nanlinearcorrelation(umbc.stemprate',giss_trend4608(:));    thecorr.ST(5) = r;
  [r,chisqr,P] = nanlinearcorrelation(era5.stemprate',giss_trend4608(:));    thecorr.ST_ERA5_GISS = r;

  iFig = iFig + 1;
  figure(iFig); plot(umbc.stemprate,era5.stemprate,'r.',umbc.stemprate,merra2.stemprate,'g.',...
                   umbc.stemprate,airsL3.stemprate,'b.',umbc.stemprate,climcapsL3.stemprate,'c.');
  axis([-1 +1 -1 +1]*0.3); line([-0.3 +0.3],[-0.3 +0.3],'color','k','linewidth',3)
  plotaxis2; hl = legend('ERA5','MERRA2','AIRSL3','CLIMCAPSL3','location','best'); title('dSKT/dt [K/yr]'); xlabel('UMBC');
 
  i100 = find(plays100 >= 100);
  for mm = 1 : length(i100)
    mmm = i100(mm);
    [r,chisqr,P] = nanlinearcorrelation(umbc.RHrate(mmm,:),era5.RHrate(mmm,:));        thecorr.RH(1,mm) = r;
    [r,chisqr,P] = nanlinearcorrelation(umbc.RHrate(mmm,:),merra2.RHrate(mmm,:));      thecorr.RH(2,mm) = r;
    junk = reshape(airsL3.RHrate,4608,100); junk = junk';
    [r,chisqr,P] = nanlinearcorrelation(umbc.RHrate(mmm,:),junk(mmm,:));         thecorr.RH(3,mm) = r;
    junk = reshape(climcapsL3.RHrate,4608,100); junk = junk';
    [r,chisqr,P] = nanlinearcorrelation(umbc.RHrate(mmm,:),junk(mmm,:));         thecorr.RH(4,mm) = r;
  end

  for mm = 1 : length(i100)
    mmm = i100(mm);
    [r,chisqr,P] = nanlinearcorrelation(umbc.waterrate(mmm,:),era5.waterrate(mmm,:));        thecorr.water(1,mm) = r;
    [r,chisqr,P] = nanlinearcorrelation(umbc.waterrate(mmm,:),merra2.waterrate(mmm,:));      thecorr.water(2,mm) = r;
    junk = reshape(airsL3.waterrate,4608,100); junk = junk';
    [r,chisqr,P] = nanlinearcorrelation(umbc.waterrate(mmm,:),junk(mmm,:));         thecorr.water(3,mm) = r;
    junk = reshape(climcapsL3.waterrate,4608,100); junk = junk';
    [r,chisqr,P] = nanlinearcorrelation(umbc.waterrate(mmm,:),junk(mmm,:));         thecorr.water(4,mm) = r;
  end

  i10 = find(plays100 >= 10);
  for mm = 1 : length(i10)
    mmm = i10(mm);
    [r,chisqr,P] = nanlinearcorrelation(umbc.ptemprate(mmm,:),era5.ptemprate(mmm,:));        thecorr.ptemp(1,mm) = r;
    [r,chisqr,P] = nanlinearcorrelation(umbc.ptemprate(mmm,:),merra2.ptemprate(mmm,:));      thecorr.ptemp(2,mm) = r;
    junk = reshape(permute(airsL3.ptemprate,[3 1 2]),100,4608); 
    [r,chisqr,P] = nanlinearcorrelation(umbc.ptemprate(mmm,:),junk(mmm,:));         thecorr.ptemp(3,mm) = r;
    junk = reshape(climcapsL3.ptemprate,4608,100); junk = junk';
    [r,chisqr,P] = nanlinearcorrelation(umbc.ptemprate(mmm,:),junk(mmm,:));         thecorr.ptemp(4,mm) = r;
  end
 
  iFig = iFig + 1;
  figure(iFig); clf;
   fprintf(1,'ST correlations vs UMBC : ERA5/MERRA2/AIRSL3/CLIMCAPSL3/GISS =  %8.6f %8.6f %8.6f %8.6f %8.6f\n',thecorr.ST);
   fprintf(1,'ST correlations         : ERA5 vs GISS =  %8.6f\n',thecorr.ST_ERA5_GISS);
   subplot(131); plot(thecorr.RH',plays100(i100),'linewidth',2); plotaxis2; hl = legend('ERA5','MERRA2','AIRSL3','CLIMCAPSL3','location','best','fontsize',8); set(gca,'ydir','reverse'); title('RH corr');
   subplot(132); plot(thecorr.water',plays100(i100),'linewidth',2); plotaxis2; hl = legend('ERA5','MERRA2','AIRSL3','CLIMCAPSL3','location','best','fontsize',8); set(gca,'ydir','reverse'); title('water corr');
   subplot(133); plot(thecorr.ptemp',plays100(i10),'linewidth',2); plotaxis2; hl = legend('ERA5','MERRA2','AIRSL3','CLIMCAPSL3','location','best','fontsize',8); set(gca,'ydir','reverse'); title('ptemp corr');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  i200 = find(plays100 >= 200,1);
  i500 = find(plays100 >= 500,1);
  i800 = find(plays100 >= 800,1);

  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.RHrate(:,:,i500);
  zyx12 = climcapsL3.RHrate(:,:,i500);
  zyx21 = reshape(merra2.RHrate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i500);
  zyx22 = reshape(era5.RHrate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i500);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = '500 mb dRH/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';
  plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';
  plotoptions.str22 = 'ERA5';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  aslmap_2x2tiledlayout(zyx11,zyx12,zyx21,zyx22,iFig,plotoptions);
  
  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.waterrate(:,:,i500);
  zyx12 = climcapsL3.waterrate(:,:,i500);
  zyx21 = reshape(merra2.waterrate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i500);
  zyx22 = reshape(era5.waterrate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i500);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.0151; plotoptions.maintitle = '500 mb dWVfrac/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';
  plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';
  plotoptions.str22 = 'ERA5';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  aslmap_2x2tiledlayout(zyx11,zyx12,zyx21,zyx22,iFig,plotoptions);
  
  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.ptemprate(:,:,i500);
  zyx12 = climcapsL3.ptemprate(:,:,i500);
  zyx21 = reshape(merra2.ptemprate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i500);
  zyx22 = reshape(era5.ptemprate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i500);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = '500 mb dT/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';
  plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';
  plotoptions.str22 = 'ERA5';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  aslmap_2x2tiledlayout(zyx11,zyx12,zyx21,zyx22,iFig,plotoptions);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % z11 = AIRS, z12 = climcaps   z21 = MERRA2     z22 = ERA5      z11x = umbc   z32 = giss

  iFig = iFig + 1;
  figure(iFig); clf;
  z32 = giss_trend4608(:); z32 = z32';
  clear plotoptions
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.str31 = 'UMBC';      plotoptions.str32 = 'GISS';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  z31 = z11x; z32 = z32;
  aslmap_3x2tiledlayout(z11,z12,z21,z22,z31,z32,iFig,plotoptions);

  figure(20); clf; aslmap(20,rlat65,rlon73,smoothn(reshape(z31,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.151); %% UMBC
  figure(21); clf; aslmap(21,rlat65,rlon73,smoothn(reshape(z32,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.151); %% GISS
  
  %%%%%%%%%% 

  clear plotoptions
  plotoptions.cx = [-1 +1]*0.151;  plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3'; plotoptions.str13 = 'MERRA2';    
  plotoptions.str21 = 'ERA5';      plotoptions.str22 = 'UMBC';        plotoptions.str23 = 'GISS';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  z31 = z11x; z32 = z32;
  aslmap_2x3tiledlayout(z11,z12,z21,z22,z31,z32,iFig,plotoptions);

  %%%%%%%%%% 

  clear plotoptions
  plotoptions.cx = [-1 +1]*0.151;  plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';       plotoptions.str12 = 'ERA5';          plotoptions.str13 = 'UMBC';    
  plotoptions.str21 = 'CLIMCAPS L3';   plotoptions.str22 = 'MERRA2';        plotoptions.str23 = 'GISS';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  z31 = z11x; z32 = z32;
  aslmap_2x3tiledlayout(z11,z22,z31,z12,z21,z32,iFig,plotoptions);

  %%%%%%%%%% 

  disp('SKT trend K/yr  | AIRS   CLIMCAPS   MERRA2      ERA5    UMBC     GISS');
  disp('----------------|-------------------------------------------------------');
  junk = [nanmean(z11)  nanmean(z12) nanmean(z21)  nanmean(z22) nanmean(z31)  nanmean(z32)];
  fprintf(1,'(ALL  )         |  %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n',junk);
  boo = find(landfrac == 0);  
  junk = [nanmean(z11(boo))  nanmean(z12(boo)) nanmean(z21(boo))  nanmean(z22(boo)) nanmean(z31(boo))  nanmean(z32(boo))];
  fprintf(1,'(OCEAN)         |  %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n',junk);
  boo = find(landfrac > 0.99);  
  junk = [nanmean(z11(boo))  nanmean(z12(boo)) nanmean(z21(boo))  nanmean(z22(boo)) nanmean(z31(boo))  nanmean(z32(boo))];
  fprintf(1,'(LAND )         |  %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n',junk);
  disp('----------------|-------------------------------------------------------');
  
  disp('SKT trend K/yr  | AIRS   CLIMCAPS   MERRA2      ERA5    UMBC     GISS');
  disp('(area weighted with cos(lat)')
  disp('----------------|--------------------------------------------------------');
  junk = [nansum(z11.*mu)/nansum(mu)  nansum(z12.*mu)/nansum(mu) nansum(z21.*mu)/nansum(mu)  ...
          nansum(z22.*mu)/nansum(mu) nansum(z31.*mu)/nansum(mu)  nansum(z32.*mu)/nansum(mu)];
  fprintf(1,'(ALL  )         |  %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n',junk);
  boo = find(landfrac == 0);  
  junk = [nansum(z11(boo).*mu(boo))/nansum(mu(boo))  nansum(z12(boo).*mu(boo))/nansum(mu(boo)) ...
          nansum(z21(boo).*mu(boo))/nansum(mu(boo))  nansum(z22(boo).*mu(boo))/nansum(mu(boo)) ...
          nansum(z31(boo).*mu(boo))/nansum(mu(boo))  nansum(z32(boo).*mu(boo))/nansum(mu(boo))];
  fprintf(1,'(OCEAN)         |  %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n',junk);
  boo = find(landfrac > 0.99);  
  junk = [nansum(z11(boo).*mu(boo))/nansum(mu(boo))  nansum(z12(boo).*mu(boo))/nansum(mu(boo)) ...
          nansum(z21(boo).*mu(boo))/nansum(mu(boo))  nansum(z22(boo).*mu(boo))/nansum(mu(boo))...
          nansum(z31(boo).*mu(boo))/nansum(mu(boo))  nansum(z32(boo).*mu(boo))/nansum(mu(boo))];
  fprintf(1,'(LAND )         |  %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n',junk);
  disp('----------------|--------------------------------------------------------');

  %%%%%%%%%%%%%%%%%%%%%%%%%

  [p.salti, p.landfrac] = usgs_deg10_dem(p.rlat, p.rlon);
  land = find(p.landfrac >= eps);              junk = zeros(4608,1); junk(land) = 1; land = junk;
  ocean = find(p.landfrac < eps);              junk = zeros(4608,1); junk(ocean) = 1; ocean = junk;
  tropics = find(abs(p.rlat) <= 30);          
  midlatsNtropics = find(abs(p.rlat) <= 60); 
  midlats = setdiff(midlatsNtropics,tropics); 
  poles = find(abs(p.rlat) > 60);             

  junk = zeros(4608,1); junk(tropics) = 1;          tropics = junk;
  junk = zeros(4608,1); junk(midlatsNtropics) = 1;  midlatsNtropics = junk;
  junk = zeros(4608,1); junk(midlats)         = 1;  midlats = junk;
  junk = zeros(4608,1); junk(poles) = 1;            poles = junk;  

  junk = [sum(land) sum(ocean) sum(tropics) sum(midlats) sum(midlatsNtropics) sum(poles)];
  fprintf(1,'breakdown : land = %4i ocean = %4i      tropics = %4i midlats = %4i midlatsNtropics = %4i poles = %4i of 4608 \n',junk);
  fprintf(1,'breakdown : land = %4.3f ocean = %4.3f    tropics = %4.3f midlats = %4.3f midlatsNtropics = %4.3f poles = %4.3f of 4608 \n',junk/4608);

  i200 = find(plays100 >= 200,1);
  i500 = find(plays100 >= 500,1);
  i800 = find(plays100 >= 800,1);

  iX = 0;
  iY = 0;

  iY = iY+1; iX = 0;
  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.RHrate(:,:,i200);
  zyx12 = climcapsL3.RHrate(:,:,i200);
  zyx21 = reshape(merra2.RHrate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i200);
  zyx22 = reshape(era5.RHrate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i200);
  zyx31 = reshape(umbc.RHrate,100,72,64);    zyx31 = permute(zyx31,[2 3 1]); zyx31 = zyx31(:,:,i200);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.251; plotoptions.maintitle = '200 mb dRH/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.str31 = 'UMBC';      plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'UMBC';      zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dRH/dt 200 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx11(:));  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx12(:));  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx21(:));  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx31(:));  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dRH/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx11(:).*tropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx12(:).*tropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx21(:).*tropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx31(:).*tropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dRH/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx11(:).*midlats);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx12(:).*midlats);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx21(:).*midlats);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx31(:).*midlats);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dRH/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx11(:).*midlatsNtropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx12(:).*midlatsNtropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx21(:).*midlatsNtropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx31(:).*midlatsNtropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dRH/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx11(:).*poles);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx12(:).*poles);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx21(:).*poles);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx31(:).*poles);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dRH/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n\n',thecorr.junk);
  
  iY = iY+1; iX = 0;
  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.waterrate(:,:,i200);
  zyx12 = climcapsL3.waterrate(:,:,i200);
  zyx21 = reshape(merra2.waterrate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i200);
  zyx22 = reshape(era5.waterrate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i200);
  zyx31 = reshape(umbc.waterrate,100,72,64);    zyx31 = permute(zyx31,[2 3 1]); zyx31 = zyx31(:,:,i200);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.0151; plotoptions.maintitle = '200 mb dWVfrac/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.str31 = 'UMBC';      plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'UMBC';      zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dWVfrac/dt 200 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx11(:));  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx12(:));  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx21(:));  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx31(:));  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dWVfrac/dt                   correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk); 
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx11(:).*tropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx12(:).*tropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx21(:).*tropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx31(:).*tropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dWVfrac/dt   TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx11(:).*midlats);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx12(:).*midlats);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx21(:).*midlats);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx31(:).*midlats);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dWVfrac/dt   MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx11(:).*midlatsNtropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx12(:).*midlatsNtropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx21(:).*midlatsNtropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx31(:).*midlatsNtropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dWVfrac/dt   MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx11(:).*poles);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx12(:).*poles);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx21(:).*poles);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx31(:).*poles);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dWVfrac/dt   POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n\n',thecorr.junk);
 
  iY = iY+1; iX = 0;
  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.ptemprate(:,:,i200);
  zyx12 = climcapsL3.ptemprate(:,:,i200);
  zyx21 = reshape(merra2.ptemprate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i200);
  zyx22 = reshape(era5.ptemprate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i200);
  zyx31 = reshape(umbc.ptemprate,100,72,64);    zyx31 = permute(zyx31,[2 3 1]); zyx31 = zyx31(:,:,i200);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = '200 mb dT/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.str31 = 'UMBC';      plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'UMBC';      zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dT/dt 200 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx11(:));  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx12(:));  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx21(:));  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx31(:));  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dT/dt                        correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx11(:).*tropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx12(:).*tropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx21(:).*tropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx31(:).*tropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dT/dt        TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx11(:).*midlats);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx12(:).*midlats);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx21(:).*midlats);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx31(:).*midlats);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dT/dt        MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx11(:).*midlatsNtropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx12(:).*midlatsNtropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx21(:).*midlatsNtropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx31(:).*midlatsNtropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dT/dt        MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx11(:).*poles);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx12(:).*poles);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx21(:).*poles);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx31(:).*poles);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'200 mb dT/dt        POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n\n',thecorr.junk);

  %%%%%%%%%%
  iY = iY+1; iX = 0;
  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.RHrate(:,:,i500);
  zyx12 = climcapsL3.RHrate(:,:,i500);
  zyx21 = reshape(merra2.RHrate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i500);
  zyx22 = reshape(era5.RHrate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i500);
  zyx31 = reshape(umbc.RHrate,100,72,64);    zyx31 = permute(zyx31,[2 3 1]); zyx31 = zyx31(:,:,i500);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.251; plotoptions.maintitle = '500 mb dRH/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.str31 = 'UMBC';      plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'UMBC';      zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dRH/dt 500 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx11(:));  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx12(:));  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx21(:));  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx31(:));  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dRH/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx11(:).*tropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx12(:).*tropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx21(:).*tropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx31(:).*tropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dRH/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx11(:).*midlats);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx12(:).*midlats);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx21(:).*midlats);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx31(:).*midlats);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dRH/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx11(:).*midlatsNtropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx12(:).*midlatsNtropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx21(:).*midlatsNtropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx31(:).*midlatsNtropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dRH/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx11(:).*poles);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx12(:).*poles);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx21(:).*poles);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx31(:).*poles);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dRH/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n\n',thecorr.junk);
  
  iY = iY+1; iX = 0;
  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.waterrate(:,:,i500);
  zyx12 = climcapsL3.waterrate(:,:,i500);
  zyx21 = reshape(merra2.waterrate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i500);
  zyx22 = reshape(era5.waterrate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i500);
  zyx31 = reshape(umbc.waterrate,100,72,64);    zyx31 = permute(zyx31,[2 3 1]); zyx31 = zyx31(:,:,i500);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.0151; plotoptions.maintitle = '500 mb dWVfrac/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.str31 = 'UMBC';      plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'UMBC';      zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dWVfrac/dt 500 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx11(:));  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx12(:));  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx21(:));  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx31(:));  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dWVfrac/dt                   correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx11(:).*tropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx12(:).*tropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx21(:).*tropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx31(:).*tropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dWVfrac/dt   TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx11(:).*midlats);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx12(:).*midlats);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx21(:).*midlats);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx31(:).*midlats);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dWVfrac/dt   MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx11(:).*midlatsNtropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx12(:).*midlatsNtropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx21(:).*midlatsNtropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx31(:).*midlatsNtropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dWVfrac/dt   MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx11(:).*poles);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx12(:).*poles);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx21(:).*poles);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx31(:).*poles);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dWVfrac/dt   POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n\n',thecorr.junk);

  iY = iY+1; iX = 0;
  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.ptemprate(:,:,i500);
  zyx12 = climcapsL3.ptemprate(:,:,i500);
  zyx21 = reshape(merra2.ptemprate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i500);
  zyx22 = reshape(era5.ptemprate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i500);
  zyx31 = reshape(umbc.ptemprate,100,72,64);    zyx31 = permute(zyx31,[2 3 1]); zyx31 = zyx31(:,:,i500);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = '500 mb dT/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.str31 = 'UMBC';      plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'UMBC';      zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dT/dt 500 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx11(:));  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx12(:));  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx21(:));  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx31(:));  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dT/dt                        correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx11(:).*tropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx12(:).*tropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx21(:).*tropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx31(:).*tropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dT/dt        TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx11(:).*midlats);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx12(:).*midlats);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx21(:).*midlats);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx31(:).*midlats);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dT/dt        MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx11(:).*midlatsNtropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx12(:).*midlatsNtropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx21(:).*midlatsNtropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx31(:).*midlatsNtropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dT/dt        MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx11(:).*poles);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx12(:).*poles);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx21(:).*poles);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx31(:).*poles);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'500 mb dT/dt        POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n\n',thecorr.junk);

  %%%%%%%%%%
  iY = iY+1; iX = 0;
  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.RHrate(:,:,i800);
  zyx12 = climcapsL3.RHrate(:,:,i800);
  zyx21 = reshape(merra2.RHrate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i800);
  zyx22 = reshape(era5.RHrate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i800);
  zyx31 = reshape(umbc.RHrate,100,72,64);    zyx31 = permute(zyx31,[2 3 1]); zyx31 = zyx31(:,:,i800);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.251; plotoptions.maintitle = '800 mb dRH/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.str31 = 'UMBC';      plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'UMBC';      zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dRH/dt 800 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx11(:));  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx12(:));  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx21(:));  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx31(:));  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dRH/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx11(:).*tropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx12(:).*tropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx21(:).*tropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx31(:).*tropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dRH/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx11(:).*midlats);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx12(:).*midlats);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx21(:).*midlats);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx31(:).*midlats);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dRH/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx11(:).*midlatsNtropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx12(:).*midlatsNtropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx21(:).*midlatsNtropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx31(:).*midlatsNtropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dRH/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx11(:).*poles);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx12(:).*poles);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx21(:).*poles);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx31(:).*poles);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dRH/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n\n',thecorr.junk);

  iY = iY+1; iX = 0;
  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.waterrate(:,:,i800);
  zyx12 = climcapsL3.waterrate(:,:,i800);
  zyx21 = reshape(merra2.waterrate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i800);
  zyx22 = reshape(era5.waterrate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i800);
  zyx31 = reshape(umbc.waterrate,100,72,64);    zyx31 = permute(zyx31,[2 3 1]); zyx31 = zyx31(:,:,i800);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.0151; plotoptions.maintitle = '800 mb dWVfrac/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.str31 = 'UMBC';      plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'UMBC';      zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dWVfrac/dt 500 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx11(:));  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx12(:));  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx21(:));  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx31(:));  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dWVfrac/dt                   correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx11(:).*tropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx12(:).*tropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx21(:).*tropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx31(:).*tropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dWVfrac/dt   TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx11(:).*midlats);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx12(:).*midlats);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx21(:).*midlats);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx31(:).*midlats);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dWVfrac/dt   MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx11(:).*midlatsNtropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx12(:).*midlatsNtropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx21(:).*midlatsNtropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx31(:).*midlatsNtropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dWVfrac/dt   MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx11(:).*poles);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx12(:).*poles);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx21(:).*poles);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx31(:).*poles);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dWVfrac/dt   POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n\n',thecorr.junk);
  
  iY = iY+1; iX = 0;
  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  zyx11 = airsL3.ptemprate(:,:,i800);
  zyx12 = climcapsL3.ptemprate(:,:,i800);
  zyx21 = reshape(merra2.ptemprate,100,72,64);  zyx21 = permute(zyx21,[2 3 1]); zyx21 = zyx21(:,:,i800);
  zyx22 = reshape(era5.ptemprate,100,72,64);    zyx22 = permute(zyx22,[2 3 1]); zyx22 = zyx22(:,:,i800);
  zyx31 = reshape(umbc.ptemprate,100,72,64);    zyx31 = permute(zyx31,[2 3 1]); zyx31 = zyx31(:,:,i800);
  maskLF = ones(size(zyx22));
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = '800 mb dT/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.str31 = 'UMBC';      plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'UMBC';      zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dT/dt 800 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx11(:));  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx12(:));  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx21(:));  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:),zyx31(:));  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dT/dt                        correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx11(:).*tropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx12(:).*tropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx21(:).*tropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*tropics,zyx31(:).*tropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dT/dt        TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx11(:).*midlats);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx12(:).*midlats);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx21(:).*midlats);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlats,zyx31(:).*midlats);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dT/dt        MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx11(:).*midlatsNtropics);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx12(:).*midlatsNtropics);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx21(:).*midlatsNtropics);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*midlatsNtropics,zyx31(:).*midlatsNtropics);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dT/dt        MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n',thecorr.junk);
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx11(:).*poles);  thecorr.junk(1) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx12(:).*poles);  thecorr.junk(2) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx21(:).*poles);  thecorr.junk(3) = r;
    [r,chisqr,P] = nanlinearcorrelation(zyx22(:).*poles,zyx31(:).*poles);  thecorr.junk(4) = r;
    iX = iX + 1; allXchi(iY,iX,:) = thecorr.junk;
    fprintf(1,'800 mb dT/dt        POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.6f %8.6f %8.6f %8.6f\n\n',thecorr.junk);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% iX = 1,2,3,4,5 === all, tropics,midlats,midlats+tropics,poles  = 5 total
  %% iY = 1,2,3  200 mb    4,5,6 500 mb    7,8,9 800 mb             = 9 total
  %% whos allXchi = 9 x 5 x 4 

    %iFig = iFig + 1;
    %figure(iFig); clf;   
    %h1=subplot(311);  bar(1:5,squeeze(allXchi(1,:,:))); ylim([-1 1]); hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','best','fontsize',6); ylabel('RH 200')
    %h2=subplot(312);  bar(1:5,squeeze(allXchi(4,:,:))); ylim([-1 1]); hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','best','fontsize',6); ylabel('RH 500')
    %h3=subplot(313);  bar(1:5,squeeze(allXchi(7,:,:))); ylim([-1 1]); hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','best','fontsize',6); ylabel('RH 800')
    %adjust(31,h1,h2,h3,'even');
    %set(gca,'xtick',[1:5],'xticklabel',names)

    names = {'ALL','T','ML','ML+T','P'};

    iFig = iFig + 1;
    figure(iFig); clf;   
    ta = tiledlayout(3,1,'TileSpacing','None', 'Padding','None');
    ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
    ymin = 0;
    tfov(1) = nexttile;  bar(1:5,squeeze(allXchi(1,:,:))); ylim([ymin 1]); ylabel('RH 200');  %hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','best','fontsize',6); ylabel('RH 200')
    tfov(2) = nexttile;  bar(1:5,squeeze(allXchi(4,:,:))); ylim([ymin 1]); ylabel('RH 500');  %hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','best','fontsize',6); ylabel('RH 500')
    tfov(3) = nexttile;  bar(1:5,squeeze(allXchi(7,:,:))); ylim([ymin 1]); ylabel('RH 800');  hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','best','fontsize',6); ylabel('RH 800')
    tfov(1).XTickLabel = '';  tfov(1).XLabel.String = [];      tfov(2).XTickLabel = '';  tfov(2).XLabel.String = [];
    set(gca,'xtick',[1:5],'xticklabel',names)
    ta.Padding = 'none';
    ta.TileSpacing = 'tight';

    iFig = iFig + 1;
    figure(iFig); clf;   
    ta = tiledlayout(3,1,'TileSpacing','None', 'Padding','None');
    ta.OuterPosition = [0.0375 0.0375 0.925 0.925];    
    ymin = -0.25;
    ymin2 = 0;
    tfov(1) = nexttile;  bar(1:5,squeeze(allXchi(2,:,:))); ylim([ymin  1]); ylabel('WV 200'); %hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','s','fontsize',5); 
    tfov(2) = nexttile;  bar(1:5,squeeze(allXchi(5,:,:))); ylim([ymin2 1]); ylabel('WV 500'); %hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','best','fontsize',6);
    tfov(3) = nexttile;  bar(1:5,squeeze(allXchi(8,:,:))); ylim([ymin2 1]); ylabel('WV 800'); %hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','ne','fontsize',5);
    tfov(1).XTickLabel = '';  tfov(1).XLabel.String = [];      tfov(2).XTickLabel = '';  tfov(2).XLabel.String = [];
    set(gca,'xtick',[1:5],'xticklabel',names)
    ta.Padding = 'none';
    ta.TileSpacing = 'tight';

    iFig = iFig + 1;
    figure(iFig); clf;   
    ta = tiledlayout(3,1,'TileSpacing','None', 'Padding','None');
    ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
    ymin = 0;
    tfov(1) = nexttile;  bar(1:5,squeeze(allXchi(3,:,:))); ylim([ymin 1]); ylabel('T 200');  %hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','best','fontsize',6); ylabel('T 200')
    tfov(2) = nexttile;  bar(1:5,squeeze(allXchi(6,:,:))); ylim([ymin 1]); ylabel('T 500');  %hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','best','fontsize',6); ylabel('T 500')
    tfov(3) = nexttile;  bar(1:5,squeeze(allXchi(9,:,:))); ylim([ymin 1]); ylabel('T 800');  %hl = legend('AIRS v7','CLIMCAPS','MERRA2','OEM','location','best','fontsize',6); ylabel('T 800')
    tfov(1).XTickLabel = '';  tfov(1).XLabel.String = [];      tfov(2).XTickLabel = '';  tfov(2).XLabel.String = [];
    set(gca,'xtick',[1:5],'xticklabel',names)
    ta.Padding = 'none';
    ta.TileSpacing = 'tight';

keyboard_nowindow
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  iY = input('Show 5 panel plots? (-1 default/+1) : ');
  if length(iY) == 0
    iY = -1;
  end

  clear plotoptions
  plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';
  plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str12 = 'CLIMCAPS';
  plotoptions.str13 = 'MERRA2';
  plotoptions.str14 = 'ERA5';
  plotoptions.str15 = 'UMBC';
  plotoptions.xstr = 'Latitude';        plotoptions.ystr = 'Pressure (mb)';
  plotoptions.yLinearOrLog = +1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits = [100 1000];  
   
  if iY > 0
    iFig = iFig + 1;
    figure(iFig); clf;   
    plotoptions.maintitle = 'RH trends [%/yr]'; 
    plotoptions.cx = [-1 +1]*0.5; 
    miaow11 = squeeze(nanmean(permute(airsL3.RHrate,[3 1 2]),2));
    miaow12 = squeeze(nanmean(permute(climcapsL3.RHrate,[3 1 2]),2));
    miaow13 = squeeze(nanmean(reshape(merra2.RHrate,100,72,64),2));
    miaow14 = squeeze(nanmean(reshape(era5.RHrate,100,72,64),2));
    miaow15 = squeeze(nanmean(reshape(umbc.RHrate,100,72,64),2));
    if iPentagonPlot < 0
      profile_plots_1x5tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow13,miaow14,miaow15,iFig,plotoptions);
    else
      plotoptions2x1x2 = plotoptions;
      plotoptions2x1x2.str11 = 'AIRS L3';
      plotoptions2x1x2.str12 = 'UMBC';
      plotoptions2x1x2.str13 = 'CLIMCAPS';
      plotoptions2x1x2.str21 = 'MERRA2';
      plotoptions2x1x2.str22 = 'ERA5';
      plotoptions2x1x2.cstr  = 'dRH/dt';
      profile_plots_2x1x2tiledlayout_tall(trend_rlat64,plays100,miaow11,miaow15,miaow12,miaow13,miaow14,iFig,plotoptions2x1x2);
    end

    iFig = iFig + 1;
    figure(iFig); clf;   
    plotoptions.maintitle = 'd(log(WV))/dt [1/yr]'; 
    plotoptions.cx = [-1 +1]*0.015; 
    miaow11 = squeeze(nanmean(permute(airsL3.waterrate,[3 1 2]),2));
    miaow12 = squeeze(nanmean(permute(climcapsL3.waterrate,[3 1 2]),2));
    miaow13 = squeeze(nanmean(reshape(merra2.waterrate,100,72,64),2));
    miaow14 = squeeze(nanmean(reshape(era5.waterrate,100,72,64),2));
    miaow15 = squeeze(nanmean(reshape(umbc.waterrate,100,72,64),2));
    if iPentagonPlot < 0
      profile_plots_1x5tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow13,miaow14,miaow15,iFig,plotoptions);
    else
      plotoptions2x1x2 = plotoptions;
      plotoptions2x1x2.str11 = 'AIRS L3';
      plotoptions2x1x2.str12 = 'UMBC';
      plotoptions2x1x2.str13 = 'CLIMCAPS';
      plotoptions2x1x2.str21 = 'MERRA2';
      plotoptions2x1x2.str22 = 'ERA5';
      plotoptions2x1x2.cstr  = 'dWVfrac/dt';
      profile_plots_2x1x2tiledlayout_tall(trend_rlat64,plays100,miaow11,miaow15,miaow12,miaow13,miaow14,iFig,plotoptions2x1x2);
    end
    
    iFig = iFig + 1;
    figure(iFig); clf;
    plotoptions.yLinearOrLog = -1;
    plotoptions.maintitle = 'Temperature Trends [K/yr]';
    plotoptions.cx = [-1 +1]*0.151;
    plotoptions.yLimits = [10 1000];  
    miaow11 = squeeze(nanmean(permute(airsL3.ptemprate,[3 1 2]),2));
    miaow12 = squeeze(nanmean(permute(climcapsL3.ptemprate,[3 1 2]),2));
    miaow13 = squeeze(nanmean(reshape(merra2.ptemprate,100,72,64),2));
    miaow14 = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
    miaow15 = squeeze(nanmean(reshape(umbc.ptemprate,100,72,64),2));
    if iPentagonPlot < 0
      profile_plots_1x5tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow13,miaow14,miaow15,iFig,plotoptions);
    else
      plotoptions2x1x2 = plotoptions;
      plotoptions2x1x2.str11 = 'AIRS L3';
      plotoptions2x1x2.str12 = 'UMBC';
      plotoptions2x1x2.str13 = 'CLIMCAPS';
      plotoptions2x1x2.str21 = 'MERRA2';
      plotoptions2x1x2.str22 = 'ERA5';
      plotoptions2x1x2.cstr  = 'dT/dt';
      profile_plots_2x1x2tiledlayout_tall(trend_rlat64,plays100,miaow11,miaow15,miaow12,miaow13,miaow14,iFig,plotoptions2x1x2);
    end
  
    %{ 
    %  quick_print_figs_compare_trends.m
    addpath /asl/matlib/plotutils
    figure(16); aslprint('QuickFigs/Strow_March2023/skt_trend.pdf');
    figure(17); aslprint('QuickFigs/Strow_March2023/rh_trend.pdf');
    figure(18); aslprint('QuickFigs/Strow_March2023/wv_trend.pdf');
    figure(19); aslprint('QuickFigs/Strow_March2023/tz_trend.pdf');
    %}
  else
    fprintf(1,'currently at iFig = %2i \n',iFig);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  iFig = iFig + 1;
  figure(iFig); clf;
  clear plotoptions
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'GISS';    
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';    
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  aslmap_2x2tiledlayout(z11,z32,z21,z22,iFig,plotoptions);

  iY = input('Show 3 panel plots? (-1 default/+1) : ');
  if length(iY) == 0
    iY = -1;
  end

  clear plotoptions
  plotoptions.plotcolors = llsmap5;
  plotoptions.str1 = 'AIRS L3';
  plotoptions.str2 = 'MERRA2';
  plotoptions.str3 = 'ERA5';
  plotoptions.xstr = 'Latitude';        plotoptions.ystr = 'Pressure (mb)';
  plotoptions.yLinearOrLog = +1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits = [100 1000];  

  if iY > 0
    iFig = iFig + 1;
    figure(iFig); clf;   
    plotoptions.maintitle = 'RH trends [%/yr]'; 
    plotoptions.cx = [-1 +1]*0.5; 
    miaow11 = squeeze(nanmean(permute(airsL3.RHrate,[3 1 2]),2));
    miaow12 = squeeze(nanmean(reshape(merra2.RHrate,100,72,64),2));
    miaow13 = squeeze(nanmean(reshape(era5.RHrate,100,72,64),2));
    profile_plots_1x3tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow13,iFig,plotoptions);
    
    iFig = iFig + 1;
    figure(iFig); clf;   
    plotoptions.maintitle = 'd(log(WV))/dt [1/yr]'; 
    plotoptions.cx = [-1 +1]*0.015; 
    miaow11 = squeeze(nanmean(permute(airsL3.waterrate,[3 1 2]),2));
    miaow12 = squeeze(nanmean(reshape(merra2.waterrate,100,72,64),2));
    miaow13 = squeeze(nanmean(reshape(era5.waterrate,100,72,64),2));
    profile_plots_1x3tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow13,iFig,plotoptions);
    
    iFig = iFig + 1;
    figure(iFig); clf;
    plotoptions.yLinearOrLog = -1;
    plotoptions.maintitle = 'Temperature Trends [K/yr]';
    plotoptions.cx = [-1 +1]*0.151;
    plotoptions.yLimits = [10 1000];  
    miaow11 = squeeze(nanmean(permute(airsL3.ptemprate,[3 1 2]),2));
    miaow12 = squeeze(nanmean(reshape(merra2.ptemprate,100,72,64),2));
    miaow13 = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
    profile_plots_1x3tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow13,iFig,plotoptions);
  
    %{ 
    %  quick_print_figs_compare_trends.m
    addpath /asl/matlib/plotutils
    figure(20); aslprint('QuickFigs/Strow_March2023/skt_trend_3panel.pdf');
    figure(21); aslprint('QuickFigs/Strow_March2023/rh_trend_3panel.pdf');
    figure(22); aslprint('QuickFigs/Strow_March2023/wv_trend_3panel.pdf');
    figure(23); aslprint('QuickFigs/Strow_March2023/tz_trend_3panel.pdf');
    %}

  else
    fprintf(1,'currently at iFig = %2i \n',iFig);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iFigiAmp = iFig;

iAmp = input('Do amplification??? (-9999 to stop [default]/ 0 to go on) : ');
if length(iAmp) == 0
  iAmp = -9999;
end

if iAmp >= -10
  disp('Enter (-7) polar    L/O')
  disp('      (-6) midlat   L/O')
  disp('      (-5) tropical L/O')
  disp('      (-4) polar land    (+4) polar ocean')
  disp('      (-3) midlat land   (+3) midlat ocean')
  disp('      (-2) tropical land (+2) tropical ocean')
  disp('      (-1) land          (+1) ocean');
  disp('      [0,default] ALL trends : ');
  iAorOorL = input('Enter region : ');
  if length(iAorOorL) == 0 
    iAorOorL = 0;
  end

  while length(intersect(iAorOorL,[-7 -6 -5 -4 -3 -2 -1 0 +1 +2 +3 +4]) == 1)
    driver_show_atmospheric_amplification
  
    disp('Enter (-7) polar    L/O')
    disp('      (-6) midlat   L/O')
    disp('      (-5) tropical L/O')
    disp('      (-4) polar land    (+4) polar ocean')
    disp('      (-3) midlat land   (+3) midlat ocean')
    disp('      (-2) tropical land (+2) tropical ocean')
    disp('      (-1) land          (+1) ocean');
    disp('      [0,default] ALL trends : ');
    iAorOorL = input('Enter region (-9999 to stop/0 default) : ');
    if length(iAorOorL) == 0 
      iAorOorL = 0;
    end
  end
else
  fprintf(1,'currently at iFig = %2i \n',iFig);
end
fprintf(1,'currently at iFig = %2i \n',iFig);

if iUMBC > 0
  iSpectra = input('Show spectral trends : (-1 /+1 yes [default]) : ');
  if length(iSpectra) == 0
    iSpectra = +1;
  end
  if iSpectra > 0
    %% umbcXspectra = load(strUMBC,'fits','rates');
    umbcXspectra = load(strUMBC,'fits');
    driver_get_the_model_trends
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPrint = -1;
if iPrint > 0
  dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/PAPER17_TRENDS/Figs/';
  figure(60); aslprint([dir0 'obs_spectralavg_'        num2str(iNumYears) '_years.pdf']);
  figure(61); aslprint([dir0 'umbc_spectralavg_'       num2str(iNumYears) '_years.pdf']);
  figure(62); aslprint([dir0 'era5_spectralavg_'       num2str(iNumYears) '_years.pdf']);
  figure(63); aslprint([dir0 'merra2_spectralavg_'     num2str(iNumYears) '_years.pdf']);
  figure(64); aslprint([dir0 'airsL3_spectralavg_'     num2str(iNumYears) '_years.pdf']);
  figure(65); aslprint([dir0 'climcapsL3_spectralavg_' num2str(iNumYears) '_years.pdf']);

%%%% this is ughugh
  dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/PAPER17_TRENDS/Figs/';

%%  figure(07); aslprint([dir0 'rh_trends_zonal_p_4panels' num2str(iNumYears) '_years.pdf']);
%%  figure(08); aslprint([dir0 'fracWV_trends_zonal_p_4panels' num2str(iNumYears) '_years.pdf']);
%%  figure(09); aslprint([dir0 'tz_trends_zonal_p_4panels' num2str(iNumYears) '_years.pdf']);

%%  figure(21); aslprint([dir0 '200mb_rhtrends_lat_lon_6panels' num2str(iNumYears) '_years.pdf']);
%%  figure(22); aslprint([dir0 '200mb_fracWVtrends_lat_lon_6panels' num2str(iNumYears) '_years.pdf']);
%%  figure(23); aslprint([dir0 '200mb_tztrends_lat_lon_6panels' num2str(iNumYears) '_years.pdf']);

%%  figure(24); aslprint([dir0 '500mb_rhtrends_lat_lon_6panels' num2str(iNumYears) '_years.pdf']);
%%  figure(25); aslprint([dir0 '500mb_fracWVtrends_lat_lon_6panels' num2str(iNumYears) '_years.pdf']);
%%  figure(26); aslprint([dir0 '500mb_tztrends_lat_lon_6panels' num2str(iNumYears) '_years.pdf']);

%%  figure(27); aslprint([dir0 '800mb_rhtrends_lat_lon_6panels' num2str(iNumYears) '_years.pdf']);
%%  figure(28); aslprint([dir0 '800mb_fracWVtrends_lat_lon_6panels' num2str(iNumYears) '_years.pdf']);
%%  figure(29); aslprint([dir0 '800mb_tztrends_lat_lon_6panels' num2str(iNumYears) '_years.pdf']);

  figure(30); aslprint([dir0 'rh_correlations_regions_models' num2str(iNumYears) '_years.pdf']);
  figure(31); aslprint([dir0 'fracWV_correlations_regions_models' num2str(iNumYears) '_years.pdf']);
  figure(32); aslprint([dir0 'tz_correlations_regions_models' num2str(iNumYears) '_years.pdf']);
  %%%% figure(30); aslprint([dir0 'rh_correlations_regions_models' num2str(iNumYears) '_years_MLS.pdf']);
  %%%% figure(31); aslprint([dir0 'fracWV_correlations_regions_models' num2str(iNumYears) '_years_MLS.pdf']);
  %%%% figure(32); aslprint([dir0 'tz_correlations_regions_models' num2str(iNumYears) '_years_MLS.pdf']);

  figure(33); aslprint([dir0 'rh_trends_zonal_p_5panels' num2str(iNumYears) '_years.pdf']);
  figure(34); aslprint([dir0 'fracWV_trends_zonal_p_5panels' num2str(iNumYears) '_years.pdf']);
  figure(35); aslprint([dir0 'tz_trends_zonal_p_5panels' num2str(iNumYears) '_years.pdf']);
  %%%% figure(33); aslprint([dir0 'rh_trends_zonal_p_5panels' num2str(iNumYears) '_years_MLS.pdf']); %% <<<<
  %%%% figure(34); aslprint([dir0 'fracWV_trends_zonal_p_5panels' num2str(iNumYears) '_years_MLS.pdf']); %% <<<<
  %%%% figure(35); aslprint([dir0 'tz_trends_zonal_p_5panels' num2str(iNumYears) '_years_MLS.pdf']); %% <<<<

  figure(20); aslprint([dir0 'sst_trends_lat_lon_6panels' num2str(iNumYears) '_years.pdf']);

  figure(21); aslprint([dir0 '200mb_rhtrends_lat_lon_5panels' num2str(iNumYears) '_years.pdf']);
  figure(22); aslprint([dir0 '200mb_fracWVtrends_lat_lon_5panels' num2str(iNumYears) '_years.pdf']);
  figure(23); aslprint([dir0 '200mb_tztrends_lat_lon_5panels' num2str(iNumYears) '_years.pdf']);

  figure(24); aslprint([dir0 '500mb_rhtrends_lat_lon_5panels' num2str(iNumYears) '_years.pdf']);
  figure(25); aslprint([dir0 '500mb_fracWVtrends_lat_lon_5panels' num2str(iNumYears) '_years.pdf']);
  figure(26); aslprint([dir0 '500mb_tztrends_lat_lon_5panels' num2str(iNumYears) '_years.pdf']);

  figure(27); aslprint([dir0 '800mb_rhtrends_lat_lon_5panels' num2str(iNumYears) '_years.pdf']);
  figure(28); aslprint([dir0 '800mb_fracWVtrends_lat_lon_5panels' num2str(iNumYears) '_years.pdf']);
  figure(29); aslprint([dir0 '800mb_tztrends_lat_lon_5panels' num2str(iNumYears) '_years.pdf']);

end
