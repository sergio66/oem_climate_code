function [era5,merra2,airsL3,climcapsL3,umbc,thecorr,amp,saverates_rlat_pres,summary_olr_mmw_stats] = driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends(strUMBC,iNumYears,iPentagonPlot,comment);

%{
  strUMBC = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithMLSL3_uncX100_50fatlayers_AIRSL3_ERA5_CMIP6_globalSSTfeedback.mat';
  strUMBC = '/asl/s1/sergio/JUNK/test9_guessstartWV_Vers1_march22_2023.mat';
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q05_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_feedback.mat';              %% not too bad at lower atm/polar!!!!
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q05_newERA5_2021jacs_startwith_MLSL3_TOA_guessWV_dRH_zero_bot_50fatlayers.mat';  %% too much oomph at gnd : use MLS L3 TOA and dRH/dt = 0 at bottom

  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q16_newERA5_2021jacs_startwith0_50fatlayers_ERA5calcs_spectraltrends.mat';       iNumYears = 20; %% compare geophysical trends from ERA5     spectral trends
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q16_newERA5_2021jacs_startwith0_50fatlayers_MERRA2calcs_spectraltrends.mat';     iNumYears = 20; %% compare geophysical trends from MERRA2   spectral trends
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q16_newERA5_2021jacs_startwith0_50fatlayers_AIRSL3calcs_spectraltrends.mat';     iNumYears = 20; %% compare geophysical trends from AIRS L3  spectral trends
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q16_newERA5_2021jacs_startwith0_50fatlayers_CLIMCAPSL3calcs_spectraltrends.mat'; iNumYears = 20; %% compare geophysical trends from CLIMCAPS spectral trends

  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset10_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 05;         %% use CarbonTracker CO2 trends
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset11_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 10;         %% use CarbonTracker CO2 trends
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset12_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 15;         %% use CarbonTracker CO2 trends

  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_MLS.mat'; iNumYears = 20;     %% use CarbonTracker CO2 trends, MLS a priori
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 20;         %% use CarbonTracker CO2 trends

  driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends(strUMBC,iNumYears);  
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iTypeIB = input('Plot the correlations of 200 mb,500 mb, 800 mb T/WV/RH as (+1/default) imagesc or as (-1) bargraphs : ');
%  if length(iTypeIB) == 0
%    iTypeIB = +1;
%  end
iTypeIB = +1;

%  disp('may as well hit +1 for the 5 panel (3 on top, 2 at bottom) plots and the correlation/mean/std dev 100-400,400-1000 mb'); 
%  iType5 = input('Show 5 panel plots? (-1/+1 default) : ');
%  if length(iType5) == 0
%    iType5 = +1;
%  end
iType5 = +1;

%  iPrint = input('Print plots for trends paper (-1 [default/+1) : ');
%  if length(iPrint) == 0
%    iPrint = -1;
%  end
iPrint = -1;

summary_olr_mmw_stats = [];

if nargin == 0

  strUMBC = [];
  iUMBC = -1;
  umbc = [];
  iNumYears = -1;
  iPentagonPlot = +1;
  comment = [];

  thecorr = [];

elseif nargin == 1
  iNumYears = 20;
  iPentagonPlot = +1;
  iUMBC = +1;
  comment = [];

elseif nargin == 2
  iPentagonPlot = +1;
  iUMBC = +1;
  comment = [];

elseif nargin == 3
  iUMBC = +1;
  comment = [];

elseif nargin == 4
  iUMBC = +1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultaxesLineWidth',1);
set(0,'DefaultaxesFontSize',16);

saverates_rlat_pres = [];
amp = [];

iWriteCorrelNumbers = +1; iBiasWRT_ERA5orUMBC = -1;   modelnames = {'THIS WORK','AIRS L3','CLIMCAPS L3','MERRA2','ERA5'};
iWriteCorrelNumbers = +1; iBiasWRT_ERA5orUMBC = +1;   modelnames = {'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK'};

fprintf(1,'iBiasWRT_ERA5orUMBC == %2i \n',iBiasWRT_ERA5orUMBC)

  %%   Name                   Size             Bytes  Class     Attributes
  %%   allXSlope              9x5x4              720  single
  %%   allX_frac_neg0pos      9x5x4             1440  double
  %%   allXchi                9x5x4              720  single
  %%   allXmean               9x5x4              720  single
  %%   allXstd                9x5x4              720  single

  %% allX_BLAH(iY,iX,:) will have size 9 x 5 x 4 
  %%     == [200/500/800 RH/WV/T] X ['ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK'] x [GLOBAL TRP MIDLAT POLAR]
  %%   first index  (iY)   1 .. 9 is  1:3=200 mb RH/WV/T   4:6=500 mb RH/WV/T   7:9=500 mb RH/WV/T
  %%   second index (iX)   1 .. 5 is  'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK' if iBiasWRT_ERA5orUMBC > 0, if iBiasWRT_ERA5orUMBC < 0 then swap UMBC, ERA5
  %%   third index  (iG)   1 .. 4 is  global,tropical,midlat,polar

  %% iY = 1,2,3  200 mb    4,5,6 500 mb    7,8,9 800 mb             = 9 total   [RH,WVfrac,T][RH,WVfrac,T][RH,WVfrac,T]
  %% iX = 1,2,3,4,5 === all, tropics,midlats,midlats+tropics,poles  = 5 total
  %% whos allXchi = 9 x 5 x 4                                       correlate ERA5 with [AIRS L3, CLIMCAPS, MERRA2, UMBC]
  %%   first index  (iY)   1 .. 9 is  1:3=200 mb RH/WV/T   4:6=500 mb RH/WV/T   7:9=500 mb RH/WV/T
  %%   second index (iX)   1 .. 5 is  'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK' if iBiasWRT_ERA5orUMBC > 0, if iBiasWRT_ERA5orUMBC < 0 then swap UMBC, ERA5
  %%   third index  (iG)   1 .. 4 is  global,tropical,midlat,polar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
addpath /asl/matlib/plotutils

load llsmap5

plays100 = load('ERA5_atm_data_2002_09_to_2022_08_trends_desc.mat','trend_plays');
plays100 = plays100.trend_plays;

%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
%% oh oh I had this backwards before Sept 28, 2023
do_XX_YY_from_X_Y

addpath /home/sergio/MATLABCODE/matlib/science/            %% for usgs_deg10_dem.m that has correct paths
[salti, landfrac] = usgs_deg10_dem(Y,X);

mu = cos(YY*pi/180);

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

clear ajunk
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
  loader = ['ajunk = load(' fin ');'];
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
  loader = ['ajunk = load(''' fin ''');'];
  eval(loader)
  fprintf(1,'loaded %2i years AIRS     geophysical trends from %s \n',iNumYearsPowWow,fin)
end

%%%%% <<<<< z11 is AIRS L3, now make it newz21 = z11 >>>>>
z11 = ajunk.thestats64x72.stemprate(:);
z11 = z11';
z11unc = ajunk.thestats64x72.stempratestd(:);
z11unc = z11unc';
newz21    = z11;
newz21unc = z11unc;

for ii = 1 : 72
  for jj = 1 : 64
    junk = squeeze(ajunk.thestats64x72.RHrate(ii,jj,:));
    airsL3.RHrate(ii,jj,:) = interp1(log10(ajunk.Qlevs),junk,log10(plays100),[],'extrap');
    junk = squeeze(ajunk.thestats64x72.waterrate(ii,jj,:));
    airsL3.waterrate(ii,jj,:) = interp1(log10(ajunk.Qlevs),junk,log10(plays100),[],'extrap');
    junk = squeeze(ajunk.thestats64x72.ptemprate(ii,jj,:));
    airsL3.ptemprate(ii,jj,:) = interp1(log10(ajunk.Tlevs),junk,log10(plays100),[],'extrap');
  end
end
airsL3.stemprate = ajunk.thestats64x72.stemprate(:)';
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

clear ajunk
fin = ['/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug' num2str(2002+iNumYears) '_' num2str(iNumYears) 'yr_desc.mat'];
iNumYearsPowWow = iNumYears;
if ~exist(fin)
  iaJunk = [5 10 15 20];
  moo = abs(iaJunk-iNumYears);
  moo = find(moo == min(moo),1);    
  fin = ['/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug' num2str(2002+iaJunk(moo)) '_' num2str(iaJunk(moo)) 'yr_desc.mat'];    
  iNumYearsPowWow = 2002+iaJunk(moo);
end
loader = ['ajunk = load(''' fin ''');'];
eval(loader)
fprintf(1,'loaded %2i years CLIMCAPS geophysical trends from %s \n',iNumYearsPowWow,fin)

%%%%% <<<<< z12 is AIRS L3, now make it newz22 = z12 >>>>>
z12 = ajunk.thestats64x72.stemprate(:);
z12 = z12';
z12unc = ajunk.thestats64x72.stempratestd(:);
z12unc = z12unc';
newz22    = z12;
newz22unc = z12unc;

for ii = 1 : 72
  for jj = 1 : 64
    junk = squeeze(ajunk.thestats64x72.RHrate(ii,jj,:));
    climcapsL3.RHrate(ii,jj,:) = interp1(log10(ajunk.Qlevs/100),100*junk,log10(plays100),[],'extrap');
    junk = squeeze(ajunk.thestats64x72.waterrate(ii,jj,:));
    climcapsL3.waterrate(ii,jj,:) = interp1(log10(ajunk.Qlevs/100),junk,log10(plays100),[],'extrap');
    junk = squeeze(ajunk.thestats64x72.ptemprate(ii,jj,:));
    climcapsL3.ptemprate(ii,jj,:) = interp1(log10(ajunk.Tlevs/100),junk,log10(plays100),[],'extrap');
  end
end
climcapsL3.stemprate = ajunk.thestats64x72.stemprate(:)';
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

clear ajunk
fin = ['MERRA2_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc.mat'];
iNumYearsPowWow = iNumYears;
if ~exist(fin)
  iaJunk = [5 10 15 20];
  moo = abs(iaJunk-iNumYears);
  moo = find(moo == min(moo),1);    
  fin = ['MERRA2_atm_data_2002_09_to_' num2str(2002+iaJunk(moo)) '_08_trends_desc.mat']; 
  iNumYearsPowWow = 2002+iaJunk(moo);
end
loader = ['ajunk = load(''' fin ''');'];
eval(loader)
fprintf(1,'loaded %2i years MERRA2   geophysical trends from %s \n',iNumYearsPowWow,fin)

%%%%% <<<<< z21 is MERRA2, now make it newz12 = z21 >>>>>
z21    = ajunk.trend_stemp;
z21unc = ajunk.trend_stemp_err;
merra2.stemprate = ajunk.trend_stemp;
merra2.RHrate   = ajunk.trend_RH;
merra2.waterrate = ajunk.trend_gas_1;
merra2.ptemprate = ajunk.trend_ptemp;
merra2.trend_mmw = ajunk.trend_mmw;
newz12    = z21;
newz12unc = z21unc;

bad = find(isnan(merra2.stemprate) | isinf(merra2.stemprate)); merra2.stemprate(bad) = 0;
bad = find(isnan(merra2.ptemprate) | isinf(merra2.ptemprate)); merra2.ptemprate(bad) = 0;
bad = find(isnan(merra2.waterrate) | isinf(merra2.waterrate)); merra2.waterrate(bad) = 0;

iFig = 14;
iFig = 1;
figure(iFig); clf
[h,ha,p,pa] = rtpread('summary_19years_all_lat_all_lon_2002_2021_monthlyMERRA2.op.rtp');
figure(13); clf; dmmw_dST_MERRA2 = dmmw_dsst_VS_lat(lfmaskA,lfmaskL,lfmaskO,mmw0,merra2.trend_mmw,merra2.stemprate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear ajunk
if iNumYears > 0
  eraX = ['ERA5_atm_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc.mat'];
  eraX = ['ERA5_atm_N_cld_data_2002_09_to_' num2str(2002+iNumYears) '_08_trends_desc.mat'];
else
  eraX = ['ERA5_atm_data_2002_09_to_' num2str(2022) '_08_trends_desc.mat'];
end
loader = ['ajunk = load(''' eraX ''');'];
fprintf(1,'loaded %2i years ERA5     geophysical trends from %s \n',iNumYears,eraX)
eval(loader);

%%%%% <<<<< z22 is ERA5, now make it newz11 = z11 >>>>>
z22    = ajunk.trend_stemp;
z22unc = ajunk.trend_stemp_err;
era5.stemprate = ajunk.trend_stemp;
era5.RHrate    = ajunk.trend_RH;
era5.waterrate = ajunk.trend_gas_1;
era5.ptemprate = ajunk.trend_ptemp;
era5.trend_mmw = ajunk.trend_mmw;
trend_rlat64   = ajunk.trend_rlat64;
newz11    = z22;
newz11unc = z22unc;

    trend_rlat64_polar    = find(abs(trend_rlat64) >= 60);
    trend_rlat64_midlat   = find(abs(trend_rlat64) <  60 & abs(trend_rlat64) >= 30);
    trend_rlat64_tropical = find(abs(trend_rlat64) <  30);
    i10  = find(plays100 >= 10 & plays100 <= 1000);
    i100 = find(plays100 >= 100 & plays100 <= 1000);
    i400 = find(plays100 >= 400 & plays100 <= 1000);
    i10_100   = find(plays100 >= 10  & plays100 <= 100);
    i100_400  = find(plays100 >= 100 & plays100 <= 400);
    i400_1000 = find(plays100 >= 400 & plays100 <= 1000);

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
fprintf(1,'(ALL  )         |  %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac == 0);  
junk = [nanmean(z11(boo))  nanmean(z12(boo)) nanmean(z21(boo))  nanmean(z22(boo))];
fprintf(1,'(OCEAN)         |  %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac > 0.99);  
junk = [nanmean(z11(boo))  nanmean(z12(boo)) nanmean(z21(boo))  nanmean(z22(boo))];
fprintf(1,'(LAND )         |  %8.3f %8.3f %8.3f %8.3f \n',junk);
disp('----------------|---------------------------------------');

disp('SKT trend K/yr  | AIRS   CLIMCAPS   MERRA2      ERA5    ');
disp('area weighted with cos(lat)')
disp('----------------|---------------------------------------');
junk = [nansum(z11.*mu)/nansum(mu)  nansum(z12.*mu)/nansum(mu) nansum(z21.*mu)/nansum(mu)  nansum(z22.*mu)/nansum(mu)];
fprintf(1,'(ALL  )         |  %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac == 0);  
junk = [nansum(z11(boo).*mu(boo))/nansum(mu(boo))  nansum(z12(boo).*mu(boo))/nansum(mu(boo)) nansum(z21(boo).*mu(boo))/nansum(mu(boo))  nansum(z22(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'(OCEAN)         |  %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac > 0.99);  
junk = [nansum(z11(boo).*mu(boo))/nansum(mu(boo))  nansum(z12(boo).*mu(boo))/nansum(mu(boo)) nansum(z21(boo).*mu(boo))/nansum(mu(boo))  nansum(z22(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'(LAND )         |  %8.3f %8.3f %8.3f %8.3f \n',junk);
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
  umbcX        = load(strUMBC,'deltaRH','deltaT','fracWV','results','resultsunc');

  strGISS = ['ChrisHTrends/giss_trends_2002_' num2str(2002+iNumYears) '.mat'];
  iNumYearsPowWow = iNumYears;
  if ~exist(strGISS)
    iaJunk = [05 10 15 818 19 20];
    moo = abs(iaJunk-iNumYears);
    moo = find(moo == min(moo),1);    
    strGISS = ['ChrisHTrends/giss_trends_2002_' num2str(2002+iaJunk(moo)) '.mat'];
    iNumYearsPowWow = 2002+iaJunk(moo);
  end

  clear ajunk agiss
  fprintf(1,'loaded %2i years %s \n',iNumYearsPowWow,strGISS);
  loader = ['agiss = load(''' strGISS ''');'];
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
  umbcX.fracWV       = umbcX.fracWV(1:100,:);
  umbcX.stemprate    = umbcX.results(:,6)';
  if isfield(umbcX,'resultsunc')
    umbcX.stemprateunc = umbcX.resultsunc(:,6)';
  else
    umbcX.stemprateunc = 0.1*umbcX.results(:,6)';
  end

  umbc.stemprate    = umbcX.stemprate;
  umbc.stemprateunc = umbcX.stemprateunc;
  umbc.ptemprate    = umbcX.deltaT;
  umbc.RHrate       = umbcX.deltaRH;
  umbc.waterrate    = umbcX.fracWV;

  clear umbcX

  z11x       = umbc.stemprate;
  z11xunc    = umbc.stemprateunc;
  newz11x    = umbc.stemprate;
  newz11xunc = umbc.stemprateunc;

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
  plotaxis2; hl = legend('MERRA2','ERA5','AIRS L3','CLIMCAPSL3','THIS WORK','location','best','fontsize',8); ylabel('d%mmw/dST'); xlabel('Latitude'); ylim([-10 +20]); title('Ocean')

  figure(14); clf
  plot(trend_rlat64,smooth(dmmw_dST_MERRA2.land,10),'g',trend_rlat64,smooth(dmmw_dST_ERA5.land,10),'r',trend_rlat64,smooth(dmmw_dST_airsL3.land,10),'b',trend_rlat64,smooth(dmmw_dST_climcapsL3.land,10),'c',...
       trend_rlat64,smooth(dmmw_dST_umbc.land,10),'k','linewidth',2); 
  plotaxis2; hl = legend('MERRA2','ERA5','AIRS L3','CLIMCAPSL3','THIS WORK','location','best','fontsize',8); ylabel('d%mmw/dST'); xlabel('Latitude'); ylim([-10 +20]); title('Land')

  figure(15); clf
  plot(trend_rlat64,smooth(dmmw_dST_MERRA2.all,10),'g',trend_rlat64,smooth(dmmw_dST_ERA5.all,10),'r',trend_rlat64,smooth(dmmw_dST_airsL3.all,10),'b',trend_rlat64,smooth(dmmw_dST_climcapsL3.all,10),'c',...
       trend_rlat64,smooth(dmmw_dST_umbc.all,10),'k','linewidth',2); 
  plotaxis2; hl = legend('MERRA2','ERA5','AIRS L3','CLIMCAPSL3','THIS WORK','location','best','fontsize',8); ylabel('d%mmw/dST'); xlabel('Latitude'); ylim([-10 +20]); title('All')

  %%%%%%%%%%%%%%%%%%%%%%%%%

  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  maskLF = ones(size(z22));
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'THIS WORK';
  plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';
  plotoptions.str22 = 'ERA5';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  aslmap_2x2tiledlayout(z11x,z12,z21,z22,iFig,plotoptions);
  
  clear plotoptions
  plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'THIS WORK';
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
  figure(iFig); clf
  miaowA = squeeze(nanmean(reshape(umbc.ptemprate,100,72,64),2));
  miaowB = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
  miaow21 = miaowA - miaowB;
  pcolor(trend_rlat64,plays100,miaow12); colormap(llsmap5); ylim([10 1000])
  title('dTz/dt : UMBC - ERA5'); shading interp; colorbar
  caxis([-1 +1]*0.05)
  set(gca,'ydir','reverse')

  %%%%%
  iFig = iFig + 1;
  figure(iFig); clf
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

  [r,chisqr,P] = nanlinearcorrelation(umbc.stemprate,era5.stemprate);           thecorr.ST(1) = r;
  [r,chisqr,P] = nanlinearcorrelation(umbc.stemprate,merra2.stemprate);         thecorr.ST(2) = r;
  [r,chisqr,P] = nanlinearcorrelation(umbc.stemprate,airsL3.stemprate);         thecorr.ST(3) = r;
  [r,chisqr,P] = nanlinearcorrelation(umbc.stemprate,climcapsL3.stemprate);     thecorr.ST(4) = r;
  [r,chisqr,P] = nanlinearcorrelation(umbc.stemprate',agiss.giss_trend4608(:)); thecorr.ST(5) = r;
  [r,chisqr,P] = nanlinearcorrelation(era5.stemprate',agiss.giss_trend4608(:)); thecorr.ST_ERA5_GISS = r;

  %% do fractional signs
  zall = [era5.stemprate; merra2.stemprate; airsL3.stemprate; climcapsL3.stemprate; umbc.stemprate; agiss.giss_trend4608(:)'];
  wall = ones(size(era5.stemprate));
  wall = cos(YY*pi/180)';
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos,frac_neg0pos_mean_std_stemp] = corrplot_weighted_mean_stddev(zall',wall',modelnames);

  iFig = iFig + 1;
  figure(iFig); clf; plot(umbc.stemprate,era5.stemprate,'r.',umbc.stemprate,merra2.stemprate,'g.',...
                   umbc.stemprate,airsL3.stemprate,'b.',umbc.stemprate,climcapsL3.stemprate,'c.');
  axis([-1 +1 -1 +1]*0.3); line([-0.3 +0.3],[-0.3 +0.3],'color','k','linewidth',3)
  plotaxis2; hl = legend('ERA5','MERRA2','AIRSL3','CLIMCAPSL3','location','best'); title('dSKT/dt [K/yr]'); 
  xlabel('THIS WORK');
 
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
   fprintf(1,'ST correlations vs UMBC : ERA5/MERRA2/AIRSL3/CLIMCAPSL3/GISS =  %8.2f %8.2f %8.2f %8.2f %8.2f\n',thecorr.ST);
   fprintf(1,'ST correlations         : ERA5 vs GISS =  %8.2f\n',thecorr.ST_ERA5_GISS);
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

  % z11 = AIRS,       z12 = climcaps      z21 = MERRA2        z22 = ERA5            z11x/z31 = umbc      z32 = giss
  % newz11 = ERA5, newz12 = MERRA2     newz21 = AIRSL3     newz22 = CLIMCAPS  newz11x/newz31 = umbc   newz32 = giss

  % z31 = z11x;    z31unc = z11xunc; 
  % z32 = z32;
  % z32    = agiss.giss_trend4608(:);     z32 = z32';
  % z32unc = agiss.giss_trend_err4608(:); z32unc = z32unc';

  newz31 = newz11x;                        newz31unc = newz11xunc; 
  newz32    = agiss.giss_trend4608(:);     newz32 = newz32';
  newz32unc = agiss.giss_trend_err4608(:); newz32unc = newz32unc';

  %%%%%%%%%%

  iFig = iFig + 1;
  figure(iFig); clf;
  clear plotoptions
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  % plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  % plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.str11 = 'ERA5';        plotoptions.str22 = 'MERRA2';
  plotoptions.str11 = 'AIRS L3';     plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str31 = 'THIS WORK';   plotoptions.str32 = 'GISS';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  aslmap_3x2tiledlayout(newz11,newz12,newz21,newz22,newz31,newz32,iFig,plotoptions);

  figure(20); clf; aslmap(20,rlat65,rlon73,smoothn(reshape(newz31,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.151); %% UMBC
  figure(21); clf; aslmap(21,rlat65,rlon73,smoothn(reshape(newz32,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.151); %% GISS
  
  %%%%%%%%%% 

  clear plotoptions
  plotoptions.cx = [-1 +1]*0.151;  plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  % plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3'; plotoptions.str13 = 'MERRA2';    
  % plotoptions.str21 = 'ERA5';      plotoptions.str22 = 'THIS WORK';   plotoptions.str23 = 'GISS';
  % z31 = z11x;    z31unc = z11xunc;  z32 = z32;
  % z31 = z11x;    z31unc = z11xunc;  z32 = z32;
  plotoptions.str11 = 'ERA5';         plotoptions.str12 = 'MERRA2';      plotoptions.str13 = 'AIRS L3';    
  plotoptions.str21 = 'CLIMCAPS L3';  plotoptions.str22 = 'THIS WORK';   plotoptions.str23 = 'GISS';

  aslmap_2x3tiledlayout(newz11,newz12,newz21,newz22,newz31,newz32,iFig,plotoptions);

  %%%%%%%%%% 

  clear plotoptions
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  plotoptions.cx = [-1 +1]*0.151;  plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  % plotoptions.str11 = 'AIRS L3';       plotoptions.str12 = 'ERA5';          plotoptions.str13 = 'THIS WORK';    
  % plotoptions.str21 = 'CLIMCAPS L3';   plotoptions.str22 = 'MERRA2';        plotoptions.str23 = 'GISS';
  % z31 = z11x;    z31unc = z11xunc; 
  % z32 = z32;
  plotoptions.str11 = 'ERA5';    plotoptions.str12 = 'THIS WORK';   plotoptions.str13 = 'AIRS L3';    
  plotoptions.str21 = 'MERRA2';  plotoptions.str22 = 'GISS';        plotoptions.str23 = 'CLIMCAPSL3';

  aslmap_2x3tiledlayout(newz11,newz31,newz21,newz12,newz32,newz22,iFig,plotoptions);

  %%%%%%%%%% 

  tropics = find(abs(YY) <= 30);          
  midlatsNtropics = find(abs(YY) <= 60); 
  midlats = setdiff(midlatsNtropics,tropics); 
  poles = find(abs(YY) > 60);             

  tropics = nan(size(YY));         junk = find(abs(YY) <= 30);                tropics(junk) = 1;
  midlatsNtropics = nan(size(YY)); junk = find(abs(YY) <= 60);                midlatsNtropics(junk) = 1;
  midlats = nan(size(YY));         junk = find(abs(YY) > 30 & abs(YY) <= 60); midlats(junk) = 1;
  poles = nan(size(YY));           junk = find(abs(YY) > 60);                 poles(junk) = 1;

  iNewOrOld = +1;
  if iNewOrOld < 0
    show_skt_trends_6models_old
  else
    show_skt_trends_6models_new
  end

  % disp('hit RET to continue to RH/Tz/Wvz plots at 200/500/800 mb '); pause
  pause(1)

  if ~exist('era5skt')
    era5skt = load('ERA5_atm_data_2002_09_to_2022_08_desc.mat','all');
    era5mmw = era5skt.all.mmw';
    era5skt = era5skt.all.stemp';
  end

  [Lera5skt, EOFsera5skt, ECera5skt, errorera5skt] = detrend_time_series_make_EOF_v1(era5skt');

  anom_1231 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/anomaly_chID_1520_Q03.mat');
  [L1231, EOFs1231, EC1231, error1231] = detrend_time_series_make_EOF_v1(anom_1231.btanom');
  for iEOF = 1 : 10
    [r,chisqr,P] = nanlinearcorrelation(EOFs1231(:,iEOF),era5.stemprate');           thecorrEOF.BT1231_ST(iEOF,1) = r;
    [r,chisqr,P] = nanlinearcorrelation(EOFs1231(:,iEOF),merra2.stemprate');         thecorrEOF.BT1231_ST(iEOF,2) = r;
    [r,chisqr,P] = nanlinearcorrelation(EOFs1231(:,iEOF),airsL3.stemprate');         thecorrEOF.BT1231_ST(iEOF,3) = r;
    [r,chisqr,P] = nanlinearcorrelation(EOFs1231(:,iEOF),climcapsL3.stemprate');     thecorrEOF.BT1231_ST(iEOF,4) = r;
    bonk = agiss.giss_trend4608; bonk = bonk(:);
    [r,chisqr,P] = nanlinearcorrelation(EOFs1231(iEOF,:),bonk');                     thecorrEOF.BT1231_ST(iEOF,5) = r;
    [r,chisqr,P] = nanlinearcorrelation(EOFs1231(:,iEOF),umbc.stemprate');           thecorrEOF.BT1231_ST(iEOF,6) = r;
  end

  iFig = iFig + 1;
  figure(iFig); clf
   aslmap(iFig,rlat65,rlon73,reshape(EOFs1231(:,1),72,64)', [-90 +90],[-180 +180]); colormap(usa2);

  clf
  ymax = 4;
  imagesc(thecorrEOF.BT1231_ST); colorbar; caxis([-1 +1]); ylabel('BT1231 EOF'); ylim([0.5 0.5+ymax]); 
  title('dSKT/dt Correlations with BT1231 EOF')
  set(gca,'ytick',[1:ymax],'yticklabel',num2str((1:ymax)'),'fontsize',10);
  set(gca,'xtick',[1:6],'xticklabel',{'ERA5','MERRA2','AIRS L3','CLIMCAPSL3','GISS','THIS WORK'},'fontsize',8);
  colormap(usa2)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  oceantropics = intersect(find(isfinite(ocean)),find(isfinite(tropics)));

  junk = [nansum(land) nansum(ocean) nansum(tropics) nansum(midlats) nansum(midlatsNtropics) nansum(poles)];
  fprintf(1,'breakdown : land = %4i ocean = %4i      tropics = %4i midlats = %4i midlatsNtropics = %4i poles = %4i of 4608 \n',junk);
  fprintf(1,'breakdown : land = %4.3f ocean = %4.3f    tropics = %4.3f midlats = %4.3f midlatsNtropics = %4.3f poles = %4.3f of 4608 \n',junk/4608);

%%%%%%%%%%%%%%%%%%%% new new
  tropics = find(abs(YY) <= 30);          
  midlatsNtropics = find(abs(YY) <= 60); 
  midlats = setdiff(midlatsNtropics,tropics); 
  poles = find(abs(YY) > 60);             

  tropics = nan(size(YY));         junk = find(abs(YY) <= 30);                tropics(junk) = 1;
  midlatsNtropics = nan(size(YY)); junk = find(abs(YY) <= 60);                midlatsNtropics(junk) = 1;
  midlats = nan(size(YY));         junk = find(abs(YY) > 30 & abs(YY) <= 60); midlats(junk) = 1;
  poles = nan(size(YY));           junk = find(abs(YY) > 60);                 poles(junk) = 1;
%%%%%%%%%%%%%%%%%%%% new new

  %%% EACH AND EVERYONE OF COMPUTATIONS IN THIS SECTION IS CORRELATING/MEAN/STDDEV/COMPARING to ERA5
  %%%   zall0 = [zyx22(:) zyx11(:) zyx12(:) zyx21(:) zyx31(:)]'; === [ERA5 AIRS CLIMCAPS MERRA UMBC]
  %%%   corrplot(zall) use the first column as the standard ===== ERA5

  i200 = find(plays100 >= 200,1);
  i500 = find(plays100 >= 500,1);
  i800 = find(plays100 >= 800,1);

  %% allX_BLAH(iY,iX,:) will have size 9 x 5 x 4 
  %%     == [200/500/800 RH/WV/T] X ['ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK'] x [GLOBAL TRP MIDLAT POLAR]
  %%   first index  (iY)   1 .. 9 is  1:3=200 mb RH/WV/T   4:6=500 mb RH/WV/T   7:9=500 mb RH/WV/T
  %%   second index (iX)   1 .. 5 is  'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK' if iBiasWRT_ERA5orUMBC > 0, if iBiasWRT_ERA5orUMBC < 0 then swap UMBC, ERA5
  %%   third index  (iG)   1 .. 4 is  global,tropical,midlat,polar

  %% iY = 1,2,3  200 mb    4,5,6 500 mb    7,8,9 800 mb             = 9 total   [RH,WVfrac,T][RH,WVfrac,T][RH,WVfrac,T]
  %% iX = 1,2,3,4,5 === all, tropics,midlats,midlats+tropics,poles  = 5 total
  %% whos allXchi = 9 x 5 x 4                                       correlate ERA5 with [AIRS L3, CLIMCAPS, MERRA2, UMBC]
  %%   first index  (iY)   1 .. 9 is  1:3=200 mb RH/WV/T   4:6=500 mb RH/WV/T   7:9=500 mb RH/WV/T
  %%   second index (iX)   1 .. 5 is  'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK' if iBiasWRT_ERA5orUMBC > 0, if iBiasWRT_ERA5orUMBC < 0 then swap UMBC, ERA5
  %%   third index  (iG)   1 .. 4 is  global,tropical,midlat,polar

%%%%%%%%%%%%%%%%%%%% new new
  iX = 0;
  iY = 0;

  iPrint = -1;

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
  plotoptions.str31 = 'THIS WORK'; plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'THIS WORK'; zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dRH/dt 200 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
  zall0 = [zyx22(:) zyx11(:) zyx12(:) zyx21(:) zyx31(:)]';  %% make ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
  if iBiasWRT_ERA5orUMBC < 0
    %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
    zall0 = zall0([5 2 3 4 1],:);
  end
  wall0 = cos(YY*pi/180)';
  zall = zall0; wall = wall0;
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dRH/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(tropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dRH/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlats)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dRH/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlatsNtropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dRH/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(poles)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dRH/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n\n',thecorr.junk); end
  
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
  plotoptions.str31 = 'THIS WORK'; plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'THIS WORK'; zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dWVfrac/dt 200 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
  zall0 = [zyx22(:) zyx11(:) zyx12(:) zyx21(:) zyx31(:)]';  %% make ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
  if iBiasWRT_ERA5orUMBC < 0
    %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
    zall0 = zall0([5 2 3 4 1],:);
  end
  wall0 = cos(YY*pi/180)';
  zall = zall0; wall = wall0;
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dWV/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(tropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dWV/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlats)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dWV/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlatsNtropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dWV/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(poles)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dWV/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n\n',thecorr.junk); end
   
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
  plotoptions.str31 = 'THIS WORK'; plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'THIS WORK'; zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dT/dt 200 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
  zall0 = [zyx22(:) zyx11(:) zyx12(:) zyx21(:) zyx31(:)]';  %% make ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
  if iBiasWRT_ERA5orUMBC < 0
    %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
    zall0 = zall0([5 2 3 4 1],:);
  end
  wall0 = cos(YY*pi/180)';
  zall = zall0; wall = wall0;
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dTz/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(tropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dTz/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlats)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dTz/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlatsNtropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dTz/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(poles)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'200 mb dTz/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n\n',thecorr.junk); end
   
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
  plotoptions.str31 = 'THIS WORK'; plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'THIS WORK'; zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dRH/dt 500 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
  zall0 = [zyx22(:) zyx11(:) zyx12(:) zyx21(:) zyx31(:)]';  %% make ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
  if iBiasWRT_ERA5orUMBC < 0
    %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
    zall0 = zall0([5 2 3 4 1],:);
  end
  wall0 = cos(YY*pi/180)';
  zall = zall0; wall = wall0;
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dRH/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(tropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dRH/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlats)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dRH/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlatsNtropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dRH/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(poles)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dRH/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n\n',thecorr.junk); end
  
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
  plotoptions.str31 = 'THIS WORK'; plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'THIS WORK'; zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dWVfrac/dt 500 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
  zall0 = [zyx22(:) zyx11(:) zyx12(:) zyx21(:) zyx31(:)]';  %% make ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
  if iBiasWRT_ERA5orUMBC < 0
    %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
    zall0 = zall0([5 2 3 4 1],:);
  end
  wall0 = cos(YY*pi/180)';
  zall = zall0; wall = wall0;
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dWV/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(tropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dWV/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlats)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dWV/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlatsNtropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dWV/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(poles)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dWV/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n\n',thecorr.junk); end

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
  plotoptions.str31 = 'THIS WORK'; plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'THIS WORK'; zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dT/dt 500 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
  zall0 = [zyx22(:) zyx11(:) zyx12(:) zyx21(:) zyx31(:)]';  %% make ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
  if iBiasWRT_ERA5orUMBC < 0
    %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
    zall0 = zall0([5 2 3 4 1],:);
  end
  wall0 = cos(YY*pi/180)';
  zall = zall0; wall = wall0;
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dTz/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(tropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dTz/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlats)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dTz/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlatsNtropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dTz/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(poles)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'500 mb dTz/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n\n',thecorr.junk); end

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
  plotoptions.str31 = 'THIS WORK'; plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'THIS WORK'; zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dRH/dt 800 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
  zall0 = [zyx22(:) zyx11(:) zyx12(:) zyx21(:) zyx31(:)]';  %% make ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
  if iBiasWRT_ERA5orUMBC < 0
    %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
    zall0 = zall0([5 2 3 4 1],:);
  end
  wall0 = cos(YY*pi/180)';
  zall = zall0; wall = wall0;
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dRH/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(tropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dRH/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlats)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dRH/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlatsNtropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dRH/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(poles)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dRH/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n\n',thecorr.junk); end

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
  plotoptions.str31 = 'THIS WORK'; plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'THIS WORK'; zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dWVfrac/dt 500 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
  zall0 = [zyx22(:) zyx11(:) zyx12(:) zyx21(:) zyx31(:)]';  %% make ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
  if iBiasWRT_ERA5orUMBC < 0
    %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
    zall0 = zall0([5 2 3 4 1],:);
  end
  wall0 = cos(YY*pi/180)';
  zall = zall0; wall = wall0;
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dWV/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(tropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dWV/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlats)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dWV/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlatsNtropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dWV/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(poles)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dWV/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n\n',thecorr.junk); end

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
  plotoptions.str31 = 'THIS WORK'; plotoptions.str32 = '    ';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  if iPentagonPlot < 0
    aslmap_3x2tiledlayout(zyx11,zyx12,zyx21,zyx22,zyx31,0*zyx31,iFig,plotoptions);
  else
    plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
    plotoptions.strzz = 'THIS WORK'; zyxzz = zyx31;
    plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
    plotoptions.cstr  = 'dT/dt 800 mb';
    aslmap_2x1x2tiledlayout(zyx11,zyx12,zyxzz,zyx21,zyx22,iFig,plotoptions);
  end
  zall0 = [zyx22(:) zyx11(:) zyx12(:) zyx21(:) zyx31(:)]';  %% make ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
  if iBiasWRT_ERA5orUMBC < 0
    %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
    zall0 = zall0([5 2 3 4 1],:);
  end
  wall0 = cos(YY*pi/180)';
  zall = zall0; wall = wall0;
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dTz/dt                       correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(tropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dTz/dt       TROPICS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlats)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dTz/dt       MIDLATS         correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(midlatsNtropics)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dTz/dt       MIDLATSNTROPICS correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n',thecorr.junk); end
  junk = find(isfinite(poles)); zall = zall0(:,junk); wall = wall0(junk);
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos] = corrplot_weighted_mean_stddev(zall',wall',modelnames); thecorr.junk = R(2:5,1);
    iX = iX + 1; allXchi(iY,iX,:) = R(2:5,1); allXmean(iY,iX,:) = m(2:5); allXstd(iY,iX,:) = s(2:5); allXSlope(iY,iX,:) = linearfit(2:5,1); allX_frac_neg0pos(iY,iX,:) = frac_neg0pos(2:5,1);
    if iPrint > 0; fprintf(1,'800 mb dTz/dt       POLES           correlations vs ERA5 : AIRSL3/CLIMCAPSL3/MERRA2/UMBC =  %8.2f %8.2f %8.2f %8.2f\n\n',thecorr.junk); end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% allX_BLAH(iY,iX,:) will have size 9 x 5 x 4
  %%   first index  (iY)   1 .. 9 is  1:3=200 mb RH/WV/T   4:6=500 mb RH/WV/T   7:9=500 mb RH/WV/T
  %%   second index (iX)   1 .. 5 is  'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK' if iBiasWRT_ERA5orUMBC > 0, if iBiasWRT_ERA5orUMBC < 0 then swap UMBC, ERA5  
  %%   third index  (iG)   1 .. 4 is  global,tropical,midlat,polar

  %% iX = 1,2,3,4,5 === all, tropics,midlats,midlats+tropics,poles  = 5 total
  %% iY = 1,2,3  200 mb    4,5,6 500 mb    7,8,9 800 mb             = 9 total   [RH,WVfrac,T][RH,WVfrac,T][RH,WVfrac,T]
  %% whos allXchi = 9 x 5 x 4                                       correlate ERA5 with [AIRS L3, CLIMCAPS, MERRA2, UMBC]
  %%   first index  (iY)   1 .. 9 is  1:3=200 mb RH/WV/T   4:6=500 mb RH/WV/T   7:9=500 mb RH/WV/T
  %%   second index (iX)   1 .. 5 is  'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK' if iBiasWRT_ERA5orUMBC > 0, if iBiasWRT_ERA5orUMBC < 0 then swap UMBC, ERA5
  %%   third index  (iG)   1 .. 4 is  global,tropical,midlat,polar

  %iFig = iFig + 1;
  %figure(iFig); clf;   
  %h1=subplot(311);  bar(1:5,squeeze(allXchi(1,:,:))); ylim([-1 1]); hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','best','fontsize',6); ylabel('RH 200')
  %h2=subplot(312);  bar(1:5,squeeze(allXchi(4,:,:))); ylim([-1 1]); hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','best','fontsize',6); ylabel('RH 500')
  %h3=subplot(313);  bar(1:5,squeeze(allXchi(7,:,:))); ylim([-1 1]); hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','best','fontsize',6); ylabel('RH 800')
  %adjust(31,h1,h2,h3,'even');
  %set(gca,'xtick',[1:5],'xticklabel',names)

  zonalnames = {'ALL','T','ML','ML+T','P'};

  %iTypeIB = input('Plot the correlations of 200 mb,500 mb, 800 mb T/WV/RH as (+1/default) imagesc or as (-1) bargraphs : ');
  %if length(iTypeIB) == 0
  %  iTypeIB = +1;
  %end

  if iTypeIB < 0
    iFig = iFig + 1;
    figure(iFig); clf;   
    ta = tiledlayout(3,1,'TileSpacing','None', 'Padding','None');
    ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
    ymin = 0;
    tfov(1) = nexttile;  bar(1:5,squeeze(allXchi(1,:,:))); ylim([ymin 1]); ylabel('RH 200');  %hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','best','fontsize',6); ylabel('RH 200')
    tfov(2) = nexttile;  bar(1:5,squeeze(allXchi(4,:,:))); ylim([ymin 1]); ylabel('RH 500');  %hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','best','fontsize',6); ylabel('RH 500')
    tfov(3) = nexttile;  bar(1:5,squeeze(allXchi(7,:,:))); ylim([ymin 1]); ylabel('RH 800');  hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','best','fontsize',6); ylabel('RH 800')
    tfov(1).XTickLabel = '';  tfov(1).XLabel.String = [];      tfov(2).XTickLabel = '';  tfov(2).XLabel.String = [];
    set(gca,'xtick',[1:5],'xticklabel',zonalnames)
    ta.Padding = 'none';
    ta.TileSpacing = 'tight';

    iFig = iFig + 1;
    figure(iFig); clf;   
    ta = tiledlayout(3,1,'TileSpacing','None', 'Padding','None');
    ta.OuterPosition = [0.0375 0.0375 0.925 0.925];    
    ymin = -0.25;
    ymin2 = 0;
    tfov(1) = nexttile;  bar(1:5,squeeze(allXchi(2,:,:))); ylim([ymin  1]); ylabel('WV 200'); %hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','s','fontsize',5); 
    tfov(2) = nexttile;  bar(1:5,squeeze(allXchi(5,:,:))); ylim([ymin2 1]); ylabel('WV 500'); %hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','best','fontsize',6);
    tfov(3) = nexttile;  bar(1:5,squeeze(allXchi(8,:,:))); ylim([ymin2 1]); ylabel('WV 800'); %hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','ne','fontsize',5);
    tfov(1).XTickLabel = '';  tfov(1).XLabel.String = [];      tfov(2).XTickLabel = '';  tfov(2).XLabel.String = [];
    set(gca,'xtick',[1:5],'xticklabel',zonalnames)
    ta.Padding = 'none';
    ta.TileSpacing = 'tight';

    iFig = iFig + 1;
    figure(iFig); clf;   
    ta = tiledlayout(3,1,'TileSpacing','None', 'Padding','None');
    ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
    ymin = 0;
    tfov(1) = nexttile;  bar(1:5,squeeze(allXchi(3,:,:))); ylim([ymin 1]); ylabel('T 200');  %hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','best','fontsize',6); ylabel('T 200')
    tfov(2) = nexttile;  bar(1:5,squeeze(allXchi(6,:,:))); ylim([ymin 1]); ylabel('T 500');  %hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','best','fontsize',6); ylabel('T 500')
    tfov(3) = nexttile;  bar(1:5,squeeze(allXchi(9,:,:))); ylim([ymin 1]); ylabel('T 800');  %hl = legend('AIRS v7','CLIMCAPS','MERRA2','THIS WORK','location','best','fontsize',6); ylabel('T 800')
    tfov(1).XTickLabel = '';  tfov(1).XLabel.String = [];      tfov(2).XTickLabel = '';  tfov(2).XLabel.String = [];
    set(gca,'xtick',[1:5],'xticklabel',zonalnames)
    ta.Padding = 'none';
    ta.TileSpacing = 'tight';
  
  else  
    ymin = 0;
    %mnames = {'AIRS v7','CLIMCAPS','MERRA2','THIS WORK'};    
    zonalnames = {'GLOBAL','TROPICAL','MIDLAT','TRP+MDL','POLAR'};
    regions = [1 2 3 5];    %% ignore region 4 = midlats+tropics
    mnames = modelnames(2:5);
    znames = zonalnames(regions);

    iFig = iFig + 1;
    figure(iFig); clf;   
      plot_correlations_200_500_800mb_RH_WV_T

    iFig = iFig + 1;
    figure(iFig); clf;   
      plot_bias_200_500_800mb_RH_WV_T

    iFig = iFig + 1;
    figure(iFig); clf;   
      plot_frac_neg0pos_200_500_800mb_RH_WV_T

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  disp('may as well hit +1 for the 5 panel (3 on top, 2 at bottom) plots and the correlation/mean/std dev 100-400,400-1000 mb'); 
%  iType5 = input('Show 5 panel plots? (-1/+1 default) : ');
%  if length(iType5) == 0
%    iType5 = +1;
%  end

  clear plotoptions
  plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';
  plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str12 = 'CLIMCAPS';
  plotoptions.str13 = 'MERRA2';
  plotoptions.str14 = 'ERA5';
  plotoptions.str15 = 'THIS WORK';
  plotoptions.xstr = 'Latitude';        plotoptions.ystr = 'Pressure (mb)';
  plotoptions.yLinearOrLog = +1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits = [100 1000];  
   
  if iType5 > 0
    i10  = find(plays100 >= 10 & plays100 <= 1000);
    i100 = find(plays100 >= 100 & plays100 <= 1000);

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
      plotoptions2x1x2.str12 = 'THIS WORK';
      plotoptions2x1x2.str13 = 'CLIMCAPS';
      plotoptions2x1x2.str21 = 'MERRA2';
      plotoptions2x1x2.str22 = 'ERA5';
      plotoptions2x1x2.cstr  = 'dRH/dt';
      profile_plots_2x1x2tiledlayout_tall(trend_rlat64,plays100,miaow11,miaow15,miaow12,miaow13,miaow14,iFig,plotoptions2x1x2);
    end

    %zall = [miaow11(:) miaow15(:) miaow12(:) miaow13(:) miaow14(:)];
    %wall = ones(size(miaow11(:)));
    miaow11x = miaow11(i100,:);      miaow15x = miaow15(i100,:);      miaow12x = miaow12(i100,:);      miaow13x = miaow13(i100,:);      miaow14x = miaow14(i100,:); 
    zall = [miaow11x(:) miaow15x(:) miaow12x(:) miaow13x(:) miaow14x(:)];  %% AIRS L3 is center of attraction [AIRSL3 UMBC CLIMCAPSL3 MERRA2 ERA5] 
    zall = zall';
    zall = zall([5 1 3 4 2],:);                                            %% ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
    if iBiasWRT_ERA5orUMBC < 0
      %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
      zall = zall([5 2 3 4 1],:);
    end
    wall = ones(size(miaow11x(:)));
    [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos,frac_neg0pos_mean_std_RH] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
    iFig = iFig + 1;
    figure(iFig); clf;   
    plot_corr_mean_std_RH
    saverates_rlat_pres.trend_rlat64        = trend_rlat64;
    saverates_rlat_pres.plays100            = plays100;
    saverates_rlat_pres.RH_z_lat.airsL3     = miaow11;
    saverates_rlat_pres.RH_z_lat.umbc       = miaow15;
    saverates_rlat_pres.RH_z_lat.climcapsL3 = miaow12;
    saverates_rlat_pres.RH_z_lat.merra2     = miaow13;
    saverates_rlat_pres.RH_z_lat.era5       = miaow14;

    %%%%%%%%%%%%%%%%%%%%%%%%%

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
      plotoptions2x1x2.str12 = 'THIS WORK';
      plotoptions2x1x2.str13 = 'CLIMCAPS';
      plotoptions2x1x2.str21 = 'MERRA2';
      plotoptions2x1x2.str22 = 'ERA5';
      plotoptions2x1x2.cstr  = 'dWVfrac/dt';
      profile_plots_2x1x2tiledlayout_tall(trend_rlat64,plays100,miaow11,miaow15,miaow12,miaow13,miaow14,iFig,plotoptions2x1x2);
    end

    %zall = [miaow11(:) miaow15(:) miaow12(:) miaow13(:) miaow14(:)];
    %wall = ones(size(miaow11(:)));
    miaow11x = miaow11(i100,:);      miaow15x = miaow15(i100,:);      miaow12x = miaow12(i100,:);      miaow13x = miaow13(i100,:);      miaow14x = miaow14(i100,:); 
    zall = [miaow11x(:) miaow15x(:) miaow12x(:) miaow13x(:) miaow14x(:)];  %% AIRS L3 is center of attraction [AIRSL3 UMBC CLIMCAPSL3 MERRA2 ERA5] 
    zall = zall';
    zall = zall([5 1 3 4 2],:);                                              %% ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
    if iBiasWRT_ERA5orUMBC < 0
      %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
      zall = zall([5 2 3 4 1],:);
    end
    wall = ones(size(miaow11x(:)));
    [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos,frac_neg0pos_mean_std_WVfrac] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
    iFig = iFig + 1;
    figure(iFig); clf;   
    plot_corr_mean_std_WVfrac

    saverates_rlat_pres.trend_rlat64            = trend_rlat64;
    saverates_rlat_pres.plays100                = plays100;
    saverates_rlat_pres.WVfrac_z_lat.airsL3     = miaow11;
    saverates_rlat_pres.WVfrac_z_lat.umbc       = miaow15;
    saverates_rlat_pres.WVfrac_z_lat.climcapsL3 = miaow12;
    saverates_rlat_pres.WVfrac_z_lat.merra2     = miaow13;
    saverates_rlat_pres.WVfrac_z_lat.era5       = miaow14;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%

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
      plotoptions2x1x2.str12 = 'THIS WORK';
      plotoptions2x1x2.str13 = 'CLIMCAPS';
      plotoptions2x1x2.str21 = 'MERRA2';
      plotoptions2x1x2.str22 = 'ERA5';
      plotoptions2x1x2.cstr  = 'dT/dt';
      profile_plots_2x1x2tiledlayout_tall(trend_rlat64,plays100,miaow11,miaow15,miaow12,miaow13,miaow14,iFig,plotoptions2x1x2);
    end

    %zall = [miaow11(:) miaow15(:) miaow12(:) miaow13(:) miaow14(:)];
    %wall = ones(size(miaow11(:)));
    miaow11x = miaow11(i10,:);      miaow15x = miaow15(i10,:);      miaow12x = miaow12(i10,:);      miaow13x = miaow13(i10,:);      miaow14x = miaow14(i10,:); 
    zall = [miaow11x(:) miaow15x(:) miaow12x(:) miaow13x(:) miaow14x(:)];  %% AIRS L3 is center of attraction [AIRSL3 UMBC CLIMCAPSL3 MERRA2 ERA5] 
    zall = zall';
    zall = zall([5 1 3 4 2],:);                                              %% ERA5 the center of attraction iBiasWRT_ERA5orUMBC > 0 [ERA5 AIRSL3 CLIMCAPSL3 MERRA2 UMBC]
    if iBiasWRT_ERA5orUMBC < 0
      %% make UMBC the center of attraction [UMBC AIRSL3 CLIMCAPSL3 MERRA2 ERA5]
      zall = zall([5 2 3 4 1],:);
    end
    wall = ones(size(miaow11x(:)));
    [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos,frac_neg0pos_mean_std_T] = corrplot_weighted_mean_stddev(zall',wall',modelnames);
    iFig = iFig + 1;
    figure(iFig); clf;   
    plot_corr_mean_std_T

    saverates_rlat_pres.trend_rlat64       = trend_rlat64;
    saverates_rlat_pres.plays100           = plays100;
    saverates_rlat_pres.T_z_lat.airsL3     = miaow11;
    saverates_rlat_pres.T_z_lat.umbc       = miaow15;
    saverates_rlat_pres.T_z_lat.climcapsL3 = miaow12;
    saverates_rlat_pres.T_z_lat.merra2     = miaow13;
    saverates_rlat_pres.T_z_lat.era5       = miaow14;

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

  %%%%%%%%%%%%%%%%%%%%%%%%%
  plot_avg_over_5_models_mean_std

  % iPrint = input('Print plots for trends paper (-1 [default/+1) : ');
  % if length(iPrint) == 0
  %   iPrint = -1;
  % end
  if strfind(strUMBC,'GULP')
    dir0 = '//home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
  elseif strfind(strUMBC,'SEQN')
    dir0 = '//home/sergio/PAPERS/SUBMITPAPERS/trends/Figs_SEQN/';
  end
  if iPrint > 0
    print_for_trendspaper
  end

  figure(40); clf; disp('in case you read in MLS a-priori retrievals .. here is dWVfrac/dt')
  miaow15 = squeeze(nanmean(reshape(umbc.waterrate,100,72,64),2));
  pcolor(trend_rlat64,plays100,miaow15); set(gca,'ydir','reverse'); ylim([100 1000]); caxis([-1 +1]*0.015);
  colormap(llsmap5); xlabel('Latitude'); ylabel('Pressure (mb)'); colorbar; shading interp

  look_at_olr_trends

  %% iWhich = 0;    %% to speed through all 11 choices with matlab nowindow
  if ~exist('iWhich')
    disp('can (+1) continue to 3 panel T(z,lat), WV(z,lat), RH(z,lat) 3x1 plots  ... and amplification or (0/default) quit or (-1) keyboard_nowindow : '); 
    iWhich = input('Enter (+1) to continue (0/default) to return (-1) to examine via keyboard_nowindow : ');
    if length(iWhich) == 0
      iWhich = 0;
    end
  end

  if iWhich == 0
    return
  elseif iWhich == -1
    keyboard_nowindow
  elseif iWhich == +1
    disp('on we go')
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  iFig = iFig + 1;
  figure(iFig); clf;
  clear plotoptions
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'GISS';    
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';    
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  aslmap_2x2tiledlayout(z11,newz32,z21,z22,iFig,plotoptions);

  iType3 = input('Show 3 panel plots? (-1 default/+1) : ');
  if length(iType3) == 0
    iType3 = -1;
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

  if iType3 > 0
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

iKB = input('separate the pentagon T(z) or WV(z) or RH(z) plots usimg saverates_rlat_pres (-1 [default]/+1) : ');
if length(iKB) == 0
  iKB = -1;
end
if iKB > 0
  split_saverates_rlat_pres_T_WV_plots
  keyboard_nowindow
end

iPrint = -1;
if iPrint > 0
  dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/PAPER17_TRENDS/Figs/';
  dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
  figure(60); aslprint([dir0 'obs_spectralavg_'        num2str(iNumYears) '_years.pdf']);
  figure(61); aslprint([dir0 'umbc_spectralavg_'       num2str(iNumYears) '_years.pdf']);
  figure(62); aslprint([dir0 'era5_spectralavg_'       num2str(iNumYears) '_years.pdf']);
  figure(63); aslprint([dir0 'merra2_spectralavg_'     num2str(iNumYears) '_years.pdf']);
  figure(64); aslprint([dir0 'airsL3_spectralavg_'     num2str(iNumYears) '_years.pdf']);
  figure(65); aslprint([dir0 'climcapsL3_spectralavg_' num2str(iNumYears) '_years.pdf']);

%%%% this is ughugh
  dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/PAPER17_TRENDS/Figs/';
  dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';

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
