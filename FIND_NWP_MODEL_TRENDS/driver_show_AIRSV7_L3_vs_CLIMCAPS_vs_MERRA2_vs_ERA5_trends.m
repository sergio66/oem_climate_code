function [era5,merra2,airsL3,climcapsL3,umbc,thecorr,amp] = driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends(strUMBC);

set(0,'DefaultaxesLineWidth',1);
set(0,'DefaultaxesFontSize',16);

amp = [];

if nargin == 0
  strUMBC = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithMLSL3_uncX100_50fatlayers_AIRSL3_ERA5_CMIP6_globalSSTfeedback.mat';
  strUMBC = '/asl/s1/sergio/JUNK/test9_guessstartWV_Vers1_march22_2023.mat';
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q05_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_feedback.mat';              %% not too bad at lower atm/polar!!!!
  strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q05_newERA5_2021jacs_startwith_MLSL3_TOA_guessWV_dRH_zero_bot_50fatlayers.mat';  ## too much oomph at gnd : use MLS L3 TOA and dRH/dt = 0 at bottom
  strUMBC = [];
  iUMBC = -1;
  umbc = [];
  thecorr = [];
else
  iUMBC = +1;
end

%% see How closely do changes in surface and column water vapor follow Clausius-Clapeyron scaling in climate-change simulations?
%% P A O’Gorman, C J Muller, https://core.ac.uk/download/pdf/4426849.pdf, or saved in PDF/change_of_RH_with_stemp_GCM_PGorman.pdf
%% P A O'Gorman and C J Muller 2010 Environ. Res. Lett. 5 025207 DOI 10.1088/1748-9326/5/2/025207

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/NANROUTINES
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools
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
iX = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iX > 0
  load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat
else
  load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug2022_20yr_desc.mat
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

load /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat

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

load MERRA2_atm_data_2002_09_to_2022_08_trends_desc.mat

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

load ERA5_atm_data_2002_09_to_2022_08_trends_desc.mat

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
  fprintf(1,'loading in %s \n',strUMBC);
  umbcX = load(strUMBC,'deltaRH','deltaT','fracWV','results');

  strGISS = 'ChrisHTrends/giss_trends_2002_2022.mat';
  fprintf(1,'loading in %s \n',strGISS);
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
  addpath /asl/matlib/plotutils
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

  plotoptions.maintitle = 'Compare UMBC/ERA5 or UMBC-ERA5'; 
  plotoptions.cx = [-1 +1]*0.151; 
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
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

  clear plotoptions
  plotoptions.cx = [-1 +1]*0.151;  plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3'; plotoptions.str13 = 'MERRA2';    
  plotoptions.str21 = 'ERA5';      plotoptions.str22 = 'UMBC';        plotoptions.str23 = 'GISS';
  plotoptions.xstr = ' ';          plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  z31 = z11x; z32 = z32;
  aslmap_2x3tiledlayout(z11,z12,z21,z22,z31,z32,iFig,plotoptions);

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
    profile_plots_1x5tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow13,miaow14,miaow15,iFig,plotoptions);
    
    iFig = iFig + 1;
    figure(iFig); clf;   
    plotoptions.maintitle = 'd(log(WV))/dt [1/yr]'; 
    plotoptions.cx = [-1 +1]*0.015; 
    miaow11 = squeeze(nanmean(permute(airsL3.waterrate,[3 1 2]),2));
    miaow12 = squeeze(nanmean(permute(climcapsL3.waterrate,[3 1 2]),2));
    miaow13 = squeeze(nanmean(reshape(merra2.waterrate,100,72,64),2));
    miaow14 = squeeze(nanmean(reshape(era5.waterrate,100,72,64),2));
    miaow15 = squeeze(nanmean(reshape(umbc.waterrate,100,72,64),2));
    profile_plots_1x5tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow13,miaow14,miaow15,iFig,plotoptions);
    
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
    profile_plots_1x5tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow13,miaow14,miaow15,iFig,plotoptions);
  
    %{ 
    %  quick_print_figs_compare_trends.m
    addpath /asl/matlib/plotutils
    figure(16); aslprint('QuickFigs/Strow_March2023/skt_trend.pdf');
    figure(17); aslprint('QuickFigs/Strow_March2023/rh_trend.pdf');
    figure(18); aslprint('QuickFigs/Strow_March2023/wv_trend.pdf');
    figure(19); aslprint('QuickFigs/Strow_March2023/tz_trend.pdf');
    %}
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
  end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iFigiAmp = iFig;

iAmp = input('Do amplification??? (-9999 to stop/ 0 [default] to go on) : ');
if length(iAmp) == 0
  iAmp = 0;
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

