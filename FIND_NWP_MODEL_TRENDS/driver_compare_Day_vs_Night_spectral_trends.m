addpath /asl/matlib/plotutils
addpath /asl/matlib/rtptools
addpath /asl/matlib/maps
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/science/
addpath /home/sergio/MATLABCODE/SHOWSTATS
addpath /home/sergio/MATLABCODE/NANROUTINES
addpath /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/

load llsmap5

[h,ha,p,pa] = rtpread('summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.op.rtp');
[salti,landfrac] =  usgs_deg10_dem(p.rlat,p.rlon);
p.landfrac = landfrac;
f = h.vchan;
fchanx = f;

do_XX_YY_from_X_Y

%% >>>>>> to make "do_correlations_of_spectra" work <<<<<<<<<<
YY = YY';
landfrac = landfrac';
%% >>>>>> to make "do_correlations_of_spectra" work <<<<<<<<<<

mu = cos(YY*pi/180);
load llsmap5

muYY2645 = (ones(2645,1) * cos(YY'*pi/180));
bad = find(p.landfrac == 0); muYYLand2645 = muYY2645;  muYYLand2645(:,bad) = NaN;
bad = find(p.landfrac == 1); muYYOcean2645 = muYY2645; muYYOcean2645(:,bad) = NaN;
cosYY = muYY2645;

oceanX = find(p.landfrac == 0);
landX  = find(p.landfrac == 1);

tropics = find(abs(YY) <= 30);          
midlatsNtropics = find(abs(YY) <= 60); 
midlats = setdiff(midlatsNtropics,tropics); 
poles = find(abs(YY) > 60);             

tropics = nan(size(YY));         junk = find(abs(YY) <= 30);                 tropics(junk) = 1;
midlatsNtropics = nan(size(YY)); junk = find(abs(YY) <= 60);                 midlatsNtropics(junk) = 1;
midlats = nan(size(YY));         junk = find(abs(YY) > 30 & abs(YY) <= 60);  midlats(junk) = 1;
poles = nan(size(YY));           junk = find(abs(YY) > 60);                  poles(junk) = 1;

nmidlats = nan(size(YY));         junk = find(YY > 30 & YY <= 60);           nmidlats(junk) = 1;
npoles = nan(size(YY));           junk = find(YY > 60);                      npoles(junk) = 1;
smidlats = nan(size(YY));         junk = find(YY >= -60 & YY <= -30);        smidlats(junk) = 1;
spoles = nan(size(YY));           junk = find(YY < -60);                     spoles(junk) = 1;

land = nan(size(p.landfrac')); land(p.landfrac' == 1) = 1; 
ocean = nan(size(p.landfrac')); ocean(p.landfrac' == 0) = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%

tropics_o = nan(size(YY));         junk = find(p.landfrac' == 0 &abs(YY) <= 30);                 tropics_o(junk) = 1;
midlatsNtropics_o = nan(size(YY)); junk = find(p.landfrac' == 0 &abs(YY) <= 60);                 midlatsNtropics_o(junk) = 1;
midlats_o = nan(size(YY));         junk = find(p.landfrac' == 0 &abs(YY) > 30 & abs(YY) <= 60);  midlats_o(junk) = 1;
poles_o = nan(size(YY));           junk = find(p.landfrac' == 0 &abs(YY) > 60);                  poles_o(junk) = 1;

nmidlats_o = nan(size(YY));         junk = find(p.landfrac' == 0 &YY > 30 & YY <= 60);           nmidlats_o(junk) = 1;
npoles_o = nan(size(YY));           junk = find(p.landfrac' == 0 &YY > 60);                      npoles_o(junk) = 1;
smidlats_o = nan(size(YY));         junk = find(p.landfrac' == 0 &YY >= -60 & YY <= -30);        smidlats_o(junk) = 1;
spoles_o = nan(size(YY));           junk = find(p.landfrac' == 0 &YY < -60);                     spoles_o(junk) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%

tropics_l = nan(size(YY));         junk = find(p.landfrac' == 1 &abs(YY) <= 30);                 tropics_l(junk) = 1;
midlatsNtropics_l = nan(size(YY)); junk = find(p.landfrac' == 1 &abs(YY) <= 60);                 midlatsNtropics_l(junk) = 1;
midlats_l = nan(size(YY));         junk = find(p.landfrac' == 1 &abs(YY) > 30 & abs(YY) <= 60);  midlats_l(junk) = 1;
poles_l = nan(size(YY));           junk = find(p.landfrac' == 1 &abs(YY) > 60);                  poles_l(junk) = 1;

nmidlats_l = nan(size(YY));         junk = find(p.landfrac' == 1 &YY > 30 & YY <= 60);           nmidlats_l(junk) = 1;
npoles_l = nan(size(YY));           junk = find(p.landfrac' == 1 &YY > 60);                      npoles_l(junk) = 1;
smidlats_l = nan(size(YY));         junk = find(p.landfrac' == 1 &YY >= -60 & YY <= -30);        smidlats_l(junk) = 1;
spoles_l = nan(size(YY));           junk = find(p.landfrac' == 1 &YY < -60);                     spoles_l(junk) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir0 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/';
if ~exist('fMERRA2')
  fMERRA2 = load([dir0 'MERRA2_SARTA_SPECTRAL_RATES/all_4608_2002_09_2022_08.mat']);
  
  fERA5_d = load([dir0 'ERA5_SARTA_SPECTRAL_RATES/all_4608_asc_2002_09_2022_08.mat']);
  fERA5_n = load([dir0 'ERA5_SARTA_SPECTRAL_RATES/all_4608_desc_2002_09_2022_08.mat']);
  
  fAIRSL3_d = load([dir0 'AIRSL3_SARTA_SPECTRAL_RATES/all_4608_asc_2002_09_2022_08.mat']);
  fAIRSL3_n = load([dir0 'AIRSL3_SARTA_SPECTRAL_RATES/all_4608_desc_2002_09_2022_08.mat']);
  
  fCLIMCAPSL3_d = load([dir0 'CLIMCAPSL3_SARTA_SPECTRAL_RATES/all_4608_asc_2002_09_2022_08.mat']);
  fCLIMCAPSL3_n = load([dir0 'CLIMCAPSL3_SARTA_SPECTRAL_RATES/all_4608_desc_2002_09_2022_08.mat']);

  get_umbc_day_night_name
  
  fUMBC_d = load(umbc_day_file,'rates','fits');
  obs_d = fUMBC_d.rates;
  fUMBC_d = fUMBC_d.fits;
  
  fUMBC_n = load(umbc_night_file,'rates','fits');
  obs_n = fUMBC_n.rates;
  fUMBC_n = fUMBC_n.fits;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii = 0;
% ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(obs_d(1520,:),72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('dBT1231/dt : AIRS L1C DAY');
% ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(obs_n(1520,:),72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('dBT1231/dt : AIRS L1C NIGHT');
% 
% ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fUMBC_d(1520,:),72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('dBT1231/dt : CHIRP\_A DAY');
% ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fUMBC_n(1520,:),72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('dBT1231/dt : CHIRP\_A NIGHT');
% 
% ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fERA5_d.trend(1520,:),72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('dBT1231/dt : ERA5 DAY');
% ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fERA5_n.trend(1520,:),72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('dBT1231/dt : ERA5 NIGHT');
% 
% ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fAIRSL3_d.trend(1520,:),72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('dBT1231/dt : AIRSL3 DAY');
% ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fAIRSL3_n.trend(1520,:),72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('dBT1231/dt : AIRSL3 NIGHT');
% 
% ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fCLIMCAPSL3_d.trend(1520,:),72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('dBT1231/dt : CLIMCAPSL3 DAY');
% ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fCLIMCAPSL3_n.trend(1520,:),72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('dBT1231/dt : CLIMCAPSL3 NIGHT');

plotoptions.cmap = llsmap5; plotoptions.cx = [-1 +1]*0.15; 
plotoptions.str11 = 'AIRS L1C \newline(D)'; plotoptions.str12 = 'CHIRP\_A \newline(D)'; plotoptions.str13 = 'AIRS L3 \newline(D)'; plotoptions.str14 = 'CLIMCAPS L3 \newline(D)';plotoptions.str15 = 'ERA5 \newline(D)';  
plotoptions.str21 = '(N)';    plotoptions.str22 = '(N)';        plotoptions.str23 = '(N)';  plotoptions.str24 = '(N)';     plotoptions.str25 = '(N)'; 
aslmap_2x5tiledlayout(obs_d(1520,:),fUMBC_d(1520,:),fAIRSL3_d.trend(1520,:),fCLIMCAPSL3_d.trend(1520,:),fERA5_d.trend(1520,:),...
                      obs_n(1520,:),fUMBC_n(1520,:),fAIRSL3_n.trend(1520,:),fCLIMCAPSL3_n.trend(1520,:),fERA5_d.trend(1520,:),...
                      1,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%

%% /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/
%% 
%%% /home/sergio/MATLABCODE/RESET_SST_colWV/GOODVERS_2020/reset_Tsurf_water.m
% Two day-and-night channels to use for determining total water
% NOTE: channel IDs must be specified in ascending order
%dayWchanlist   = [ 1282  1291 ];  % 1226.679 cm-1 and 1231.330 cm-1
% The corresponding minimum passing dbtc*ddbtx for these 2 channels
%dayWminpass = 0.20;

%% tile_fits_quantiles.m uses 1228 and 1231 for SKT hmmm
plotoptions.cmap = llsmap5; plotoptions.cx = [-1 +1]*0.02; 
i1231 = 1520;
i1227 = find(fchanx >= 1227,1); i1227 = 1513;
i1227 = find(fchanx >= 1226,1); i1227 = 1511;

smX = 3;
figure(2); clf
  moo0 = -obs_d(i1227,:)+obs_d(i1231,:);                              moo0 = nanmean(reshape(moo0,72,64),1);
  moo1 = -fUMBC_d(i1227,:)+fUMBC_d(i1231,:);                          moo1 = nanmean(reshape(moo1,72,64),1);
  moo2 = -fAIRSL3_d.trend(i1227,:)+fAIRSL3_d.trend(i1231,:);          moo2 = nanmean(reshape(moo2,72,64),1);
  moo3 = -fCLIMCAPSL3_d.trend(i1227,:)+fCLIMCAPSL3_d.trend(i1231,:);  moo3 = nanmean(reshape(moo3,72,64),1);
  moo4 = -fERA5_d.trend(i1227,:)+fERA5_d.trend(i1231,:);              moo4 = nanmean(reshape(moo4,72,64),1);
  moo5 = -fMERRA2.trend(i1227,:)+fMERRA2.trend(i1231,:);              moo5 = nanmean(reshape(moo5,72,64),1);
plot(rlat,smooth(moo0,smX),'color',[1 1 1]*0.6,'linewidth',4); hold on
plot(rlat,smooth(moo1,smX),'k',rlat,smooth(moo2,smX),'b',rlat,smooth(moo3,smX),'g',rlat,smooth(moo4,smX),'r',rlat,smooth(moo5,smX),'m','linewidth',4); 
hold off; plotaxis2; 
%title('BT1226-BT1231 = colWV'); 
hl = legend('AIRS L1C','CHIRP\_A','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',10);
xlim([-1 +1]*90); xlabel('Latitude'); ylabel('\delta BT1226-BT1231 trend [K/yr]')

aslmap_2x5tiledlayout(-obs_d(i1227,:)+obs_d(i1231,:),-fUMBC_d(i1227,:)+fUMBC_d(i1231,:),-fAIRSL3_d.trend(i1227,:)+fAIRSL3_d.trend(i1231,:),...
                      -fCLIMCAPSL3_d.trend(i1227,:)+fCLIMCAPSL3_d.trend(i1231,:),-fERA5_d.trend(i1227,:)+fERA5_d.trend(i1231,:),...
                      -obs_n(i1227,:)+obs_n(i1231,:),-fUMBC_n(i1227,:)+fUMBC_n(i1231,:),-fAIRSL3_n.trend(i1227,:)+fAIRSL3_n.trend(i1231,:),...
                      -fCLIMCAPSL3_n.trend(i1227,:)+fCLIMCAPSL3_n.trend(i1231,:),-fERA5_n.trend(i1227,:)+fERA5_n.trend(i1231,:),...
                      3,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%

plotoptions.cmap = llsmap5; plotoptions.cx = [-1 +1]*0.15/2; 
i1419 = find(fchanx >= 1419,1);
aslmap_2x5tiledlayout(obs_d(i1419,:),fUMBC_d(i1419,:),fAIRSL3_d.trend(i1419,:),fCLIMCAPSL3_d.trend(i1419,:),fERA5_d.trend(i1419,:),...
                      obs_n(i1419,:),fUMBC_n(i1419,:),fAIRSL3_n.trend(i1419,:),fCLIMCAPSL3_n.trend(i1419,:),fERA5_d.trend(i1419,:),...
                      4,plotoptions);

plotoptions.cmap = llsmap5; plotoptions.cx = [-1 +1]*0.15/2; 
i1519 = find(fchanx >= 1519,1);
aslmap_2x5tiledlayout(obs_d(i1519,:),fUMBC_d(i1519,:),fAIRSL3_d.trend(i1519,:),fCLIMCAPSL3_d.trend(i1519,:),fERA5_d.trend(i1519,:),...
                      obs_n(i1519,:),fUMBC_n(i1519,:),fAIRSL3_n.trend(i1519,:),fCLIMCAPSL3_n.trend(i1519,:),fERA5_d.trend(i1519,:),...
                      5,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFrac = input('Show Day (-1) or Night (+1/default) or average (0) or difference (-10) or stop (-9999): ');
if length(iFrac) == 0
  iFrac = 1;
end

while iFrac > -11
  
  if iFrac == +1
    nightfrac = 1;
    dayfrac = 0;
  elseif iFrac == -1
    nightfrac = 0;
    dayfrac = 1;
  elseif iFrac == 0
    nightfrac = 0.5;
    dayfrac = 0.5;
  elseif iFrac == -10
    nightfrac = -1.0;
    dayfrac = +1.0;
  end
  
  obsrates.rates                       = nightfrac * obs_n               + dayfrac * obs_d;
  umbcL3.umbcL3_spectral_rates         = nightfrac * fUMBC_n             + dayfrac * fUMBC_d;
  airsL3.airsL3_spectral_rates         = nightfrac * fAIRSL3_n.trend     + dayfrac * fAIRSL3_d.trend;
  climcapsL3.climcapsL3_spectral_rates = nightfrac * fCLIMCAPSL3_n.trend + dayfrac * fCLIMCAPSL3_d.trend;
  era5.era5_spectral_rates             = nightfrac * fERA5_n.trend       + dayfrac * fERA5_d.trend;
  if iFrac > -10
    merra2.merra2_spectral_rates         = fMERRA2.trend;
  else
    merra2.merra2_spectral_rates         = fMERRA2.trend * 0; 
  end

%  obsXM    = nansum(obsX .* muYY2645,2) ./ nansum(muYY2645,2);
%  umbcXM   = nansum(obsX .* muYY2645,2) ./ nansum(muYY2645,2);
%  aL3XM    = nansum(aL3X .* muYY2645,2) ./ nansum(muYY2645,2);
%  cL3XM    = nansum(cL3X .* muYY2645,2) ./ nansum(muYY2645,2);
%  eraXM    = nansum(eraX .* muYY2645,2) ./ nansum(muYY2645,2);
%  merraXM  = nansum(merraX .* muYY2645,2) ./ nansum(muYY2645,2);

%  obsrates.rates = obsX;
%  era5.era5_spectral_rates = eraX;
%  merra2.merra2_spectral_rates = merraX;
%  airsL3.airsL3_spectral_rates = aL3X;
%  climcapsL3.climcapsL3_spectral_rates = cL3X;
%  umbcL3.umbcL3_spectral_rates = umbcX;
  
  monitor_memory_whos;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  iOffSet = 40;    %% see driver_get_the_model_trends.m --> plot_spectral_get_the_model_trends2.m
  iOffSet = -1;    %% do not do this if already saved off the (D+N)/2 for the paper. besides you can see it in Fig 4 from spectral_comparisons_ocean_vs_land_obs_chirp_nwp_l3_v0
  if iOffSet == 40
    do_correlations_of_spectra
  end

  spectral_comparisons_ocean_vs_land_obs_chirp_nwp_l3_v0             %%% does all           : A, T, ML, P        1 plot
    spectral_comparisons_ocean_vs_land_obs_chirp_nwp_l3_v0_land      %%% does all           : A, T, ML, P        1 plot
    spectral_comparisons_ocean_vs_land_obs_chirp_nwp_l3_v0_ocean     %%% does all           : A, T, ML, P        1 plot
  %spectral_comparisons_ocean_vs_land_obs_chirp_nwp_l3_v1            %%% does Ocean vs Land : T, NML,SML,NP,SP   6 different plots

  disp(' ')
  disp(' ')
  disp(' ')
  disp(' ')
  iFrac = input('Show Day (-1) or Night (+1/default) or average (0) or difference (-10) or stop (-9999): ');
  if length(iFrac) == 0
    iFrac = 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
print_Day_vs_Night_spectral_trendspaper
%}

