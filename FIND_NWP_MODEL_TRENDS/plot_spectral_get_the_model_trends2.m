disp(' ')
disp(' ... starting plot_spectral_get_the_model_trends2 ...')

addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS/
load llsmap5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncYY = nansum(cosYY,2);
iOffSet = 40;
fprintf(1,'Figs %2i : %2i showing Biases : All lats (Obs-Cal) cosine averaged, N/S Polar \n',iOffSet+1,iOffSet+3)
figure(iOffSet+1); clf;
plot(fchanx,nansum(cosYY.*obsrates.rates - cosYY.*obsrates.rates,2)./ncYY,'k',fchanx,nansum(cosYY.*obsrates.rates - cosYY.*umbcL3.umbcL3_spectral_rates,2)./ncYY,'g',...
     fchanx,nansum(cosYY.*obsrates.rates - cosYY.*airsL3.airsL3_spectral_rates,2)./ncYY,'b',fchanx,nansum(cosYY.*obsrates.rates - cosYY.*climcapsL3.climcapsL3_spectral_rates,2)./ncYY,'c',......
     fchanx,nansum(cosYY.*obsrates.rates - cosYY.*era5.era5_spectral_rates,2)./ncYY,'r',fchanx,nansum(cosYY.*obsrates.rates - cosYY.*merra2.merra2_spectral_rates,2)./ncYY,'m','linewidth',2);
 xlim([640 1640]); ylim([-1 +1]*0.03)
 plotaxis2;
 title('All lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Bias dBT/dt K/yr')

cosYYSPolar = cos(YY*pi/180)';
cosYYSPolar = ones(size(YY))';
cosYYSPolar(YY > -60) = 0;
cosYYSPolar = ones(2645,1) * cosYYSPolar;
ncYY = nansum(cosYYSPolar,2);
figure(iOffSet+2); clf;
plot(fchanx,nansum(cosYYSPolar.*obsrates.rates - cosYYSPolar.*obsrates.rates,2)./ncYY,'k',fchanx,nansum(cosYYSPolar.*obsrates.rates - cosYYSPolar.*umbcL3.umbcL3_spectral_rates,2)./ncYY,'g',...
     fchanx,nansum(cosYYSPolar.*obsrates.rates - cosYYSPolar.*airsL3.airsL3_spectral_rates,2)./ncYY,'b',fchanx,nansum(cosYYSPolar.*obsrates.rates - cosYYSPolar.*climcapsL3.climcapsL3_spectral_rates,2)./ncYY,'c',...
     fchanx,nansum(cosYYSPolar.*obsrates.rates - cosYYSPolar.*era5.era5_spectral_rates,2)./ncYY,'r',fchanx,nansum(cosYYSPolar.*obsrates.rates - cosYYSPolar.*merra2.merra2_spectral_rates,2)./ncYY,'m','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01]*4)
 plotaxis2;
 title('SPolar lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Bias dBT/dt K/yr')

cosYYNPolar = cos(YY*pi/180)';
cosYYNPolar = ones(size(YY))';
cosYYNPolar(YY < +60) = 0;
cosYYNPolar = ones(2645,1) * cosYYNPolar;
ncYY = nansum(cosYYNPolar,2);
figure(iOffSet+3); clf;
plot(fchanx,nansum(cosYYNPolar.*obsrates.rates - cosYYNPolar.*obsrates.rates,2)./ncYY,'k',fchanx,nansum(cosYYNPolar.*obsrates.rates - cosYYNPolar.*umbcL3.umbcL3_spectral_rates,2)./ncYY,'g',...
     fchanx,nansum(cosYYNPolar.*obsrates.rates - cosYYNPolar.*airsL3.airsL3_spectral_rates,2)./ncYY,'b',fchanx,nansum(cosYYNPolar.*obsrates.rates - cosYYNPolar.*climcapsL3.climcapsL3_spectral_rates,2)./ncYY,'c',...
     fchanx,nansum(cosYYNPolar.*obsrates.rates - cosYYNPolar.*era5.era5_spectral_rates,2)./ncYY,'r',fchanx,nansum(cosYYNPolar.*obsrates.rates - cosYYNPolar.*merra2.merra2_spectral_rates,2)./ncYY,'m','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01]*4)
 plotaxis2;
 title('NPolar lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Bias dBT/dt K/yr')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncYY = nansum(cosYY,2);
figure(iOffSet+4); clf;
fprintf(1,'Figs %2i : %2i showing StdDev : All lats (Obs-Cal) cosine averaged, N/S Polar \n',iOffSet+4,iOffSet+6)
plot(fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*obsrates.rates'),'k',fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*umbcL3.umbcL3_spectral_rates'),'g',...
     fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*airsL3.airsL3_spectral_rates'),'b',fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*climcapsL3.climcapsL3_spectral_rates'),'c',......
     fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*era5.era5_spectral_rates'),'r',fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*merra2.merra2_spectral_rates'),'m','linewidth',2);
 xlim([640 1640]); ylim([0 +0.01]*4)
  plotaxis2;
 title('All lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Std Dev dBT/dt K/yr')

cosYYSPolar = cos(YY*pi/180)';
cosYYSPolar = ones(size(YY))';
cosYYSPolar(YY > -60) = 0;
cosYYSPolar = ones(2645,1) * cosYYSPolar;
ncYY = nansum(cosYYSPolar,2);
figure(iOffSet+5); clf;
plot(fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*obsrates.rates'),'k',fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*umbcL3.umbcL3_spectral_rates'),'g',...
     fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*airsL3.airsL3_spectral_rates'),'b',fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*climcapsL3.climcapsL3_spectral_rates'),'c',...
     fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*era5.era5_spectral_rates'),'r',fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*merra2.merra2_spectral_rates'),'m','linewidth',2);
 xlim([640 1640]); ylim([0 +0.01]*4)
 plotaxis2;
 title('SPolar lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Std Dev dBT/dt K/yr')

cosYYNPolar = cos(YY*pi/180)';
cosYYNPolar = ones(size(YY))';
cosYYNPolar(YY < +60) = 0;
cosYYNPolar = ones(2645,1) * cosYYNPolar;
ncYY = nansum(cosYYNPolar,2);
figure(iOffSet+6); clf;
plot(fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*obsrates.rates'),'k',fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*umbcL3.umbcL3_spectral_rates'),'g',...
     fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*airsL3.airsL3_spectral_rates'),'b',fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*climcapsL3.climcapsL3_spectral_rates'),'c',...
     fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*era5.era5_spectral_rates'),'r',fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*merra2.merra2_spectral_rates'),'m','linewidth',2);
 xlim([640 1640]); ylim([0 +0.01]*4)
 plotaxis2;
 title('NPolar lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Std Dev dBT/dt K/yr')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncYY = nansum(cosYY,2);
figure(iOffSet+7); clf;
fprintf(1,'Figs %2i : %2i showing All lats (Obs and Cal) cosine averaged, N/S Polar \n',iOffSet+7,iOffSet+9)
plot(fchanx,nansum(cosYY.*obsrates.rates,2)./ncYY,'k',fchanx,nansum(cosYY.*umbcL3.umbcL3_spectral_rates,2)./ncYY,'g',...
     fchanx,nansum(cosYY.*airsL3.airsL3_spectral_rates,2)./ncYY,'b',fchanx,nansum(cosYY.*climcapsL3.climcapsL3_spectral_rates,2)./ncYY,'c',......
     fchanx,nansum(cosYY.*era5.era5_spectral_rates,2)./ncYY,'r',fchanx,nansum(cosYY.*merra2.merra2_spectral_rates,2)./ncYY,'m','linewidth',2);
 xlim([640 1640]); ylim([-1 +1]*0.04)
 plotaxis2;
 title('All lats : (Obs or Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYSPolar = cos(YY*pi/180)';
cosYYSPolar = ones(size(YY))';
cosYYSPolar(YY > -60) = 0;
cosYYSPolar = ones(2645,1) * cosYYSPolar;
ncYY = nansum(cosYYSPolar,2);
figure(iOffSet+8); clf;
plot(fchanx,nansum(cosYYSPolar.*obsrates.rates,2)./ncYY,'k',fchanx,nansum(cosYYSPolar.*umbcL3.umbcL3_spectral_rates,2)./ncYY,'g',...
     fchanx,nansum(cosYYSPolar.*airsL3.airsL3_spectral_rates,2)./ncYY,'b',fchanx,nansum(cosYYSPolar.*climcapsL3.climcapsL3_spectral_rates,2)./ncYY,'c',...
     fchanx,nansum(cosYYSPolar.*era5.era5_spectral_rates,2)./ncYY,'r',fchanx,nansum(cosYYSPolar.*merra2.merra2_spectral_rates,2)./ncYY,'m','linewidth',2);
 xlim([640 1640]); ylim([-1 +1]*0.04)
 plotaxis2;
 title('SPolar lats : (Obs or Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYNPolar = cos(YY*pi/180)';
cosYYNPolar = ones(size(YY))';
cosYYNPolar(YY < +60) = 0;
cosYYNPolar = ones(2645,1) * cosYYNPolar;
ncYY = nansum(cosYYNPolar,2);
figure(iOffSet+9); clf;
plot(fchanx,nansum(cosYYNPolar.*obsrates.rates,2)./ncYY,'k',fchanx,nansum(cosYYNPolar.*umbcL3.umbcL3_spectral_rates,2)./ncYY,'g',...
     fchanx,nansum(cosYYNPolar.*airsL3.airsL3_spectral_rates,2)./ncYY,'b',fchanx,nansum(cosYYNPolar.*climcapsL3.climcapsL3_spectral_rates,2)./ncYY,'c',...
     fchanx,nansum(cosYYNPolar.*era5.era5_spectral_rates,2)./ncYY,'r',fchanx,nansum(cosYYNPolar.*merra2.merra2_spectral_rates,2)./ncYY,'m','linewidth',2);
 xlim([640 1640]); ylim([-1 +1]*0.04)
 plotaxis2;
 title('NPolar lats : (Obs or Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('showing *** ALL *** 4608 lats/lons Figs 50-56 so you may see a lot a lot of crazy structure!!!!')
fprintf(1,'iOffSet = %3i \n',iOffSet)
fprintf(1,'iOffSet+10 = %3i \n',iOffSet+10)
figure(iOffSet+10); clf;
  pcolor(fchanx,YY,obsrates.rates');                       shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('Obs Rates \newline all72 lonbins/latbin'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+11); clf;
  pcolor(fchanx,YY,umbcL3.umbcL3_spectral_rates');         shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('UMBC Rates \newline all72 lonbins/latbin'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+12); clf;
  pcolor(fchanx,YY,era5.era5_spectral_rates');             shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('ERA5 Rates \newline all72 lonbins/latbin'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+13); clf;
  pcolor(fchanx,YY,merra2.merra2_spectral_rates');         shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('MERRA2 Rates \newline all72 lonbins/latbin'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+14); clf;
  pcolor(fchanx,YY,airsL3.airsL3_spectral_rates');         shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('AIRS L3 Rates \newline all72 lonbins/latbin'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+15); clf;
  pcolor(fchanx,YY,climcapsL3.climcapsL3_spectral_rates'); shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('CLIMCAPS L3 Rates \newline all72 lonbins/latbin'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])

keyboard_nowindow
iX = 1 : 2645;
iX = find(fchanx <= 1608);
z11 = squeeze(nanmean(reshape(obsrates.rates',72,64,2645),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates',72,64,2645),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates',72,64,2645),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates',72,64,2645),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates',72,64,2645),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates',72,64,2645),1));
  plotoptions6.str11 = 'AIRS L1C';   plotoptions6.str12 = 'UMBC';
  plotoptions6.str21 = 'ERA52';      plotoptions6.str22 = 'MERRA2';
  plotoptions6.str31 = 'AIRS L3';   plotoptions6.str32 = 'CLIMCAPS';
  plotoptions6.ystr = 'Latitude';   plotoptions6.xstr = 'Wavenumber / [cm-1]';
  plotoptions6.maintitle = 'dBT/dt / [K/yr]';
  plotoptions6.cx = [-1 +1]*0.25; plotoptions6.cmap = llsmap5; plotoptions6.yReverseDir = -1; plotoptions6.yLinearOrLog = +1;
  tiled_3x2layout(z11(:,iX),z12(:,iX),z21(:,iX),z22(:,iX),z31(:,iX),z32(:,iX),iOffSet+16,plotoptions6,fchanx(iX),unique(YY));

%{
!!! takes too long!!!
iX = 2645;
z11 = obsrates.rates(1:iX,:);
z12 = umbcL3.umbcL3_spectral_rates(1:iX,:);;
z21 = era5.era5_spectral_rates(1:iX,:);;
z22 = merra2.merra2_spectral_rates(1:iX,:);;
z31 = airsL3.airsL3_spectral_rates(1:iX,:);;
z32 = climcapsL3.climcapsL3_spectral_rates(1:iX,:);;
  plotoptions6.str11 = 'AIRS L1C';   plotoptions6.str12 = 'UMBC';
  plotoptions6.str21 = 'ERA52';      plotoptions6.str22 = 'MERRA2';
  plotoptions6.str31 = 'AIRS L3';   plotoptions6.str32 = 'CLIMCAPS';
  plotoptions6.ystr = 'Latitude';   plotoptions6.xstr = 'Wavenumber / [cm-1]';
  plotoptions6.maintitle = 'dBT/dt / [K/yr]';
  plotoptions6.cx = [-1 +1]*0.25; plotoptions6.cmap = llsmap5; plotoptions6.yReverseDir = -1; plotoptions6.yLinearOrLog = +1;
  tiled_3x2layout(z11',z12',z21',z22',z31',z32',iOffSet+16,plotoptions6,fchanx(1:iX),YY);
%}

%zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
%figure(iOffSet+17); clf
%%varnames = {'OBS','UMBC','ERA5','MERRA2','AIRS L3','CLIMCAPS'};
%%[R,Pvalue] = corrplot(zall',Varnames = varnames)
%%hfig = gcf;
%%haxes = findobj(hfig, 'Type', 'Axes');
%%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('showing *** zonal averages *** in Figs 60-66 : zonally averaged 64 latbins so cleaner/less structure than the previous five')
fprintf(1,'iOffSet = %3i \n',iOffSet)
fprintf(1,'iOffSet+20 = %3i \n',iOffSet+20)
figure(iOffSet+20); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(obsrates.rates',72,64,2645),1)));                       
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('Obs Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+21); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates',72,64,2645),1)));         
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('UMBC Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+22); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(era5.era5_spectral_rates',72,64,2645),1)));             
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('ERA5 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+23); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(merra2.merra2_spectral_rates',72,64,2645),1)));         
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('MERRA2 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+24); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates',72,64,2645),1)));         
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('AIRS L3 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+25); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates',72,64,2645),1))); 
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('CLIMCAPS L3 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])

savetherates = struct;
savetherates.rlatYY   = reshape(YY,72,64);
savetherates.rlonXX   = XX;
savetherates.obs_rates  = obsrates.rates;
savetherates.umbc_rates = umbcL3.umbcL3_spectral_rates;
savetherates.era5_rates = era5.era5_spectral_rates;
savetherates.fchanx     = fchanx;
savetherates.latavg     = nanmean(reshape(YY,72,64),1);
savetherates.lfavg      = nanmean(reshape(landfrac,72,64),1);
savetherates.obs        = squeeze(nanmean(reshape(obsrates.rates',72,64,2645),1));
savetherates.umbc       = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates',72,64,2645),1));
savetherates.era5       = squeeze(nanmean(reshape(era5.era5_spectral_rates',72,64,2645),1));
savetherates.merra2     = squeeze(nanmean(reshape(merra2.merra2_spectral_rates',72,64,2645),1));
savetherates.alrsL3     = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates',72,64,2645),1));
savetherates.climcapsL3 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates',72,64,2645),1));
savetherates.comment1   = 'strUMBC = ''/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2D.mat'';   iNumYears = 20;';
savetherates.comment2   = 'driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends(strUMBC,iNumYears);';
%% save spectral_rate_avgs_umbc_obs_era5_merra2_airsL3_climcapsL3.mat savetherates

data_fig5c.comment  = 'see ~/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/plot_spectral_get_the_model_trends2.m called by ~/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends(str,iNum)';
data_fig5c.comment2 = 'strUMBC = ''/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat''; iNumYears = 20;';
data_fig5c.obs   = savetherates.obs;
data_fig5c.era5  = savetherates.era5;
data_fig5c.rlat  = savetherates.latavg;
data_fig5c.vchan = savetherates.fchanx;
%% save /umbc/xfs2/strow/asl/s1/sergio/home/git/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig5c.mat data_fig5c

%%%%%%%%%%%%%%%%%%%%%%%%%
oldmat = load('spectral_rate_avgs_umbc_obs_era5_merra2_airsL3_climcapsL3.mat');  %% saved Jan 20, 2025

figure(1); clf; pcolor(data_fig5c.vchan,data_fig5c.rlat,data_fig5c.obs);  shading interp; ylabel('Latitude (deg)'); xlabel('Wavenumber cm^{-1}'); xlim([650 1620]); caxis([-1 +1]*0.1); colormap(usa2); colorbar
figure(3); clf; pcolor(data_fig5c.vchan,data_fig5c.rlat,data_fig5c.era5); shading interp; ylabel('Latitude (deg)'); xlabel('Wavenumber cm^{-1}'); xlim([650 1620]); caxis([-1 +1]*0.1); colormap(usa2); colorbar

figure(4); clf; pcolor(oldmat.savetherates.fchanx,oldmat.savetherates.latavg,oldmat.savetherates.obs);  shading interp; ylabel('Latitude (deg)'); xlabel('Wavenumber cm^{-1}'); xlim([650 1620]); caxis([-1 +1]*0.1); colormap(usa2); colorbar
figure(5); clf; pcolor(oldmat.savetherates.fchanx,oldmat.savetherates.latavg,oldmat.savetherates.era5); shading interp; ylabel('Latitude (deg)'); xlabel('Wavenumber cm^{-1}'); xlim([650 1620]); caxis([-1 +1]*0.1); colormap(usa2); colorbar

figure(6); clf; pcolor(oldmat.savetherates.fchanx,oldmat.savetherates.latavg,oldmat.savetherates.obs-data_fig5c.obs);   shading interp; ylabel('Latitude (deg)'); xlabel('Wavenumber cm^{-1}'); xlim([650 1620]); caxis([-1 +1]*0.1); colormap(usa2); colorbar
figure(7); clf; pcolor(oldmat.savetherates.fchanx,oldmat.savetherates.latavg,oldmat.savetherates.era5-data_fig5c.era5); shading interp; ylabel('Latitude (deg)'); xlabel('Wavenumber cm^{-1}'); xlim([650 1620]); caxis([-1 +1]*0.1); colormap(usa2); colorbar

figure(8); clf; pcolor(oldmat.savetherates.fchanx,oldmat.savetherates.latavg,smoothn(oldmat.savetherates.obs,1));  shading interp; ylabel('Latitude (deg)'); xlabel('Wavenumber cm^{-1}'); xlim([650 1620]); caxis([-1 +1]*0.1); colormap(usa2); colorbar
figure(9); clf; pcolor(oldmat.savetherates.fchanx,oldmat.savetherates.latavg,smoothn(oldmat.savetherates.era5,1)); shading interp; ylabel('Latitude (deg)'); xlabel('Wavenumber cm^{-1}'); xlim([650 1620]); caxis([-1 +1]*0.1); colormap(usa2); colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%

disp('keyboard_nowindow so you can make fig5 for trends paper')
keyboard_nowindow

quick_tropical_perts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk = input('do the correlations of spectra : takes a while to make the plots but some plots are in eg trends paper ... not needed for OLR calcs .... (-1 [default]/+1) : ');
if length(junk) == 0
  junk = -1;
end

if junk > 0
  do_correlations_of_spectra
end

disp(' ... ending plot_spectral_get_the_model_trends2 ...')
disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure(iOffSet+30); clf;
%  wah = obsrates.rates; wah = reshape(wah,2645,72,64); wah = squeeze(nanmean(wah,2));
%  pcolor(fchanx,nanmean(reshape(YY,72,64),1),wah);
%  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('Obs Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
%keyboard_nowindow
