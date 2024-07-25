figure(1); clf;
  iFig = 1;
  newz11 = (fUMBC_day.ptemptrend+fUMBC_night.ptemptrend)*0.5;  newz12 = (fAIRSL3_day.ptemptrend+fAIRSL3_night.ptemptrend)*0.5; newz13 = (fCLIMCAPSL3_day.ptemptrend+fCLIMCAPSL3_night.ptemptrend)*0.5;
                                       newz21 = fMERRA2.ptemptrend;  newz22 = (fERA5_day.ptemptrend+fERA5_night.ptemptrend)*0.5; 
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2));     newz13 = squeeze(nanmean(reshape(newz13,100,72,64),2)); 
              newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2));     

  plotoptions.xstr = ' '; plotoptions.ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dT/dt'; plotoptions.cmap = llsmap5;
  plotoptions.barstr = 'dT/dt [K/yr]';
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits      = [10 1000];
    plotoptions.str1 = 'AIRS\_RT';
    plotoptions.str2 = 'ERA5';
    plotoptions.barstr = 'dT/dt [K/yr]';
    profile_plots_1x2tiledlayout(rlat,plays100,newz11,newz22,iFig,plotoptions);
clear plotoptions

figure(2); clf;
  iFig = 2;

  figure(iFig); clf;
  newz11 = (fUMBC_day.gas_1trend+fUMBC_night.gas_1trend)*0.5;  newz12 = (fAIRSL3_day.gas_1trend+fAIRSL3_night.gas_1trend)*0.5; newz13 = (fCLIMCAPSL3_day.gas_1trend+fCLIMCAPSL3_night.gas_1trend)*0.5;
                                       newz21 = fMERRA2.gas_1trend;  newz22 = (fERA5_day.gas_1trend+fERA5_night.gas_1trend)*0.5; 
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2));     newz13 = squeeze(nanmean(reshape(newz13,100,72,64),2)); 
              newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2));     
  plotoptions.xstr = ' '; plotoptions.ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151/10; plotoptions.maintitle = 'dfracWV/dt'; plotoptions.cmap = llsmap5;
  plotoptions.barstr = 'dfracWV/dt [1/yr]';
  plotoptions.yLinearOrLog = +1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits      = [100 1000];
    plotoptions.str1 = 'AIRS\_RT';
    plotoptions.str2 = 'ERA5';
    plotoptions.barstr = 'dfracWV/dt [K/yr]';
    profile_plots_1x2tiledlayout(rlat,plays100,newz11,newz22,iFig,plotoptions);
clear plotoptions

%{
dir0 = '/home/sergio/PAPERS/CONFERENCES/ClimateNCAR_Mar2024_Poster/QuickFigs';
if iv3 < 6
  figure(1); aslprint_asis([dir0 '/umbc_era5_dT_dt.pdf']);
  figure(2); aslprint_asis([dir0 '/umbc_era5_dWV_dt.pdf']);
else
  figure(1); aslprint_asis([dir0 '/umbc_era5_dT_dt_v6.pdf']);
  figure(2); aslprint_asis([dir0 '/umbc_era5_dWV_dt_v6.pdf']);
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3); clf;
  iFig = 3;
  newz11 = (fUMBC_day.ptemptrend+fUMBC_night.ptemptrend)*0.5;  newz12 = (fAIRSL3_day.ptemptrend+fAIRSL3_night.ptemptrend)*0.5; newz13 = (fCLIMCAPSL3_day.ptemptrend+fCLIMCAPSL3_night.ptemptrend)*0.5;
                                       newz21 = fMERRA2.ptemptrend;  newz22 = (fERA5_day.ptemptrend+fERA5_night.ptemptrend)*0.5; 
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2));     newz13 = squeeze(nanmean(reshape(newz13,100,72,64),2)); 
              newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2));     

  plotoptions.xstr = ' '; plotoptions.ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dT/dt'; plotoptions.cmap = llsmap5;
  plotoptions.barstr = 'dT/dt [K/yr]';
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits      = [10 1000];
    plotoptions.str1 = 'AIRS\_RT';
    plotoptions.str2 = 'ERA5';
    plotoptions.str3 = 'AIRS v7 L3';
    plotoptions.barstr = 'dT/dt [K/yr]';
    profile_plots_1x3tiledlayout(rlat,plays100,newz11,newz22,newz12,iFig,plotoptions);
clear plotoptions

figure(4); clf;
  iFig = 4;

  figure(iFig); clf;
  newz11 = (fUMBC_day.gas_1trend+fUMBC_night.gas_1trend)*0.5;  newz12 = (fAIRSL3_day.gas_1trend+fAIRSL3_night.gas_1trend)*0.5; newz13 = (fCLIMCAPSL3_day.gas_1trend+fCLIMCAPSL3_night.gas_1trend)*0.5;
                                       newz21 = fMERRA2.gas_1trend;  newz22 = (fERA5_day.gas_1trend+fERA5_night.gas_1trend)*0.5; 
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2));     newz13 = squeeze(nanmean(reshape(newz13,100,72,64),2)); 
              newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2));     
  plotoptions.xstr = ' '; plotoptions.ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151/10; plotoptions.maintitle = 'dfracWV/dt'; plotoptions.cmap = llsmap5;
  plotoptions.barstr = 'dfracWV/dt [1/yr]';
  plotoptions.yLinearOrLog = +1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits      = [100 1000];
    plotoptions.str1 = 'AIRS\_RT';
    plotoptions.str2 = 'ERA5';
    plotoptions.str3 = 'AIRS v7 L3';
    plotoptions.barstr = 'dfracWV/dt [K/yr]';
    profile_plots_1x3tiledlayout(rlat,plays100,newz11,newz22,newz12,iFig,plotoptions);
clear plotoptions

figure(5); clf; pcolor(rlat,plays100,newz11); caxis([-1 +1]*0.015); colormap(llsmap5); colorbar('location','eastoutside'); xlabel('Latitude'); ylabel('Pressure (mb)');  ylim([100 1000])
  shading interp; set(gca,'ydir','reverse'); set(gca,'fontsize',12)
%% dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs_DN_Temp';; sergioprintfig([dir0 '/wvfrac_trends_20years_MLS_UA']);

%{
dir0 = '/home/sergio/PAPERS/CONFERENCES/ClimateNCAR_Mar2024_Poster/QuickFigs';
if iV3 < 6
  figure(3); aslprint_asis([dir0 '/umbc_era5_airsL3_dT_dt.pdf']);
  figure(4); aslprint_asis([dir0 '/umbc_era5_airsL3_dWV_dt.pdf']);
else
  figure(3); aslprint_asis([dir0 '/umbc_era5_airsL3_dT_dt_v6.pdf']);
  figure(4); aslprint_asis([dir0 '/umbc_era5_airsL3_dWV_dt_v6.pdf']);
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); sizefig; ; clf; aslmap(6,rlat65,rlon73,smoothn(reshape(fUMBC_night.results(:,6),72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT NIGHT');
  caxis([-1 +1]*0.151); colormap(llsmap5);

figure(7); clf; 
junkmmw0 = mmwater_rtp(h,p);
junkmmwN = load(umbc_night_file,'pert');
junkmmwN = mmwater_rtp(h,junkmmwN.pert);
junkmmwD = load(umbc_day_file,'pert');
junkmmwD = mmwater_rtp(h,junkmmwD.pert);
plot(rlat,nanmean(reshape(junkmmwN-junkmmw0,72,64),1),'b',rlat,nanmean(reshape(junkmmwD-junkmmw0,72,64),1),'g',rlat,0.5*nanmean(reshape((junkmmwN-junkmmw0)+(junkmmwD-junkmmw0),72,64),1),'k','linewidth',2); 
plotaxis2; title('UMBC night COL WV trends'); hl = legend('night','day','average','location','best','fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% this was done with editing /asl/s1/sergio/JUNK/gitjunk5/oem_climate_code/AIRS_gridded_STM_May2021_trendsonlyCLR/set_CO2_CH4_N2O_ESRL.m so it had same <<<<<< till Dec 31, 2023 MAJOR DIFF >>>>>> as the one in this dir
gitt2 = load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_gitjunk5_dec29_2023_v2.mat');  

%% this is original commit in /asl/s1/sergio/JUNK/gitjunk5/oem_climate_code/AIRS_gridded_STM_May2021_trendsonlyCLR/set_CO2_CH4_N2O_ESRL.m     Dec 29, 2023
gitt = load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_gitjunk5_dec29_2023.mat');

%% this is currently in sergio/home/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/set_CO2_CH4_N2O_ESRL.m, March 11, 2024
current = load('/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss.mat');

mmw0       = mmwater_rtp(gitt2.h,gitt2.p);
mmwgit2    = mmwater_rtp(gitt2.h,gitt2.pert);
mmwgit     = mmwater_rtp(gitt2.h,gitt.pert);
mmwcurrent = mmwater_rtp(gitt2.h,current.pert);

imagesc(gitt.xbWV'-current.xbWV'); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*1e-3);
imagesc(gitt.xbWV');               shading interp; colorbar; colormap(usa2); caxis([-1 +1]*1e-3);
imagesc(current.xbWV');            shading interp; colorbar; colormap(usa2); caxis([-1 +1]*1e-3);
imagesc(gitt.xbWV');               shading interp; colorbar; colormap(usa2); caxis([-1 +1]*1e-3);

imagesc(current.xbWV');            shading interp; colorbar; colormap(usa2); caxis([-1 +1]*1e-3);
plot(current.xbWV(30,:))
plot(current.xbWV(:,30))
plot(1:100,current.xbWV(1:100,30))
plot(1:100,current.xbWV(1:100,30),1:100,gitt.xbWV(1:100,30))
ind = 1:100;      plot(ind,current.xbWV(ind,30),ind,gitt.xbWV(ind,30))
ind = 1:100:4608; plot(ind,current.xbWV(ind,30),ind,gitt.xbWV(ind,30))
ind = 1:100:4608; plot(ind,current.xbWV(ind,30),ind,gitt.xbWV(ind,30)); plotaxis2;
ind = 1:100:4608; plot(ind,current.xbWV(ind,30),'bx-',ind,gitt.xbWV(ind,30),'ro-'); plotaxis2;
ind = 1:50:4608;  plot(ind,current.xbWV(ind,30),'bx-',ind,gitt.xbWV(ind,30),'ro-'); plotaxis2;
plot(rlat,nanmean(reshape(mmwgit2-mmwcurrent,72,64),1))
plot(rlat,nanmean(reshape(mmwgit2-mmw0,72,64),1))
plot(rlat,nanmean(reshape(mmwgit2-mmw0,72,64),1),rlat,nanmean(reshape(mmwgit-mmw0,72,64),1))
plot(rlat,nanmean(reshape(mmwgit2-mmw0,72,64),1),rlat,nanmean(reshape(mmwgit-mmw0,72,64),1)); plotaxis2;
plot(rlat,nanmean(reshape(mmwgit2-mmw0,72,64),1),rlat,nanmean(reshape(mmwcurrent-mmw0,72,64),1)); plotaxis2;
plot(rlat,nanmean(reshape(mmwgit2-mmw0,72,64),1),rlat,nanmean(reshape(mmwcurrent-mmw0,72,64),1),'linewidth',2); plotaxis2;
%}


