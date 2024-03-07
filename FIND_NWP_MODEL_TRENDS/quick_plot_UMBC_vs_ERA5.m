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
figure(5); sizefig; ; clf; aslmap(5,rlat65,rlon73,smoothn(reshape(fUMBC_night.results(:,6),72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT NIGHT');
  caxis([-1 +1]*0.151); colormap(llsmap5);
