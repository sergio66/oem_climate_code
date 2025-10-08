figure(65); figure(68); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(165); clf
data_fig6.era5_plays = era5_plays;
data_fig6.rlat       = rlat;
data_fig6.retr_WV = squeeze(nanmean(reshape(fracWV,101,72,64),2)); data_fig6.retr_WV = data_fig6.retr_WV(1:100,:);
data_fig6.era5_WV = squeeze(nanmean(reshape(era5_wvrate,100,72,64),2)); data_fig6.era5_WV = data_fig6.era5_WV(1:100,:);
clear plotoptions
  plotoptions.cx = [-1 +1]*0.0101;
  plotoptions.cmap = llsmap5;
  %plotoptions.maintitle = ' ';
  %plotoptions.maintitle = 'dfracWV(z,lat)/dt'; 
  plotoptions.yLimits = [100 1000]; plotoptions.xLimits = [-85 +85];
  plotoptions.yReverseDir = +1;
  plotoptions.yLinearOrLog = +1;
  plotoptions.ystr = 'Pressure [mb]';
  plotoptions.xstr = 'Latitude [deg]';
  plotoptions.smooth = iSmooth;
profile_plots_2x1tiledlayout(data_fig6.rlat,data_fig6.era5_plays,data_fig6.era5_WV,data_fig6.retr_WV,165,plotoptions);
plotoptionsWV = plotoptions;

figure(168); clf
data_fig6.retr_T = squeeze(nanmean(reshape(deltaT,101,72,64),2)); data_fig6.retr_T = data_fig6.retr_T(1:100,:);
data_fig6.era5_T = squeeze(nanmean(reshape(era5_tzrate,100,72,64),2)); data_fig6.era5_T = data_fig6.era5_T(1:100,:);
clear plotoptions
  plotoptions.cx = [-1 +1]*0.101;
  plotoptions.cmap = llsmap5;
  %plotoptions.maintitle = ' ';
  %plotoptions.maintitle = 'dT(z,lat)/dt'; 
  plotoptions.yLimits = [10 1000]; plotoptions.xLimits = [-85 +85];
  plotoptions.yReverseDir = +1;
  plotoptions.yLinearOrLog = -1;
  plotoptions.ystr = 'Pressure [mb]';
  plotoptions.xstr = 'Latitude [deg]';
  plotoptions.smooth = iSmooth;
profile_plots_2x1tiledlayout(data_fig6.rlat,data_fig6.era5_plays,data_fig6.era5_T,data_fig6.retr_T,168,plotoptions);
plotoptionsT = plotoptions;

data_fig6.comment = 'see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_trendspaper_figure15_uncertainty_ztest.m iVers = 9 --> quick_compare_era5_retrieval.m';  
%% save /home/sergio/MATLABCODE/oem_pkg_run/MATFILES_for_JGR_trends_paper/fig6.mat data_fig6

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(169); close
figure(169); 
profile_figure6_wiley_2x1tiledlayout(data_fig6.rlat,data_fig6.era5_plays, data_fig6.era5_T,data_fig6.retr_T, data_fig6.era5_WV,data_fig6.retr_WV, 169 ,plotoptionsT,plotoptionsWV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
