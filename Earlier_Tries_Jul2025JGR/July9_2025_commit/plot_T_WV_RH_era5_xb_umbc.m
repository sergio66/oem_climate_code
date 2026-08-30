%% plot UMBC (retrieval) vs xb (initialization, typically 0) vs ERA5 trends
clear plotoptions;
plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dT/dt'; plotoptions.cmap = llsmap5;
plotoptions.str1 = 'UMBC';   
plotoptions.str2 = 'XB=INIT'; 
plotoptions.str3 = 'ERA5';   
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = -1;
plotoptions.yReverseDir = +1;
plotoptions.yLimits = [10 1000];
plotoptions.yLimits = [1 1000];
plotoptions.surface = spres_avg;

plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dT/dt'; plotoptions.cmap = llsmap5;
plotoptions.yLimits = [10 1000];
plotoptions.yLinearOrLog = -1;
z1 = resultsT';          z1 = reshape(z1,length(pjunk20),72,64); z1 = squeeze(nanmean(z1,2));
z2 = xbT';               z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
z3 = temprate_ak0_era5'; z3 = reshape(z3,length(pjunk20),72,64); z3 = squeeze(nanmean(z3,2));
iFig = 52; figure(iFig); clf; profile_plots_1x3tiledlayout(rlat,pjunk20,z1,z2,z3,iFig,plotoptions);

plotoptions.cx = [-1 +1]*0.015; plotoptions.maintitle = 'dWVfrac/dt'; plotoptions.cmap = llsmap5;
plotoptions.yLimits = [100 1000];
plotoptions.yLinearOrLog = +1;
z1 = resultsWV';          z1 = reshape(z1,length(pjunk20),72,64); z1 = squeeze(nanmean(z1,2));
z2 = xbWV';               z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
z3 = waterrate_ak0_era5'; z3 = reshape(z3,length(pjunk20),72,64); z3 = squeeze(nanmean(z3,2));
iFig = 53; figure(iFig); clf; profile_plots_1x3tiledlayout(rlat,pjunk20,z1,z2,z3,iFig,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot ratio : UMBC (retrieval)/ERA5 trends vs xb (initialization, typically 0)/ERA5 trends vs ERA5 trends/ERA5 trends
clear plotoptions;
plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dT/dt'; plotoptions.cmap = llsmap5;
plotoptions.str1 = 'UMBC vs ERA5';   
plotoptions.str2 = 'XB(INIT) vs ERA5'; 
plotoptions.str3 = 'ERA5 vs ERA5 = 1';   
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = -1;
plotoptions.yReverseDir = +1;
plotoptions.yLimits = [10 1000];
plotoptions.yLimits = [1 1000];
plotoptions.surface = spres_avg;

plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dT(ERA5)/dt - dTX/dt'; plotoptions.cmap = llsmap5;
plotoptions.yLimits = [10 1000];
plotoptions.yLinearOrLog = -1;
z1 = resultsT';          z1 = reshape(z1,length(pjunk20),72,64); z1 = squeeze(nanmean(z1,2));
z2 = xbT';               z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
z3 = temprate_ak0_era5'; z3 = reshape(z3,length(pjunk20),72,64); z3 = squeeze(nanmean(z3,2));
z1 = z3 - z1;
z2 = z3 - z2;
z3 = z3 - z3;
iFig = 54; figure(iFig); clf; profile_plots_1x3tiledlayout(rlat,pjunk20,z1,z2,z3,iFig,plotoptions);

plotoptions.cx = 1 + ([-1 +1]*2); plotoptions.maintitle = '1 - (dWVXfrac/dt / dWVXfrac(REA5)/dt)'; plotoptions.cmap = llsmap5;
plotoptions.cx = [-1 +1];         plotoptions.maintitle = '1 - (dWVXfrac/dt / dWVXfrac(REA5)/dt)'; plotoptions.cmap = llsmap5;
plotoptions.cx = [-1 +1];         plotoptions.maintitle = '(dWVXfrac/dt / dWVXfrac(REA5)/dt) - 1'; plotoptions.cmap = llsmap5;
plotoptions.yLimits = [100 1000];
plotoptions.yLinearOrLog = +1;
z1 = resultsWV';          z1 = reshape(z1,length(pjunk20),72,64); z1 = squeeze(nanmean(z1,2));
z2 = xbWV';               z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
z3 = waterrate_ak0_era5'; z3 = reshape(z3,length(pjunk20),72,64); z3 = squeeze(nanmean(z3,2));
z1 = 1 - z1./z3;
z2 = 1 - z2./z3;
z3 = 1 - z3./z3;
iFig = 55; figure(iFig); clf; profile_plots_1x3tiledlayout(rlat,pjunk20,-z1,-z2,-z3,iFig,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

