iDoDiff = +1;
if topts.ocb_set == 1
  iDoDiff = input('this is synthetic cal retrieval, so enter (+1/Default) to plot UMBC-ERA5*AK (-1) to plot ERA5 x AK : ');
  if length(iDoDiff) == 0
    iDoDiff = 1;
  end
end

clear plotoptions;
plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dT/dt'; plotoptions.cmap = llsmap5;
plotoptions.str1 = 'UMBC';   
if topts.ocb_set == 0
  plotoptions.str2 = 'ERA5 x AK'; 
elseif topts.ocb_set == 1 & iDoDiff > 0
  plotoptions.str2 = 'UMBC-ERA5*AK'; 
elseif topts.ocb_set == 1 & iDoDiff < 0  
  plotoptions.str2 = 'ERA5 x AK'; 
end
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
if topts.ocb_set == 0 
  z2 = temprate_akF_era5'; z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
else
  if topts.ocb_set == 1 & iDoDiff == 1
    z2 = resultsT' - temprate_ak0_era5'; z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
    z2 = resultsT' - temprate_akF_era5'; z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
    disp('WARNING : sice this is retrieval using SYNTHETIC rates, fig 49/50 : center panel is not AK * ERA5rate, but is UMBC - AK * ERA5rate')
  else
    z2 = temprate_akF_era5'; z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
  end
end
z3 = temprate_ak0_era5'; z3 = reshape(z3,length(pjunk20),72,64); z3 = squeeze(nanmean(z3,2));
iFig = 49; figure(iFig); clf; profile_plots_1x3tiledlayout(rlat,pjunk20,z1,z2,z3,iFig,plotoptions);

plotoptions.cx = [-1 +1]*0.015; plotoptions.maintitle = 'dWVfrac/dt'; plotoptions.cmap = llsmap5;
plotoptions.yLimits = [100 1000];
plotoptions.yLinearOrLog = +1;
z1 = resultsWV';          z1 = reshape(z1,length(pjunk20),72,64); z1 = squeeze(nanmean(z1,2));
if topts.ocb_set == 0 
  z2 = waterrate_akF_era5'; z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
else
  if topts.ocb_set == 1 & iDoDiff == 1
    z2 = 1-resultsWV'./(waterrate_ak0_era5'+eps); z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
    z2 = resultsWV' - waterrate_ak0_era5'; z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
    z2 = resultsWV' - waterrate_akF_era5'; z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
    disp('WARNING : sice this is retrieval using SYNTHETIC rates, fig 49/50 : center panel is not AK * ERA5rate, but is UMBC - AK * ERA5rate')
  else
    z2 = waterrate_akF_era5'; z2 = reshape(z2,length(pjunk20),72,64); z2 = squeeze(nanmean(z2,2));
  end
end
z3 = waterrate_ak0_era5'; z3 = reshape(z3,length(pjunk20),72,64); z3 = squeeze(nanmean(z3,2));
iFig = 50; figure(iFig); clf; profile_plots_1x3tiledlayout(rlat,pjunk20,z1,z2,z3,iFig,plotoptions);

plotoptions.cx = [-1 +1]*0.5; plotoptions.maintitle = 'dRH/dt'; plotoptions.cmap = llsmap5;
plotoptions.yLimits = [100 1000];
plotoptions.yLinearOrLog = +1;
plotoptions.str2 = 'ERA5';   
if isfield(plotoptions,'str3')
  plotoptions = rmfield(plotoptions,'str3');
end
z1 = deltaRHlat'; 
z2 = era5.trend_RH; z2 = reshape(z2,100,72,64); z2 = squeeze(nanmean(z2,2));
iFig = 51; figure(iFig); clf; profile_plots_1x2tiledlayout(rlat,plays,z1,z2,iFig,plotoptions);
