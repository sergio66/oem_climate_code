%% from compare_SKT_trends_Day_vs_Night.m

clear newz*
iFig = 22;
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
%{
  newz11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  newz12 = (fAIRSL3_day.thestats64x72.stemprate+fAIRSL3_night.thestats64x72.stemprate)*0.5; 
    newz13 = (fCLIMCAPSL3_day.thestats64x72.stemprate+fCLIMCAPSL3_night.thestats64x72.stemprate)*0.5;
  newz21 = fGISS.giss_trend4608;                                   newz22 = (fERA5_day.trend_stemp+fERA5_night.trend_stemp)*0.5;                             newz23 = fMERRA2.trend_stemp;
    newz11 = reshape(newz11,72,64);   newz22 = reshape(newz22,72,64);   newz23 = reshape(newz23,72,64);  
%}

  save skt_trends_strow.mat newz11 newz12 newz13 newz21 newz22 newz23 iFig plotoptions

  figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);
