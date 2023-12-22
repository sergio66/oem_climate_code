%% see driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m
% z11 = AIRS,       z12 = climcaps      z21 = MERRA2        z22 = ERA5            z11x/z31 = umbc      z32 = giss
% newz11 = ERA5, newz12 = MERRA2     newz21 = AIRSL3     newz22 = CLIMCAPS  newz11x/newz31 = umbc   newz32 = giss

% z31 = z11x;    z31unc = z11xunc;
% z32 = z32;
% z32    = agiss.giss_trend4608(:);     z32 = z32';
% z32unc = agiss.giss_trend_err4608(:); z32unc = z32unc';

% xnewz11 = newz11 + newz11unc;
% xnewz12 = newz12 + newz12unc;
% xnewz21 = newz21 + newz21unc;
% xnewz22 = newz22 + newz22unc;
% xnewz31 = newz11x + newz11xunc;
% xnewz32 = newz32 + newz32unc;

iShowSKTvsColWVunits = +1;

disp(' >->->->->->-')
clear newz*
disp(' DO NOT BELIEVE UNC    DO NOT BELIEVE UNC    D/N')
newz11x = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  
newz12x = (fAIRSL3_day.thestats64x72.stemprate+fAIRSL3_night.thestats64x72.stemprate)*0.5; 
newz13x = (fCLIMCAPSL3_day.thestats64x72.stemprate+fCLIMCAPSL3_night.thestats64x72.stemprate)*0.5;
newz21x = fGISS.giss_trend4608;                                   
newz22x = (fERA5_day.trend_stemp+fERA5_night.trend_stemp)*0.5;                             
newz23x = fMERRA2.trend_stemp;
newz11x = reshape(newz11x,72,64);   newz22x = reshape(newz22x,72,64);   newz23x = reshape(newz23x,72,64);  
newz11 = newz22x(:); newz12 = newz23x(:); newz21 = newz12x(:); newz22 = newz13x(:); newz31 = newz11x(:); newz32 = newz21x(:); 
%% ERA5     MERRA2   AIRSL3   CLIMCAPSL3   UMBC  GISS
clear newz*x; newz11x = newz31; newz11xunc = newz11x*0.1;
newz11unc = newz11*0.1; newz12unc = newz11*0.1; newz21unc = newz21*0.1; newz22unc = newz11*0.1; newz31unc = newz11*0.1; newz32unc = newz11*0.1;
iShowSKTUnc = -1;
show_skt_trends_6models_new
disp(' >->->->->->-')

  %% correlations against UMBC, except sixth (last) is against ERA5
  [r,chisqr,P] = nanlinearcorrelation(newz31,newz11);    thecorr.ST(1) = r;
  [r,chisqr,P] = nanlinearcorrelation(newz31,newz12);    thecorr.ST(2) = r;
  [r,chisqr,P] = nanlinearcorrelation(newz31,newz21);    thecorr.ST(3) = r;
  [r,chisqr,P] = nanlinearcorrelation(newz31,newz22);    thecorr.ST(4) = r;
  [r,chisqr,P] = nanlinearcorrelation(newz31,newz32);   thecorr.ST(5) = r;
  [r,chisqr,P] = nanlinearcorrelation(newz11,newz32);   thecorr.ST_ERA5_GISS = r;
  disp('correlations of UMBC with ERA5/MERRA2/AIRS/CLIMCAPS/GISS')
  thecorr.ST

  modelnames = {'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK'};
  %% do fractional signs
  zall = [newz11 newz12 newz21 newz22 newz31 newz32]';
  wall = ones(size(newz31));
  wall = cos(YY*pi/180)';
  %corr(zall')
  [R,Pvalue,m,s,m0,s0,linearfit,frac_neg0pos,frac_neg0pos_mean_std_stemp] = corrplot_weighted_mean_stddev(zall',wall',modelnames);

disp(' >->->->->->-')

clear newz*
disp(' DO NOT BELIEVE UNC    DO NOT BELIEVE UNC    D')
newz11x = fUMBC_day.results(:,6);
newz12x = fAIRSL3_day.thestats64x72.stemprate;
newz13x = fCLIMCAPSL3_day.thestats64x72.stemprate;
newz21x = fGISS.giss_trend4608;                                   
newz22x = fERA5_day.trend_stemp;
newz23x = fMERRA2.trend_stemp;
newz11x = reshape(newz11x,72,64);   newz22x = reshape(newz22x,72,64);   newz23x = reshape(newz23x,72,64);  
newz11 = newz22x(:); newz12 = newz23x(:); newz21 = newz12x(:); newz22 = newz13x(:); newz31 = newz11x(:); newz32 = newz21x(:); 
clear newz*x; newz11x = newz31; newz11xunc = newz11x*0.1;
newz11unc = newz11*0.1; newz12unc = newz11*0.1; newz21unc = newz21*0.1; newz22unc = newz11*0.1; newz31unc = newz11*0.1; newz32unc = newz11*0.1;
iShowSKTUnc = -1;
show_skt_trends_6models_new
disp(' >->->->->->-')

clear newz*
disp(' DO NOT BELIEVE UNC    DO NOT BELIEVE UNC    N')
newz11x = fUMBC_night.results(:,6);
newz12x = fAIRSL3_night.thestats64x72.stemprate;
newz13x = fCLIMCAPSL3_night.thestats64x72.stemprate;
newz21x = fGISS.giss_trend4608;                                   
newz22x = fERA5_night.trend_stemp;
newz23x = fMERRA2.trend_stemp;
newz11x = reshape(newz11x,72,64);   newz22x = reshape(newz22x,72,64);   newz23x = reshape(newz23x,72,64);  
newz11 = newz22x(:); newz12 = newz23x(:); newz21 = newz12x(:); newz22 = newz13x(:); newz31 = newz11x(:); newz32 = newz21x(:); 
clear newz*x; newz11x = newz31; newz11xunc = newz11x*0.1;
newz11unc = newz11*0.1; newz12unc = newz11*0.1; newz21unc = newz21*0.1; newz22unc = newz11*0.1; newz31unc = newz11*0.1; newz32unc = newz11*0.1;
iShowSKTUnc = -1;
show_skt_trends_6models_new
disp(' >->->->->->-')
