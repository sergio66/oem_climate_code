omi = load('omi_tcwv_trends.mat');
[omiY,omiX] = meshgrid(omi.omi_colwv.lat,omi.omi_colwv.lon);
[m,s,m0,s0] = weighted_mean_stddev(omi.omi_colwv.trend,cos(omiY*pi/180));                         fprintf(1,'cosine weighted OMI      mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);

boo = find(abs(omiY) <= 30);
[m,s,m0,s0] = weighted_mean_stddev(omi.omi_colwv.trend(boo),ones(size(boo)));                     fprintf(1,'tropical        OMI mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);

iShowSKTvsColWVunits = -1;

disp(' COL WV TRENDS!!! mm/yr!!!! >->->->->->-')
clear newz*
disp('D/N')
newz11x = (fUMBC_day.mmwtrend+fUMBC_night.mmwtrend)*0.5;  
newz12x = (fAIRSL3_day.mmwtrend+fAIRSL3_night.mmwtrend)*0.5; 
newz13x = (fCLIMCAPSL3_day.mmwtrend+fCLIMCAPSL3_night.mmwtrend)*0.5;
newz21x = 0*fGISS.giss_trend4608;                                   
newz22x = (fERA5_day.mmwtrend+fERA5_night.mmwtrend)*0.5;                             
newz23x = fMERRA2.mmwtrend;
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
disp('D')
newz11x = fUMBC_day.mmwtrend;
newz12x = fAIRSL3_day.mmwtrend;
newz13x = fCLIMCAPSL3_day.mmwtrend;
newz21x = 0*fGISS.giss_trend4608;                                   
newz22x = fERA5_day.mmwtrend;
newz23x = fMERRA2.mmwtrend;
newz11x = reshape(newz11x,72,64);   newz22x = reshape(newz22x,72,64);   newz23x = reshape(newz23x,72,64);  
newz11 = newz22x(:); newz12 = newz23x(:); newz21 = newz12x(:); newz22 = newz13x(:); newz31 = newz11x(:); newz32 = newz21x(:); 
clear newz*x; newz11x = newz31; newz11xunc = newz11x*0.1;
newz11unc = newz11*0.1; newz12unc = newz11*0.1; newz21unc = newz21*0.1; newz22unc = newz11*0.1; newz31unc = newz11*0.1; newz32unc = newz11*0.1;
iShowSKTUnc = -1;
show_skt_trends_6models_new
disp(' >->->->->->-')

clear newz*
disp('N')
newz11x = fUMBC_night.mmwtrend;
newz12x = fAIRSL3_night.mmwtrend;
newz13x = fCLIMCAPSL3_night.mmwtrend;
newz21x = 0*fGISS.giss_trend4608;                                   
newz22x = fERA5_night.mmwtrend;
newz23x = fMERRA2.mmwtrend;
newz11x = reshape(newz11x,72,64);   newz22x = reshape(newz22x,72,64);   newz23x = reshape(newz23x,72,64);  
newz11 = newz22x(:); newz12 = newz23x(:); newz21 = newz12x(:); newz22 = newz13x(:); newz31 = newz11x(:); newz32 = newz21x(:); 
clear newz*x; newz11x = newz31; newz11xunc = newz11x*0.1;
newz11unc = newz11*0.1; newz12unc = newz11*0.1; newz21unc = newz21*0.1; newz22unc = newz11*0.1; newz31unc = newz11*0.1; newz32unc = newz11*0.1;
iShowSKTUnc = -1;
show_skt_trends_6models_new
disp(' >->->->->->-')
