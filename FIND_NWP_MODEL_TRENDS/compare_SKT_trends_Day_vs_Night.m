ii = 0; 
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fUMBC_day.results(:,6),72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT DAY');
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fUMBC_night.results(:,6),72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT NIGHT');

ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fERA5_day.trend_stemp,72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : ERA5 DAY');
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fERA5_night.trend_stemp,72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : ERA5 NIGHT');

ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(fAIRSL3_day.thestats64x72.stemprate',1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS L3 DAY');
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(fAIRSL3_night.thestats64x72.stemprate',1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS L3 NIGHT');

ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(fCLIMCAPSL3_day.thestats64x72.stemprate',1),[-90 +90],[-180 +180]); title('dSKT/dt : CLIMCAPS L3 DAY');
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(fCLIMCAPSL3_night.thestats64x72.stemprate',1),[-90 +90],[-180 +180]); title('dSKT/dt : CLIMCAPS L3 NIGHT');

ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fMERRA2.trend_stemp,72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : MERRA2 DAY/NIGHT');
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(fGISS.giss_trend4608',1),[-90 +90],[-180 +180]); title('dSKT/dt : GISS2 DAY/NIGHT');

ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5,72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT DAY/NIGHT');
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn((fAIRSL3_day.thestats64x72.stemprate'+fAIRSL3_night.thestats64x72.stemprate')*0.5,1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS L3 DAY/NIGHT');
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn((fCLIMCAPSL3_day.thestats64x72.stemprate'+fCLIMCAPSL3_night.thestats64x72.stemprate')*0.5,1),[-90 +90],[-180 +180]); title('dSKT/dt : CLIMCAPS L3 DAY/NIGHT');
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fERA5_day.trend_stemp+fERA5_night.trend_stemp)*0.5,72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : ERA5 DAY/NIGHT');

ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fUMBC_day.results(:,6)-fUMBC_night.results(:,6))*1.0,72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS\_RT DAY - NIGHT');
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn((fAIRSL3_day.thestats64x72.stemprate'-fAIRSL3_night.thestats64x72.stemprate')*1.0,1),[-90 +90],[-180 +180]); title('dSKT/dt : AIRS L3 DAY - NIGHT');
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn((fCLIMCAPSL3_day.thestats64x72.stemprate'-fCLIMCAPSL3_night.thestats64x72.stemprate')*1.0,1),[-90 +90],[-180 +180]); title('dSKT/dt : CLIMCAPS L3 DAY - NIGHT');
ii = ii + 1; figure(ii); sizefig; ; clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fERA5_day.trend_stemp-fERA5_night.trend_stemp)*1.0,72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : ERA5 DAY - NIGHT');

iiMax = ii;
for ii = 1 : iiMax; figure(ii); caxis([-1 +1]*0.15); colormap(llsmap5); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newz*
iFig = 20; 
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = '(D) AIRS\_RT';     plotoptions.str21 = '(N)';
  plotoptions.str12 = '(D) AIRS L3';      plotoptions.str22 = '(N)';
  plotoptions.str13 = '(D) CLIMCAPS L3';  plotoptions.str23 = '(N)';
  plotoptions.str14 = '(D) ERA5';         plotoptions.str24 = '(N)';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
  newz11 = fUMBC_day.results(:,6);                  newz21 = fUMBC_night.results(:,6);                           
  newz12 = fAIRSL3_day.thestats64x72.stemprate;     newz22 = fAIRSL3_night.thestats64x72.stemprate; 
  newz13 = fCLIMCAPSL3_day.thestats64x72.stemprate; newz23 = fCLIMCAPSL3_night.thestats64x72.stemprate;
  newz14 = fERA5_day.trend_stemp;                   newz24 = fERA5_night.trend_stemp; 
    newz11 = reshape(newz11,72,64);   newz21 = reshape(newz21,72,64);
    newz14 = reshape(newz14,72,64);   newz24 = reshape(newz24,72,64);
  figure(iFig); sizefig; ; clf; aslmap_2x4tiledlayout(newz11,newz12,newz13,newz14,newz21,newz22,newz23,newz24,iFig,plotoptions);

clear newz*
iFig = 20; 
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = '(D) AIRS\_RT';     plotoptions.str21 = '(N)';   plotoptions.str31 = '(D-N)';
  plotoptions.str12 = '(D) AIRS L3';      plotoptions.str22 = '(N)';   plotoptions.str32 = '(D-N)';
  plotoptions.str13 = '(D) CLIMCAPS L3';  plotoptions.str23 = '(N)';   plotoptions.str33 = '(D-N)';
  plotoptions.str14 = '(D) ERA5';         plotoptions.str24 = '(N)';   plotoptions.str34 = '(D-N)';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
  newz11 = fUMBC_day.results(:,6);                  newz21 = fUMBC_night.results(:,6);                       newz31 = newz11 - newz21;
  newz12 = fAIRSL3_day.thestats64x72.stemprate;     newz22 = fAIRSL3_night.thestats64x72.stemprate;          newz32 = newz12 - newz22;
  newz13 = fCLIMCAPSL3_day.thestats64x72.stemprate; newz23 = fCLIMCAPSL3_night.thestats64x72.stemprate;      newz33 = newz13 - newz23;
  newz14 = fERA5_day.trend_stemp;                   newz24 = fERA5_night.trend_stemp;                        newz34 = newz14 - newz24;
    newz11 = reshape(newz11,72,64);   newz21 = reshape(newz21,72,64);   newz31 = reshape(newz31,72,64);   
    newz14 = reshape(newz14,72,64);   newz24 = reshape(newz24,72,64);   newz34 = reshape(newz34,72,64);   
  figure(iFig); sizefig; ; clf; aslmap_3x4tiledlayout(newz11,newz12,newz13,newz14,newz21,newz22,newz23,newz24,newz31,newz32,newz33,newz34,iFig,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newz*
iFig = 21; 
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'ERA5';        
  plotoptions.str21 = 'AIRS L3';      plotoptions.str22 = 'CLIMCAPS L3'; 
  newz11 = (fUMBC_day.results(:,6)-fUMBC_night.results(:,6))*1.0;                           newz12 = (fERA5_day.trend_stemp-fERA5_night.trend_stemp)*1.0; 
  newz21 = (fAIRSL3_day.thestats64x72.stemprate-fAIRSL3_night.thestats64x72.stemprate)*1.0; newz22 = (fCLIMCAPSL3_day.thestats64x72.stemprate-fCLIMCAPSL3_night.thestats64x72.stemprate)*1.0;
    newz11 = reshape(newz11,72,64);   newz12 = reshape(newz12,72,64);
  figure(iFig); sizefig; ; clf; aslmap_2x2tiledlayout(newz11,newz12,newz21,newz22,iFig,plotoptions);

clear newz*
iFig = 22;
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = 'dSKT/dt [K/yr]';
  newz11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  newz12 = (fAIRSL3_day.thestats64x72.stemprate+fAIRSL3_night.thestats64x72.stemprate)*0.5; 
    newz13 = (fCLIMCAPSL3_day.thestats64x72.stemprate+fCLIMCAPSL3_night.thestats64x72.stemprate)*0.5;
  newz21 = fGISS.giss_trend4608;                                   newz22 = (fERA5_day.trend_stemp+fERA5_night.trend_stemp)*0.5;                             newz23 = fMERRA2.trend_stemp;
    newz11 = reshape(newz11,72,64);   newz22 = reshape(newz22,72,64);   newz23 = reshape(newz23,72,64);  
  figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);

junk = [newz11(:) newz12(:) newz13(:) newz21(:) newz22(:) newz23(:)];
moo0 = nanmean(junk,2);
moo1 = abs(nanmean(junk,2));
moo2 = nanstd(junk,[],2);
moo3 = zeros(size(moo1)); moo3(moo1 > moo2) = 1;
moo4 = moo1./moo2; moo4(moo4 < 1) = 1; moo4 = log10(moo4);

iFig = 23; figure(iFig); sizefig; ; clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(nanmean(junk,2),72,64)',1),[-90 +90],[-180 +180]);   cx = [-1 +1]*0.151; caxis(cx); colormap(llsmap5);
iFig = 24; figure(iFig); sizefig; ; clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(nanstd(junk,[],2),72,64)',1),[-90 +90],[-180 +180]); cx = [0 +1]*0.031;  caxis(cx); colormap(jet);

iFig = 25; figure(iFig); sizefig; ; clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(moo3,72,64)',1),[-90 +90],[-180 +180]); cx = [0 +1];  caxis(cx); colormap(jet);
iFig = 25; figure(iFig); sizefig; ; clf; aslmap(iFig,rlat65,rlon73,reshape(moo3,72,64)',[-90 +90],[-180 +180]);            cx = [0 +1]; caxis(cx); colormap(usa2x);

iFig = 26; figure(iFig); sizefig; ; clf; aslmap(iFig,rlat65,rlon73,reshape(moo4,72,64)',[-90 +90],[-180 +180]);        cx = [0 +2]; caxis(cx); colormap(usa2x);
iFig = 26; figure(iFig); sizefig; ; clf; aslmap(iFig,rlat65,rlon73,reshape(10.^(moo4),72,64)',[-90 +90],[-180 +180]);  cx = (10.^[0 +1]); caxis(cx); colormap(jett);

mooUMBC_num   = 0.5*(fUMBC_day.results(:,6) + fUMBC_night.results(:,6));
mooUMBC_denom = sqrt(fUMBC_day.resultsunc(:,6).^2 + fUMBC_night.resultsunc(:,6).^2);
mooUMBC = abs(mooUMBC_num)./mooUMBC_denom;
iFig = 27; figure(iFig); sizefig; ; clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(mooUMBC_num,72,64)',1),[-90 +90],[-180 +180]);   cx = ([-1 +1]*0.15); caxis(cx); colormap(llsmap5);
iFig = 28; figure(iFig); sizefig; ; clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(mooUMBC_denom,72,64)',1),[-90 +90],[-180 +180]); cx = ([0 +1]*0.15/3); caxis(cx); colormap(jet);
iFig = 29; figure(iFig); sizefig; ; clf; aslmap(iFig,rlat65,rlon73,reshape(10.^(mooUMBC),72,64)',[-90 +90],[-180 +180]);            cx = (10.^[0 +1]); caxis(cx); colormap(jett);
iFig = 29; figure(iFig); sizefig; ; clf; aslmap(iFig,rlat65,rlon73,reshape(mooUMBC,72,64)',[-90 +90],[-180 +180]);                  cx = ([0 +5]); caxis(cx); colormap(jett);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(90); clf
  ta = tiledlayout(2,1,'TileSpacing','None', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tfov(1) = nexttile;
    aslmapSergio(rlat65,rlon73,smoothn(reshape(moo0,72,64)',1), [-90 +90],[-180 +180]); colormap(tfov(1),llsmap5); caxis([-1 +1]*0.151); colorbar
    text(-2.5,1.15,'(a)')
  tfov(2) = nexttile;
    aslmapSergio(rlat65,rlon73,smoothn(reshape(moo2,72,64)',1), [-90 +90],[-180 +180]); colormap(tfov(2),jet); caxis([ 0 +1]/5*0.15); colorbar
    text(-2.5,1.15,'(b)')
colormap(tfov(1),llsmap5);
colormap(tfov(2),jet);
ta.Padding = 'compact';
ta.TileSpacing = 'compact';

