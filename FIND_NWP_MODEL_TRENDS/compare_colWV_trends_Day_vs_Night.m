ii = 0; 

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fUMBC_day.mmwtrend,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : AIRS\_RT DAY');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fUMBC_night.mmwtrend,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : AIRS\_RT NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fERA5_day.mmwtrend,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : ERA5 DAY');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fERA5_night.mmwtrend,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : ERA5 NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fAIRSL3_day.mmwtrend,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : AIRS L3 DAY');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fAIRSL3_night.mmwtrend,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : AIRS L3 NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fCLIMCAPSL3_day.mmwtrend,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : CLIMCAPS L3 DAY');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fCLIMCAPSL3_night.mmwtrend,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : CLIMCAPS L3 NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fMERRA2.mmwtrend,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : MERRA2 DAY/NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fUMBC_day.mmwtrend+fUMBC_night.mmwtrend)*0.5,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : AIRS\_RT DAY/NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fAIRSL3_day.mmwtrend+fAIRSL3_night.mmwtrend)*0.5,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : AIRS L3 DAY/NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fCLIMCAPSL3_day.mmwtrend+fCLIMCAPSL3_night.mmwtrend)*0.5,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : CLIMCAPS L3 DAY/NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fERA5_day.mmwtrend+fERA5_night.mmwtrend)*0.5,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : ERA5 DAY/NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fUMBC_day.mmwtrend-fUMBC_night.mmwtrend)*1.0,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : AIRS\_RT DAY - NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fAIRSL3_day.mmwtrend-fAIRSL3_night.mmwtrend)*1.0,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : AIRS L3 DAY - NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fCLIMCAPSL3_day.mmwtrend-fCLIMCAPSL3_night.mmwtrend)*1.0,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : CLIMCAPS L3 DAY - NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fERA5_day.mmwtrend-fERA5_night.mmwtrend)*1.0,72,64)',1),[-90 +90],[-180 +180]); title('dcolWV/dt : ERA5 DAY - NIGHT');

iiMax = ii;
for ii = 1 : iiMax; figure(ii); caxis([-1 +1]*0.15); colormap(llsmap5); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newz*
iFig = 30; 
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dcolWV/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = '(D) AIRS\_RT';     plotoptions.str21 = '(N)';
  plotoptions.str12 = '(D) AIRS L3';      plotoptions.str22 = '(N)';
  plotoptions.str13 = '(D) CLIMCAPS L3';  plotoptions.str23 = '(N)';
  plotoptions.str14 = '(D) ERA5';         plotoptions.str24 = '(N)';
  plotoptions.barstr = 'dcolWV/dt [mm/yr]';
  newz11 = fUMBC_day.mmwtrend;         newz21 = fUMBC_night.mmwtrend;                           
  newz12 = fAIRSL3_day.mmwtrend;       newz22 = fAIRSL3_night.mmwtrend; 
  newz13 = fCLIMCAPSL3_day.mmwtrend;   newz23 = fCLIMCAPSL3_night.mmwtrend;
  newz14 = fERA5_day.mmwtrend;         newz24 = fERA5_night.mmwtrend; 
    newz11 = reshape(newz11,72,64);   newz21 = reshape(newz21,72,64);
    newz12 = reshape(newz12,72,64);   newz22 = reshape(newz22,72,64);
    newz13 = reshape(newz13,72,64);   newz23 = reshape(newz23,72,64);
    newz14 = reshape(newz14,72,64);   newz24 = reshape(newz24,72,64);
  figure(iFig); clf; aslmap_2x4tiledlayout(newz11,newz12,newz13,newz14,newz21,newz22,newz23,newz24,iFig,plotoptions);

clear newz*
iFig = 30; 
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = '(D) AIRS\_RT';     plotoptions.str21 = '(N)';   plotoptions.str31 = '(D-N)';
  plotoptions.str12 = '(D) AIRS L3';      plotoptions.str22 = '(N)';   plotoptions.str32 = '(D-N)';
  plotoptions.str13 = '(D) CLIMCAPS L3';  plotoptions.str23 = '(N)';   plotoptions.str33 = '(D-N)';
  plotoptions.str14 = '(D) ERA5';         plotoptions.str24 = '(N)';   plotoptions.str34 = '(D-N)';
  plotoptions.barstr = 'dcolWV/dt [mm/yr]';
  newz11 = fUMBC_day.mmwtrend;                 newz21 = fUMBC_night.mmwtrend;            newz31 = newz11 - newz21;
  newz12 = fAIRSL3_day.mmwtrend;               newz22 = fAIRSL3_night.mmwtrend;          newz32 = newz12 - newz22;
  newz13 = fCLIMCAPSL3_day.mmwtrend;           newz23 = fCLIMCAPSL3_night.mmwtrend;      newz33 = newz13 - newz23;
  newz14 = fERA5_day.mmwtrend;                 newz24 = fERA5_night.mmwtrend;            newz34 = newz14 - newz24;
    newz11 = reshape(newz11,72,64);   newz21 = reshape(newz21,72,64);   newz31 = reshape(newz31,72,64);   
    newz14 = reshape(newz14,72,64);   newz24 = reshape(newz24,72,64);   newz34 = reshape(newz34,72,64);   
  figure(iFig); clf; aslmap_3x4tiledlayout(newz11,newz12,newz13,newz14,newz21,newz22,newz23,newz24,newz31,newz32,newz33,newz34,iFig,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newz*
iFig = 31; 
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dcolWV/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'ERA5';        
  plotoptions.str21 = 'AIRS L3';      plotoptions.str22 = 'CLIMCAPS L3'; 
  plotoptions.barstr = 'dcolWV/dt [mm/yr]';
  newz11 = (fUMBC_day.mmwtrend-fUMBC_night.mmwtrend)*1.0;                           newz12 = (fERA5_day.mmwtrend-fERA5_night.mmwtrend)*1.0; 
  newz21 = (fAIRSL3_day.mmwtrend-fAIRSL3_night.mmwtrend)*1.0; newz22 = (fCLIMCAPSL3_day.mmwtrend-fCLIMCAPSL3_night.mmwtrend)*1.0;
    newz11 = reshape(newz11,72,64);   newz12 = reshape(newz12,72,64);
  figure(iFig); clf; aslmap_2x2tiledlayout(newz11,newz12,newz21,newz22,iFig,plotoptions);

clear newz*
iFig = 32;
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dcolWV/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'N/A';          plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = 'dcolWV/dt [mm/yr]';
  newz11 = (fUMBC_day.mmwtrend+fUMBC_night.mmwtrend)*0.5;    newz12 = (fAIRSL3_day.mmwtrend+fAIRSL3_night.mmwtrend)*0.5;   newz13 = (fCLIMCAPSL3_day.mmwtrend+fCLIMCAPSL3_night.mmwtrend)*0.5;
  newz21 = fGISS.giss_trend4608*NaN;                         newz22 = (fERA5_day.mmwtrend+fERA5_night.mmwtrend)*0.5;       newz23 = fMERRA2.mmwtrend;
    newz11 = reshape(newz11,72,64);   newz22 = reshape(newz22,72,64);   newz23 = reshape(newz23,72,64);  
  figure(iFig); clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);
  junk = [newz11(:) newz12(:) newz13(:) newz21(:) newz22(:) newz23(:)];

clear newz*
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dcolWV/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.strzz = 'AIRS\_RT';  
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  plotoptions.barstr = 'dcolWV/dt [mm/yr]';
  newzAA = (fUMBC_day.mmwtrend+fUMBC_night.mmwtrend)*0.5;    newz11 = (fAIRSL3_day.mmwtrend+fAIRSL3_night.mmwtrend)*0.5; newz12 = (fCLIMCAPSL3_day.mmwtrend+fCLIMCAPSL3_night.mmwtrend)*0.5;
  newz22 = (fERA5_day.mmwtrend+fERA5_night.mmwtrend)*0.5;    newz21 = fMERRA2.mmwtrend;
    newz11 = reshape(newz11,72,64);   newz12 = reshape(newz12,72,64);  newz21 = reshape(newz21,72,64);   newz22 = reshape(newz22,72,64);  newzAA = reshape(newzAA,72,64);    
  aslmap_2x1x2tiledlayout(newz11,newz12,newzAA,newz21,newz22,iFig,plotoptions);
  junk = [newz11(:) newz12(:) newzAA(:) newz21(:) newz22(:)];

clear newz*
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dcolWV/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';   plotoptions.str12 = 'AIRS L3';   plotoptions.str13 = 'CLIMCAPS L3';  
  plotoptions.str21 = 'ERA5';      plotoptions.str22 = 'MERRA2';
  plotoptions.barstr = 'dcolWV/dt [mm/yr]';
  newz11 = (fUMBC_day.mmwtrend+fUMBC_night.mmwtrend)*0.5;    newz12 = (fAIRSL3_day.mmwtrend+fAIRSL3_night.mmwtrend)*0.5; newz13 = (fCLIMCAPSL3_day.mmwtrend+fCLIMCAPSL3_night.mmwtrend)*0.5;
  newz21 = (fERA5_day.mmwtrend+fERA5_night.mmwtrend)*0.5;    newz22 = fMERRA2.mmwtrend;
    newz11 = reshape(newz11,72,64);   newz12 = reshape(newz12,72,64);  newz21 = reshape(newz21,72,64);   newz22 = reshape(newz22,72,64);  newz13 = reshape(newz13,72,64);    
  aslmap_2x1x2tiledlayout_tall(newz11,newz12,newz13,newz21,newz22,iFig,plotoptions);  
  junk = [newz11(:) newz12(:) newz13(:) newz21(:) newz22(:)];

moo0 = nanmean(junk,2);
moo1 = abs(nanmean(junk,2));
moo2 = nanstd(junk,[],2);
moo3 = zeros(size(moo1)); moo3(moo1 > moo2) = 1;
moo4 = moo1./moo2; moo4(moo4 < 1) = 1; moo4 = log10(moo4);

iFig = 33; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(nanmean(junk,2),72,64)',1),[-90 +90],[-180 +180]);   cx = [-1 +1]*0.151; caxis(cx); colormap(llsmap5);
iFig = 34; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(nanstd(junk,[],2),72,64)',1),[-90 +90],[-180 +180]); cx = [0 +1]*0.071;  caxis(cx); colormap(jet);

iFig = 35; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(moo3,72,64)',1),[-90 +90],[-180 +180]); cx = [0 +1];  caxis(cx); colormap(jet);
iFig = 35; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,reshape(moo3,72,64)',[-90 +90],[-180 +180]); cx = [0 +1]; caxis(cx); colormap(usa2x);

iFig = 36; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,reshape(moo4,72,64)',[-90 +90],[-180 +180]);        cx = [0 +2]; caxis(cx); colormap(usa2x);
iFig = 36; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,reshape(10.^(moo4),72,64)',[-90 +90],[-180 +180]);  cx = (10.^[0 +1]); caxis(cx); colormap(jett);

mooUMBC_num   = 0.5*(fUMBC_day.mmwtrend + fUMBC_night.mmwtrend);
mooUMBC_denom = sqrt(fUMBC_day.mmwtrendunc.^2 + fUMBC_night.mmwtrendunc.^2);
mooUMBC = abs(mooUMBC_num)./mooUMBC_denom;
iFig = 37; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(mooUMBC_num,72,64)',1),[-90 +90],[-180 +180]);   cx = ([-1 +1]*0.15); caxis(cx); colormap(llsmap5);
iFig = 38; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(mooUMBC_denom,72,64)',1),[-90 +90],[-180 +180]); cx = ([0 +1]*0.15/2); caxis(cx); colormap(jet);
iFig = 39; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,reshape(10.^(mooUMBC),72,64)',[-90 +90],[-180 +180]);            cx = (10.^[0 +1]); caxis(cx); colormap(jett);
iFig = 39; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,reshape(mooUMBC,72,64)',[-90 +90],[-180 +180]);                  cx = ([0 +5]); caxis(cx); colormap(jett);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(91); clf
  ta = tiledlayout(2,1,'TileSpacing','None', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tfov(1) = nexttile;
    aslmapSergio(rlat65,rlon73,smoothn(reshape(moo0,72,64)',1), [-90 +90],[-180 +180]); colormap(tfov(1),llsmap5); caxis([-1 +1]*0.151); colorbar
    text(-2.5,1.15,'(a)')
  tfov(2) = nexttile;
    aslmapSergio(rlat65,rlon73,smoothn(reshape(moo2,72,64)',1), [-90 +90],[-180 +180]); colormap(tfov(2),jet); caxis([ 0 +1]*0.071); colorbar
    text(-2.5,1.15,'(b)')
colormap(tfov(1),llsmap5);
colormap(tfov(2),jet);
ta.Padding = 'compact';
ta.TileSpacing = 'compact';

figure(92); clf
plot(rlat,nanmean(newz11,1),'k',rlat,nanmean(newz12,1),'b',rlat,nanmean(newz13,1),'gx-',rlat,nanmean(newz21,1),'r',rlat,nanmean(newz22,1),'m','linewidth',4);
plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('d mmw/dt [mm/yr]');

fD = 2.0; fN = 0.0;
newz11 = (fD*fUMBC_day.mmwtrend+fN*fUMBC_night.mmwtrend)*0.5;    newz12 = (fD*fAIRSL3_day.mmwtrend+fN*fAIRSL3_night.mmwtrend)*0.5; newz13 = (fD*fCLIMCAPSL3_day.mmwtrend+fN*fCLIMCAPSL3_night.mmwtrend)*0.5;
newz21 = (fD*fERA5_day.mmwtrend+fN*fERA5_night.mmwtrend)*0.5;    newz22 = fMERRA2.mmwtrend;
    newz11 = reshape(newz11,72,64);   newz12 = reshape(newz12,72,64);  newz21 = reshape(newz21,72,64);   newz22 = reshape(newz22,72,64);  newz13 = reshape(newz13,72,64);    
figure(97); clf
plot(rlat,nanmean(newz11,1),'k',rlat,nanmean(newz12,1),'b',rlat,nanmean(newz13,1),'gx-',rlat,nanmean(newz21,1),'r',rlat,nanmean(newz22,1),'m','linewidth',4);
plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('DAY d mmw/dt [mm/yr]');

fD = 0.0; fN = 2.0;
newz11 = (fD*fUMBC_day.mmwtrend+fN*fUMBC_night.mmwtrend)*0.5;    newz12 = (fD*fAIRSL3_day.mmwtrend+fN*fAIRSL3_night.mmwtrend)*0.5; newz13 = (fD*fCLIMCAPSL3_day.mmwtrend+fN*fCLIMCAPSL3_night.mmwtrend)*0.5;
newz21 = (fD*fERA5_day.mmwtrend+fN*fERA5_night.mmwtrend)*0.5;    newz22 = fMERRA2.mmwtrend;
    newz11 = reshape(newz11,72,64);   newz12 = reshape(newz12,72,64);  newz21 = reshape(newz21,72,64);   newz22 = reshape(newz22,72,64);  newz13 = reshape(newz13,72,64);    
figure(98); clf
plot(rlat,nanmean(newz11,1),'k',rlat,nanmean(newz12,1),'b',rlat,nanmean(newz13,1),'gx-',rlat,nanmean(newz21,1),'r',rlat,nanmean(newz22,1),'m','linewidth',4);
plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('NIGHT d mmw/dt [mm/yr]');
