ii = 0; 

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fUMBC_day.RHsurftrend,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : CHIRP\_A DAY');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fUMBC_night.RHsurftrend,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : CHIRP\_A NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fERA5_day.RHsurftrend,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : ERA5 DAY');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fERA5_night.RHsurftrend,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : ERA5 NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fAIRSL3_day.RHsurftrend,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : AIRS L3 DAY');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fAIRSL3_night.RHsurftrend,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : AIRS L3 NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fCLIMCAPSL3_day.RHsurftrend,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : CLIMCAPS L3 DAY');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fCLIMCAPSL3_night.RHsurftrend,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : CLIMCAPS L3 NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(fMERRA2.RHsurftrend,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : MERRA2 DAY/NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fUMBC_day.RHsurftrend+fUMBC_night.RHsurftrend)*0.5,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : CHIRP\_A DAY/NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fAIRSL3_day.RHsurftrend+fAIRSL3_night.RHsurftrend)*0.5,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : AIRS L3 DAY/NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fCLIMCAPSL3_day.RHsurftrend+fCLIMCAPSL3_night.RHsurftrend)*0.5,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : CLIMCAPS L3 DAY/NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fERA5_day.RHsurftrend+fERA5_night.RHsurftrend)*0.5,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : ERA5 DAY/NIGHT');

ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fUMBC_day.RHsurftrend-fUMBC_night.RHsurftrend)*1.0,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : CHIRP\_A DAY - NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fAIRSL3_day.RHsurftrend-fAIRSL3_night.RHsurftrend)*1.0,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : AIRS L3 DAY - NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fCLIMCAPSL3_day.RHsurftrend-fCLIMCAPSL3_night.RHsurftrend)*1.0,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : CLIMCAPS L3 DAY - NIGHT');
ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((fERA5_day.RHsurftrend-fERA5_night.RHsurftrend)*1.0,72,64)',1),[-90 +90],[-180 +180]); title('dRHsurf/dt : ERA5 DAY - NIGHT');

iiMax = ii;
for ii = 1 : iiMax; figure(ii); caxis([-1 +1]*0.15); colormap(llsmap5); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newz*
iFig = 70; 
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dRHsurf/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = '(D) CHIRP\_A';     plotoptions.str21 = '(N)';
  plotoptions.str12 = '(D) AIRS L3';      plotoptions.str22 = '(N)';
  plotoptions.str13 = '(D) CLIMCAPS L3';  plotoptions.str23 = '(N)';
  plotoptions.str14 = '(D) ERA5';         plotoptions.str24 = '(N)';

  newz11 = fUMBC_day.RHsurftrend;         newz21 = fUMBC_night.RHsurftrend;                           
  newz12 = fAIRSL3_day.RHsurftrend;       newz22 = fAIRSL3_night.RHsurftrend; 
  newz13 = fCLIMCAPSL3_day.RHsurftrend;   newz23 = fCLIMCAPSL3_night.RHsurftrend;
  newz14 = fERA5_day.RHsurftrend;         newz24 = fERA5_night.RHsurftrend; 
    newz11 = reshape(newz11,72,64);   newz21 = reshape(newz21,72,64);
    newz12 = reshape(newz12,72,64);   newz22 = reshape(newz22,72,64);
    newz13 = reshape(newz13,72,64);   newz23 = reshape(newz23,72,64);
    newz14 = reshape(newz14,72,64);   newz24 = reshape(newz24,72,64);
  figure(iFig); clf; aslmap_2x4tiledlayout(newz11,newz12,newz13,newz14,newz21,newz22,newz23,newz24,iFig,plotoptions);

clear newz*
iFig = 70; 
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = '(D) CHIRP\_A';     plotoptions.str21 = '(N)';   plotoptions.str31 = '(D-N)';
  plotoptions.str12 = '(D) AIRS L3';      plotoptions.str22 = '(N)';   plotoptions.str32 = '(D-N)';
  plotoptions.str13 = '(D) CLIMCAPS L3';  plotoptions.str23 = '(N)';   plotoptions.str33 = '(D-N)';
  plotoptions.str14 = '(D) ERA5';         plotoptions.str24 = '(N)';   plotoptions.str34 = '(D-N)';

  newz11 = fUMBC_day.RHsurftrend;                 newz21 = fUMBC_night.RHsurftrend;            newz31 = newz11 - newz21;
  newz12 = fAIRSL3_day.RHsurftrend;               newz22 = fAIRSL3_night.RHsurftrend;          newz32 = newz12 - newz22;
  newz13 = fCLIMCAPSL3_day.RHsurftrend;           newz23 = fCLIMCAPSL3_night.RHsurftrend;      newz33 = newz13 - newz23;
  newz14 = fERA5_day.RHsurftrend;                 newz24 = fERA5_night.RHsurftrend;            newz34 = newz14 - newz24;
    newz11 = reshape(newz11,72,64);   newz21 = reshape(newz21,72,64);   newz31 = reshape(newz31,72,64);   
    newz14 = reshape(newz14,72,64);   newz24 = reshape(newz24,72,64);   newz34 = reshape(newz34,72,64);   
  figure(iFig); clf; aslmap_3x4tiledlayout(newz11,newz12,newz13,newz14,newz21,newz22,newz23,newz24,newz31,newz32,newz33,newz34,iFig,plotoptions);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newz*
iFig = 71; 
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dRHsurf/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'CHIRP\_A';     plotoptions.str12 = 'ERA5';        
  plotoptions.str21 = 'AIRS L3';      plotoptions.str22 = 'CLIMCAPS L3'; 
  newz11 = (fUMBC_day.RHsurftrend-fUMBC_night.RHsurftrend)*1.0;                           newz12 = (fERA5_day.RHsurftrend-fERA5_night.RHsurftrend)*1.0; 
  newz21 = (fAIRSL3_day.RHsurftrend-fAIRSL3_night.RHsurftrend)*1.0; newz22 = (fCLIMCAPSL3_day.RHsurftrend-fCLIMCAPSL3_night.RHsurftrend)*1.0;
    newz11 = reshape(newz11,72,64);   newz12 = reshape(newz12,72,64);
  figure(iFig); clf; aslmap_2x2tiledlayout(newz11,newz12,newz21,newz22,iFig,plotoptions);

clear newz*
iFig = 72;
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dRHsurf/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'CHIRP\_A';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'N/A';          plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  newz11 = (fUMBC_day.RHsurftrend+fUMBC_night.RHsurftrend)*0.5;    newz12 = (fAIRSL3_day.RHsurftrend+fAIRSL3_night.RHsurftrend)*0.5;   newz13 = (fCLIMCAPSL3_day.RHsurftrend+fCLIMCAPSL3_night.RHsurftrend)*0.5;
  newz21 = fGISS.giss_trend4608*NaN;                         newz22 = (fERA5_day.RHsurftrend+fERA5_night.RHsurftrend)*0.5;       newz23 = fMERRA2.RHsurftrend;
    newz11 = reshape(newz11,72,64);   newz22 = reshape(newz22,72,64);   newz23 = reshape(newz23,72,64);  
  figure(iFig); clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);
  junk = [newz11(:) newz12(:) newz13(:) newz21(:) newz22(:) newz23(:)];

clear newz*
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dRHsurf/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS L3';   plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.strzz = 'CHIRP\_A';  
  plotoptions.str21 = 'MERRA2';    plotoptions.str22 = 'ERA5';
  newzAA = (fUMBC_day.RHsurftrend+fUMBC_night.RHsurftrend)*0.5;    newz11 = (fAIRSL3_day.RHsurftrend+fAIRSL3_night.RHsurftrend)*0.5; newz12 = (fCLIMCAPSL3_day.RHsurftrend+fCLIMCAPSL3_night.RHsurftrend)*0.5;
  newz22 = (fERA5_day.RHsurftrend+fERA5_night.RHsurftrend)*0.5;    newz21 = fMERRA2.RHsurftrend;
    newz11 = reshape(newz11,72,64);   newz12 = reshape(newz12,72,64);  newz21 = reshape(newz21,72,64);   newz22 = reshape(newz22,72,64);  newzAA = reshape(newzAA,72,64);    
  aslmap_2x1x2tiledlayout(newz11,newz12,newzAA,newz21,newz22,iFig,plotoptions);
  junk = [newz11(:) newz12(:) newzAA(:) newz21(:) newz22(:)];

clear newz*
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dRHsurf/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'CHIRP\_A';   plotoptions.str12 = 'AIRS L3';   plotoptions.str13 = 'CLIMCAPS L3';  
  plotoptions.str21 = 'ERA5';      plotoptions.str22 = 'MERRA2';
  newz11 = (fUMBC_day.RHsurftrend+fUMBC_night.RHsurftrend)*0.5;    newz12 = (fAIRSL3_day.RHsurftrend+fAIRSL3_night.RHsurftrend)*0.5; newz13 = (fCLIMCAPSL3_day.RHsurftrend+fCLIMCAPSL3_night.RHsurftrend)*0.5;
  newz21 = (fERA5_day.RHsurftrend+fERA5_night.RHsurftrend)*0.5;    newz22 = fMERRA2.RHsurftrend;
    newz11 = reshape(newz11,72,64);   newz12 = reshape(newz12,72,64);  newz21 = reshape(newz21,72,64);   newz22 = reshape(newz22,72,64);  newz13 = reshape(newz13,72,64);    
  aslmap_2x1x2tiledlayout_tall(newz11,newz12,newz13,newz21,newz22,iFig,plotoptions);  
  junk = [newz11(:) newz12(:) newz13(:) newz21(:) newz22(:)];

moo0 = nanmean(junk,2);
moo1 = abs(nanmean(junk,2));
moo2 = nanstd(junk,[],2);
moo3 = zeros(size(moo1)); moo3(moo1 > moo2) = 1;
moo4 = moo1./moo2; moo4(moo4 < 1) = 1; moo4 = log10(moo4);

iFig = 73; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(nanmean(junk,2),72,64)',1),[-90 +90],[-180 +180]);   cx = [-1 +1]*0.151; caxis(cx); colormap(llsmap5);
iFig = 74; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(nanstd(junk,[],2),72,64)',1),[-90 +90],[-180 +180]); cx = [0 +1]*0.071;  caxis(cx); colormap(jet);

iFig = 75; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(moo3,72,64)',1),[-90 +90],[-180 +180]); cx = [0 +1];  caxis(cx); colormap(jet);
iFig = 75; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,reshape(moo3,72,64)',[-90 +90],[-180 +180]); cx = [0 +1]; caxis(cx); colormap(usa2x);

iFig = 76; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,reshape(moo4,72,64)',[-90 +90],[-180 +180]);        cx = [0 +2]; caxis(cx); colormap(usa2x);
iFig = 76; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,reshape(10.^(moo4),72,64)',[-90 +90],[-180 +180]);  cx = (10.^[0 +1]); caxis(cx); colormap(jett);

mooUMBC_num   = 0.5*(fUMBC_day.RHsurftrend + fUMBC_night.RHsurftrend);
mooUMBC_denom = sqrt(fUMBC_day.RHsurftrendunc.^2 + fUMBC_night.RHsurftrendunc.^2);
mooUMBC = abs(mooUMBC_num)./mooUMBC_denom;
iFig = 77; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(mooUMBC_num,72,64)',1),[-90 +90],[-180 +180]);   cx = ([-1 +1]*0.15); caxis(cx); colormap(llsmap5);
iFig = 78; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(mooUMBC_denom,72,64)',1),[-90 +90],[-180 +180]); cx = ([0 +1]*0.15/2); caxis(cx); colormap(jet);
iFig = 79; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,reshape(10.^(mooUMBC),72,64)',[-90 +90],[-180 +180]);            cx = (10.^[0 +1]); caxis(cx); colormap(jett);
iFig = 79; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,reshape(mooUMBC,72,64)',[-90 +90],[-180 +180]);                  cx = ([0 +5]); caxis(cx); colormap(jett);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(96); clf
  ta = tiledlayout(2,1,'TileSpacing','None', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tfov(1) = nexttile;
    aslmapSergio(rlat65,rlon73,smoothn(reshape(moo0,72,64)',1), [-90 +90],[-180 +180]); colormap(tfov(1),llsmap5); caxis([-1 +1]*0.151); colorbar
    text(-2.5,1.15,'(a)')
  tfov(2) = nexttile;
    aslmapSergio(rlat65,rlon73,smoothn(reshape(moo2,72,64)',1), [-90 +90],[-180 +180]); colormap(tfov(2),jet); caxis([ 0 +1]*0.151); colorbar
    text(-2.5,1.15,'(b)')
colormap(tfov(1),llsmap5);
colormap(tfov(2),jet);
ta.Padding = 'compact';
ta.TileSpacing = 'compact';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OCEAN + LAND
newzST11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  newzST12 = (fAIRSL3_day.thestats64x72.stemprate+fAIRSL3_night.thestats64x72.stemprate)*0.5; newzST13 = (fCLIMCAPSL3_day.thestats64x72.stemprate+fCLIMCAPSL3_night.thestats64x72.stemprate)*0.5;
newzST21 = (fERA5_day.trend_stemp+fERA5_night.trend_stemp)*0.5;    newzST22 = fMERRA2.trend_stemp;
newzST11 = reshape(newzST11,72,64);   newzST12 = reshape(newzST12,72,64);  newzST21 = reshape(newzST21,72,64);   newzST22 = reshape(newzST22,72,64);  newzST13 = reshape(newzST13,72,64);    

figure(80); clf
plot(rlat,nanmean(newz11,1),'k',rlat,nanmean(newz12,1),'b',rlat,nanmean(newz13,1),'gx-',rlat,nanmean(newz21,1),'r',rlat,nanmean(newz22,1),'m','linewidth',2);
plotaxis2; hl = legend('CHIRP\_A','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('dRHsurf/dt');

figure(81); clf
plot(rlat,nanmean(newzST11,1),'k',rlat,nanmean(newzST12,1),'b',rlat,nanmean(newzST13,1),'gx-',rlat,nanmean(newzST21,1),'r',rlat,nanmean(newzST22,1),'m','linewidth',2);
plotaxis2; hl = legend('CHIRP\_A','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('dST/dt');

figure(82); clf
junkA = newz11(:); junkB = newzST11(:); junk11 = junkA./junkB; junk11(abs(junkB) < 1e-4) = NaN; junk11 = reshape(junk11,72,64); % plot(rlat,nanmean(junk11,1));
junkA = newz12(:); junkB = newzST12(:); junk12 = junkA./junkB; junk12(abs(junkB) < 1e-4) = NaN; junk12 = reshape(junk12,72,64); % plot(rlat,nanmean(junk12,1));
junkA = newz13(:); junkB = newzST13(:); junk13 = junkA./junkB; junk13(abs(junkB) < 1e-4) = NaN; junk13 = reshape(junk13,72,64); % plot(rlat,nanmean(junk13,1));
junkA = newz21(:); junkB = newzST21(:); junk21 = junkA./junkB; junk21(abs(junkB) < 1e-4) = NaN; junk21 = reshape(junk21,72,64); % plot(rlat,nanmean(junk21,1));
junkA = newz22(:); junkB = newzST22(:); junk22 = junkA./junkB; junk22(abs(junkB) < 1e-4) = NaN; junk22 = reshape(junk22,72,64); % plot(rlat,nanmean(junk22,1));
plot(rlat,nanmean(junk11,1),'k',rlat,nanmean(junk12,1),'b',rlat,nanmean(junk13,1),'gx-',rlat,nanmean(junk21,1),'r',rlat,nanmean(junk22,1),'m','linewidth',2);
plot(rlat,smooth(nanmean(junk11,1),10),'k',rlat,smooth(nanmean(junk12,1),10),'b',rlat,smooth(nanmean(junk13,1),10),'gx-',rlat,smooth(nanmean(junk21,1),10),'r',rlat,smooth(nanmean(junk22,1),10),'m','linewidth',2);
ylim([-1 +1]*10);
plotaxis2; hl = legend('CHIRP\_A','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('dRHsurf/dST Land+Ocean'); ylim([-1 +1]*10);

figure(82); clf
plot(rlat,nanmean(newz11,1)./nanmean(newzST11,1),'k',rlat,nanmean(newz12,1)./nanmean(newzST12,1),'b',rlat,nanmean(newz13,1)./nanmean(newzST13,1),'gx-',...
     rlat,nanmean(newz21,1)./nanmean(newzST21,1),'r',rlat,nanmean(newz22,1)./nanmean(newzST22,1),'m','linewidth',2);
plot(rlat,smooth(nanmean(newz11,1)./nanmean(newzST11,1),10),'k',rlat,smooth(nanmean(newz12,1)./nanmean(newzST12,1),10),'b',rlat,smooth(nanmean(newz13,1)./nanmean(newzST13,1),10),'gx-',...
     rlat,smooth(nanmean(newz21,1)./nanmean(newzST21,1),10),'r',rlat,smooth(nanmean(newz22,1)./nanmean(newzST22,1),10),'m','linewidth',2);
ylim([-1 +1]*20);
plotaxis2; hl = legend('CHIRP\_A','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('dRHsurf/dST Ocean'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OCEAN ONLY
newz11 = (fUMBC_day.RHsurftrend+fUMBC_night.RHsurftrend)*0.5;    newz12 = (fAIRSL3_day.RHsurftrend+fAIRSL3_night.RHsurftrend)*0.5; newz13 = (fCLIMCAPSL3_day.RHsurftrend+fCLIMCAPSL3_night.RHsurftrend)*0.5;
newz21 = (fERA5_day.RHsurftrend+fERA5_night.RHsurftrend)*0.5;    newz22 = fMERRA2.RHsurftrend;
  newz11 = reshape(newz11,72,64);   newz12 = reshape(newz12,72,64);  newz21 = reshape(newz21,72,64);   newz22 = reshape(newz22,72,64);  newz13 = reshape(newz13,72,64);    
newzST11 = (fUMBC_day.results(:,6)+fUMBC_night.results(:,6))*0.5;  newzST12 = (fAIRSL3_day.thestats64x72.stemprate+fAIRSL3_night.thestats64x72.stemprate)*0.5; newzST13 = (fCLIMCAPSL3_day.thestats64x72.stemprate+fCLIMCAPSL3_night.thestats64x72.stemprate)*0.5;
newzST21 = (fERA5_day.trend_stemp+fERA5_night.trend_stemp)*0.5;    newzST22 = fMERRA2.trend_stemp;
newzST11 = reshape(newzST11,72,64);   newzST12 = reshape(newzST12,72,64);  newzST21 = reshape(newzST21,72,64);   newzST22 = reshape(newzST22,72,64);  newzST13 = reshape(newzST13,72,64);    

newz11(landfrac > eps) = NaN;   newz12(landfrac > eps) = NaN;   newz13(landfrac > eps) = NaN;   newz21(landfrac > eps) = NaN;   newz22(landfrac > eps) = NaN;
newzST11(landfrac > eps) = NaN; newzST12(landfrac > eps) = NaN; newzST13(landfrac > eps) = NaN; newzST21(landfrac > eps) = NaN; newzST22(landfrac > eps) = NaN;

%figure(86);
landfrac = reshape(p.landfrac,72,64);
moo_landfrac = ones(size(landfrac)); moo_landfrac(landfrac < eps) = 0; moo_oceanfrac = 1 - moo_landfrac;
%plot(rlat,nansum(moo_oceanfrac,1));
%plot(rlat,nanmean(moo_oceanfrac,1));
%plot(rlat,nanmean(moo_oceanfrac,1).*cos(rlat'*pi/180));
xrlat  = nanmean(moo_oceanfrac,1).*cos(rlat'*pi/180);
x2rlat = nanmean(moo_oceanfrac,1).*sin(rlat'*pi/180);
%plot(rlat,xrlat,rlat,cos(rlat*pi/180))
%plot(rlat,x2rlat,rlat,sin(rlat*pi/180))
%plot(xrlat,nanmean(newzST11,1),'k',xrlat,nanmean(newzST12,1),'b',xrlat,nanmean(newzST13,1),'gx-',xrlat,nanmean(newzST21,1),'r',xrlat,nanmean(newzST22,1),'m','linewidth',2);
%plot(xrlat.*sign(rlat'),nanmean(newzST11,1),'k',xrlat.*sign(rlat'),nanmean(newzST12,1),'b',xrlat.*sign(rlat'),nanmean(newzST13,1),'gx-',xrlat.*sign(rlat'),nanmean(newzST21,1),'r',xrlat.*sign(rlat'),nanmean(newzST22,1),'m','linewidth',2);
%plot(x2rlat,nanmean(newzST11,1),'k',x2rlat,nanmean(newzST12,1),'b',x2rlat,nanmean(newzST13,1),'gx-',x2rlat,nanmean(newzST21,1),'r',x2rlat,nanmean(newzST22,1),'m','linewidth',2);
%plotaxis2; hl = legend('CHIRP\_A','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
%xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('dST/dt');

figure(83); clf
plot(rlat,nanmean(newz11,1),'k',rlat,nanmean(newz12,1),'b',rlat,nanmean(newz13,1),'gx-',rlat,nanmean(newz21,1),'r',rlat,nanmean(newz22,1),'m','linewidth',2);
plotaxis2; hl = legend('CHIRP\_A','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('dRHsurf/dt Ocean');

figure(84); clf
plot(rlat,nanmean(newzST11,1),'k',rlat,nanmean(newzST12,1),'b',rlat,nanmean(newzST13,1),'gx-',rlat,nanmean(newzST21,1),'r',rlat,nanmean(newzST22,1),'m','linewidth',2);
plotaxis2; hl = legend('CHIRP\_A','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('dST/dt Ocean');

figure(85); clf
junkA = newz11(:); junkB = newzST11(:); junk11 = junkA./junkB; junk11(abs(junkB) < 1e-4 | p.landfrac' > eps) = NaN; junk11 = reshape(junk11,72,64); % plot(rlat,nanmean(junk11,1));
junkA = newz12(:); junkB = newzST12(:); junk12 = junkA./junkB; junk12(abs(junkB) < 1e-4 | p.landfrac' > eps) = NaN; junk12 = reshape(junk12,72,64); % plot(rlat,nanmean(junk12,1));
junkA = newz13(:); junkB = newzST13(:); junk13 = junkA./junkB; junk13(abs(junkB) < 1e-4 | p.landfrac' > eps) = NaN; junk13 = reshape(junk13,72,64); % plot(rlat,nanmean(junk13,1));
junkA = newz21(:); junkB = newzST21(:); junk21 = junkA./junkB; junk21(abs(junkB) < 1e-4 | p.landfrac' > eps) = NaN; junk21 = reshape(junk21,72,64); % plot(rlat,nanmean(junk21,1));
junkA = newz22(:); junkB = newzST22(:); junk22 = junkA./junkB; junk22(abs(junkB) < 1e-4 | p.landfrac' > eps) = NaN; junk22 = reshape(junk22,72,64); % plot(rlat,nanmean(junk22,1));
plot(rlat,nanmean(junk11,1),'k',rlat,nanmean(junk12,1),'b',rlat,nanmean(junk13,1),'gx-',rlat,nanmean(junk21,1),'r',rlat,nanmean(junk22,1),'m','linewidth',2);
plot(rlat,smooth(nanmean(junk11,1),10),'k',rlat,smooth(nanmean(junk12,1),10),'b',rlat,smooth(nanmean(junk13,1),10),'gx-',rlat,smooth(nanmean(junk21,1),10),'r',rlat,smooth(nanmean(junk22,1),10),'m','linewidth',2);
ylim([-1 +1]*20);
plotaxis2; hl = legend('CHIRP\_A','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('dRHsurf/dST Ocean'); 

figure(85); clf
plot(rlat,nanmean(newz11,1)./nanmean(newzST11,1),'k',rlat,nanmean(newz12,1)./nanmean(newzST12,1),'b',rlat,nanmean(newz13,1)./nanmean(newzST13,1),'gx-',...
     rlat,nanmean(newz21,1)./nanmean(newzST21,1),'r',rlat,nanmean(newz22,1)./nanmean(newzST22,1),'m','linewidth',2);
plot(rlat,smooth(nanmean(newz11,1)./nanmean(newzST11,1),10),'k',rlat,smooth(nanmean(newz12,1)./nanmean(newzST12,1),10),'b',rlat,smooth(nanmean(newz13,1)./nanmean(newzST13,1),10),'gx-',...
     rlat,smooth(nanmean(newz21,1)./nanmean(newzST21,1),10),'r',rlat,smooth(nanmean(newz22,1)./nanmean(newzST22,1),10),'m','linewidth',2);
ylim([-1 +1]*10);
plotaxis2; hl = legend('CHIRP\_A','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlim([-1 +1]*90); xlabel('Latitude [deg]'); ylabel('dRHsurf/dST Ocean'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear junk*
