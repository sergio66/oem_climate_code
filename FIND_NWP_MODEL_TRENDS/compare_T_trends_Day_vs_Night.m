ii = 0; 

ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fUMBC_day.ptemptrend,100,72,64),2)),1)); title('dT/dt : AIRS\_RT DAY');
ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fUMBC_night.ptemptrend,100,72,64),2)),1)); title('dT/dt : AIRS\_RT NIGHT');

ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fERA5_day.ptemptrend,100,72,64),2)),1)); title('dT/dt : ERA5 DAY');
ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fERA5_night.ptemptrend,100,72,64),2)),1)); title('dT/dt : ERA5 NIGHT');

ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fAIRSL3_day.ptemptrend,100,72,64),2)),1)); title('dT/dt : AIRSL3 DAY');
ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fAIRSL3_night.ptemptrend,100,72,64),2)),1)); title('dT/dt : AIRSL3 NIGHT');

ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fCLIMCAPSL3_day.ptemptrend,100,72,64),2)),1)); title('dT/dt : CLIMCAPSL3 DAY');
ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fCLIMCAPSL3_night.ptemptrend,100,72,64),2)),1)); title('dT/dt : CLIMCAPSL3 NIGHT');

ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fMERRA2.ptemptrend,100,72,64),2)),1)); title('dT/dt : MERRA2 DAY');

ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fUMBC_day.ptemptrend + fUMBC_night.ptemptrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/2); title('dT/dt : AIRS\_RT DAY/NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fERA5_day.ptemptrend + fERA5_night.ptemptrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/2); title('dT/dt : ERA5 DAY/NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fAIRSL3_day.ptemptrend + fAIRSL3_night.ptemptrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/2); title('dT/dt : AIRSL3 DAY/NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fCLIMCAPSL3_day.ptemptrend + fCLIMCAPSL3_night.ptemptrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/2); title('dT/dt : CLIMCAPSL3 DAY/NIGHT');

ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fUMBC_day.ptemptrend - fUMBC_night.ptemptrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/1); title('dT/dt : AIRS\_RT DAY - NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fERA5_day.ptemptrend - fERA5_night.ptemptrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/1); title('dT/dt : ERA5 DAY - NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fAIRSL3_day.ptemptrend - fAIRSL3_night.ptemptrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/1); title('dT/dt : AIRSL3 DAY - NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fCLIMCAPSL3_day.ptemptrend - fCLIMCAPSL3_night.ptemptrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/1); title('dT/dt : CLIMCAPSL3 DAY - NIGHT');

iiMax = ii;
for ii = 1 : iiMax; figure(ii); caxis([-1 +1]*0.15); colormap(llsmap5); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newz*
iFig = 40;
  figure(iFig); clf;
  clear plotoptions
  plotoptions.xstr = ' '; plotoptions.ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dT/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = '(D) AIRS\_RT';     plotoptions.str21 = '(N)';
  plotoptions.str12 = '(D) AIRS L3';      plotoptions.str22 = '(N)';
  plotoptions.str13 = '(D) CLIMCAPS L3';  plotoptions.str23 = '(N)';
  plotoptions.str14 = '(D) ERA5';         plotoptions.str24 = '(N)';
  plotoptions.barstr = 'dT/dt [K/yr]';
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits      = [10 1000];

  newz11 = fUMBC_day.ptemptrend;        newz21 = fUMBC_night.ptemptrend;                           
  newz12 = fAIRSL3_day.ptemptrend;      newz22 = fAIRSL3_night.ptemptrend; 
  newz13 = fCLIMCAPSL3_day.ptemptrend;  newz23 = fCLIMCAPSL3_night.ptemptrend;
  newz14 = fERA5_day.ptemptrend;        newz24 = fERA5_night.ptemptrend; 
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2)); 
    newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2)); 
    newz13 = squeeze(nanmean(reshape(newz13,100,72,64),2));     newz23 = squeeze(nanmean(reshape(newz23,100,72,64),2)); 
    newz14 = squeeze(nanmean(reshape(newz14,100,72,64),2));     newz24 = squeeze(nanmean(reshape(newz24,100,72,64),2)); 
  profile_plots_2x4tiledlayout(rlat,plays100,newz11,newz12,newz13,newz14,newz21,newz22,newz23,newz24,iFig,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newz*
iFig = 41; 
  figure(iFig); clf;
  clear plotoptions
  plotoptions.xstr = ' '; plotoptions.ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dT/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'ERA5';        
  plotoptions.str21 = 'AIRS L3';      plotoptions.str22 = 'CLIMCAPS L3'; 
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits      = [10 1000];
  plotoptions.barstr = 'dT/dt [K/yr]';
  newz11 = (fUMBC_day.ptemptrend-fUMBC_night.ptemptrend)*1.0;        newz12 = (fERA5_day.ptemptrend-fERA5_night.ptemptrend)*1.0; 
  newz21 = (fAIRSL3_day.ptemptrend-fAIRSL3_night.ptemptrend)*1.0;    newz22 = (fCLIMCAPSL3_day.ptemptrend-fCLIMCAPSL3_night.ptemptrend)*1.0;
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2)); 
    newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2)); 
  profile_plots_2x2tiledlayout(rlat,plays100,newz11,newz12,newz21,newz22,iFig,plotoptions);

clear newz*
iFig = 42;
  figure(iFig); clf;
  clear plotoptions
  plotoptions.xstr = ' '; plotoptions.ystr = 'P [mb]';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dT/dt'; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'N/A';          plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits      = [10 1000];
  plotoptions.barstr = 'dT/dt [K/yr]';
  newz11 = (fUMBC_day.ptemptrend+fUMBC_night.ptemptrend)*0.5;  newz12 = (fAIRSL3_day.ptemptrend+fAIRSL3_night.ptemptrend)*0.5; newz13 = (fCLIMCAPSL3_day.ptemptrend+fCLIMCAPSL3_night.ptemptrend)*0.5;
  newz21 = newz11*0;                                           newz22 = (fERA5_day.ptemptrend+fERA5_night.ptemptrend)*0.5;     newz23 = fMERRA2.ptemptrend;
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2));     newz13 = squeeze(nanmean(reshape(newz13,100,72,64),2)); 
    newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2));     newz23 = squeeze(nanmean(reshape(newz23,100,72,64),2)); 
  profile_plots_2x3tiledlayout(rlat,plays100,newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);

clear newz*
  figure(iFig); clf;
  newz11 = (fUMBC_day.ptemptrend+fUMBC_night.ptemptrend)*0.5;  newz12 = (fAIRSL3_day.ptemptrend+fAIRSL3_night.ptemptrend)*0.5; newz13 = (fCLIMCAPSL3_day.ptemptrend+fCLIMCAPSL3_night.ptemptrend)*0.5;
                                       newz21 = fMERRA2.ptemptrend;  newz22 = (fERA5_day.ptemptrend+fERA5_night.ptemptrend)*0.5; 
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2));     newz13 = squeeze(nanmean(reshape(newz13,100,72,64),2)); 
              newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2));     
    plotoptions2x1x2 = plotoptions;
    plotoptions2x1x2.str11 = 'AIRS\_RT';
    plotoptions2x1x2.str12 = 'AIRS L3';
    plotoptions2x1x2.str13 = 'CLIMCAPS';
    plotoptions2x1x2.str21 = 'MERRA2';
    plotoptions2x1x2.str22 = 'ERA5';
    plotoptions.barstr = 'dT/dt [K/yr]';
plotoptions2x1x2.cx = [-1 +1]*0.151;
plotoptions2x1x2.cx = [-1 +1]*0.101;
      profile_plots_2x1x2tiledlayout_tall(rlat,plays100,newz11,newz12,newz13,newz21,newz22,iFig,plotoptions2x1x2);

moo = find(plays100 > 10,1);
cjunkUMBC = newz11(moo:97,:); cjunkL3 = newz12(moo:97,:); cjunkC3 = newz13(moo:97,:); cjunkMERRA2 = newz21(moo:97,:); cjunkERA5 = newz22(moo:97,:);
da5x5correl = ...
[nanlinearcorrelation(cjunkUMBC(:),cjunkUMBC(:))  nanlinearcorrelation(cjunkUMBC(:),cjunkL3(:))   nanlinearcorrelation(cjunkUMBC(:),cjunkC3(:))   nanlinearcorrelation(cjunkUMBC(:),cjunkMERRA2(:))   nanlinearcorrelation(cjunkUMBC(:),cjunkERA5(:)); ...
nanlinearcorrelation(cjunkL3(:),cjunkUMBC(:))     nanlinearcorrelation(cjunkL3(:),cjunkL3(:))     nanlinearcorrelation(cjunkL3(:),cjunkC3(:))     nanlinearcorrelation(cjunkL3(:),cjunkMERRA2(:))     nanlinearcorrelation(cjunkL3(:),cjunkERA5(:)); ...
nanlinearcorrelation(cjunkC3(:),cjunkUMBC(:))     nanlinearcorrelation(cjunkC3(:),cjunkL3(:))     nanlinearcorrelation(cjunkC3(:),cjunkC3(:))     nanlinearcorrelation(cjunkC3(:),cjunkMERRA2(:))     nanlinearcorrelation(cjunkC3(:),cjunkERA5(:)); ...
nanlinearcorrelation(cjunkMERRA2(:),cjunkUMBC(:)) nanlinearcorrelation(cjunkMERRA2(:),cjunkL3(:)) nanlinearcorrelation(cjunkMERRA2(:),cjunkC3(:)) nanlinearcorrelation(cjunkMERRA2(:),cjunkMERRA2(:)) nanlinearcorrelation(cjunkMERRA2(:),cjunkERA5(:));...
nanlinearcorrelation(cjunkERA5(:),cjunkUMBC(:))   nanlinearcorrelation(cjunkERA5(:),cjunkL3(:))   nanlinearcorrelation(cjunkERA5(:),cjunkC3(:))   nanlinearcorrelation(cjunkERA5(:),cjunkMERRA2(:))   nanlinearcorrelation(cjunkERA5(:),cjunkERA5(:));];
printarray(da5x5correl,'nanlinearcorrelation T trends : UMBC  AIRSL3  CLIMCAPSL3  MERRA2  ERA5');

[corrcoef(cjunkUMBC(:),cjunkUMBC(:))  corrcoef(cjunkUMBC(:),cjunkL3(:))   corrcoef(cjunkUMBC(:),cjunkC3(:))   corrcoef(cjunkUMBC(:),cjunkMERRA2(:))   corrcoef(cjunkUMBC(:),cjunkERA5(:)); ...
corrcoef(cjunkL3(:),cjunkUMBC(:))     corrcoef(cjunkL3(:),cjunkL3(:))     corrcoef(cjunkL3(:),cjunkC3(:))     corrcoef(cjunkL3(:),cjunkMERRA2(:))     corrcoef(cjunkL3(:),cjunkERA5(:)); ...
corrcoef(cjunkC3(:),cjunkUMBC(:))     corrcoef(cjunkC3(:),cjunkL3(:))     corrcoef(cjunkC3(:),cjunkC3(:))     corrcoef(cjunkC3(:),cjunkMERRA2(:))     corrcoef(cjunkC3(:),cjunkERA5(:)); ...
corrcoef(cjunkMERRA2(:),cjunkUMBC(:)) corrcoef(cjunkMERRA2(:),cjunkL3(:)) corrcoef(cjunkMERRA2(:),cjunkC3(:)) corrcoef(cjunkMERRA2(:),cjunkMERRA2(:)) corrcoef(cjunkMERRA2(:),cjunkERA5(:));...
corrcoef(cjunkERA5(:),cjunkUMBC(:))   corrcoef(cjunkERA5(:),cjunkL3(:))   corrcoef(cjunkERA5(:),cjunkC3(:))   corrcoef(cjunkERA5(:),cjunkMERRA2(:))   corrcoef(cjunkERA5(:),cjunkERA5(:));];
printarray(da5x5correl,'corrcoef             T trends : UMBC  AIRSL3  CLIMCAPSL3  MERRA2  ERA5');

%junk = [newz11(:) newz12(:) newz13(:) newz21(:) newz22(:) newz23(:)];
%junk = reshape(junk,100,64,6);
junk = [newz11(:) newz12(:) newz13(:) newz21(:) newz22(:)];
junk = reshape(junk,100,64,5);
moo0 = nanmean(junk,3);
moo1 = abs(nanmean(junk,3));
moo2 = nanstd(junk,[],3);
moo3 = zeros(size(moo1)); moo3(moo1 > moo2) = 1;
moo4 = moo1./moo2; moo4(moo4 < 1) = 1; moo4 = log10(moo4);
 
iFig = 43; figure(iFig); clf; pcolor(rlat,plays100,smoothn(nanmean(junk,3),1));   cx = [-1 +1]*0.151; caxis(cx); colormap(llsmap5);
iFig = 44; figure(iFig); clf; pcolor(rlat,plays100,smoothn(nanstd(junk,[],3),1)); cx = [0 +1]*0.031;  caxis(cx); colormap(jet);

iFig = 45; figure(iFig); clf; pcolor(rlat,plays100,smoothn(moo3,1));         cx = [0 +1];  caxis(cx); colormap(jet);
iFig = 45; figure(iFig); clf; pcolor(rlat,plays100,moo3);                    cx = [0 +1]; caxis(cx); colormap(usa2x);

iFig = 46; figure(iFig); clf; pcolor(rlat,plays100,smoothn(moo4,1));         cx = [0 +2];       caxis(cx); colormap(usa2x);
iFig = 46; figure(iFig); clf; pcolor(rlat,plays100,10.^moo4);                cx = (10.^[0 +1]); caxis(cx); colormap(jett);

mooUMBC_num   = 0.5*(fUMBC_day.ptemptrend + fUMBC_night.ptemptrend);                 mooUMBC_num   = squeeze(nanmean(reshape(mooUMBC_num,100,72,64),2)); 
mooUMBC_denom = sqrt(fUMBC_day.ptemptrendunc.^2 + fUMBC_night.ptemptrendunc.^2);     mooUMBC_denom = squeeze(nanmean(reshape(mooUMBC_denom,100,72,64),2)); 
mooUMBC = abs(mooUMBC_num)./mooUMBC_denom;
iFig = 47; figure(iFig); clf; pcolor(rlat,plays100,smoothn(mooUMBC_num,1));   cx = ([-1 +1]*0.15);  caxis(cx); colormap(llsmap5);
iFig = 48; figure(iFig); clf; pcolor(rlat,plays100,smoothn(mooUMBC_denom,1)); cx = ([0 +1]*0.15/2); caxis(cx); colormap(jet);
iFig = 49; figure(iFig); clf; pcolor(rlat,plays100,10.^(mooUMBC));            cx = (10.^[0 +1]);    caxis(cx); colormap(jett);
iFig = 49; figure(iFig); clf; pcolor(rlat,plays100,mooUMBC);                  cx = ([0 +1]);        caxis(cx); colormap(jett);

for ii = 43 : 49; figure(ii); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); colorbar; shading interp; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(93); clf
  ta = tiledlayout(2,1,'TileSpacing','None', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tfov(1) = nexttile;
    pcolor(rlat,plays100,smoothn(moo0,1));
    colorbar; shading interp; set(gca,'ydir','reverse','yscale','log'); ylim([10 1000]); caxis([-1 +1]*0.15); colormap(llsmap5);
    text(-80,+200,'(a)');  ylabel('Pressure [mb]');
  tfov(2) = nexttile;
    pcolor(rlat,plays100,smoothn(moo2,1));
    colorbar; shading interp; set(gca,'ydir','reverse','yscale','log'); ylim([10 1000]); caxis([0 +1]/5*0.15); colormap(jet);
    text(-80,+200,'(b)'); ylabel('Pressure [mb]'); xlabel('Latitude [deg]');
colormap(tfov(1),llsmap5);
colormap(tfov(2),jet);
ta.Padding = 'compact';
ta.TileSpacing = 'compact';
    tfov(1).XTickLabel = '';  tfov(1).XLabel.String = [];
