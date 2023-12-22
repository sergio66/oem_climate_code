ii = 0; 

ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fUMBC_day.RHtrend,100,72,64),2)),1)); title('dRH/dt : CHIRP\_A DAY');
ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fUMBC_night.RHtrend,100,72,64),2)),1)); title('dRH/dt : CHIRP\_A NIGHT');

ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fERA5_day.RHtrend,100,72,64),2)),1)); title('dRH/dt : ERA5 DAY');
ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fERA5_night.RHtrend,100,72,64),2)),1)); title('dRH/dt : ERA5 NIGHT');

ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fAIRSL3_day.RHtrend,100,72,64),2)),1)); title('dRH/dt : AIRSL3 DAY');
ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fAIRSL3_night.RHtrend,100,72,64),2)),1)); title('dRH/dt : AIRSL3 NIGHT');

ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fCLIMCAPSL3_day.RHtrend,100,72,64),2)),1)); title('dRH/dt : CLIMCAPSL3 DAY');
ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fCLIMCAPSL3_night.RHtrend,100,72,64),2)),1)); title('dRH/dt : CLIMCAPSL3 NIGHT');

ii = ii + 1; figure(ii); clf; pcolor(rlat,plays100,smoothn(squeeze(nanmean(reshape(fMERRA2.RHtrend,100,72,64),2)),1)); title('dRH/dt : MERRA2 DAY');

ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fUMBC_day.RHtrend + fUMBC_night.RHtrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/2); title('dRH/dt : CHIRP\_A DAY/NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fERA5_day.RHtrend + fERA5_night.RHtrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/2); title('dRH/dt : ERA5 DAY/NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fAIRSL3_day.RHtrend + fAIRSL3_night.RHtrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/2); title('dRH/dt : AIRSL3 DAY/NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fCLIMCAPSL3_day.RHtrend + fCLIMCAPSL3_night.RHtrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/2); title('dRH/dt : CLIMCAPSL3 DAY/NIGHT');

ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fUMBC_day.RHtrend - fUMBC_night.RHtrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/1); title('dRH/dt : CHIRP\_A DAY - NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fERA5_day.RHtrend - fERA5_night.RHtrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/1); title('dRH/dt : ERA5 DAY - NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fAIRSL3_day.RHtrend - fAIRSL3_night.RHtrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/1); title('dRH/dt : AIRSL3 DAY - NIGHT');
ii = ii + 1; figure(ii); clf; junk = smoothn(squeeze(nanmean(reshape(fCLIMCAPSL3_day.RHtrend - fCLIMCAPSL3_night.RHtrend,100,72,64),2)),1); pcolor(rlat,plays100,junk/1); title('dRH/dt : CLIMCAPSL3 DAY - NIGHT');

iiMax = ii;
for ii = 1 : iiMax; figure(ii); caxis([-1 +1]*0.50); colormap(llsmap5); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); colorbar; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newz*
iFig = 60;
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.501; plotoptions.maintitle = 'dRH/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = '(D) CHIRP\_A';     plotoptions.str21 = '(N)';
  plotoptions.str12 = '(D) AIRS L3';      plotoptions.str22 = '(N)';
  plotoptions.str13 = '(D) CLIMCAPS L3';  plotoptions.str23 = '(N)';
  plotoptions.str14 = '(D) ERA5';         plotoptions.str24 = '(N)';

  plotoptions.yLinearOrLog = +1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits      = [100 1000];

  newz11 = fUMBC_day.RHtrend;        newz21 = fUMBC_night.RHtrend;                           
  newz12 = fAIRSL3_day.RHtrend;      newz22 = fAIRSL3_night.RHtrend; 
  newz13 = fCLIMCAPSL3_day.RHtrend;  newz23 = fCLIMCAPSL3_night.RHtrend;
  newz14 = fERA5_day.RHtrend;        newz24 = fERA5_night.RHtrend; 
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2)); 
    newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2)); 
    newz13 = squeeze(nanmean(reshape(newz13,100,72,64),2));     newz23 = squeeze(nanmean(reshape(newz23,100,72,64),2)); 
    newz14 = squeeze(nanmean(reshape(newz14,100,72,64),2));     newz24 = squeeze(nanmean(reshape(newz24,100,72,64),2)); 
  profile_plots_2x4tiledlayout(rlat,plays100,newz11,newz12,newz13,newz14,newz21,newz22,newz23,newz24,iFig,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear newz*
iFig = 61; 
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.501; plotoptions.maintitle = 'dRH/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'CHIRP\_A';     plotoptions.str12 = 'ERA5';        
  plotoptions.str21 = 'AIRS L3';      plotoptions.str22 = 'CLIMCAPS L3'; 
  plotoptions.yLinearOrLog = +1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits      = [100 1000];
  newz11 = (fUMBC_day.RHtrend-fUMBC_night.RHtrend)*1.0;        newz12 = (fERA5_day.RHtrend-fERA5_night.RHtrend)*1.0; 
  newz21 = (fAIRSL3_day.RHtrend-fAIRSL3_night.RHtrend)*1.0;    newz22 = (fCLIMCAPSL3_day.RHtrend-fCLIMCAPSL3_night.RHtrend)*1.0;
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2)); 
    newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2)); 
  profile_plots_2x2tiledlayout(rlat,plays100,newz11,newz12,newz21,newz22,iFig,plotoptions);

clear newz*
iFig = 62;
  figure(iFig); clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.501; plotoptions.maintitle = 'dRH/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'CHIRP\_A';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'N/A';          plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.yLinearOrLog = +1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits      = [100 1000];
  newz11 = (fUMBC_day.RHtrend+fUMBC_night.RHtrend)*0.5;  newz12 = (fAIRSL3_day.RHtrend+fAIRSL3_night.RHtrend)*0.5; newz13 = (fCLIMCAPSL3_day.RHtrend+fCLIMCAPSL3_night.RHtrend)*0.5;
  newz21 = newz11*0;                                           newz22 = (fERA5_day.RHtrend+fERA5_night.RHtrend)*0.5;     newz23 = fMERRA2.RHtrend;
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2));     newz13 = squeeze(nanmean(reshape(newz13,100,72,64),2)); 
    newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2));     newz23 = squeeze(nanmean(reshape(newz23,100,72,64),2)); 
  profile_plots_2x3tiledlayout(rlat,plays100,newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);

clear newz*
  figure(iFig); clf;
  newz11 = (fUMBC_day.RHtrend+fUMBC_night.RHtrend)*0.5;  newz12 = (fAIRSL3_day.RHtrend+fAIRSL3_night.RHtrend)*0.5; newz13 = (fCLIMCAPSL3_day.RHtrend+fCLIMCAPSL3_night.RHtrend)*0.5;
                                       newz21 = fMERRA2.RHtrend;  newz22 = (fERA5_day.RHtrend+fERA5_night.RHtrend)*0.5; 
    newz11 = squeeze(nanmean(reshape(newz11,100,72,64),2));     newz12 = squeeze(nanmean(reshape(newz12,100,72,64),2));     newz13 = squeeze(nanmean(reshape(newz13,100,72,64),2)); 
              newz21 = squeeze(nanmean(reshape(newz21,100,72,64),2));     newz22 = squeeze(nanmean(reshape(newz22,100,72,64),2));     
    plotoptions2x1x2 = plotoptions;
    plotoptions2x1x2.str11 = 'CHIRP\_A';
    plotoptions2x1x2.str12 = 'AIRS L3';
    plotoptions2x1x2.str13 = 'CLIMCAPS';
    plotoptions2x1x2.str21 = 'MERRA2';
    plotoptions2x1x2.str22 = 'ERA5';
    plotoptions2x1x2.cstr  = 'dWVfrac/dt';
      profile_plots_2x1x2tiledlayout_tall(rlat,plays100,newz11,newz12,newz13,newz21,newz22,iFig,plotoptions2x1x2);

%junk = [newz11(:) newz12(:) newz13(:) newz21(:) newz22(:) newz23(:)];
%junk = reshape(junk,100,64,6);
junk = [newz11(:) newz12(:) newz13(:) newz21(:) newz22(:)];
junk = reshape(junk,100,64,5);
moo0 = nanmean(junk,3);
moo1 = abs(nanmean(junk,3));
moo2 = nanstd(junk,[],3);
moo3 = zeros(size(moo1)); moo3(moo1 > moo2) = 1;
moo4 = moo1./moo2; moo4(moo4 < 1) = 1; moo4 = log10(moo4);
 
iFig = 63; figure(iFig); clf; pcolor(rlat,plays100,smoothn(nanmean(junk,3),1));   cx = [-1 +1]*0.501; caxis(cx); colormap(llsmap5);
iFig = 64; figure(iFig); clf; pcolor(rlat,plays100,smoothn(nanstd(junk,[],3),1)); cx = [0 +1]*0.301;  caxis(cx); colormap(jet);

iFig = 65; figure(iFig); clf; pcolor(rlat,plays100,smoothn(moo3,1));         cx = [0 +1];  caxis(cx); colormap(jet);
iFig = 65; figure(iFig); clf; pcolor(rlat,plays100,moo3);                    cx = [0 +1]; caxis(cx); colormap(usa2x);

iFig = 66; figure(iFig); clf; pcolor(rlat,plays100,smoothn(moo4,1));         cx = [0 +2];       caxis(cx); colormap(usa2x);
iFig = 66; figure(iFig); clf; pcolor(rlat,plays100,10.^moo4);                cx = (10.^[0 +1]); caxis(cx); colormap(jett);

mooUMBC_num   = 0.5*(fUMBC_day.RHtrend + fUMBC_night.RHtrend);                 mooUMBC_num   = squeeze(nanmean(reshape(mooUMBC_num,100,72,64),2)); 
mooUMBC_denom = sqrt(fUMBC_day.RHtrendunc.^2 + fUMBC_night.RHtrendunc.^2);     mooUMBC_denom = squeeze(nanmean(reshape(mooUMBC_denom,100,72,64),2)); 
%mooUMBC_denom = mooUMBC_num/10;
mooUMBC = abs(mooUMBC_num)./mooUMBC_denom;
iFig = 67; figure(iFig); clf; pcolor(rlat,plays100,smoothn(mooUMBC_num,1));   cx = ([-1 +1]*0.501);  caxis(cx); colormap(llsmap5);
iFig = 68; figure(iFig); clf; pcolor(rlat,plays100,smoothn(mooUMBC_denom,1)); cx = ([0 +1]*0.501/2/10); caxis(cx); colormap(jet);
iFig = 69; figure(iFig); clf; pcolor(rlat,plays100,10.^(mooUMBC));            cx = (10.^[0 +1]);        caxis(cx); colormap(jett);
iFig = 69; figure(iFig); clf; pcolor(rlat,plays100,mooUMBC);                  cx = ([0 +1]);            caxis(cx); colormap(jett);

for ii = 63 : 69; figure(ii); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); colorbar; shading interp; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(95); clf
  ta = tiledlayout(2,1,'TileSpacing','None', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  tfov(1) = nexttile;
    pcolor(rlat,plays100,smoothn(moo0,1));
    colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([-1 +1]*0.501); colormap(llsmap5);
    text(-80,+200,'(a)');  ylabel('Pressure [mb]');
  tfov(2) = nexttile;
    pcolor(rlat,plays100,smoothn(moo2,1));
    colorbar; shading interp; set(gca,'ydir','reverse','yscale','linear'); ylim([100 1000]); caxis([0 +1]*0.31); colormap(jet);
    text(-80,+200,'(b)'); ylabel('Pressure [mb]'); xlabel('Latitude [deg]');
colormap(tfov(1),llsmap5);
colormap(tfov(2),jet);
ta.Padding = 'compact';
ta.TileSpacing = 'compact';
    tfov(1).XTickLabel = '';  tfov(1).XLabel.String = [];

