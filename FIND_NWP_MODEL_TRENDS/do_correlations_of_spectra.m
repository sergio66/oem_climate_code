iCheckMath = -1;
if iCheckMath > 0
  disp('Checking obs-merra calcs')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('All 2645 chans')
  sah = (obsrates.rates - merra2.merra2_spectral_rates); 
  wah = cosYY.*sah;
  gah = reshape(wah,2645,72,64); gah = squeeze(nanmean(gah,2));
  
  figure(66); clf; plot(fchanx,nanmean(wah,2),'rx-',fchanx,nansum(wah,2)./nansum(cosYY,2),'m',fchanx,nanmean(gah,2),'k'); 
  plotaxis2; ylim([-1 +1]*0.03); 
  hl = legend('4608 incorrect mean','4608 correct mean','64 incorrect mean','location','best','fontsize',10);
  
  bah = nansum(wah,2)./nansum(cosYY,2);
  [m,s,m0,s0] = weighted_mean_stddev(sah,cosYY(iX,:));   %%% SAH IS UNWEIGHTED!!!! COSYY ARE THE WGTS
  fprintf(1,' mean(bah(:)) = %8.6f K      m s m0 s0 = %8.6f +/- %8.6f K, %8.6f +/- %8.6f K \n',mean(bah(:)),m,s,m0,s0);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('MW/MW chans only')
  iX = find(fchanx >=1605,1);
  iX = 1 : iX;
  sah = (obsrates.rates - merra2.merra2_spectral_rates); 
  wah = cosYY.*(obsrates.rates - merra2.merra2_spectral_rates); 
  gah = reshape(wah,2645,72,64); gah = squeeze(nanmean(gah,2)); 
  wah = wah(iX,:); gah = gah(iX,:); sah = sah(iX,:);
  figure(67); clf; plot(fchanx(iX),nanmean(wah,2),'rx-',fchanx(iX),nansum(wah,2)./nansum(cosYY(iX,:),2),'m',fchanx(iX),nanmean(gah,2),'k'); 
  plotaxis2; ylim([-1 +1]*0.03); 
  hl = legend('4608 incorrect mean','4608 correct mean','64 incorrect mean','location','best','fontsize',10);
  
  bah = nansum(wah,2)./nansum(cosYY(iX,:),2);
  [m,s,m0,s0] = weighted_mean_stddev(sah,cosYY(iX,:));   %%% SAH IS UNWEIGHTED!!!! COSYY ARE THE WGTS
  fprintf(1,' mean(bah(:)) = %8.6f K      m s m0 s0 = %8.6f +/- %8.6f K, %8.6f +/- %8.6f K \n',mean(bah(:)),m,s,m0,s0);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% [R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall,wall,varnames);
%%% [R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall,wall,varnames);
%%% [R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall,wall,varnames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varnames = {'OBS','UMBC','ERA5','MERRA2','AIRS L3','CLIMCAPS'};
varnames = {'OBS','This work','ERA5','MERRA2','AIRS L3','CLIMCAPS'};
varnames = {'OBS','AIRS\_RT','ERA5','MERRA2','AIRS L3','CLIMCAPS'};

plotoptions6.str11 = 'AIRS L1C';   plotoptions6.str12 = 'AIRS\_RT';
plotoptions6.str21 = 'ERA52';      plotoptions6.str22 = 'MERRA2';
plotoptions6.str31 = 'AIRS L3';   plotoptions6.str32 = 'CLIMCAPS';
plotoptions6.ystr = 'Latitude';   plotoptions6.xstr = 'Wavenumber / [cm-1]';
plotoptions6.maintitle = 'dBT/dt / [K/yr]';
plotoptions6.cx = [-1 +1]*0.25; plotoptions6.cmap = llsmap5; plotoptions6.yReverseDir = -1; plotoptions6.yLinearOrLog = +1;

YYY = nanmean(reshape(YY,72,64),1);

dtropical = find(abs(YYY) <= 30);
dmidlat   = find(abs(YYY) > 30 & abs(YYY) <= 60);
dpolar    = find(abs(YYY) > 60);

dnottropical = setdiff(1:64,dtropical);
dnotmidlat   = setdiff(1:64,dmidlat);
dnotpolar    = setdiff(1:64,dpolar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('R = 6 x 4 matrix of R^2 with the  rows being :   64-1640cm-1    640-800 cm-1   800-960cm-1   1370-1620cm-1')
disp('and the columns being')
disp('OBS')
disp('AIRS\_RT')
disp('ERA5')
disp('MERRA2')    
disp('AIRS L3 v7')
disp('CLIMCAPS L3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iX = find(fchanx >= 1608,1);
ijunk = length(1:iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(1:iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(1:iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(1:iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z11 = z11(dtropical,:);
z12 = z12(dtropical,:);
z21 = z21(dtropical,:);
z22 = z22(dtropical,:);
z31 = z31(dtropical,:);
z32 = z32(dtropical,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(1:iX),YYY(dtropical));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
% addpath /home/sergio/MATLABCODE/NANROUTINES/
% R(1) = 1.0;
% for ii = 2 : 6
%   [R(ii),chisqr] = nanlinearcorrelation(zall(1,:),zall(ii,:));
% end
% R is the same as the first colum or row of these : 
%corrcoef(zall')
%corr(zall')
figure(iOffSet+27); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rtropical(:,1) = R(:,1); Mtropical(:,1) = m'; Stropical(:,1) = s';
title('Correlation Matrix : ALL')
hfig = gcf;
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 790,1);
ijunk = length(1:iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(1:iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(1:iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(1:iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z11 = z11(dtropical,:);
z12 = z12(dtropical,:);
z21 = z21(dtropical,:);
z22 = z22(dtropical,:);
z31 = z31(dtropical,:);
z32 = z32(dtropical,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(1:iX),YYY(dtropical));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+28); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rtropical(:,2) = R(:,1); Mtropical(:,2) = m'; Stropical(:,2) = s';
title('Correlation Matrix : 15 um')
hfig = gcf;
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 800 & fchanx <= 962);
ijunk = length(iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z11 = z11(dtropical,:);
z12 = z12(dtropical,:);
z21 = z21(dtropical,:);
z22 = z22(dtropical,:);
z31 = z31(dtropical,:);
z32 = z32(dtropical,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(iX),YYY(dtropical));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+29); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rtropical(:,3) = R(:,1); Mtropical(:,3) = m'; Stropical(:,3) = s';
hfig = gcf;
title('Correlation Matrix : Window Region 10-12um')
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 1370 & fchanx <= 1620);
ijunk = length(iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z11 = z11(dtropical,:);
z12 = z12(dtropical,:);
z21 = z21(dtropical,:);
z22 = z22(dtropical,:);
z31 = z31(dtropical,:);
z32 = z32(dtropical,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(iX),YYY(dtropical));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+30); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rtropical(:,4) = R(:,1); Mtropical(:,4) = m'; Stropical(:,4) = s';
hfig = gcf;
title('Correlation Matrix : WV')
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

Rtropical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iX = find(fchanx >= 1608,1);
ijunk = length(1:iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(1:iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(1:iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(1:iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z11 = z11(dmidlat,:);
z12 = z12(dmidlat,:);
z21 = z21(dmidlat,:);
z22 = z22(dmidlat,:);
z31 = z31(dmidlat,:);
z32 = z32(dmidlat,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(1:iX),YYY(dmidlat));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+27); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rmidlat(:,1) = R(:,1); Mmidlat(:,1) = m'; Smidlat(:,1) = s';
title('Correlation Matrix : ALL')
hfig = gcf;
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 790,1);
ijunk = length(1:iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(1:iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(1:iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(1:iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z11 = z11(dmidlat,:);
z12 = z12(dmidlat,:);
z21 = z21(dmidlat,:);
z22 = z22(dmidlat,:);
z31 = z31(dmidlat,:);
z32 = z32(dmidlat,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(1:iX),YYY(dmidlat));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+28); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rmidlat(:,2) = R(:,1); Mmidlat(:,2) = m'; Smidlat(:,2) = s';
title('Correlation Matrix : 15 um')
hfig = gcf;
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 800 & fchanx <= 962);
ijunk = length(iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z11 = z11(dmidlat,:);
z12 = z12(dmidlat,:);
z21 = z21(dmidlat,:);
z22 = z22(dmidlat,:);
z31 = z31(dmidlat,:);
z32 = z32(dmidlat,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(iX),YYY(dmidlat));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+29); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rmidlat(:,3) = R(:,1); Mmidlat(:,3) = m'; Smidlat(:,3) = s';
hfig = gcf;
title('Correlation Matrix : Window Region 10-12um')
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 1370 & fchanx <= 1620);
ijunk = length(iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z11 = z11(dmidlat,:);
z12 = z12(dmidlat,:);
z21 = z21(dmidlat,:);
z22 = z22(dmidlat,:);
z31 = z31(dmidlat,:);
z32 = z32(dmidlat,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(iX),YYY(dmidlat));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+30); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rmidlat(:,4) = R(:,1); Mmidlat(:,4) = m'; Smidlat(:,4) = s';
hfig = gcf;
title('Correlation Matrix : WV')
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

Rmidlat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iX = find(fchanx >= 1608,1);
ijunk = length(1:iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(1:iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(1:iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(1:iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z11 = z11(dpolar,:);
z12 = z12(dpolar,:);
z21 = z21(dpolar,:);
z22 = z22(dpolar,:);
z31 = z31(dpolar,:);
z32 = z32(dpolar,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(1:iX),YYY(dpolar));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+27); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rpolar(:,1) = R(:,1); Mpolar(:,1) = m'; Spolar(:,1) = s';
title('Correlation Matrix : ALL')
hfig = gcf;
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 790,1);
ijunk = length(1:iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(1:iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(1:iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(1:iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z11 = z11(dpolar,:);
z12 = z12(dpolar,:);
z21 = z21(dpolar,:);
z22 = z22(dpolar,:);
z31 = z31(dpolar,:);
z32 = z32(dpolar,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(1:iX),YYY(dpolar));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+28); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rpolar(:,2) = R(:,1); Mpolar(:,2) = m'; Spolar(:,2) = s';
title('Correlation Matrix : 15 um')
hfig = gcf;
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 800 & fchanx <= 962);
ijunk = length(iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z11 = z11(dpolar,:);
z12 = z12(dpolar,:);
z21 = z21(dpolar,:);
z22 = z22(dpolar,:);
z31 = z31(dpolar,:);
z32 = z32(dpolar,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(iX),YYY(dpolar));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+29); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rpolar(:,3) = R(:,1); Mpolar(:,3) = m'; Spolar(:,3) = s';
hfig = gcf;
title('Correlation Matrix : Window Region 10-12um')
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 1370 & fchanx <= 1620);
ijunk = length(iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z11 = z11(dpolar,:);
z12 = z12(dpolar,:);
z21 = z21(dpolar,:);
z22 = z22(dpolar,:);
z31 = z31(dpolar,:);
z32 = z32(dpolar,:);
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(iX),YYY(dpolar));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+30); clf
wall = ones(size(z11(:)));
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rpolar(:,4) = R(:,1); Mpolar(:,4) = m'; Spolar(:,4) = s';
hfig = gcf;
title('Correlation Matrix : WV')
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

Rpolar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iX = find(fchanx >= 1608,1);
ijunk = length(1:iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(1:iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(1:iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(1:iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(1:iX),nanmean(reshape(YY,72,64),1));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+27); clf
wall = ones(size(z11(:)));
wall = cos(YYY'*pi/180)*ones(1,ijunk); wall = wall(:);
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rglobal(:,1) = R(:,1); Mglobal(:,1) = m'; Sglobal(:,1) = s';
title('Correlation Matrix : ALL')
hfig = gcf;
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 790,1);
ijunk = length(1:iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(1:iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(1:iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(1:iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(1:iX,:)',72,64,ijunk),1));
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(1:iX),nanmean(reshape(YY,72,64),1));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+28); clf
wall = ones(size(z11(:)));
wall = cos(YYY'*pi/180)*ones(1,ijunk); wall = wall(:);
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rglobal(:,2) = R(:,1); Mglobal(:,2) = m'; Sglobal(:,2) = s';
title('Correlation Matrix : 15 um')
hfig = gcf;
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 800 & fchanx <= 962);
ijunk = length(iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(iX,:)',72,64,ijunk),1));
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(iX),nanmean(reshape(YY,72,64),1));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+29); clf
wall = ones(size(z11(:)));
wall = cos(YYY'*pi/180)*ones(1,ijunk); wall = wall(:);
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rglobal(:,3) = R(:,1); Mglobal(:,3) = m'; Sglobal(:,3) = s';
hfig = gcf;
title('Correlation Matrix : Window Region 10-12um')
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

iX = find(fchanx >= 1370 & fchanx <= 1620);
ijunk = length(iX);
z11 = squeeze(nanmean(reshape(obsrates.rates(iX,:)',72,64,ijunk),1));
z12 = squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates(iX,:)',72,64,ijunk),1));
z21 = squeeze(nanmean(reshape(era5.era5_spectral_rates(iX,:)',72,64,ijunk),1));
z22 = squeeze(nanmean(reshape(merra2.merra2_spectral_rates(iX,:)',72,64,ijunk),1));
z31 = squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates(iX,:)',72,64,ijunk),1));
z32 = squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates(iX,:)',72,64,ijunk),1));
  tiled_3x2layout(z11,z12,z21,z22,z31,z32,iOffSet+26,plotoptions6,fchanx(iX),nanmean(reshape(YY,72,64),1));
zall = [z11(:) z12(:) z21(:) z22(:) z31(:) z32(:)]';
%corrcoef(zall')
%corr(zall')
figure(iOffSet+30); clf
wall = ones(size(z11(:)));
wall = cos(YYY'*pi/180)*ones(1,ijunk); wall = wall(:);
[R,Pvalue,m,s,m0,s0] = corrplot_weighted_mean_stddev(zall',wall',varnames);
Rglobal(:,4) = R(:,1); Mglobal(:,4) = m'; Sglobal(:,4) = s';
hfig = gcf;
title('Correlation Matrix : WV')
%haxes = findobj(hfig, 'Type', 'Axes');
%arrayfun(@(ax) xlim(ax, [-1 +1]*0.25), haxes);

Rglobal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Rx = 6 x 4 matrix of R^2 with the  rows being :   GLOBAL  TROPICAL MIDLAT POLAR')

disp('and the columns being')
disp('OBS')
disp('AIRS\_RT')
disp('ERA5')
disp('MERRA2')    
disp('AIRS L3 v7')
disp('CLIMCAPS L3')

disp('first are regression results : Rglobal Rtropical Rmidlat Rpolar')
Rallchans  = [Rglobal(:,1) Rtropical(:,1) Rmidlat(:,1) Rpolar(:,1)]
R15umchans = [Rglobal(:,2) Rtropical(:,2) Rmidlat(:,2) Rpolar(:,2)]
Rwinchans  = [Rglobal(:,3) Rtropical(:,3) Rmidlat(:,3) Rpolar(:,3)]
RWVchans   = [Rglobal(:,4) Rtropical(:,4) Rmidlat(:,4) Rpolar(:,4)]

%% put negative sign because ../FIND_NWP_MODEL_TRENDS/plot_spectral_get_the_model_trends2.m does (blah(:,ii)-blah(:,1))
%% and here blah(:,1) ==== obs  
disp('next are mean O-C bias : Mglobal Mtropical Mmidlat Mpolar')
Mallchans  = -[Mglobal(:,1) Mtropical(:,1) Mmidlat(:,1) Mpolar(:,1)]
M15umchans = -[Mglobal(:,2) Mtropical(:,2) Mmidlat(:,2) Mpolar(:,2)]
Mwinchans  = -[Mglobal(:,3) Mtropical(:,3) Mmidlat(:,3) Mpolar(:,3)]
MWVchans   = -[Mglobal(:,4) Mtropical(:,4) Mmidlat(:,4) Mpolar(:,4)]

disp('then are O-C stddev : Sglobal Stropical Smidlat Spolar')
Sallchans  = [Sglobal(:,1) Stropical(:,1) Smidlat(:,1) Spolar(:,1)]
S15umchans = [Sglobal(:,2) Stropical(:,2) Smidlat(:,2) Spolar(:,2)]
Swinchans  = [Sglobal(:,3) Stropical(:,3) Smidlat(:,3) Spolar(:,3)]
SWVchans   = [Sglobal(:,4) Stropical(:,4) Smidlat(:,4) Spolar(:,4)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('now go run save_matrix_plot_correlations_of_spectraV2D/E')
disp('now go run save_matrix_plot_correlations_of_spectraV2D/E')
disp('now go run save_matrix_plot_correlations_of_spectraV2D/E')

%{
comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/do_correlations_of_spectra.m';
save /home/sergio/PAPERS/SUBMITPAPERS/trends/Matlab/obsVScal_stats_R_mean_std.mat R*chans M*chans S*chans

[Rallchans(2:5,:) R15umchans(2:5,:) Rwinchans(2:5,:) RWVchans(2:5,:)] 
[Rallchans(2:5,:); R15umchans(2:5,:); Rwinchans(2:5,:); RWVchans(2:5,:)] 

save_matrix_plot_correlations_of_spectra
save_matrix_plot_correlations_of_spectraV2

save_matrix_plot_correlations_of_spectraV2A
save_matrix_plot_correlations_of_spectraV2B
save_matrix_plot_correlations_of_spectraV2C
save_matrix_plot_correlations_of_spectraV2D
save_matrix_plot_correlations_of_spectraV2E
%} 

save_matrix_plot_correlations_of_spectraV2E
%{
addpath /asl/matlib/plotutils
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
figure(40); aslprint([dir0 'land_ocean_spectra_all_tropics_midlats_polar.pdf']);
figure(41); aslprint([dir0 'ocean_spectra_all_tropics_midlats_polar.pdf']);
figure(42); aslprint([dir0 'land_spectra_all_tropics_midlats_polar.pdf']);
figure(70); aslprint([dir0 'land_ocean_spectra_matrix_correlation.pdf']);
figure(71); aslprint([dir0 'land_ocean_spectra_matrix_bias.pdf']);
figure(72); aslprint([dir0 'land_ocean_spectra_matrix_stddev.pdf']);
%}

