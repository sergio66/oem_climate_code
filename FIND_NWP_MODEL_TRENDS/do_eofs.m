[Lera5skt, EOFsera5skt, ECera5skt, errorera5skt] = detrend_time_series_make_EOF_v1(era5skt');
[Lera5mmw, EOFsera5mmw, ECera5mmw, errorera5mmw] = detrend_time_series_make_EOF_v1(era5mmw');

anom_1231 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/anomaly_chID_1520_Q03.mat');
  [L1231, EOFs1231, EC1231, error1231] = detrend_time_series_make_EOF_v1(anom_1231.btanom');
anom_1227 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/anomaly_chID_1511_Q03.mat');
  [L1227, EOFs1227, EC1227, error1227] = detrend_time_series_make_EOF_v1(anom_1227.btanom');
anom_1419 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/anomaly_chID_1861_Q03.mat');
  [L1419, EOFs1419, EC1419, error1419] = detrend_time_series_make_EOF_v1(anom_1419.btanom');
anom_1521 = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/anomaly_chID_2025_Q03.mat');
  [L1521, EOFs1521, EC1521, error1521] = detrend_time_series_make_EOF_v1(anom_1521.btanom');

EOFs1 = EOFs1231;
EOFs1 = EOFs1231 - EOFs1227;
EOFs1 = EOFs1419;
EOFs1 = EOFs1521;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(:,iEOF),era5.stemprate');           thecorrEOF.era5skt_ST(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(:,iEOF),merra2.stemprate');         thecorrEOF.era5skt_ST(iEOF,2) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(:,iEOF),airsL3.stemprate');         thecorrEOF.era5skt_ST(iEOF,3) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(:,iEOF),climcapsL3.stemprate');     thecorrEOF.era5skt_ST(iEOF,4) = r;
  bonk = agiss.giss_trend4608; bonk = bonk(:);
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(iEOF,:),bonk');                     thecorrEOF.era5skt_ST(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(:,iEOF),umbc.stemprate');           thecorrEOF.era5skt_ST(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(:,iEOF),era5.stemprate');           thecorrEOF.BT1231_ST(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(:,iEOF),merra2.stemprate');         thecorrEOF.BT1231_ST(iEOF,2) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(:,iEOF),airsL3.stemprate');         thecorrEOF.BT1231_ST(iEOF,3) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(:,iEOF),climcapsL3.stemprate');     thecorrEOF.BT1231_ST(iEOF,4) = r;
  bonk = agiss.giss_trend4608; bonk = bonk(:);
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(iEOF,:),bonk');                     thecorrEOF.BT1231_ST(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(:,iEOF),umbc.stemprate');           thecorrEOF.BT1231_ST(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(:,iEOF),olr.trend_mmw');            thecorrEOF.era5mmw_MMW(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(:,iEOF),merra2_mmwtrend');          thecorrEOF.era5mmw_MMW(iEOF,2) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(:,iEOF),airsL3.stemprate');        thecorrEOF.era5mmw_MMW(iEOF,3) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(:,iEOF),climcapsL3.stemprate');    thecorrEOF.era5mmw_MMW(iEOF,4) = r;
  %bonk = agiss.giss_trend4608; bonk = bonk(:);
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(iEOF,:),bonk');                    thecorrEOF.era5mmw_MMW(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(:,iEOF),umbc_mmwrate');             thecorrEOF.era5mmw_MMW(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(:,iEOF),olr.trend_mmw');            thecorrEOF.BT1519_MMW(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(:,iEOF),merra2_mmwtrend');          thecorrEOF.BT1519_MMW(iEOF,2) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(:,iEOF),airsL3.stemprate');        thecorrEOF.BT1519_MMW(iEOF,3) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(:,iEOF),climcapsL3.stemprate');    thecorrEOF.BT1519_MMW(iEOF,4) = r;
  %bonk = agiss.giss_trend4608; bonk = bonk(:);
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(iEOF,:),bonk');                    thecorrEOF.BT1519_MMW(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(:,iEOF),umbc_mmwrate');             thecorrEOF.BT1519_MMW(iEOF,6) = r;
end

disp('GLOBAL : ')
fprintf(1,'correlations of [EOF1 of BT1231]   with SKT trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.BT1231_ST(1,:))
fprintf(1,'correlations of [EOF1 of ERA5 SKT] with SKT trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.era5skt_ST(1,:))
fprintf(1,'correlations of [EOF1 of BT1519]   with MMW trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.BT1519_MMW(1,:))
fprintf(1,'correlations of [EOF1 of ERA5 MMW] with MMW trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.era5mmw_MMW(1,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(tropics),iEOF),era5.stemprate(isfinite(tropics))');           thecorrEOF.tropical_era5skt_ST(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(tropics),iEOF),merra2.stemprate(isfinite(tropics))');         thecorrEOF.tropical_era5skt_ST(iEOF,2) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(tropics),iEOF),airsL3.stemprate(isfinite(tropics))');         thecorrEOF.tropical_era5skt_ST(iEOF,3) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(tropics),iEOF),climcapsL3.stemprate(isfinite(tropics))');     thecorrEOF.tropical_era5skt_ST(iEOF,4) = r;
  bonk = agiss.giss_trend4608; bonk = bonk(:)';
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(tropics),iEOF),bonk(isfinite(tropics))');                     thecorrEOF.tropical_era5skt_ST(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(tropics),iEOF),umbc.stemprate(isfinite(tropics))');           thecorrEOF.tropical_era5skt_ST(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(tropics),iEOF),era5.stemprate(isfinite(tropics))');           thecorrEOF.tropical_BT1231_ST(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(tropics),iEOF),merra2.stemprate(isfinite(tropics))');         thecorrEOF.tropical_BT1231_ST(iEOF,2) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(tropics),iEOF),airsL3.stemprate(isfinite(tropics))');         thecorrEOF.tropical_BT1231_ST(iEOF,3) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(tropics),iEOF),climcapsL3.stemprate(isfinite(tropics))');     thecorrEOF.tropical_BT1231_ST(iEOF,4) = r;
  bonk = agiss.giss_trend4608; bonk = bonk(:)';
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(tropics),iEOF),bonk(isfinite(tropics))');                     thecorrEOF.tropical_BT1231_ST(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(tropics),iEOF),umbc.stemprate(isfinite(tropics))');           thecorrEOF.tropical_BT1231_ST(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(tropics),iEOF),olr.trend_mmw(isfinite(tropics))');            thecorrEOF.tropical_era5mmw_MMW(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(tropics),iEOF),merra2_mmwtrend(isfinite(tropics))');          thecorrEOF.tropical_era5mmw_MMW(iEOF,2) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(tropics),iEOF),airsL3.stemprate(isfinite(tropics))');        thecorrEOF.tropical_era5mmw_MMW(iEOF,3) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(tropics),iEOF),climcapsL3.stemprate(isfinite(tropics))');    thecorrEOF.tropical_era5mmw_MMW(iEOF,4) = r;
  %bonk = agiss.giss_trend4608; bonk = bonk(:)';
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(tropics),iEOF),bonk(isfinite(tropics))');                    thecorrEOF.tropical_era5mmw_MMW(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(tropics),iEOF),umbc_mmwrate(isfinite(tropics))');             thecorrEOF.tropical_era5mmw_MMW(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(tropics),iEOF),olr.trend_mmw(isfinite(tropics))');            thecorrEOF.tropical_BT1519_MMW(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(tropics),iEOF),merra2_mmwtrend(isfinite(tropics))');          thecorrEOF.tropical_BT1519_MMW(iEOF,2) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(tropics),iEOF),airsL3.stemprate(isfinite(tropics))');        thecorrEOF.tropical_BT1519_MMW(iEOF,3) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(tropics),iEOF),climcapsL3.stemprate(isfinite(tropics))');    thecorrEOF.tropical_BT1519_MMW(iEOF,4) = r;
  %bonk = agiss.giss_trend4608; bonk = bonk(:)';
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(tropics),iEOF),bonk(isfinite(tropics))');                    thecorrEOF.tropical_BT1519_MMW(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(tropics),iEOF),umbc_mmwrate(isfinite(tropics))');             thecorrEOF.tropical_BT1519_MMW(iEOF,6) = r;
end

disp('TROPICAL : ')
fprintf(1,'correlations of [EOF1 of BT1231]   with SKT trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.tropical_BT1231_ST(1,:))
fprintf(1,'correlations of [EOF1 of ERA5 SKT] with SKT trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.tropical_era5skt_ST(1,:))
fprintf(1,'correlations of [EOF1 of BT1519]   with MMW trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.tropical_BT1519_MMW(1,:))
fprintf(1,'correlations of [EOF1 of ERA5 MMW] with MMW trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.tropical_era5mmw_MMW(1,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(oceantropics),iEOF),era5.stemprate(isfinite(oceantropics))');           thecorrEOF.oceantropical_era5skt_ST(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(oceantropics),iEOF),merra2.stemprate(isfinite(oceantropics))');         thecorrEOF.oceantropical_era5skt_ST(iEOF,2) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(oceantropics),iEOF),airsL3.stemprate(isfinite(oceantropics))');         thecorrEOF.oceantropical_era5skt_ST(iEOF,3) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(oceantropics),iEOF),climcapsL3.stemprate(isfinite(oceantropics))');     thecorrEOF.oceantropical_era5skt_ST(iEOF,4) = r;
  bonk = agiss.giss_trend4608; bonk = bonk(:)';
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(oceantropics),iEOF),bonk(isfinite(oceantropics))');                     thecorrEOF.oceantropical_era5skt_ST(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(isfinite(oceantropics),iEOF),umbc.stemprate(isfinite(oceantropics))');           thecorrEOF.oceantropical_era5skt_ST(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(oceantropics),iEOF),era5.stemprate(isfinite(oceantropics))');           thecorrEOF.oceantropical_BT1231_ST(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(oceantropics),iEOF),merra2.stemprate(isfinite(oceantropics))');         thecorrEOF.oceantropical_BT1231_ST(iEOF,2) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(oceantropics),iEOF),airsL3.stemprate(isfinite(oceantropics))');         thecorrEOF.oceantropical_BT1231_ST(iEOF,3) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(oceantropics),iEOF),climcapsL3.stemprate(isfinite(oceantropics))');     thecorrEOF.oceantropical_BT1231_ST(iEOF,4) = r;
  bonk = agiss.giss_trend4608; bonk = bonk(:)';
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(oceantropics),iEOF),bonk(isfinite(oceantropics))');                     thecorrEOF.oceantropical_BT1231_ST(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(isfinite(oceantropics),iEOF),umbc.stemprate(isfinite(oceantropics))');           thecorrEOF.oceantropical_BT1231_ST(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(oceantropics),iEOF),olr.trend_mmw(isfinite(oceantropics))');            thecorrEOF.oceantropical_era5mmw_MMW(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(oceantropics),iEOF),merra2_mmwtrend(isfinite(oceantropics))');          thecorrEOF.oceantropical_era5mmw_MMW(iEOF,2) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(oceantropics),iEOF),airsL3.stemprate(isfinite(oceantropics))');        thecorrEOF.oceantropical_era5mmw_MMW(iEOF,3) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(oceantropics),iEOF),climcapsL3.stemprate(isfinite(oceantropics))');    thecorrEOF.oceantropical_era5mmw_MMW(iEOF,4) = r;
  %bonk = agiss.giss_trend4608; bonk = bonk(:)';
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(oceantropics),iEOF),bonk(isfinite(oceantropics))');                    thecorrEOF.oceantropical_era5mmw_MMW(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(isfinite(oceantropics),iEOF),umbc_mmwrate(isfinite(oceantropics))');             thecorrEOF.oceantropical_era5mmw_MMW(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(oceantropics),iEOF),olr.trend_mmw(isfinite(oceantropics))');            thecorrEOF.oceantropical_BT1519_MMW(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(oceantropics),iEOF),merra2_mmwtrend(isfinite(oceantropics))');          thecorrEOF.oceantropical_BT1519_MMW(iEOF,2) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(oceantropics),iEOF),airsL3.stemprate(isfinite(oceantropics))');        thecorrEOF.oceantropical_BT1519_MMW(iEOF,3) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(oceantropics),iEOF),climcapsL3.stemprate(isfinite(oceantropics))');    thecorrEOF.oceantropical_BT1519_MMW(iEOF,4) = r;
  %bonk = agiss.giss_trend4608; bonk = bonk(:)';
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(oceantropics),iEOF),bonk(isfinite(oceantropics))');                    thecorrEOF.oceantropical_BT1519_MMW(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(isfinite(oceantropics),iEOF),umbc_mmwrate(isfinite(oceantropics))');             thecorrEOF.oceantropical_BT1519_MMW(iEOF,6) = r;
end

disp('TROPICAL OCEAN : ')
fprintf(1,'correlations of [EOF1 of BT1231]   with SKT trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.oceantropical_BT1231_ST(1,:))
fprintf(1,'correlations of [EOF1 of ERA5 SKT] with SKT trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.oceantropical_era5skt_ST(1,:))
fprintf(1,'correlations of [EOF1 of BT1519]   with MMW trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.oceantropical_BT1519_MMW(1,:))
fprintf(1,'correlations of [EOF1 of ERA5 MMW] with MMW trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.oceantropical_era5mmw_MMW(1,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oocean = find(ocean == 1);

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(oocean,iEOF),era5.stemprate(oocean)');           thecorrEOF.ocean_era5skt_ST(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(oocean,iEOF),merra2.stemprate(oocean)');         thecorrEOF.ocean_era5skt_ST(iEOF,2) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(oocean,iEOF),airsL3.stemprate(oocean)');         thecorrEOF.ocean_era5skt_ST(iEOF,3) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(oocean,iEOF),climcapsL3.stemprate(oocean)');     thecorrEOF.ocean_era5skt_ST(iEOF,4) = r;
  bonk = agiss.giss_trend4608; bonk = bonk(:)';
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(oocean,iEOF),bonk(oocean)');                     thecorrEOF.ocean_era5skt_ST(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5skt(oocean,iEOF),umbc.stemprate(oocean)');           thecorrEOF.ocean_era5skt_ST(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(oocean,iEOF),era5.stemprate(oocean)');           thecorrEOF.ocean_BT1231_ST(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(oocean,iEOF),merra2.stemprate(oocean)');         thecorrEOF.ocean_BT1231_ST(iEOF,2) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(oocean,iEOF),airsL3.stemprate(oocean)');         thecorrEOF.ocean_BT1231_ST(iEOF,3) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(oocean,iEOF),climcapsL3.stemprate(oocean)');     thecorrEOF.ocean_BT1231_ST(iEOF,4) = r;
  bonk = agiss.giss_trend4608; bonk = bonk(:)';
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(oocean,iEOF),bonk(oocean)');                     thecorrEOF.ocean_BT1231_ST(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1231(oocean,iEOF),umbc.stemprate(oocean)');           thecorrEOF.ocean_BT1231_ST(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(oocean,iEOF),olr.trend_mmw(oocean)');            thecorrEOF.ocean_era5mmw_MMW(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(oocean,iEOF),merra2_mmwtrend(oocean)');          thecorrEOF.ocean_era5mmw_MMW(iEOF,2) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(oocean,iEOF),airsL3.stemprate(oocean)');        thecorrEOF.ocean_era5mmw_MMW(iEOF,3) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(oocean,iEOF),climcapsL3.stemprate(oocean)');    thecorrEOF.ocean_era5mmw_MMW(iEOF,4) = r;
  %bonk = agiss.giss_trend4608; bonk = bonk(:)';
  %[r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(oocean,iEOF),bonk(oocean)');                    thecorrEOF.ocean_era5mmw_MMW(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFsera5mmw(oocean,iEOF),umbc_mmwrate(oocean)');             thecorrEOF.ocean_era5mmw_MMW(iEOF,6) = r;
end

for iEOF = 1 : 10
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(oocean,iEOF),olr.trend_mmw(oocean)');            thecorrEOF.ocean_BT1519_MMW(iEOF,1) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(oocean,iEOF),merra2_mmwtrend(oocean)');          thecorrEOF.ocean_BT1519_MMW(iEOF,2) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(oocean,iEOF),airsL3.stemprate(oocean)');        thecorrEOF.ocean_BT1519_MMW(iEOF,3) = r;
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(oocean,iEOF),climcapsL3.stemprate(oocean)');    thecorrEOF.ocean_BT1519_MMW(iEOF,4) = r;
  %bonk = agiss.giss_trend4608; bonk = bonk(:)';
  %[r,chisqr,P] = nanlinearcorrelation(EOFs1(oocean,iEOF),bonk(oocean)');                    thecorrEOF.ocean_BT1519_MMW(iEOF,5) = r;
  [r,chisqr,P] = nanlinearcorrelation(EOFs1(oocean,iEOF),umbc_mmwrate(oocean)');             thecorrEOF.ocean_BT1519_MMW(iEOF,6) = r;
end

disp('OCEAN : ')
fprintf(1,'correlations of [EOF1 of BT1231]   with SKT trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.ocean_BT1231_ST(1,:))
fprintf(1,'correlations of [EOF1 of ERA5 SKT] with SKT trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.ocean_era5skt_ST(1,:))
fprintf(1,'correlations of [EOF1 of BT1519]   with MMW trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.ocean_BT1519_MMW(1,:))
fprintf(1,'correlations of [EOF1 of ERA5 MMW] with MMW trend = %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  \n',thecorrEOF.ocean_era5mmw_MMW(1,:))
