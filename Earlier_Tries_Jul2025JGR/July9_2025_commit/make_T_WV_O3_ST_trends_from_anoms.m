%%see ~/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/tile_fits_quantiles_anomalies.m
fdirpre      = '../DATAObsStats_StartSept2002_CORRECT_LatLon/';   %% symbolic link to ./DATAObsStats_StartSept2002_CORRECT_LatLon -> /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_CORRECT_LatLon
fdirpre      = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/';
lati = 32; loni = 36; i16daysSteps = 502;
fn_summary = sprintf('LatBin%1$02d/LonBin%2$02d/iQAX_3_summarystats_LatBin%1$02d_LonBin%2$02d_timesetps_001_%3$03d_V1.mat',lati,loni,i16daysSteps);
fn_summary = fullfile(fdirpre,fn_summary);
boo_summary = load(fn_summary,'tai93_desc','timestep_notfound');
addpath /asl/matlib/time
mtime = tai2dtime(airs2tai(boo_summary.tai93_desc));
dtime = datenum(mtime);
dtime = dtime(setdiff(1:length(dtime),boo_summary.timestep_notfound));

iNumLay6 = 6;
waha = reshape(results,iNumAnomTimeSteps,iNumAnomTiles,6);
for jjunk = 1 : iNumAnomTiles  %% 1 = global avg, 2 : 29 = latins     or 1 = global avg, 2 = tropics, 3 : 30 = latins
  for ijunk = 1 : iNumLay6
    P = nanpolyfit(dtime,squeeze(waha(:,jjunk,ijunk)),1);
    trendScalar(jjunk,ijunk) = P(1)*365.25;
  end
end

waha = reshape(resultsT,iNumAnomTimeSteps,iNumAnomTiles,iNumLay);
for jjunk = 1 : iNumAnomTiles  %% 1 = global avg, 2 : 29 = latins
  for ijunk = 1 : iNumLay
    P = nanpolyfit(dtime,squeeze(waha(:,jjunk,ijunk)),1);
    trendTz(jjunk,ijunk) = P(1)*365.25;
  end
end

waha = reshape(resultsWV,iNumAnomTimeSteps,iNumAnomTiles,iNumLay);
for jjunk = 1 : iNumAnomTiles  %% 1 = global avg, 2 : 29 = latins
  for ijunk = 1 : iNumLay
    P = nanpolyfit(dtime,squeeze(waha(:,jjunk,ijunk)),1);
    trendWV(jjunk,ijunk) = P(1)*365.25;
  end
end

waha = reshape(resultsO3,iNumAnomTimeSteps,iNumAnomTiles,iNumLay);
for jjunk = 1 : iNumAnomTiles  %% 1 = global avg, 2 : 29 = latins
  for ijunk = 1 : iNumLay
    P = nanpolyfit(dtime,squeeze(waha(:,jjunk,ijunk)),1);
    trendO3(jjunk,ijunk) = P(1)*365.25;
  end
end

iFig = iFig + 1; figure(iFig); clf; junk = trendScalar(iStartOffset:iNumAnomTiles,6); plot(rlat,junk); title('ST trends'); 
  plotaxis2; ylabel('dSKT/dt [Kyr]'); xlabel('Latitude')
iFig = iFig + 1; figure(iFig); clf; junk = trendTz(iStartOffset:iNumAnomTiles,:); pcolor(rlat,pavg,junk'); 
  xlabel('Latitude'); ylabel('Pressure [mb]');
  shading interp;  colorbar; colormap(llsmap5); caxis([-1 +1]*0.15); title('UMBC T trends'); set(gca,'ydir','reverse'); ylim([10 1000]); set(gca,'yscale','log')
iFig = iFig + 1; figure(iFig); clf; junk = trendWV(iStartOffset:iNumAnomTiles,:); pcolor(rlat,pavg,junk'); 
  shading interp;  colorbar; colormap(llsmap5); caxis([-1 +1]*0.015); title('UMBC WV trends'); set(gca,'ydir','reverse'); ylim([10 1000])
  xlabel('Latitude'); ylabel('Pressure [mb]');
iFig = iFig + 1; figure(iFig); clf; junk = trendO3(iStartOffset:iNumAnomTiles,:); pcolor(rlat,pavg,junk'); 
  shading interp;  colorbar; colormap(llsmap5); caxis([-1 +1]*0.015); title('UMBC O3 trends'); set(gca,'ydir','reverse'); ylim([10 1000]); set(gca,'yscale','log')
  xlabel('Latitude'); ylabel('Pressure [mb]');
