addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /asl/matlib/aslutil/
addpath /asl/matlib/maps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data are available for each month since January 1850, on a 5 degree grid.

numtime = (2020-1850+1)*12;
airsS   = (2002-1850+0)*12 + 9;  %% 2002/09
airsE   = (2020-1850+0)*12 + 8;  %% 2020/08
fprintf(1,'number of AIRS years to process = %2i \n',(airsE-airsS+1)/12)

for ii = 1 : 100
  fname = ['/asl/models/hadcrut4/HadCRUT.4.6.0.0.anomalies.' num2str(ii) '.nc'];
  junk = read_netcdf_lls(fname);
  data(ii,:,:,:) = junk.temperature_anomaly;
end

lat = junk.latitude;
lon = junk.longitude;
[Lon,Lat] = meshgrid(lon,lat); Lon = Lon(:); Lat = Lat(:);

figure(1); colormap(usa2);
junkS = squeeze(junk.temperature_anomaly(:,:,airsS));
junkE = squeeze(junk.temperature_anomaly(:,:,airsE));
  junkS = junkS'; junkS = junkS(:); scatter_coast(Lon,Lat,50,junkS); 
  caxis([-5 +5]); pause(0.1)

iCnt = 0;
for yy = 2002 : 2020
  if yy == 2002
    mS = 09; mE = 12;
  elseif yy == 2020
    mS = 01; mE = 08;
  else
    mS = 01; mE = 12;
  end
  for mm = mS : mE
    iCnt = iCnt + 1;
    rtime(iCnt) = utc2taiSergio(yy,mm,15,12.0);
  end
end
[yy,mm,dd,hh] = tai2utcSergio(rtime);
daysSince2002 = change2days(yy,mm,dd,2002);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [B, err, stats]=Math_tsfit_lin_robust_NanInfDump(x0,y0,n,iPosOrAll)
hadCRU4_results.Lon = Lon;
hadCRU4_results.Lat = Lat;
hadCRU4_results.lon = lon;
hadCRU4_results.lat = lat;

newdata = data(:,:,:,airsS:airsE);
newdata = nanmean(newdata,1);
alldata = squeeze(newdata);
for jj = 1 : 36
  for ii = 1 : 72
    ydata = squeeze(alldata(ii,jj,:));
    [B, err, stats] = Math_tsfit_lin_robust_NanInfDump(daysSince2002,ydata,4,-1);
    hadCRU4_results.trend_t(ii,jj) = B(2);
    hadCRU4_results.trend_t_err(ii,jj) = err(2);
  end
end
hadCRU4_results.trend_t = hadCRU4_results.trend_t'; hadCRU4_results.trend_t = hadCRU4_results.trend_t(:);
figure(1); scatter_coast(Lon,Lat,50,hadCRU4_results.trend_t); title('dST/dt K/yr'); caxis([-0.15 +0.15]);

comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_find_hadCRUT4_trends.m';
save hadCRU4_trends.mat  hadCRU4_results comment

