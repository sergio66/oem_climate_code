addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /asl/matlib/aslutil/
addpath /asl/matlib/maps

%% see driver_find_GISS_AIRSL3_ERA5_skt_trends.m

do_XX_YY_from_X_Y

iDoGiss = +1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iNumCyclesFit = 4; %% standard fit
iNumCyclesFit = 1; %% giss data is already an anomaly

fname = '/asl/models/gistemp4/gistemp1200_GHCNv4_ERSSTv5.nc';
fname = '/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/GISTEMP/F77/gistemp1200_ERSST.nc';
  
figure(1); clf; 
  
giss = read_netcdf_lls(fname);
giss.lat = double(giss.lat);
giss.lon = double(giss.lon);
giss.time = double(giss.time);
giss.time_bnds = double(giss.time_bnds);

%% want to do doy since 01/01/1800
yS = 2002; doyS = change2days(yS,09,01,1800);
yE = 2022; doyE = change2days(yE,08,31,1800);

a = giss;
  
yy = 2002;
mm = 08;

iaNumYears = [20];
for tt = 1 : iaNumYears*12

  mm = mm + 1;
  if mm > 12
    yy = yy + 1;
    mm = 1;
  end
  ysave(tt) = yy;
  msave(tt) = mm;

  doy = change2days(yy,mm,01,1800);
  boo = find(a.time >= doy,1);
  time_giss(tt) = boo;

  for jj = 1 : 64
    for ii = 1 : 72
      xlim1 = rlon73(ii);
      xlim2 = rlon73(ii+1);
      ylim1 = rlat65(jj);
      ylim2 = rlat65(jj+1);

      booX = find(a.lon >= xlim1 & a.lon < xlim2);
      booY = find(a.lat >= ylim1 & a.lat < ylim2);
      wah = squeeze(a.tempanomaly(booX,booY,boo));
      wah = wah(:);
      mean_stempIJ_giss(tt,ii,jj) = nanmean(wah);    
      mean_rlonIJ_giss(ii,jj)  = nanmean(a.lon(booX));
      mean_rlatIJ_giss(ii,jj)  = nanmean(a.lat(booY));    

      booX = find(a.lon >= rlon(ii),1);
      booY = find(a.lat >= rlat(jj),1);
      wah = squeeze(a.tempanomaly(booX,booY,boo));
      wah = wah(:);
      cntr_stempIJ_giss(tt,ii,jj) = wah;    
      cntr_rlonIJ_giss(ii,jj)  = a.lon(booX);
      cntr_rlatIJ_giss(ii,jj)  = a.lat(booY);    

    end
  end
end

comment = '/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_giss_surfaceT.m';
save L3_SKT_TIMESERIES_2002_09_to_2022_08/giss_skt_mean.mat mean_rlonIJ_giss mean_rlatIJ_giss mean_stempIJ_giss cntr_rlonIJ_giss cntr_rlatIJ_giss cntr_stempIJ_giss time_giss comment ysave msave
figure(1); clf; plot(time_giss)
figure(2); clf; scatter_coast(mean_rlonIJ_giss,mean_rlatIJ_giss,50,squeeze(nanmean(mean_stempIJ_giss,1))); colormap jet
figure(3); clf; scatter_coast(mean_rlonIJ_giss,mean_rlatIJ_giss,50,squeeze(nanmean(cntr_stempIJ_giss,1))); colormap jet
figure(4); clf; scatter_coast(mean_rlonIJ_giss,mean_rlatIJ_giss,50,squeeze(nanmean(cntr_stempIJ_giss-mean_stempIJ_giss,1))); colormap(usa2)

