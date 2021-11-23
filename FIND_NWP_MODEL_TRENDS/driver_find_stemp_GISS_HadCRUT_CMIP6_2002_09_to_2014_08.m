addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /asl/matlib/aslutil/
addpath /asl/matlib/maps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = '/asl/models/gistemp4/gistemp1200_GHCNv4_ERSSTv5.nc';
giss = read_netcdf_lls(fname);

%% want to do doy since 01/01/1800
%% want to overlap with CMIP6
doyS = change2days(2002,09,01,1800);
doyE = change2days(2014,08,31,1800);

plot(giss.time)
  line([0 length(giss.time)],[doyS doyS]);
  line([0 length(giss.time)],[doyE doyE]);
oo = find(giss.time >= doyS & giss.time <= doyE);

[Y,X] = meshgrid(giss.lat,giss.lon);

warning off
for ii = 1 : 180
  if mod(ii,50) == 0
    fprintf(1,'+')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end
  for jj = 1 : 90
    data = squeeze(giss.tempanomaly(ii,jj,:));
    aha = find(isfinite(data(oo)));
    if length(aha) > 20
      [B stats] = Math_tsfit_lin_robust(double(giss.time(oo(aha))),double(data(oo(aha))),4);
      giss_trend(ii,jj) = B(2);  
      giss_trend_err(ii,jj) = stats.se(2);  
    else
      giss_trend(ii,jj) = NaN;  
      giss_trend_err(ii,jj) = NaN;
    end
  end
end
warning on

simplemap(Y(:),X(:),giss_trend(:))
colormap(usa2);
caxis([-0.1 +0.1])

aslmap(1,-90:2:+90,-180:2:+180,smoothn(giss_trend',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('d/dt GISS K/yr'); 
caxis([-0.1 +0.1]);

%% so now interp2 to the Howard tiles
%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2; 
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y4608,X4608] = meshgrid(rlat,rlon);
giss_trend4608     = interp2(Y,X,giss_trend,Y4608,X4608);
giss_trend_err4608 = interp2(Y,X,giss_trend_err,Y4608,X4608);

aslmap(1,rlat65,rlon73,smoothn(giss_trend4608',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('2002/09-2014/08 d/dt GISS K/yr'); 
caxis([-0.1 +0.1]);

comment = 'see driver_find_stemp_GISS_HadCRUT_CMIP6_2002_09_to_2014_08.m';
save giss_hadcrut_cmip6_stemp_trends_2002_09_to_2014_08_CMIP6.mat giss_* X Y comment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear data

numtime = (2020-1850+1)*12;
airsS   = (2002-1850+0)*12 + 9;  %% 2002/09
airsE   = (2014-1850+0)*12 + 8;  %% 2014/08
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
for yy = 2002 : 2014
  if yy == 2002
    mS = 09; mE = 12;
  elseif yy == 2014
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

%hadCRU4_results.trend_t = hadCRU4_results.trend_t'; hadCRU4_results.trend_t = hadCRU4_results.trend_t(:);
rlat37 = -090 : 5 : +090;
rlat73 = -180 : 5 : +180;

aslmap(2,rlat37,rlon73,smoothn(hadCRU4_results.trend_t',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('2002/09-2014/08 d/dt HadCRUT K/yr'); 
caxis([-0.1 +0.1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_CMIP6/clust_compute_cmip6_profile_rtpfiles.m

cmip6  = read_netcdf_lls('/asl/s1/rkramer/CMIP6_modelmean_Sergio/tas_CMIP6_historical_modelmean.nc');
cmip6x = read_netcdf_lls('/asl/s1/rkramer/CMIP6_modelmean_Sergio/ta_CMIP6_historical_modelmean.nc');

numtime = (2015-1850)*12;   %% yes, the 1980 timesteps Ryan has saved

allyears = [];
allmonths = [];
for ii = 1850:2014
  junkyear = ii*ones(1,12);
  junkmonth = [1 2 3 4 5 6 7 8 9 10 11 12];
  allyears = [allyears junkyear];
  allmonths = [allmonths junkmonth];
end

i2002 = (2002-1850)*12 + 9;  %% 2002/09
i2014 = (2014-1850)*12 + 8;  %% 2014/09
[allyears(i2002) allmonths(i2002) allyears(i2014) allmonths(i2014)]
oo = i2002:i2014;

%{
for ii = i2002:i2014
  data = squeeze(cmip6.tas(:,:,ii));
  aslmap(1,rlat65,rlon73,smoothn(data',1), [-90 +90],[-180 +180]);  colormap(jet);  title(['STEMP ' num2str(allyears(ii)) '/'  num2str(allmonths(ii))]);   caxis([200 300]); pause(0.1)
end
for ii = i2002:i2014
  data = squeeze(cmip6x.ta(:,:,3,ii));
  aslmap(1,rlat65,rlon73,smoothn(data',1), [-90 +90],[-180 +180]);  colormap(jet);  title(['AIR TEMP 850 mb ' num2str(allyears(ii)) '/'  num2str(allmonths(ii))]);   caxis([200 300]); pause(0.1)
end
%}

warning off
for ii = 1 : 72
  if mod(ii,50) == 0
    fprintf(1,'+')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end
  for jj = 1 : 64
    data = squeeze(cmip6.tas(ii,jj,:));
    aha = find(isfinite(data(oo)));
    if length(aha) > 20
      [B stats] = Math_tsfit_lin_robust(double(cmip6.time(oo(aha))/24),double(data(oo(aha))),4);
      cmip6_trend(ii,jj) = B(2);  
      cmip6_trend_err(ii,jj) = stats.se(2);  
    else
      cmip6_trend(ii,jj) = NaN;  
      cmip6_trend_err(ii,jj) = NaN;
    end
  end
end
warning on

warning off
for ii = 1 : 72
  if mod(ii,50) == 0
    fprintf(1,'+')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end
  for jj = 1 : 64
    for kk = 1:19
      data = squeeze(cmip6x.ta(ii,jj,kk,:));
      aha = find(isfinite(data(oo)));
      if length(aha) > 20
        [B stats] = Math_tsfit_lin_robust(double(cmip6.time(oo(aha))/24),double(data(oo(aha))),4);
        cmip6ta_trend(kk,ii,jj) = B(2);  
        cmip6ta_trend_err(kk,ii,jj) = stats.se(2);  
      else
        cmip6ta_trend(kk,ii,jj) = NaN;  
        cmip6ta_trend_err(kk,ii,jj) = NaN;
      end
    end
  end
end
warning on

aslmap(1,-90:2:+90,-180:2:+180,smoothn(giss_trend',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('2002/09-2014/08 d/dt GISS K/yr'); 
caxis([-0.1 +0.1]);

aslmap(2,rlat37,rlon73,smoothn(hadCRU4_results.trend_t',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('2002/09-2014/08 d/dt HadCRUT K/yr'); 
caxis([-0.1 +0.1]);

aslmap(3,rlat65,rlon73,smoothn(cmip6_trend',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('2002/09-2014/08 d/dt CMIP6 K/yr'); 
caxis([-0.1 +0.1]);

for kk=1:2:19; 
  aslmap(4,rlat65,rlon73,smoothn(squeeze(cmip6ta_trend(kk,:,:))',1), [-90 +90],[-180 +180]);  colormap(usa2);  title([num2str(cmip6x.plev(kk)/100) ' mb']); caxis([-0.1 +0.1]); pause;
end
figure(4); clf; pcolor(rlat,cmip6x.plev/100,squeeze(nanmean(cmip6ta_trend,2))); shading interp; colormap(usa2); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); caxis([-0.15 +0.15])

save giss_hadcrut_cmip6_stemp_trends_2002_09_to_2014_08_CMIP6.mat giss_* cmip6*_trend* hadCRU4_results X Y comment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
miaw = load('CMIP6_atm_data_2002_09_to_2014_08.mat');
for ii = 1 : 144
  scatter_coast(miaw.all.rlon,miaw.all.rlat,50,miaw.all.stemp(ii,:)); 
  title(['STEMP ' num2str(allyears(ii+i2002)) '/'  num2str(allmonths(ii+i2002))]); colormap(jet);   caxis([200 300]); pause(0.1)
end
for ii = 1 : 144
  scatter_coast(miaw.all.rlon,miaw.all.rlat,50,squeeze(miaw.all.ptemp(ii,95,:))); 
  title(['AIR TEMP LAY 95 ' num2str(allyears(ii+i2002)) '/'  num2str(allmonths(ii+i2002))]); colormap(jet);   caxis([200 300]); pause(0.1)
end

doy = change2days(miaw.all.yy,miaw.all.mm,miaw.all.dd,2002);
oo = 1 : 144;
warning off
for ii = 1 : 4608
  if mod(ii,1000) == 0
    fprintf(1,'+')
  elseif mod(ii,100) == 0
    fprintf(1,'.')
  end
  for kk = 1:100
    data = squeeze(miaw.all.ptemp(:,kk,ii));
    aha = find(isfinite(data(oo)));
    if length(aha) > 20
      [B stats] = Math_tsfit_lin_robust(doy(oo(aha)),double(data(oo(aha))),4);
      cmip6rtp_trend(kk,ii) = B(2);  
      cmip6rtp_trend_err(kk,ii) = stats.se(2);  
    else
      cmip6rtp_trend(kk,ii) = NaN;  
      cmip6rtp_trend_err(kk,ii) = NaN;
    end
  end
end
warning on

boo = load('/home/sergio/MATLABCODE/airslevels.dat');
  pjunkN = boo(1:100)-boo(2:101);
  pjunkD = log(boo(1:100)./boo(2:101));
  pavgLAY = pjunkN./pjunkD;
  pavgLAY = flipud(pavgLAY);

figure(4); clf; pcolor(rlat,cmip6x.plev/100,squeeze(nanmean(cmip6ta_trend,2))); shading interp; colormap(usa2); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); caxis([-0.15 +0.15])
  title('dT/dt from Ryan levels'); xlabel('Latitude'); ylabel('P(mb)'); colorbar
figure(5); clf; pcolor(rlat,pavgLAY(1:100),squeeze(nanmean(reshape(cmip6rtp_trend,100,72,64),2))); shading interp; colormap(usa2); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); caxis([-0.15 +0.15])
  title('dT/dt from SARTA layers'); xlabel('Latitude'); ylabel('P(mb)'); colorbar
%}
%{
era5 = load('ERA5_atm_data_2002_09_to_2021_07_desc.mat');
figure(6)

for ii = 1 : 144
  scatter_coast(era5.all.rlon,era5.all.rlat,50,era5.all.stemp(ii,:)); 
  title(['STEMP ' num2str(allyears(ii+i2002)) '/'  num2str(allmonths(ii+i2002))]); colormap(jet);   caxis([200 300]); pause(0.1)
end
for ii = 1 : 144
  scatter_coast(era5.all.rlon,era5.all.rlat,50,squeeze(era5.all.ptemp(ii,95,:))); 
  title(['AIR TEMP LAY 95 ' num2str(allyears(ii+i2002)) '/'  num2str(allmonths(ii+i2002))]); colormap(jet);   caxis([200 300]); pause(0.1)
end

doy = change2days(era5.all.yy,era5.all.mm,era5.all.dd,2002);
oo = 1 : 144;
warning off
for ii = 1 : 4608
  if mod(ii,1000) == 0
    fprintf(1,'+')
  elseif mod(ii,100) == 0
    fprintf(1,'.')
  end
  for kk = 1:100
    data = squeeze(era5.all.ptemp(:,kk,ii));
    aha = find(isfinite(data(oo)));
    if length(aha) > 20
      [B stats] = Math_tsfit_lin_robust(doy(oo(aha)),double(data(oo(aha))),4);
      era5_trend(kk,ii) = B(2);  
      era5_trend_err(kk,ii) = stats.se(2);  
    else
      era5_trend(kk,ii) = NaN;  
      era5_trend_err(kk,ii) = NaN;
    end
  end
end
warning on

boo = load('/home/sergio/MATLABCODE/airslevels.dat');
  pjunkN = boo(1:100)-boo(2:101);
  pjunkD = log(boo(1:100)./boo(2:101));
  pavgLAY = pjunkN./pjunkD;
  pavgLAY = flipud(pavgLAY);

load llsmap5

figure(4); clf; pcolor(rlat,cmip6x.plev/100,squeeze(nanmean(cmip6ta_trend,2))); shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); caxis([-0.15 +0.15])
  title('CMIP6 dT/dt from Ryan levels'); xlabel('Latitude'); ylabel('P(mb)'); colorbar
figure(5); clf; pcolor(rlat,pavgLAY(1:100),squeeze(nanmean(reshape(cmip6rtp_trend,100,72,64),2))); shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); caxis([-0.15 +0.15])
  title('CMIP6 dT/dt from SARTA layers'); xlabel('Latitude'); ylabel('P(mb)'); colorbar
figure(6); clf; pcolor(rlat,pavgLAY(1:100),squeeze(nanmean(reshape(era5_trend,100,72,64),2))); shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); caxis([-0.15 +0.15])
  title('ERA5 dT/dt from SARTA layers'); xlabel('Latitude'); ylabel('P(mb)'); colorbar

%}

