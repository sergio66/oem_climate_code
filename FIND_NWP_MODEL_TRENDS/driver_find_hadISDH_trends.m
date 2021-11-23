%{
https://www.metoffice.gov.uk/hadobs/hadisdh/downloadblend1112020.html

Met Office Hadley Centre observations datasets

HadISDH.blend data download - version 1.1.1.2020f


The current version of HadISDH.blend is 1.1.1.2020f. For previous versions please contact the dataset maintainers.

These data are available under a free and unrestricted Open Government Licence UK.

NetCDF gridded files contain fields on 5° by 5° resolution from -177.5° W/87.5° N to 177.5° E/-87.5° S.

NetCDF station files contain monthly time series from January 1973 to December 2020.

Missing data are represented by -1e30.

Variables are described below.

HadISDH.blend*.1.1.1.2020f grids NetCDF	Land and Marine surface grids	
HadISDH.blendq.1.1.1.2020f_FLATgridHOMBClocalSHIPboth5by5_anoms8110.nc, 
HadISDH.blendRH.1.1.1.2020f_FLATgridHOMBClocalSHIPboth5by5_anoms8110.nc, 
HadISDH.blende.1.1.1.2020f_FLATgridHOMBClocalSHIPboth5by5_anoms8110.nc, 
HadISDH.blendTd.1.1.1.2020f_FLATgridHOMDPDBClocalSHIPboth5by5_anoms8110.nc, 
HadISDH.blendTw.1.1.1.2020f_FLATgridHOMBClocalSHIPboth5by5_anoms8110.nc, 
HadISDH.blendT.1.1.1.2020f_FLATgridHOMBClocalSHIPboth5by5_anoms8110.nc, 
HadISDH.blendDPD.1.1.1.2020f_FLATgridHOMBClocalSHIPboth5by5_anoms8110.nc,
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
https://www.metoffice.gov.uk/hadobs/hadisdh/

Gridded products are available for six humidity variables in addition to temperature:


Specific humidity (q), expressed in g kg-1. The ratio of the mass of
water vapour to the mass of moist air.

Relative humidity (RH), expressed as a percentage (%rh). The amount of
water vapour in the air compared to how much water could potentially
be held as a vapour at that temperature.

Dew point temperature (Td), expressed in °C. The temperature at which
the air becomes saturated at that current level of water vapour,
measured by artificially cooling a surface until water condenses onto
it.

Wet bulb temperature (Tw), expressed in °C. The amount of evaporative
cooling of a thermometer in a moistened wick. Air that is not
saturated will evaporate water from the wick, cooling the 'wet bulb'
thermometer.

Vapour pressure (e), expressed in hPa. The partial pressure exerted by
water vapour alone.

Dew point depression (DPD), expressed in °C. The amount the air has to
be cooled by to reach its dew point temperature.

Temperature (T), expressed in °C. The temperature measured by the dry
bulb thermometer.
%}

% The dataset begins in January 1973 and is updated annually till Dec 2020

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/FIND_TRENDS/
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies

q = read_netcdf_lls('/asl/models/hadISDH/2020/HadISDH.blendq.1.1.1.2020f_FLATgridHOMBClocalSHIPboth5by5_anoms8110.nc');
RH = read_netcdf_lls('/asl/models/hadISDH/2020/HadISDH.blendRH.1.1.1.2020f_FLATgridHOMBClocalSHIPboth5by5_anoms8110.nc');
%vp = read_netcdf_lls('/asl/models/hadISDH/2020/HadISDH.blende.1.1.1.2020f_FLATgridHOMBClocalSHIPboth5by5_anoms8110.nc');
tw = read_netcdf_lls('/asl/models/hadISDH/2020/HadISDH.blendTw.1.1.1.2020f_FLATgridHOMBClocalSHIPboth5by5_anoms8110.nc');
t = read_netcdf_lls('/asl/models/hadISDH/2020/HadISDH.blendT.1.1.1.2020f_FLATgridHOMBClocalSHIPboth5by5_anoms8110.nc');

numtime = (2020-1973+1)*12;
airsS   = (2002-1973+0)*12 + 9;  %% 2002/09
airsE   = (2020-1973+0)*12 + 8;  %% 2020/08
fprintf(1,'number of AIRS years to process = %2i \n',(airsE-airsS+1)/12)

lat = t.latitude;
lon = t.longitude;
[Lon,Lat] = meshgrid(lon,lat); Lon = Lon(:); Lat = Lat(:);

figure(1); colormap(usa2);

%for ttime = 1 : numtime
for ttime = airsS : 12 : airsE
  timeYY = floor(ttime/12);
  timeMM = ttime-(timeYY)*12;
  if timeMM == 0
    %timeYY = timeYY + 1;
    timeMM = 12;
  end
  timeYY = 1973 + timeYY;
  data = squeeze(t.t_anoms(:,:,ttime)); data = data'; data = data(:); scatter_coast(Lon,Lat,50,data); title([num2str(timeYY) '/' num2str(timeMM)]); 
  caxis([-5 +5]); pause(0.1)
end

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
hadISDH_results.Lon = Lon;
hadISDH_results.Lat = Lat;
hadISDH_results.lon = lon;
hadISDH_results.lat = lat;

alldata = t.t_anoms(:,:,airsS:airsE);
for jj = 1 : 36
  for ii = 1 : 72
    ydata = squeeze(alldata(ii,jj,:));
    [B, err, stats] = Math_tsfit_lin_robust_NanInfDump(daysSince2002,ydata,4,-1);
    hadISDH_results.trend_t(ii,jj) = B(2);
    hadISDH_results.trend_t_err(ii,jj) = err(2);
  end
end
hadISDH_results.trend_t = hadISDH_results.trend_t'; hadISDH_results.trend_t = hadISDH_results.trend_t(:);
figure(1); scatter_coast(Lon,Lat,50,hadISDH_results.trend_t); title('dST/dt K/yr'); caxis([-0.15 +0.15]);

alldata = tw.tw_anoms(:,:,airsS:airsE);
for jj = 1 : 36
  for ii = 1 : 72
    ydata = squeeze(alldata(ii,jj,:));
    [B, err, stats] = Math_tsfit_lin_robust_NanInfDump(daysSince2002,ydata,4,-1);
    hadISDH_results.trend_twet(ii,jj) = B(2);
    hadISDH_results.trend_twet_err(ii,jj) = err(2);
  end
end
hadISDH_results.trend_twet = hadISDH_results.trend_twet'; hadISDH_results.trend_twet = hadISDH_results.trend_twet(:);
figure(2); scatter_coast(Lon,Lat,50,hadISDH_results.trend_twet); title('dSTwet/dt K/yr');  caxis([-0.15 +0.15]);

alldata = q.q_anoms(:,:,airsS:airsE);
for jj = 1 : 36
  for ii = 1 : 72
    ydata = squeeze(alldata(ii,jj,:));
    ydata = ydata/nanmean(ydata+eps);
    [B, err, stats] = Math_tsfit_lin_robust_NanInfDump(daysSince2002,ydata,4,-1);
    hadISDH_results.trend_fracQ(ii,jj) = B(2);
    hadISDH_results.trend_fracQ_err(ii,jj) = err(2);
  end
end
hadISDH_results.trend_fracQ = hadISDH_results.trend_fracQ'; hadISDH_results.trend_fracQ = hadISDH_results.trend_fracQ(:);
figure(3); scatter_coast(Lon,Lat,50,hadISDH_results.trend_fracQ); title('d(frac q) /dt /yr'); caxis([-1 +1]);

alldata = q.q_anoms(:,:,airsS:airsE);
for jj = 1 : 36
  for ii = 1 : 72
    ydata = squeeze(alldata(ii,jj,:));
    [B, err, stats] = Math_tsfit_lin_robust_NanInfDump(daysSince2002,ydata,4,-1);
    hadISDH_results.trend_Q(ii,jj) = B(2);
    hadISDH_results.trend_Q_err(ii,jj) = err(2);
  end
end
hadISDH_results.trend_Q = hadISDH_results.trend_Q'; hadISDH_results.trend_Q = hadISDH_results.trend_Q(:);
figure(3); scatter_coast(Lon,Lat,50,hadISDH_results.trend_Q); title('d(frac q) /dt /yr'); caxis([-1 +1]);

alldata = RH.rh_anoms(:,:,airsS:airsE);
for jj = 1 : 36
  for ii = 1 : 72
    ydata = squeeze(alldata(ii,jj,:));
    [B, err, stats] = Math_tsfit_lin_robust_NanInfDump(daysSince2002,ydata,4,-1);
    hadISDH_results.trend_rh(ii,jj) = B(2);
    hadISDH_results.trend_rh_err(ii,jj) = err(2);
  end
end
hadISDH_results.trend_rh = hadISDH_results.trend_rh'; hadISDH_results.trend_rh = hadISDH_results.trend_rh(:);
figure(4); scatter_coast(Lon,Lat,50,hadISDH_results.trend_rh); title('d(RH) /dt /yr'); caxis([-0.1 +0.1]);

comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/hadISDH.m';
save hadISDH_trends.mat  hadISDH_results comment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_HadSurf_trends
