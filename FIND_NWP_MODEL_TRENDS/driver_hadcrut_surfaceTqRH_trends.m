addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies

if ~exist('had')
  hadq = read_netcdf_lls('HadISDH.blend/2021/HadISDH.blendq.1.3.0.2021f_FLATgridHOMBClocalSHIPboth5by5_anoms9120.nc');
end

if ~exist('hadrh')
  hadrh = read_netcdf_lls('HadISDH.blend/2021/HadISDH.blendRH.1.3.0.2021f_FLATgridHOMBClocalSHIPboth5by5_anoms9120.nc');
end

if ~exist('hadt')
  hadt = read_netcdf_lls('HadISDH.blend/2021/HadISDH.blendT.1.3.0.2021f_FLATgridHOMBClocalSHIPboth5by5_anoms9120.nc');
end

[hadq.Y,hadq.X]   = meshgrid(hadq.latitude,hadq.longitude);
[hadrh.Y,hadrh.X] = meshgrid(hadrh.latitude,hadrh.longitude);
[hadt.Y,hadt.X] = meshgrid(hadt.latitude,hadt.longitude);

figure(1); pcolor(hadq.q_abs(:,:,100)'); shading flat; colorbar; title('SH g/g')
figure(2); pcolor_coast(hadq.X,hadq.Y,hadq.q_abs(:,:,100)); shading flat; colorbar; title('SH g/g')

figure(3); pcolor(hadrh.rh_abs(:,:,100)'); shading flat; colorbar; title('RH')
figure(4); pcolor_coast(hadrh.X,hadrh.Y,hadrh.rh_abs(:,:,100)); shading flat; colorbar; title('RH')

figure(5); pcolor(hadt.t_abs(:,:,100)'); shading flat; colorbar; title('T')
figure(6); pcolor_coast(hadt.X,hadt.Y,hadt.t_abs(:,:,100)); shading flat; colorbar; title('T')

for ii = 1 : 6
  figure(ii); colormap(jet);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% first month is 1973/01
%% last  month is 2021/12

ii = 0;
for yy = 1973:2021
  for mm = 1 : 12
    ii = ii + 1;
    iaY(ii) = yy;
    iaM(ii) = mm;
    iaD(ii) = 15;
  end
end

boo1 = find(iaY == 2002 & iaM == 09);
boo2 = find(iaY == 2021 & iaM == 08);
boo = boo1:boo2;

daysSince2002 = change2days(iaY(boo),iaM(boo),iaD(boo),2002);

for jj = 1 :36
  if mod(jj,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end
  for ii = 1 : 72
    data = squeeze(hadq.q_abs(ii,jj,boo));
    moo = find(isfinite(data));
    if length(moo) > 10
      [B,err] = Math_tsfit_lin_robust(daysSince2002(moo),data(moo)/mean(data(moo)),4);
      trend.q(ii,jj) = B(2);
      trend.q_unc(ii,jj) = err.se(2);
    else
      trend.q(ii,jj) = NaN;
      trend.q_unc(ii,jj) = NaN;
    end
  end
end
fprintf(1,'\n');

for jj = 1 :36
  if mod(jj,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end
  for ii = 1 : 72
    data = squeeze(hadrh.rh_abs(ii,jj,boo));
    moo = find(isfinite(data));
    if length(moo) > 10
      [B,err] = Math_tsfit_lin_robust(daysSince2002(moo),data(moo),4);
      trend.rh(ii,jj) = B(2);
      trend.rh_unc(ii,jj) = err.se(2);
    else
      trend.rh(ii,jj) = NaN;
      trend.rh_unc(ii,jj) = NaN;
    end
  end
end
fprintf(1,'\n');

for jj = 1 :36
  if mod(jj,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end
  for ii = 1 : 72
    data = squeeze(hadt.t_abs(ii,jj,boo));
    moo = find(isfinite(data));
    if length(moo) > 10
      [B,err] = Math_tsfit_lin_robust(daysSince2002(moo),data(moo),4);
      trend.t(ii,jj) = B(2);
      trend.t_unc(ii,jj) = err.se(2);
    else
      trend.t(ii,jj) = NaN;
      trend.t_unc(ii,jj) = NaN;
    end
  end
end
fprintf(1,'\n');

figure(1); pcolor(hadq.X,hadq.Y,trend.q); shading flat; colorbar; caxis([-1 +1]*0.025); title('SH trend (fraction/yr)')
figure(2); pcolor_coast(hadq.X,hadq.Y,trend.q); shading flat; colorbar; caxis([-1 +1]*0.025); title('SH trend (fraction/yr)')

figure(3); pcolor(hadrh.X,hadrh.Y,trend.rh); shading flat; colorbar; caxis([-1 +1]*2); title('RH trend (percent/yr)')
figure(4); pcolor_coast(hadrh.X,hadrh.Y,trend.rh); shading flat; colorbar; caxis([-1 +1]*2); title('RH trend (percent/yr)')

figure(5); pcolor(hadt.X,hadt.Y,trend.t); shading flat; colorbar; caxis([-1 +1]*0.25); title('T trend (K/yr)')
figure(6); pcolor_coast(hadt.X,hadt.Y,trend.t); shading flat; colorbar; caxis([-1 +1]*0.25); title('T trend (K/yr)')

for ii = 1 : 6
  figure(ii); colormap(usa2);
end

trend.lat = hadrh.Y;
trend.lon = hadrh.X;
trend.longitude = hadrh.longitude;
trend.latitude = hadrh.latitude;
trend.comment = 'see driver_hadcrut_surfaceTqRH_trends.m';

save hadcrut_surfaceTqRH_trends_2002_2021.mat trend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now turn this into our 72x64 tiles
    load ../AIRS_gridded_STM_May2021_trendsonlyCLR/latB64.mat
    rlat65 = latB2; rlon73 = -180 : 5 : +180;
    rlon = -180 : 5 : +180;  rlat = latB2; 
    rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
    rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

trend.t_72x64 = nan(72,64);
trend.t_unc_72x64 = nan(72,64);
trend.q_72x64 = nan(72,64);
trend.q_unc_72x64 = nan(72,64);
trend.rh_72x64 = nan(72,64);
trend.rh_unc_72x64 = nan(72,64);
for ii = 1: 72
  for jj = 1 : 64
    booX = find(trend.longitude >= rlon73(ii) & trend.longitude < rlon73(ii+1));
    booY = find(trend.latitude >= rlat65(jj) & trend.latitude < rlat65(jj+1));

    scaleY = 1;    scaleX = 0;
    scaleY = 1;    scaleX = 1;
    booX = find(trend.longitude >= rlon73(ii)-5*scaleX & trend.longitude < rlon73(ii+1)+5*scaleX);
    booY = find(trend.latitude >= rlat65(jj)-3*scaleY & trend.latitude < rlat65(jj+1)+3*scaleY);

    if length(booX) > 0 & length(booY) > 0
      junk = trend.t(booX,booY);     trend.t_72x64(ii,jj) = nanmean(junk(:));
      junk = trend.t_unc(booX,booY); trend.t_unc_72x64(ii,jj) = nanmean(junk(:));

      junk = trend.q(booX,booY);     trend.q_72x64(ii,jj) = nanmean(junk(:));
      junk = trend.q_unc(booX,booY); trend.q_unc_72x64(ii,jj) = nanmean(junk(:));

      junk = trend.rh(booX,booY);     trend.rh_72x64(ii,jj) = nanmean(junk(:));
      junk = trend.rh_unc(booX,booY); trend.rh_unc_72x64(ii,jj) = nanmean(junk(:));
    end
  end
end

trend.rlat65 = rlat65; trend.rlat = rlat;
trend.rlon73 = rlon73; trend.rlon = rlon;
[trend.Y72x64,trend.X72x64] = meshgrid(trend.rlat,trend.rlon);
save hadcrut_surfaceTqRH_trends_2002_2021.mat trend

figure(1); pcolor_coast(trend.X72x64,trend.Y72x64,trend.q_72x64,+1); shading flat; colorbar; caxis([-1 +1]*0.025); title('SH trend (fraction/yr)')
figure(2); pcolor_coast(trend.X72x64,trend.Y72x64,trend.q_72x64,-1); shading flat; colorbar; caxis([-1 +1]*0.025); title('SH trend (fraction/yr)')

figure(3); pcolor_coast(trend.X72x64,trend.Y72x64,trend.rh_72x64,+1); shading flat; colorbar; caxis([-1 +1]*2); title('RH trend (percent/yr)')
figure(4); pcolor_coast(trend.X72x64,trend.Y72x64,trend.rh_72x64,-1); shading flat; colorbar; caxis([-1 +1]*2); title('RH trend (percent/yr)')

figure(5); pcolor_coast(trend.X72x64,trend.Y72x64,trend.t_72x64,+1); shading flat; colorbar; caxis([-1 +1]*0.25); title('T trend (K/yr)')
figure(6); pcolor_coast(trend.X72x64,trend.Y72x64,trend.t_72x64,-1); shading flat; colorbar; caxis([-1 +1]*0.25); title('T trend (K/yr)')

figure(7); pcolor(hadrh.X,hadrh.Y,trend.rh); shading flat; colorbar; caxis([-1 +1]*2); title('RH trend (percent/yr)')
figure(8); pcolor_coast(hadrh.X,hadrh.Y,trend.rh); shading flat; colorbar; caxis([-1 +1]*2); title('RH trend (percent/yr)')

figure(9); plot(nanmean(trend.rh_72x64,1),trend.rlat,nanmean(trend.rh,1),trend.latitude)

for ii = 1 : 8
  figure(ii); colormap(usa2);
end
