addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
%% Q	kg/kg	30	A	Specific humidity
%% ozone maybe mol/mol
% [Twb,Teq,epott,relhum] = wetbulb(cesm_stemp-273.15,cesm_spres*100,squeeze(cesm_gas_1(:,58,:,:)),0);
% Airs_TwetSurf     = real(get_wet_bulb_temperature(Airs_STemp,Airs_RHSurf));
% figure(1);  pcolor(flipud(squeeze(nanmean(double(Airs_RHSurf),1)))); shading interp; colorbar; caxis([0 120]); title('mean RH Surf'); colormap(jet)
% figure(2);  pcolor(flipud(squeeze(nanmean(double(Airs_STemp),1))));  shading interp; colorbar; title('mean TSurf');    caxis([200 320]); colormap(jet)
% figure(3);  pcolor(flipud(squeeze(nanmean(double(Airs_TwetSurf),1))));   shading interp; colorbar; title('mean TwetSurf'); caxis([200 320]); colormap(jet);

[yy,mm,dd,hh] = matime2yymmddhh(cesm_Date);
days = change2days(yy,mm,dd,2002);
days = unique(days,'stable');

yy   = yy(woo);
mm   = mm(woo);
dd   = dd(woo);
hh   = hh(woo);
days = days(woo);

figure(1)
for ii = 1 : 1
  boo = squeeze(cesm_stemp(woo(ii),:,:));
  zoo = find(isnan(boo)); boo(zoo) = 250;
  simplemap(Lat,wrapTo180(Lon),boo); caxis([200 320]); colorbar
  %simplemap(boo); caxis([200 350]); colorbar
  pause(1)
end

figure(2)
[salti,landfrac] =  usgs_deg10_dem(Lat,Lon);
simplemap(Lat,Lon,landfrac);
simplemap(Lat,Lon,landfrac,2);

ocean = find(landfrac <= 0.001);
ocean = find(landfrac >= 0.0);
for ii = 1 : 12 : length(days)
  stemp = squeeze(cesm_stemp(ii,:,:));
  figure(2);
  simplemap(Lat(ocean),Lon(ocean),stemp(ocean),2);
  title([num2str(yy(ii)) '/' num2str(mm(ii))])
  caxis([220 310]); colorbar
  pause(1)
end

ii = find(yy == min(2012,StopY) & mm == 7);
  stemp = squeeze(cesm_stemp(ii,:,:));
  figure(3);
  simplemap(Lat(ocean),Lon(ocean),stemp(ocean),2);
  title([num2str(yy(ii)) '/' num2str(mm(ii))])
  caxis([220 310]); colorbar

iDebug = false;
bo = cesm_stemp;
[tmax,aa,bb] = size(bo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : length(latbins)-1
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1) & ...
              landfrac < 0.001);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1));
  for tt = 1 : tmax
    zz = cesm_stemp(tt,:,:);
    xx = zz(ix); xx = xx(:);
    good = find(xx > 0);
    save_stemp(ii,tt) = nanmean(xx(good));
    if iDebug
      figure(8)
      simplemap(Lat(ix(good)),Lon(ix(good)),zz(ix(good)));
      caxis([220 310]); colorbar
      title(num2str(latbins(ii))); colorbar; pause(1);
    end
  end
end
figure(1); pcolor(double(save_stemp)); colorbar; shading interp;  title('Stemp'); pause(1)

%{
for ii = 1 : length(latbins)-1
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1) & ...
              landfrac < 0.001);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1));
  for tt = 1 : tmax
    zz = cesm_OLR(tt,:,:);
    xx = zz(ix); xx = xx(:);
    good = find(xx > 0);
    save_olr(ii,tt) = nanmean(xx(good));
    if iDebug
      figure(8)
      simplemap(Lat(ix(good)),Lon(ix(good)),zz(ix(good)));
      caxis([220 310]); colorbar
      title(num2str(latbins(ii))); colorbar; pause(1);
    end
  end
end
figure(2); pcolor(double(save_olr)); colorbar; shading interp;  title('OLR'); pause(1)

for ii = 1 : length(latbins)-1
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1) & ...
              landfrac < 0.001);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1));
  for tt = 1 : tmax
    zz = cesm_ClrOLR(tt,:,:);
    xx = zz(ix); xx = xx(:);
    good = find(xx > 0);
    save_clrolr(ii,tt) = nanmean(xx(good));
    if iDebug
      figure(8)
      simplemap(Lat(ix(good)),Lon(ix(good)),zz(ix(good)));
      caxis([220 310]); colorbar
      title(num2str(latbins(ii))); colorbar; pause(1);
    end
  end
end
figure(3); pcolor(double(save_clrolr)); colorbar; shading interp;  title('CLR OLR'); pause(1)

for ii = 1 : length(latbins)-1
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1) & ...
              landfrac < 0.001);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1));
  for tt = 1 : tmax
    zz = cesm_TwetSurf(tt,:,:);
    xx = zz(ix); xx = xx(:);
    good = find(xx > 0);
    save_TWetSurf(ii,tt) = nanmean(xx(good));
    if iDebug
      figure(8)
      simplemap(Lat(ix(good)),Lon(ix(good)),zz(ix(good)));
      caxis([220 310]); colorbar
      title(num2str(latbins(ii))); colorbar; pause(1);
    end
  end
end
figure(4); pcolor(double(save_TWetSurf)); colorbar; shading interp;  title('TWetSurf'); pause(1)

for ii = 1 : length(latbins)-1
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1) & ...
              landfrac < 0.001);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1));
  for tt = 1 : tmax
    zz = cesm_RHSurf(tt,:,:);
    xx = zz(ix); xx = xx(:);
    good = find(xx > 0);
    save_RHSurf(ii,tt) = nanmean(xx(good));
    if iDebug
      figure(8)
      simplemap(Lat(ix(good)),Lon(ix(good)),zz(ix(good)));
      caxis([220 310]); colorbar
      title(num2str(latbins(ii))); colorbar; pause(1);
    end
  end
end
figure(5); pcolor(double(save_RHSurf)); colorbar; shading interp;  title('RH Surf'); pause(1)
%}

for ii = 1 : length(latbins)-1
  fprintf(1,'Q latbin = %3i \n',ii);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1) & ...
              landfrac < 0.001);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1));
  boo = cesm_gas_1;
  donk = size(boo);
  for jj = 1 : donk(2)
    for tt = 1 : tmax
      zz = boo(tt,jj,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save_Q(ii,jj,tt) = nanmean(xx(good));
      if iDebug
        figure(8)
        simplemap(Lat(ix),Lon(ix),zz(ix));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; pause(1);
      end
    end
  end
end
figure(4); pcolor(double(squeeze(save_Q(:,5,:)))); colorbar; shading interp;  title('WV(5)'); pause(1)

for ii = 1 : length(latbins)-1
  fprintf(1,'RH latbin = %3i \n',ii);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1) & ...
              landfrac < 0.001);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1));
  boo = cesm_rh;
  donk = size(boo);
  for jj = 1 : donk(2)
    for tt = 1 : tmax
      zz = boo(tt,jj,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save_RH(ii,jj,tt) = nanmean(xx(good));
      if iDebug
        figure(8)
        simplemap(Lat(ix),Lon(ix),zz(ix));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; pause(1);
      end
    end
  end
end
figure(5); pcolor(double(squeeze(save_RH(:,5,:)))); colorbar; shading interp;  title('RH(5)'); pause(1)

for ii = 1 : length(latbins)-1
  fprintf(1,'T latbin = %3i \n',ii);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1) & ...
              landfrac < 0.001);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1));
  boo = cesm_ptemp;
  donk = size(boo);
  for jj = 1 : donk(2)
    for tt = 1 : tmax
      zz = boo(tt,jj,:,:); zz = squeeze(zz);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save_T(ii,jj,tt) = nanmean(xx(good));
      if iDebug
        figure(8)
        simplemap(Lat(ix),Lon(ix),zz(ix));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; pause(1);
      end
    end
  end
end
figure(6); pcolor(double(squeeze(save_T(:,5,:)))); colorbar; shading interp;  title('T(5)'); pause(1)

for ii = 1 : length(latbins)-1
  fprintf(1,'O3 latbin = %3i \n',ii);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1) & ...
              landfrac < 0.001);
  ix = find(Lat >= latbins(ii) & Lat < latbins(ii+1));
  boo = cesm_gas_3;
  donk = size(boo);
  for jj = 1 : donk(2)
    for tt = 1 : tmax
      zz = boo(tt,jj,:,:); zz = squeeze(zz);
      xx = zz(ix); xx = xx(:);
       good = find(xx > 0);
      save_O3(ii,jj,tt) = nanmean(xx(good));
      if iDebug
        figure(8)
        simplemap(Lat(ix),Lon(ix),zz(ix));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; pause(1);
      end
    end
  end
end
figure(7); pcolor(double(squeeze(save_O3(:,5,:)))); colorbar; shading interp;  title('O3(5)'); pause(1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

latbins = equal_area_spherical_bands(20);    
save_lat = 0.5*(latbins(1:end-1)+latbins(2:end));
Tlevs = nanmean(cesm_plev,1);
Qlevs = nanmean(cesm_plev,1);
saver = ['save /asl/s1/sergio/CESM3/cesm_v7_zonal_rates_' savestr_version '.mat save_* days Tlevs Qlevs latbins'];

eval(saver)

