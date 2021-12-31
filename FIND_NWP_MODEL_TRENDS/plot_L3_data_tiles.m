%Airs_Temp    = (Airs_Temp_A + Airs_Temp_D)/2;
%Airs_STemp   = (Airs_STemp_A + Airs_STemp_D)/2;
%Airs_H2OVap  = (Airs_H2OVap_A + Airs_H2OVap_D)/2;
%Airs_Ozone   = (Airs_Ozone_A + Airs_Ozone_D)/2;
%Airs_OLR     = (Airs_OLR_A + Airs_OLR_D)/2;
%Airs_Clr_OLR = (Airs_ClrOLR_A + Airs_ClrOLR_D)/2;

if iDorA > 0
  Airs_Temp    = Airs_Temp_D;
  Airs_STemp   = Airs_STemp_D;
  Airs_H2OVap  = Airs_H2OVap_D;
  Airs_RHSurf  = Airs_RHSurf_D;
  Airs_RH      = Airs_RH_D;
  Airs_Ozone   = Airs_Ozone_D;
  Airs_CO      = Airs_CO_D;
  Airs_CH4     = Airs_CH4_D;
  Airs_OLR     = Airs_OLR_D;
  Airs_ClrOLR = Airs_ClrOLR_D;
  
  Airs_LiqWater = Airs_LiqWater_D;
  Airs_IceT     = Airs_IceT_D;
  Airs_IceSze   = Airs_IceSze_D;
  Airs_IceOD    = Airs_IceOD_D;
  Airs_CldPres  = Airs_CldPres_D;
  Airs_CldFrac  = Airs_CldFrac_D;
else
  Airs_Temp    = Airs_Temp_A;
  Airs_STemp   = Airs_STemp_A;
  Airs_H2OVap  = Airs_H2OVap_A;
  Airs_RHSurf  = Airs_RHSurf_A;
  Airs_RH      = Airs_RH_A;
  Airs_Ozone   = Airs_Ozone_A;
  Airs_CO      = Airs_CO_A;
  Airs_CH4     = Airs_CH4_A;
  Airs_OLR     = Airs_OLR_A;
  Airs_ClrOLR = Airs_ClrOLR_A;
  
  Airs_LiqWater = Airs_LiqWater_A;
  Airs_IceT     = Airs_IceT_A;
  Airs_IceSze   = Airs_IceSze_A;
  Airs_IceOD    = Airs_IceOD_A;
  Airs_CldPres  = Airs_CldPres_A;
  Airs_CldFrac  = Airs_CldFrac_A;
end
  
Airs_TwetSurf     = real(get_wet_bulb_temperature(Airs_STemp,Airs_RHSurf));
figure(1);  pcolor(flipud(squeeze(nanmean(double(Airs_RHSurf),1)))); shading interp; colorbar; caxis([0 120]); title('mean RH Surf'); colormap(jet)
figure(2);  pcolor(flipud(squeeze(nanmean(double(Airs_STemp),1))));  shading interp; colorbar; title('mean TSurf');    caxis([200 320]); colormap(jet)
figure(3);  pcolor(flipud(squeeze(nanmean(double(Airs_TwetSurf),1))));   shading interp; colorbar; title('mean TwetSurf'); caxis([200 320]); colormap(jet);

[yy,mm,dd,hh] = matime2yymmddhh(Airs_Date);
days = change2days(yy,mm,dd,2002);
days = unique(days);

yy   = yy(woo);
mm   = mm(woo);
dd   = dd(woo);
hh   = hh(woo);
days = days(woo);

figure(1)
for ii = 1 : 1
  simplemap(Airs_Lat,Airs_Lon,squeeze(Airs_STemp(ii,:,:))); caxis([200 350]); colorbar; colormap(jet); shading interp; pause(1)
end

figure(2)
[salti,landfrac] =  usgs_deg10_dem(Airs_Lat,Airs_Lon);
simplemap(Airs_Lat,Airs_Lon,landfrac);
simplemap(Airs_Lat,Airs_Lon,landfrac,2);

ocean = find(landfrac <= 0.001);
ocean = find(landfrac >= 0.0);
for ii = 1 : 12 : length(days)
  stemp = squeeze(Airs_STemp(ii,:,:));
  figure(2);
  simplemap(Airs_Lat(ocean),Airs_Lon(ocean),stemp(ocean),2);
  title([num2str(yy(ii)) '/' num2str(mm(ii))])
  caxis([220 310]); colorbar
  pause(1)
end

ii = find(yy == min(2012,StopY) & mm == 7);
  stemp = squeeze(Airs_STemp(ii,:,:));
  figure(3);
  simplemap(Airs_Lat(ocean),Airs_Lon(ocean),stemp(ocean),2);
  title([num2str(yy(ii)) '/' num2str(mm(ii))])
  caxis([220 310]); colorbar

iDebug = false;
boo = Airs_STemp;
[tmax,aa,bb] = size(boo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for tt = 1 : tmax
      zz = Airs_STemp(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_stemp(jj,ii,tt) = nanmean(xx(good));
      if iDebug
        figure(3)
        simplemap(Airs_Lat(ix(good)),Airs_Lon(ix(good)),zz(ix(good)));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
      end
    end
  end
end
figure(1); pcolor(mean(double(save64x72_stemp),3)); colorbar; colormap(jet); shading interp; pause(1)

for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for tt = 1 : tmax
      zz = Airs_OLR(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_olr(jj,ii,tt) = nanmean(xx(good));
      if iDebug
        figure(3)
        simplemap(Airs_Lat(ix(good)),Airs_Lon(ix(good)),zz(ix(good)));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
      end
    end
  end
end
figure(2); pcolor(mean(double(save64x72_olr),3)); colorbar; colormap(jet); shading interp; pause(1)

for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for tt = 1 : tmax
      zz = Airs_ClrOLR(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_clrolr(jj,ii,tt) = nanmean(xx(good));
      if iDebug
        figure(3)
        simplemap(Airs_Lat(ix(good)),Airs_Lon(ix(good)),zz(ix(good)));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
      end
    end
  end
end
figure(3); pcolor(mean(double(save64x72_clrolr),3)); colorbar; colormap(jet); shading interp; pause(1)

for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for tt = 1 : tmax
      zz = Airs_TwetSurf(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_TWetSurf(jj,ii,tt) = nanmean(xx(good));
      if iDebug
        figure(3)
        simplemap(Airs_Lat(ix(good)),Airs_Lon(ix(good)),zz(ix(good)));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
      end
    end
  end
end
figure(1); pcolor(mean(double(save64x72_TWetSurf),3)); colorbar; colormap(jet); shading interp; pause(1)

for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for tt = 1 : tmax
      zz = Airs_RHSurf(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_RHSurf(jj,ii,tt) = nanmean(xx(good));
      if iDebug
        figure(3)
        simplemap(Airs_Lat(ix(good)),Airs_Lon(ix(good)),zz(ix(good)));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
      end
    end
  end
end
figure(1); pcolor(mean(double(save64x72_RHSurf),3)); colorbar; colormap(jet); shading interp; pause(1)

boo = Airs_RH;
donk = size(boo);
for jj = 1 : length(rlat)-1
  fprintf(1,'RH latbin = %3i \n',jj);
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for ll = 1 : donk(2)
      for tt = 1 : tmax
        zz = boo(tt,ll,:,:);
        xx = zz(ix); xx = xx(:);
        good = find(xx > 0);
        save64x72_RH(jj,ii,ll,tt) = nanmean(xx(good));
        if iDebug
          figure(3)
          simplemap(Airs_Lat(ix),Airs_Lon(ix),zz(ix));
          caxis([220 310]); colorbar
          title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
        end
      end
    end
  end
end
figure(2); pcolor(double(squeeze(mean(save64x72_RH(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)

boo = Airs_H2OVap;
donk = size(boo);
for jj = 1 : length(rlat)-1
  fprintf(1,'WV latbin = %3i \n',jj);
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for ll = 1 : donk(2)
      for tt = 1 : tmax
        zz = boo(tt,ll,:,:);
        xx = zz(ix); xx = xx(:);
        good = find(xx > 0);
        save64x72_Q(jj,ii,ll,tt) = nanmean(xx(good));
        if iDebug
          figure(3)
          simplemap(Airs_Lat(ix),Airs_Lon(ix),zz(ix));
          caxis([220 310]); colorbar
          title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
        end
      end
    end
  end
end
figure(2); pcolor(double(squeeze(mean(save64x72_Q(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)

boo = Airs_Temp;
donk = size(boo);
for jj = 1 : length(rlat)-1
  fprintf(1,'T latbin = %3i \n',jj);
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for ll = 1 : donk(2)
      for tt = 1 : tmax
        zz = boo(tt,ll,:,:); zz = squeeze(zz);
        xx = zz(ix); xx = xx(:);
        good = find(xx > 0);
        save64x72_T(jj,ii,ll,tt) = nanmean(xx(good));
        if iDebug
          figure(3)
          simplemap(Airs_Lat(ix),Airs_Lon(ix),zz(ix));
          caxis([220 310]); colorbar
          title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
        end
      end
    end
  end
end
figure(3); pcolor(double(squeeze(mean(save64x72_T(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)

boo = Airs_Ozone;
donk = size(boo);
for jj = 1 : length(rlat)-1
  fprintf(1,'O3 latbin = %3i \n',jj);
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for ll = 1 : donk(2)
      for tt = 1 : tmax
        zz = boo(tt,ll,:,:); zz = squeeze(zz);
        xx = zz(ix); xx = xx(:);
        good = find(xx > 0);
        save64x72_O3(jj,ii,ll,tt) = nanmean(xx(good));
        if iDebug
          figure(3)
          simplemap(Airs_Lat(ix),Airs_Lon(ix),zz(ix));
          caxis([220 310]); colorbar
          title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
        end
      end
    end
  end
end
figure(4); pcolor(double(squeeze(mean(save64x72_O3(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)

boo = Airs_CH4;
donk = size(boo);
for jj = 1 : length(rlat)-1
  fprintf(1,'CH4 latbin = %3i \n',jj);
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for ll = 1 : donk(2)
      for tt = 1 : tmax
        zz = boo(tt,ll,:,:); zz = squeeze(zz);
        xx = zz(ix); xx = xx(:);
        good = find(xx > 0);
        save64x72_CH4(jj,ii,ll,tt) = nanmean(xx(good));
        if iDebug
          figure(3)
          simplemap(Airs_Lat(ix),Airs_Lon(ix),zz(ix));
          caxis([220 310]); colorbar
          title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
        end
      end
    end
  end
end
figure(6); pcolor(double(squeeze(mean(save64x72_CH4(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)

boo = Airs_CO;
donk = size(boo);
for jj = 1 : length(rlat)-1
  fprintf(1,'CO latbin = %3i \n',jj);
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for ll = 1 : donk(2)
      for tt = 1 : tmax
        zz = boo(tt,ll,:,:); zz = squeeze(zz);
        xx = zz(ix); xx = xx(:);
        good = find(xx > 0);
        save64x72_CO(jj,ii,ll,tt) = nanmean(xx(good));
        if iDebug
          figure(3)
          simplemap(Airs_Lat(ix),Airs_Lon(ix),zz(ix));
          caxis([220 310]); colorbar
          title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
        end
      end
    end
  end
end
figure(5); pcolor(double(squeeze(mean(save64x72_CO(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now doing scalars eg iceOD')

for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for tt = 1 : tmax
      zz = Airs_IceT(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_iceT(jj,ii,tt) = nanmean(xx(good));
    end
  end
end

for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for tt = 1 : tmax
      zz = Airs_IceSze(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_icesze(jj,ii,tt) = nanmean(xx(good));
    end
  end
end

for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for tt = 1 : tmax
      zz = Airs_IceOD(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_ice_od(jj,ii,tt) = nanmean(xx(good));
    end
  end
end

for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for tt = 1 : tmax
      zz = Airs_LiqWater(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_liq_water(jj,ii,tt) = nanmean(xx(good));
    end
  end
end

figure(6); pcolor(mean(double(save64x72_ice_od),3)); colorbar; colormap(jet); shading interp; pause(1)
figure(7); pcolor(mean(double(save64x72_liq_water),3)); colorbar; colormap(jet); shading interp; pause(1)

boo = Airs_CldPres;
donk = size(boo);
for jj = 1 : length(rlat)-1
  fprintf(1,'CldPress latbin = %3i \n',jj);
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for ll = 1 : donk(2)
      for tt = 1 : tmax
        zz = boo(tt,ll,:,:); zz = squeeze(zz);
        xx = zz(ix); xx = xx(:);
        good = find(xx > 0);
        save64x72_cld_pres(jj,ii,ll,tt) = nanmean(xx(good));
      end
    end
  end
end

boo = Airs_CldFrac;
donk = size(boo);
for jj = 1 : length(rlat)-1
  fprintf(1,'CldFrac latbin = %3i \n',jj);
  for ii = 1 : length(rlon)-1
    ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));
    for ll = 1 : donk(2)
      for tt = 1 : tmax
        zz = boo(tt,ll,:,:); zz = squeeze(zz);
        xx = zz(ix); xx = xx(:);
        good = find(xx > 0);
        save64x72_cld_frac(jj,ii,ll,tt) = nanmean(xx(good));
      end
    end
  end
end

Tlevs = Airs_PT;
Qlevs = Airs_PQ;
if iDorA > 0
  saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_' savestr_version 'desc.mat save64x72_olr save64x72_clrolr save64x72_O3 save64x72_CH4 save64x72_CO save64x72_Q save64x72_T save64x72_stemp days Tlevs Qlevs'];
  saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_' savestr_version '_desc.mat save64x72_* days Tlevs Qlevs save_l*64*72'];
else
  saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_' savestr_version 'asc.mat save64x72_olr save64x72_clrolr save64x72_O3 save64x72_CH4 save64x72_CO save64x72_Q save64x72_T save64x72_stemp days Tlevs Qlevs'];
  saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_' savestr_version '_asc.mat save64x72_* days Tlevs Qlevs save_l*64*72'];
end
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iDebug
  for ii = 1 : 24
    figure(4); pcolor(double(squeeze(mean(save64x72_T(:,:,ii,:),2)))); colorbar; 
    title(['T ' num2str(ii) ' press (mb) = ' num2str(Airs_PT(ii)) ' mb']); 
    xlabel('time in months since 2002/08');   ylabel('latbin');
    pause(1)
  end

  for ii = 1 : 24
    figure(4); pcolor(double(squeeze(mean(save64x72_O3(:,:,ii,:),2)))); colorbar; 
    title(['O3 ' num2str(ii) ' press (mb) = ' num2str(Airs_PT(ii)) ' mb']); 
    xlabel('time in months since 2002/08');   ylabel('latbin');
    pause(1)
  end

  for ii = 1 : 12
    figure(4); pcolor(double(squeeze(mean(save64x72_Q(:,:,ii,:),2)))); colorbar; 
    title(['WV ' num2str(ii) ' press (mb) = ' num2str(Airs_PQ(ii)) ' mb']); 
    xlabel('time in months since 2002/08');   ylabel('latbin');
    pause(1)
  end
end

%% no need to flip
figure(1); plot(squeeze(save64x72_T(18-10:3:18+10,36,:,56)),(Tlevs)); set(gca,'ydir','reverse')
figure(2); plot(squeeze(save64x72_Q(18-10:3:18+10,36,:,56)),(Qlevs)); set(gca,'ydir','reverse')
figure(3); pcolor(squeeze(double(save64x72_stemp(:,36,:)))); colorbar; title('stemp')
figure(4); pcolor(squeeze(double(save64x72_olr(:,36,:)))); colorbar; title('OLR')
figure(5); pcolor(squeeze(double(save64x72_clrolr(:,36,:)))); colorbar; title('Clr OLR')
figure(6); pcolor(squeeze(double(save64x72_ice_od(:,36,:)))); colorbar; title('Ice OD')
figure(7); pcolor(squeeze(double(save64x72_liq_water(:,36,:)))); colorbar; title('liq water')
for ii = 1 : 7
  figure(ii); colormap jet
end


