%% this is like eg ~/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/GRIB/ffill_era5_monthly.m

set_A_or_D_toneeded
  
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

rlat64 = meanvaluebin(rlat65);
rlon72 = meanvaluebin(rlon73);
[X,Y] = ndgrid(rlat64,rlon72);

z.s_longitude = Airs_Lon;
z.s_latitude  = Airs_Lat;

%ii=1; miaow = griddedInterpolant(Airs_Lat,Airs_Lon,squeeze(Airs_STemp(ii,:,:)),'linear');
[iimax,~,~] = size(Airs_STemp);
ii=1; miaow = interpn(Airs_Lat,Airs_Lon,squeeze(Airs_STemp(ii,:,:)),X,Y);
pcolor(Y,X,miaow); colorbar; shading flat; caxis([200 300]); colormap jet;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
tic
disp('doing ST slow method ')
for jj = 1 : length(rlat)-1
  if mod(jj,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end
  for ii = 1 : length(rlon)-1
    find_set_L3_data_tiles_ii_jj
    for tt = 1 : tmax
      zz = Airs_STemp(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      xsave64x72_stemp(jj,ii,tt) = nanmean(xx(good));
    end
  end
end
toc
xsave64x72_stemp(xsave64x72_stemp < 180) = NaN;
figure(1); pcolor(mean(double(xsave64x72_stemp),3)); colorbar; colormap(jet); shading interp; pause(1)
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
disp('doing ST')
for tt = 1 : tmax
  zz = squeeze(Airs_STemp(tt,:,:));
  save64x72_stemp(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
end
toc
save64x72_stemp(save64x72_stemp < 180) = NaN;
figure(2); pcolor(mean(double(save64x72_stemp),3)); colorbar; colormap(jet); shading interp; pause(1)

disp('doing OLR')
for tt = 1 : tmax
  zz = squeeze(Airs_OLR(tt,:,:));
  save64x72_olr(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
end
save64x72_olr(save64x72_olr < 0) = NaN;
figure(2); pcolor(mean(double(save64x72_olr),3)); colorbar; colormap(jet); shading interp; pause(1)

disp('doing CLR OLR')
for tt = 1 : tmax
  zz = squeeze(Airs_ClrOLR(tt,:,:));
  save64x72_clrolr(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
end
save64x72_clrolr(save64x72_clrolr < 0) = NaN;
figure(3); pcolor(mean(double(save64x72_clrolr),3)); colorbar; colormap(jet); shading interp; pause(1)

disp('doing TwetSurf')
for tt = 1 : tmax
  zz = squeeze(Airs_TwetSurf(tt,:,:));
  save64x72_TWetSurf(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
end
save64x72_TWetSurf(save64x72_TWetSurf < 150) = NaN;
figure(1); pcolor(mean(double(save64x72_TWetSurf),3)); colorbar; colormap(jet); shading interp; pause(1)

%% size(Airs_STemp)
%%    240   180   360
disp('doing RH Surf')
for tt = 1 : tmax
  zz = squeeze(Airs_RHSurf(tt,:,:));
  save64x72_RHSurf(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
end
save64x72_RHSurf(save64x72_RHSurf < 0) = NaN;
figure(1); pcolor(mean(double(save64x72_RHSurf),3)); colorbar; colormap(jet); shading interp; pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% size(Airs_RH)
%%   240    12   180   360
boo = Airs_RH;
donk = size(boo);
save64x72_RH = zeros(length(rlat)-1,length(rlon)-1,donk(2),tmax);
disp('RH(z)')
for tt = 1 : tmax
  if mod(tt,100) == 0
    fprintf(1,'+')
  elseif mod(tt,10) == 0
    fprintf(1,'.')
  end
  for ll = 1 : donk(2)
    zz = squeeze(boo(tt,ll,:,:));
    save64x72_RH(:,:,ll,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
  end
end
save64x72_RH(save64x72_RH < 0) = NaN;
figure(2); pcolor(double(squeeze(nanmean(save64x72_RH(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)

boo = Airs_H2OVap;
donk = size(boo);
save64x72_Q = zeros(length(rlat)-1,length(rlon)-1,donk(2),tmax);
disp('Q(z)')
for tt = 1 : tmax
  if mod(tt,100) == 0
    fprintf(1,'+')
  elseif mod(tt,10) == 0
    fprintf(1,'.')
  end
  for ll = 1 : donk(2)
    zz = squeeze(boo(tt,ll,:,:));
    save64x72_Q(:,:,ll,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
  end
end
save64x72_Q(save64x72_Q < 0) = NaN;
figure(2); pcolor(double(squeeze(mean(save64x72_Q(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)

boo = Airs_Temp;
donk = size(boo);
save64x72_T = zeros(length(rlat)-1,length(rlon)-1,donk(2),tmax);
disp('T(z)')
for tt = 1 : tmax
  if mod(tt,100) == 0
    fprintf(1,'+')
  elseif mod(tt,10) == 0
    fprintf(1,'.')
  end
  for ll = 1 : donk(2)
    zz = squeeze(boo(tt,ll,:,:));
    save64x72_T(:,:,ll,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
  end
end
save64x72_T(save64x72_T < 150) = NaN;
figure(3); pcolor(double(squeeze(mean(save64x72_T(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iL3orCLIMCAPS == +1
  boo = Airs_Ozone;
  donk = size(boo);
  save64x72_O3 = zeros(length(rlat)-1,length(rlon)-1,donk(2),tmax);
  disp('O3(z)')
  for tt = 1 : tmax
    if mod(tt,100) == 0
      fprintf(1,'+')
    elseif mod(tt,10) == 0
      fprintf(1,'.')
    end
    for ll = 1 : donk(2)
      zz = squeeze(boo(tt,ll,:,:));
      save64x72_O3(:,:,ll,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
    end
  end
  save64x72_O3(save64x72_O3 < 0) = NaN;
  figure(4); pcolor(double(squeeze(mean(save64x72_O3(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)
else
  boo = Airs_Ozone;
  save64x72_O3 = zeros(length(rlat)-1,length(rlon)-1,tmax);
  disp('O3(col)')
  for tt = 1 : tmax
    zz = squeeze(Airs_RHSurf(tt,:,:));
    save64x72_O3(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
  end
  save64x72_O3(save64x72_O3 < 0) = NaN;
  figure(4); pcolor(double(squeeze(mean(save64x72_O3(:,:,:),3)))); colorbar; colormap(jet); shading interp; pause(1)
end

if iL3orCLIMCAPS == +1
  boo = Airs_CH4;
  donk = size(boo);
  save64x72_CH4 = zeros(length(rlat)-1,length(rlon)-1,donk(2),tmax);
  disp('CH4(z)')
  for tt = 1 : tmax
    if mod(tt,100) == 0
      fprintf(1,'+')
    elseif mod(tt,10) == 0
      fprintf(1,'.')
    end
    for ll = 1 : donk(2)
      zz = squeeze(boo(tt,ll,:,:));
      save64x72_CH4(:,:,ll,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
    end
  end
  save64x72_CH4(save64x72_CH4 < 0) = NaN;
  figure(4); pcolor(double(squeeze(mean(save64x72_CH4(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)
else
  boo = Airs_CH4;
  save64x72_CH4 = zeros(length(rlat)-1,length(rlon)-1,tmax);
  disp('CH4(col)')
  for tt = 1 : tmax
    zz = squeeze(Airs_RHSurf(tt,:,:));
    save64x72_CH4(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
  end
  save64x72_CH4(save64x72_CH4 < 0) = NaN;
  figure(6); pcolor(double(squeeze(mean(save64x72_CH4(:,:,:),3)))); colorbar; colormap(jet); shading interp; pause(1)
end

if iL3orCLIMCAPS == +1
  boo = Airs_CO;
  donk = size(boo);
  donk = size(boo);
  save64x72_CO = zeros(length(rlat)-1,length(rlon)-1,donk(2),tmax);
  disp('CO(z)')
  for tt = 1 : tmax
    if mod(tt,100) == 0
      fprintf(1,'+')
    elseif mod(tt,10) == 0
      fprintf(1,'.')
    end
    for ll = 1 : donk(2)
      zz = squeeze(boo(tt,ll,:,:));
      save64x72_CO(:,:,ll,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
    end
  end
  save64x72_CO(save64x72_CO < 0) = NaN;
  figure(4); pcolor(double(squeeze(mean(save64x72_CO(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)
else
  boo = Airs_CO;
  save64x72_CO = zeros(length(rlat)-1,length(rlon)-1,tmax);
  disp('CO(col)')
  for tt = 1 : tmax
    zz = squeeze(Airs_RHSurf(tt,:,:));
    save64x72_CO(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
  end
  save64x72_CO(save64x72_CO < 0) = NaN;
  figure(6); pcolor(double(squeeze(mean(save64x72_CO(:,:,:),3)))); colorbar; colormap(jet); shading interp; pause(1)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now doing scalars eg iceOD')

for tt = 1 : tmax
  zz = squeeze(Airs_IceT(tt,:,:));
  save64x72_iceT(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
end
save64x72_iceT(save64x72_iceT < 180) = NaN;

for tt = 1 : tmax
  zz = squeeze(Airs_IceSze(tt,:,:));
  save64x72_icesze(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
end
save64x72_icesze(save64x72_icesze < 0) = NaN;

for tt = 1 : tmax
  zz = squeeze(Airs_IceOD(tt,:,:));
  save64x72_ice_od(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
end
save64x72_ice_od(save64x72_ice_od < 0) = NaN;

for tt = 1 : tmax
  zz = squeeze(Airs_LiqWater(tt,:,:));
  save64x72_liq_water(:,:,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
end
save64x72_liq_water(save64x72_liq_water < 0) = NaN;

figure(6); pcolor(mean(double(save64x72_ice_od),3)); colorbar; colormap(jet); shading interp; pause(1)
figure(7); pcolor(mean(double(save64x72_liq_water),3)); colorbar; colormap(jet); shading interp; pause(1)

bonk = size(Airs_CldPres)
for tt = 1 : tmax
  for ll =  1 : bonk(2)
    zz = squeeze(Airs_CldPres(tt,ll,:,:));
    save64x72_cld_pres(:,:,ll,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
  end
end
save64x72_cld_pres(save64x72_cld_pres < 0) = NaN;

bonk = size(Airs_CldFrac)
for tt = 1 : tmax
  for ll =  1 : bonk(2)
    zz = squeeze(Airs_CldFrac(tt,ll,:,:));
    save64x72_cld_frac(:,:,ll,tt) = interpn(Airs_Lat,Airs_Lon,zz,X,Y);
  end
end
save64x72_cld_frac(save64x72_cld_frac < 0) = NaN;

Tlevs = Airs_PT;
Qlevs = Airs_PQ;
if iL3orCLIMCAPS == +1
  if iDorA > 0
    saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_' savestr_version 'desc.mat save64x72_olr save64x72_clrolr save64x72_O3 save64x72_CH4 save64x72_CO save64x72_Q save64x72_T save64x72_stemp days Tlevs Qlevs'];
    saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_' savestr_version '_desc.mat save64x72_* days Tlevs Qlevs save_l*64*72'];
  else
    saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_' savestr_version 'asc.mat save64x72_olr save64x72_clrolr save64x72_O3 save64x72_CH4 save64x72_CO save64x72_Q save64x72_T save64x72_stemp days Tlevs Qlevs'];
    saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_' savestr_version '_asc.mat save64x72_* days Tlevs Qlevs save_l*64*72'];
  end
elseif iL3orCLIMCAPS == -1
  if iDorA > 0
    saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_' savestr_version 'desc.mat save64x72_olr save64x72_clrolr save64x72_O3 save64x72_CH4 save64x72_CO save64x72_Q save64x72_T save64x72_stemp days Tlevs Qlevs'];
    saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_' savestr_version '_desc.mat save64x72_* days Tlevs Qlevs save_l*64*72'];
  else
    saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_' savestr_version 'asc.mat save64x72_olr save64x72_clrolr save64x72_O3 save64x72_CH4 save64x72_CO save64x72_Q save64x72_T save64x72_stemp days Tlevs Qlevs'];
    saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_' savestr_version '_asc.mat save64x72_* days Tlevs Qlevs save_l*64*72'];
  end
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


