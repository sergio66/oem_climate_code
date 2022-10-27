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
  simplemap(Lat,Lon,squeeze(cesm_stemp(ii,:,:))); caxis([200 350]); colorbar; colormap(jet); shading interp; pause(1)
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
boo = cesm_stemp;
[tmax,aa,bb] = size(boo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save64x72_stempJUNK = [];

save_lat64x72 = 0.5*(rlat(1:end-1)+rlat(2:end));
save_lon64x72 = 0.5*(rlon(1:end-1)+rlon(2:end));

Tlevs = nanmean(cesm_plev,11);
Qlevs = nanmean(cesm_plev,1);
Tlevs = nanmean(squeeze(cesm_plev(:,:,96,144)),1);
Qlevs = nanmean(squeeze(cesm_plev(:,:,96,144)),1);

saver = ['save /asl/s1/sergio/CESM3/cesm_64x72_rates_' savestr_version '.mat save64x72_olr save64x72_clrolr save64x72_O3 save64x72_CH4 save64x72_CO save64x72_Q save64x72_T save64x72_stemp days Tlevs Qlevs'];
saver = ['save /asl/s1/sergio/CESM3/cesm_64x72_rates_' savestr_version '.mat save64x72_* days Tlevs Qlevs save_l*64*72'];
eval(saver)

clear save64x72_stempJUNK

pause(0.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% to use code "find_set_L3_data_tiles_ii_jj"
Airs_Lat = Lat;
Airs_Lon = wrapTo180(Lon);

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
drlon = 5; 
rlon = -180 : drlon : +180;      %% 73 edges, so 72 bins
rlat = latB2;                    %% 65 edges, so 64 bins

if ~exist('save64x72_stemp')
  disp('doing ST')
  for jj = 1 : length(rlat)-1
    for ii = 1 : length(rlon)-1
      find_set_L3_data_tiles_ii_jj
      for tt = 1 : tmax
        zz = cesm_stemp(tt,:,:);
        xx = zz(ix); xx = xx(:);
        good = find(xx > 0);
        save64x72_stemp(jj,ii,tt) = nanmean(xx(good));
        if iDebug
          figure(3)
          simplemap(Lat(ix(good)),Lon(ix(good)),zz(ix(good)));
          caxis([220 310]); colorbar
          title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
        end
      end
    end
  end
else
  disp('save64x72_stemp already exists, skipping ....')
end

figure(1); pcolor(mean(double(save64x72_stemp),3)); colorbar; colormap(jet); shading interp; pause(1)
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('save64x72_RH')
  boo = cesm_rh;
  donk = size(boo);
  save64x72_RH = zeros(length(rlat)-1,length(rlon)-1,donk(2),tmax);
  for jj = 1 : length(rlat)-1
    fprintf(1,'RH latbin = %3i \n',jj);
    for ii = 1 : length(rlon)-1
      find_set_L3_data_tiles_ii_jj
      for ll = 1 : donk(2)
        for tt = 1 : tmax
          zz = boo(tt,ll,:,:);
          xx = zz(ix); xx = xx(:);
          good = find(xx > 0);
          save64x72_RH(jj,ii,ll,tt) = nanmean(xx(good));
          if iDebug
            figure(3)
            simplemap(Lat(ix),Lon(ix),zz(ix));
            caxis([220 310]); colorbar
            title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
          end
        end
      end
    end
  end
else
  disp('save64x72_RH already exists, skipping ....')
end

figure(2); pcolor(double(squeeze(mean(save64x72_RH(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('save64x72_Q')
  boo = cesm_gas_1;
  donk = size(boo);
  save64x72_Q = zeros(length(rlat)-1,length(rlon)-1,donk(2),tmax);
  for jj = 1 : length(rlat)-1
    fprintf(1,'WV latbin = %3i \n',jj);
    for ii = 1 : length(rlon)-1
      find_set_L3_data_tiles_ii_jj
      for ll = 1 : donk(2)
        for tt = 1 : tmax
          zz = boo(tt,ll,:,:);
          xx = zz(ix); xx = xx(:);
          good = find(xx > 0);
          save64x72_Q(jj,ii,ll,tt) = nanmean(xx(good));
          if iDebug
            figure(3)
            simplemap(Lat(ix),Lon(ix),zz(ix));
            caxis([220 310]); colorbar
            title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
          end
        end
      end
    end
  end
else
  disp('save64x72_Q already exists, skipping ....')
end

figure(2); pcolor(double(squeeze(mean(save64x72_Q(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('save64x72_T')
  boo = cesm_ptemp;
  donk = size(boo);
  save64x72_T = zeros(length(rlat)-1,length(rlon)-1,donk(2),tmax);
  for jj = 1 : length(rlat)-1
    fprintf(1,'T latbin = %3i \n',jj);
    for ii = 1 : length(rlon)-1
      find_set_L3_data_tiles_ii_jj
      for ll = 1 : donk(2)
        for tt = 1 : tmax
          zz = boo(tt,ll,:,:); zz = squeeze(zz);
          xx = zz(ix); xx = xx(:);
          good = find(xx > 0);
          save64x72_T(jj,ii,ll,tt) = nanmean(xx(good));
          if iDebug
            figure(3)
            simplemap(Lat(ix),Lon(ix),zz(ix));
            caxis([220 310]); colorbar
            title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
          end
        end
      end
    end
  end
else
  disp('save64x72_T already exists, skipping ....')
end

figure(3); pcolor(double(squeeze(mean(save64x72_T(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('save64x72_O3')
  boo = cesm_gas_3;
  donk = size(boo);
  save64x72_O3 = zeros(length(rlat)-1,length(rlon)-1,donk(2),tmax);
  for jj = 1 : length(rlat)-1
    fprintf(1,'O3 latbin = %3i \n',jj);
    for ii = 1 : length(rlon)-1
      ix = find(Lat >= rlat(jj) & Lat < rlat(jj+1) & Lon >= rlon(ii) & Lon < rlon(ii+1) & landfrac < 0.001);
      ix = find(Lat >= rlat(jj) & Lat < rlat(jj+1) & Lon >= rlon(ii) & Lon < rlon(ii+1));
      for ll = 1 : donk(2)
        for tt = 1 : tmax
          zz = boo(tt,ll,:,:); zz = squeeze(zz);
          xx = zz(ix); xx = xx(:);
          good = find(xx > 0);
          save64x72_O3(jj,ii,ll,tt) = nanmean(xx(good));
          if iDebug
            figure(3)
            simplemap(Lat(ix),Lon(ix),zz(ix));
            caxis([220 310]); colorbar
            title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
          end
        end
      end
    end
  end
else
  disp('save64x72_O3 already exists, skipping ....')
end

figure(4); pcolor(double(squeeze(mean(save64x72_O3(:,:,5,:),4)))); colorbar; colormap(jet); shading interp; pause(1)
eval(saver)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
disp('doing OLR')
for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    find_set_L3_data_tiles_ii_jj
    for tt = 1 : tmax
      zz = cesm_OLR(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_olr(jj,ii,tt) = nanmean(xx(good));
      if iDebug
        figure(3)
        simplemap(Lat(ix(good)),Lon(ix(good)),zz(ix(good)));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
      end
    end
  end
end
figure(2); pcolor(mean(double(save64x72_olr),3)); colorbar; colormap(jet); shading interp; pause(1)

disp('doing CLR OLR')
for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    find_set_L3_data_tiles_ii_jj
    for tt = 1 : tmax
      zz = cesm_ClrOLR(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_clrolr(jj,ii,tt) = nanmean(xx(good));
      if iDebug
        figure(3)
        simplemap(Lat(ix(good)),Lon(ix(good)),zz(ix(good)));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
      end
    end
  end
end
figure(3); pcolor(mean(double(save64x72_clrolr),3)); colorbar; colormap(jet); shading interp; pause(1)

disp('doing TwetSurf')
for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    find_set_L3_data_tiles_ii_jj
    for tt = 1 : tmax
      zz = cesm_TwetSurf(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_TWetSurf(jj,ii,tt) = nanmean(xx(good));
      if iDebug
        figure(3)
        simplemap(Lat(ix(good)),Lon(ix(good)),zz(ix(good)));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
      end
    end
  end
end
figure(1); pcolor(mean(double(save64x72_TWetSurf),3)); colorbar; colormap(jet); shading interp; pause(1)

disp('doing RH Surf')
for jj = 1 : length(rlat)-1
  for ii = 1 : length(rlon)-1
    find_set_L3_data_tiles_ii_jj
    for tt = 1 : tmax
      zz = cesm_RHSurf(tt,:,:);
      xx = zz(ix); xx = xx(:);
      good = find(xx > 0);
      save64x72_RHSurf(jj,ii,tt) = nanmean(xx(good));
      if iDebug
        figure(3)
        simplemap(Lat(ix(good)),Lon(ix(good)),zz(ix(good)));
        caxis([220 310]); colorbar
        title(num2str(latbins(ii))); colorbar; colormap(jet); shading interp; pause(1);
      end
    end
  end
end
figure(1); pcolor(mean(double(save64x72_RHSurf),3)); colorbar; colormap(jet); shading interp; pause(1)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iDebug
  for ii = 1 : 24
    figure(4); pcolor(double(squeeze(mean(save64x72_T(:,:,ii,:),2)))); colorbar; 
    title(['T ' num2str(ii) ' press (mb) = ' num2str(cesm_PT(ii)) ' mb']); 
    xlabel('time in months since 2002/08');   ylabel('latbin');
    pause(1)
  end

  for ii = 1 : 24
    figure(4); pcolor(double(squeeze(mean(save64x72_O3(:,:,ii,:),2)))); colorbar; 
    title(['O3 ' num2str(ii) ' press (mb) = ' num2str(cesm_PT(ii)) ' mb']); 
    xlabel('time in months since 2002/08');   ylabel('latbin');
    pause(1)
  end

  for ii = 1 : 12
    figure(4); pcolor(double(squeeze(mean(save64x72_Q(:,:,ii,:),2)))); colorbar; 
    title(['WV ' num2str(ii) ' press (mb) = ' num2str(cesm_Q(ii)) ' mb']); 
    xlabel('time in months since 2002/08');   ylabel('latbin');
    pause(1)
  end
end



