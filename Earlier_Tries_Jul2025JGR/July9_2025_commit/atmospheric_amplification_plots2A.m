%% this is same as atmospheric_amplification_plots2.m except we use the "nwp_spectral_trends_cmip6_era5_airsL3_umbc" variable instead of eg cmip6,era5,airsL3 vars
%% and vice versa

epsx = 1e-3;
epsx = 1e-2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); clf
  kapow = reshape(columnSST.umbc,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (deltaT(1:97,:)')./(columnSST.umbc(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('UMBC dT/dSKT K/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(7); clf
  kapow = reshape(columnSST.airsL3,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp(1:97,:)')./(columnSST.airsL3(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('AIRSL3 dT/dSKT K/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(8); clf
  kapow = reshape(columnSST.cmip6,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp(1:97,:)')./(columnSST.cmip6(1:97,:)');
  boo = (cmip6.trend_ptemp(1:97,:)')./(columnSST.cmip6(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('CMIP6 dT/dSKT K/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(9); clf
  kapow = reshape(columnSST.era5,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp(1:97,:)')./(columnSST.era5(1:97,:)');
  boo = (era5.trend_ptemp(1:97,:)')./(columnSST.era5(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA5 dT/dSKT K/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([10 1000]); shading interp; end
for ii = 6 : 9; figure(ii); caxis([-1 +1]*5); ylim([0 1000]); shading interp; set(gca,'yscale','linear'); end
for ii = 6 : 9; figure(ii); caxis([-1 +1]*5); ylim([10 1000]); shading interp; set(gca,'yscale','log'); end
disp('showed dT/dSKT, <ret> to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); clf
  kapow = reshape(columnSST.umbc,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (fracWV(1:97,:)')./(columnSST.umbc(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('UMBC dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(7); clf
  kapow = reshape(columnSST.airsL3,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1(1:97,:)')./(columnSST.airsL3(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('AIRSL3 dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(8); clf
  kapow = reshape(columnSST.cmip6,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.gas_1(1:97,:)')./(columnSST.cmip6(1:97,:)');
  boo = (cmip6.trend_gas_1(1:97,:)')./(columnSST.cmip6(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('CMIP6 dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(9); clf
  kapow = reshape(columnSST.era5,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.gas_1(1:97,:)')./(columnSST.era5(1:97,:)');
  boo = (era5.trend_gas_1(1:97,:)')./(columnSST.era5(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA5 dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([10 1000]); shading interp; end
for ii = 6 : 9; figure(ii); caxis([-1 +1]/5); ylim([0 1000]); shading interp; set(gca,'yscale','linear'); end
for ii = 6 : 9; figure(ii); caxis([-1 +1]/5); ylim([10 1000]); shading interp; set(gca,'yscale','log'); end
disp('showed dfracWV/dSKT, <ret> to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); clf
  kapow = reshape(columnSST.umbc,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.umbc.final(1:97,:)'-compute_deltaRH.umbc.orig(1:97,:)')./(columnSST.umbc(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('UMBC dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(7); clf
  kapow = reshape(columnSST.airsL3,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.airsL3.final(1:97,:)'-compute_deltaRH.airsL3.orig(1:97,:)')./(columnSST.airsL3(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('AIRSL3 dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(8); clf
  kapow = reshape(columnSST.cmip6,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.cmip6.final(1:97,:)'-compute_deltaRH.cmip6.orig(1:97,:)')./(columnSST.cmip6(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('CMIP6 dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(9); clf
  kapow = reshape(columnSST.era5,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.era5.final(1:97,:)'-compute_deltaRH.era5.orig(1:97,:)')./(columnSST.era5(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA5 dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([10 1000]); shading interp; end
for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([0 1000]); shading interp; set(gca,'yscale','linear'); end
for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([10 1000]); shading interp; set(gca,'yscale','log'); end
disp('showed dRH/dSKT, <ret> to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
