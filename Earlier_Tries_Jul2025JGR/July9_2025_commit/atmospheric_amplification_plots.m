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
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('UMBC dT/dSKT K/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(7); clf
  kapow = reshape(columnSST.airsL3,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_era_era5_airsL3_umbc.airsL3_100_layertrends.ptemp(1:97,:)')./(columnSST.airsL3(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('AIRSL3 dT/dSKT K/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(8); clf
  kapow = reshape(columnSST.era,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_era_era5_airsL3_umbc.era_100_layertrends.ptemp(1:97,:)')./(columnSST.era(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA dT/dSKT K/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(9); clf
  kapow = reshape(columnSST.era5,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_era_era5_airsL3_umbc.era5_100_layertrends.ptemp(1:97,:)')./(columnSST.era5(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA5 dT/dSKT K/K')
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
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('UMBC dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(7); clf
  kapow = reshape(columnSST.airsL3,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_era_era5_airsL3_umbc.airsL3_100_layertrends.gas_1(1:97,:)')./(columnSST.airsL3(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('AIRSL3 dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(8); clf
  kapow = reshape(columnSST.era,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_era_era5_airsL3_umbc.era_100_layertrends.gas_1(1:97,:)')./(columnSST.era(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(9); clf
  kapow = reshape(columnSST.era5,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_era_era5_airsL3_umbc.era5_100_layertrends.gas_1(1:97,:)')./(columnSST.era5(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA5 dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([10 1000]); shading interp; end
for ii = 6 : 9; figure(ii); caxis([-1 +1]/5); ylim([0 1000]); shading interp; set(gca,'yscale','linear'); end
for ii = 6 : 9; figure(ii); caxis([-1 +1]/5); ylim([10 1000]); shading interp; set(gca,'yscale','log'); end
disp('showed dT/dSKT, <ret> to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); clf
  kapow = reshape(columnSST.umbc,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.umbc.final(1:97,:)'-compute_deltaRH.umbc.orig(1:97,:)')./(columnSST.umbc(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('UMBC dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(7); clf
  kapow = reshape(columnSST.airsL3,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.airsL3.final(1:97,:)'-compute_deltaRH.airsL3.orig(1:97,:)')./(columnSST.airsL3(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('AIRSL3 dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(8); clf
  kapow = reshape(columnSST.era,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.era.final(1:97,:)'-compute_deltaRH.era.orig(1:97,:)')./(columnSST.era(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(9); clf
  kapow = reshape(columnSST.era5,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.era5.final(1:97,:)'-compute_deltaRH.era5.orig(1:97,:)')./(columnSST.era5(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,squeeze(nanmean(boo,2)));; colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA5 dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([10 1000]); shading interp; end
for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([0 1000]); shading interp; set(gca,'yscale','linear'); end
for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([10 1000]); shading interp; set(gca,'yscale','log'); end
disp('showed dRH/dSKT, <ret> to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
