%% this is same as atmospheric_amplification_plots2.m except we use the "nwp_spectral_trends_cmip6_era5_airsL3_umbc" variable instead of eg cmip6,era5,airsL3 vars
%% and vice versa

%% and here we do smoothing
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
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('UMBC dT/dSKT K/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(7); clf
  kapow = reshape(columnSST.airsL3,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp(1:97,:)')./(columnSST.airsL3(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('AIRSL3 dT/dSKT K/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(8); clf
  kapow = reshape(columnSST.cmip6,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp(1:97,:)')./(columnSST.cmip6(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('CMIP6 dT/dSKT K/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(9); clf
  kapow = reshape(columnSST.era5,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp(1:97,:)')./(columnSST.era5(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA5 dT/dSKT K/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([50 1000]); shading interp; end
for ii = 6 : 9; figure(ii); caxis([-1 +1]*5); ylim([0 1000]); shading interp; set(gca,'yscale','linear'); end
for ii = 6 : 9; figure(ii); caxis([-1 +1]*5); ylim([50 1000]); shading interp; set(gca,'yscale','log'); end
disp('showed dT/dSKT, <ret> to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); clf
  kapow = reshape(columnSST.umbc,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (fracWV(1:97,:)')./(columnSST.umbc(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('UMBC dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(7); clf
  kapow = reshape(columnSST.airsL3,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1(1:97,:)')./(columnSST.airsL3(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('AIRSL3 dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(8); clf
  kapow = reshape(columnSST.cmip6,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.gas_1(1:97,:)')./(columnSST.cmip6(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('CMIP6 dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(9); clf
  kapow = reshape(columnSST.era5,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.gas_1(1:97,:)')./(columnSST.era5(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA5 dfracWV/dSKT 1/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([50 1000]); shading interp; end
for ii = 6 : 9; figure(ii); caxis([-1 +1]/5); ylim([0 1000]); shading interp; set(gca,'yscale','linear'); end
for ii = 6 : 9; figure(ii); caxis([-1 +1]/5); ylim([50 1000]); shading interp; set(gca,'yscale','log'); end
disp('showed dfracWV/dSKT, <ret> to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6); clf
  kapow = reshape(columnSST.umbc,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.umbc.final(1:97,:)'-compute_deltaRH.umbc.orig(1:97,:)')./(columnSST.umbc(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('UMBC dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(7); clf
  kapow = reshape(columnSST.airsL3,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.airsL3.final(1:97,:)'-compute_deltaRH.airsL3.orig(1:97,:)')./(columnSST.airsL3(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('AIRSL3 dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(8); clf
  kapow = reshape(columnSST.cmip6,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.cmip6.final(1:97,:)'-compute_deltaRH.cmip6.orig(1:97,:)')./(columnSST.cmip6(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('CMIP6 dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

figure(9); clf
  kapow = reshape(columnSST.era5,101,72,64);
  kapow = kapow(1:97,:,:);
  kamask = ones(size(kapow));
  kamask(abs(kapow) < epsx) = NaN;
  boo = (compute_deltaRH.era5.final(1:97,:)'-compute_deltaRH.era5.orig(1:97,:)')./(columnSST.era5(1:97,:)');
  boo = boo'; boo = reshape(boo,97,72,64); boo = boo.*kamask;
  pcolor(rlat,playsjunk,smoothn(squeeze(nanmean(boo,2)),1)); colorbar; caxis([-20 +20]); shading flat; xlabel('Latitude'); ylabel('Longitude'); title('ERA5 dRH/dSKT percent/K')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(llsmap5)

for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([50 1000]); shading interp; end
for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([0 1000]); shading interp; set(gca,'yscale','linear'); end
for ii = 6 : 9; figure(ii); caxis([-10 +10]); ylim([50 1000]); shading interp; set(gca,'yscale','log'); end
disp('showed dRH/dSKT, <ret> to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
boo_lps0 = lps0.lapse_othermethod;
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
pjunk = p; pjunk.ptemp(1:101,:) = pjunk.ptemp(1:101,:) + deltaT; [~,junk] = mmwater_rtp(h,pjunk); boo_lpsUMBC = junk.lapse_othermethod;
pjunk = p; pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp; [~,junk] = mmwater_rtp(h,pjunk); boo_lpsAIRSL3 = junk.lapse_othermethod;
pjunk = p; pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp;  [~,junk] = mmwater_rtp(h,pjunk); boo_lpsERA5   = junk.lapse_othermethod;
pjunk = p; pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp;   [~,junk] = mmwater_rtp(h,pjunk); boo_lpsCMIP6  = junk.lapse_othermethod;

figure(6); clf
  aslmap(6,rlat65,rlon73,100*maskLFmatr.*smoothn((reshape(boo_lpsUMBC-boo_lps0,72,64)'),1),[-90 +90],[-180 +180]);   colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/century UMBC');  
figure(7); clf
  aslmap(7,rlat65,rlon73,100*maskLFmatr.*smoothn((reshape(boo_lpsAIRSL3-boo_lps0,72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/century AIRSL3');  
figure(8); clf
  aslmap(8,rlat65,rlon73,100*maskLFmatr.*smoothn((reshape(boo_lpsCMIP6-boo_lps0,72,64)'),1),[-90 +90],[-180 +180]);  colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/century CMIP6');  
figure(9); clf
  aslmap(9,rlat65,rlon73,100*maskLFmatr.*smoothn((reshape(boo_lpsERA5-boo_lps0,72,64)'),1),[-90 +90],[-180 +180]);   colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/century ERA5');  
disp('showed dlapse rate/dt, <ret> to continue'); pause

figure(6); clf
  junk = smoothn(reshape(((boo_lpsUMBC-boo_lps0)./results(:,6)'),72,64)',1);  aslmap(6,rlat65,rlon73,maskLFmatr.*junk,[-90 +90],[-180 +180]);   
  colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/K UMBC');  
figure(7); clf
  junk = smoothn(reshape(((boo_lpsAIRSL3-boo_lps0)./reshape(airsL3.thestats64x72.stemprate,1,72*64)),72,64)',1);  aslmap(7,rlat65,rlon73,maskLFmatr.*junk,[-90 +90],[-180 +180]);   
  colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/K AIRSL3');  
figure(8); clf
  junk = smoothn(reshape(((boo_lpsCMIP6-boo_lps0)./cmip6.trend_stemp),72,64)',1);  aslmap(8,rlat65,rlon73,maskLFmatr.*junk,[-90 +90],[-180 +180]);   
  colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/K CMIP6');  
figure(9); clf
  junk = smoothn(reshape(((boo_lpsERA5-boo_lps0)./era5.trend_stemp),72,64)',1);  aslmap(9,rlat65,rlon73,maskLFmatr.*junk,[-90 +90],[-180 +180]);   
  colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/K ERA5');  
disp('showed dlapse rate/dSKT, <ret> to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

