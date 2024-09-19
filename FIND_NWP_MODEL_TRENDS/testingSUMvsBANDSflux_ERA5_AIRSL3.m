%% see driver_read64x72flux72_save64zonalflux_ERA5_AIRSL3.m

%{
for ii = 1 : 6; figure(ii); close; end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('checking AIRSL3 latbin 32')

%% see cluster_do_the_fits_airsL3_ratesv7_tiles_radiances.m
wah = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_latbin_32.mat');
aha = wah.thestatsradtrend64x72.anomflux;
hahaSum = squeeze(sum(aha(:,2:14,:),2));  %% sum over bands
haha0   = squeeze(aha(:,1,:));            %% this should be all

figure(1); pcolor(hahaSum); shading interp; colorbar; colormap(jet);  title('AIRS L3 : all')
  xlael('Time'); ylabel('Lonbin')
figure(2); pcolor(hahaSum-haha0); shading interp; colorbar; colormap(usa2);  title('AIRS L3 : all-sum(bands)')
  xlael('Time'); ylabel('Lonbin')
figure(3); plot(1:262,nansum(haha0,1),'r.-',1:262,nansum(hahaSum,1),'b'); 
  title('sum over 72 lonbins for AIRS L3'); legend('All chans','Sum Over 14 Bands','location','best','fontsize',10);

gah3 = mean(squeeze(mean(aha,1)),2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('checking ERA2 latbin 32')

%% see driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
wah = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5/reconstruct_era5_spectra_geo_rlat32_2002_09_2024_08.mat');
aha = wah.thesave.xanomflux;
hahaSum = squeeze(sum(aha(:,2:14,:),2));  %% sum over bands
haha0   = squeeze(aha(:,1,:));            %% this should be all

figure(4); pcolor(hahaSum); shading interp; colorbar; colormap(jet); title('ERA5 : all')
  xlael('Time'); ylabel('Lonbin')
figure(5); pcolor(hahaSum-haha0); shading interp; colorbar; colormap(usa2); title('ERA5 : all-sum(bands)')
  xlael('Time'); ylabel('Lonbin')
figure(6); plot(1:264,nansum(haha0,1),'r.-',1:264,nansum(hahaSum,1),'b')
  title('sum over 72 lonbins for ERA5'); legend('All chans','Sum Over 14 Bands','location','best','fontsize',10);

figure(1); caxis([-1 +1]*2.5)
figure(4); caxis([-1 +1]*2.5)
figure(3); ylim([-1 +1]*60)
figure(6); ylim([-1 +1]*60)

gah5 = mean(squeeze(mean(aha,1)),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rrtmX = [100 meanvaluebin(wah.thesave.RRTM_bands)];
figure(7); plot(rrtmX,gah3,rrtmX,gah5); plotaxis2; title('sum(anom) : (b) AIRS L3 (r) ERA5');
  xlabel('RRTM band'); ylabel('Mean Anom Flux');
