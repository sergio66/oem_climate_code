%% see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/contrast_05_10_15_20.m
load fig5b.mat

%% D (night) trends 
figure(1); clf; pcolor(data_fig5b.vchan,data_fig5b.rlat,data_fig5b.l1c_trendD); shading flat; colorbar; colormap(usa2); xlabel('Wavenumber cm^{-1}'); ylabel('Latitude'); 
  xlim([640 1640]); caxis([-1 +1]*0.1)
  text(1700,95,'K/yr');
  set(gca,'fontsize',14);

%% D (night) unc
figure(2); clf; pcolor(data_fig5b.vchan,data_fig5b.rlat,data_fig5b.l1c_uncD);   shading flat; colorbar; colormap(jet); xlabel('Wavenumber cm^{-1}'); ylabel('Latitude');
  xlim([640 1640]);
  caxis([0 0.06])
  text(1700,95,'K/yr');
  set(gca,'fontsize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see ~/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends --> ~/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/plot_spectral_get_the_model_trends2.m
load fig5c.mat

%% D (night) trends ERA5
figure(3); clf; pcolor(data_fig5c.vchan,data_fig5c.rlat,data_fig5c.era5); shading flat; colorbar; colormap(usa2); xlabel('Wavenumber cm^{-1}'); ylabel('Latitude'); 
  xlim([640 1640]); caxis([-1 +1]*0.1)
  text(1700,95,'K/yr');
  set(gca,'fontsize',14);

%% D (night) unc
figure(4); clf; pcolor(data_fig5b.vchan,data_fig5b.rlat,data_fig5b.l1c_trendD - data_fig5c.obs);   shading flat; colorbar; colormap(jet); xlabel('Wavenumber cm^{-1}'); ylabel('Latitude');
  xlim([640 1640]);
  xlim([640 1640]); caxis([-1 +1]*0.1/100)
  text(1700,95,'K/yr');
  set(gca,'fontsize',14);
  title('Night L1C : fig5b - fig5c')

%{
figure(1); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/new_fig5a_L1Cspectraltrends_bigfont');
figure(2); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/new_fig5b_L1Cspectraltrends_unc_bigfont');
figure(3); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/new_fig5c_ERA5spectraltrends_bigfont');
%}
