%% see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/contrast_05_10_15_20.m

%% see /home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/gm_redofigs.sc
%% see /umbc/xfs2/strow/asl/s1/sergio/home/git/oem_pkg_run/MATFILES_for_JGR_trends_paper/make_jgr_fig5_colorbarhorizontal.m

load fig5b.mat

addpath /home/sergio/MATLABCODE/COLORMAPS

figure(300); close
junk = figure(300);

junk.Position = [100 100 560*3 420];

pos300 = gcf_location_position(300);

ta = tiledlayout(1,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
%% D (night) trends
 pcolor(data_fig5b.vchan,data_fig5b.rlat,data_fig5b.l1c_trendD); shading flat; colorbar('location','northoutside'); colormap(tafov(1),usa2); xlabel('Wavenumber cm^{-1}'); ylabel('Latitude'); 
  xlim([640 1640]); caxis([-1 +1]*0.1)
  text(500,100,'K/yr','fontsize',14);
  set(gca,'fontsize',14);
shading flat;

tafov(2) = nexttile;
%% D (night) unc
pcolor(data_fig5b.vchan,data_fig5b.rlat,data_fig5b.l1c_uncD);   shading flat; colorbar('location','northoutside'); colormap(tafov(2),jet); xlabel('Wavenumber cm^{-1}'); ylabel('Latitude');
  xlim([640 1640]);
  caxis([0 0.06])
  %text(1700,100,'K/yr');
  set(gca,'fontsize',14);
shading flat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see ~/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends --> ~/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/plot_spectral_get_the_model_trends2.m
load fig5c.mat

%% D (night) trends ERA5
tafov(3) = nexttile;
 pcolor(data_fig5c.vchan,data_fig5c.rlat,data_fig5c.era5); shading flat; colorbar('location','northoutside'); colormap(tafov(3),usa2); xlabel('Wavenumber cm^{-1}'); ylabel('Latitude'); 
  xlim([640 1640]); caxis([-1 +1]*0.1)
  %text(1700,100,'K/yr');
  set(gca,'fontsize',14);
shading flat;

%  % Get rid of all extra space I can
%  ta.Padding = 'none';
%  ta.TileSpacing = 'none';
%  % Get rid of all extra space I can
%  ta.Padding = 'tight';
%  ta.TileSpacing = 'tight';

%%%
%{
 sergioprintfig(['junk_fig5_sept2025'],300,-1,-1);
%}
