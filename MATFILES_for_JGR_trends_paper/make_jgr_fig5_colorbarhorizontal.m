%% see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/contrast_05_10_15_20.m
load fig5b.mat

%% D (night) trends 
figure(1); clf; pcolor(data_fig5b.vchan,data_fig5b.rlat,data_fig5b.l1c_trendD); shading flat; colorbar('horizontal'); colormap(usa2); xlabel('Wavenumber cm^{-1}'); ylabel('Latitude'); 
  xlim([640 1640]); caxis([-1 +1]*0.1)
  text(550,-150,'K/yr');
  set(gca,'fontsize',14);
  
%% D (night) unc
figure(2); clf; pcolor(data_fig5b.vchan,data_fig5b.rlat,data_fig5b.l1c_uncD);   shading flat; colorbar('horizontal'); colormap(jet); xlabel('Wavenumber cm^{-1}'); %ylabel('Latitude');
  xlim([640 1640]);
  caxis([0 0.06])
  text(550,-150,'K/yr');
  set(gca,'fontsize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see ~/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends --> ~/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/plot_spectral_get_the_model_trends2.m
load fig5c.mat

%% D (night) trends ERA5
figure(3); clf; pcolor(data_fig5c.vchan,data_fig5c.rlat,data_fig5c.era5); shading flat; colorbar('horizontal'); colormap(usa2); xlabel('Wavenumber cm^{-1}'); %ylabel('Latitude'); 
  xlim([640 1640]); caxis([-1 +1]*0.1)
  text(550,-150,'K/yr');
  set(gca,'fontsize',14);

%% D (night) unc
figure(4); clf; pcolor(data_fig5b.vchan,data_fig5b.rlat,data_fig5b.l1c_trendD - data_fig5c.obs);   shading flat; colorbar('horizontal'); colormap(jet); xlabel('Wavenumber cm^{-1}'); ylabel('Latitude');
  xlim([640 1640]);
  xlim([640 1640]); caxis([-1 +1]*0.1/100)
  text(550,-150,'K/yr');
  set(gca,'fontsize',14);
  title('Night L1C : fig5b - fig5c')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load llsmap5

figure(5); clf
ta = tiledlayout(1,3,'TileSpacing','None', 'Padding','None');
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  %ta.OuterPosition = [0.0375 0.0375 0.875 0.875];
  %ta.OuterPosition = [0.0375 0.0375 0.825 0.825];

tafov(1) = nexttile;
%% D (night) trends 
pcolor(data_fig5b.vchan,data_fig5b.rlat,data_fig5b.l1c_trendD); shading flat; 
  ylabel('Latitude (deg)'); 
  xlim([640 1640]); 
  set(gca,'fontsize',10);
daspect(tafov(1),[5 1 1]);  % <---- move to after-plot
pbaspect(tafov(1),[5 1 1]); % <---- move to after-plot
  text(400,+117,'K/yr');
  cb = colorbar('horizontal','northoutside'); colormap(tafov(1),usa2); caxis([-1 +1]*0.1); %xlabel('Wavenumber cm^{-1}'); 
  %cb = colorbar; cb.Layout.Tile = 'south';
%  pos = cb.Position; pos(4) = 0.95 * pos(4); cb.Position = pos;
cb;

tafov(2) = nexttile;
%% D (night) unc
pcolor(data_fig5b.vchan,data_fig5b.rlat,data_fig5b.l1c_uncD);   shading flat; xlabel('Wavenumber (cm^{-1})'); %ylabel('Latitude');
  xlim([640 1640]);
  %text(550,-150,'K/yr');
  set(gca,'fontsize',10);
daspect(tafov(2),[5 1 1]);  % <---- move to after-plot
pbaspect(tafov(2),[5 1 1]); % <---- move to after-plot
  %cb = colorbar; cb.Layout.Tile = 'south';
  cb = colorbar('horizontal','northoutside'); colormap(tafov(2),jet); caxis([0 0.06])
%  pos = cb.Position; pos(4) = 0.95 * pos(4); cb.Position = pos;
cb;

tafov(3) = nexttile;
%% D (night) trends ERA5
pcolor(data_fig5c.vchan,data_fig5c.rlat,data_fig5c.era5); shading flat; %xlabel('Wavenumber cm^{-1}'); %ylabel('Latitude'); 
  xlim([640 1640]);
  %text(550,-150,'K/yr');
  set(gca,'fontsize',10);
daspect(tafov(3),[5 1 1]);  % <---- move to after-plot
pbaspect(tafov(3),[5 1 1]); % <---- move to after-plot
  %cb = colorbar; cb.Layout.Tile = 'south';
  cb = colorbar('horizontal','northoutside'); colormap(tafov(3),usa2); caxis([-1 +1]*0.1)
%  pos = cb.Position; pos(4) = 0.95 * pos(4); cb.Position = pos;
cb;

  ta.Padding = 'tight';
  ta.TileSpacing = 'tight';

  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';

  %ta.Padding = 'loose';
  %ta.TileSpacing = 'loose';

% Remove all ytick labels except for 1st column
for ii = [2 3]
   tafov(ii).YTickLabel = '';
   tafov(ii).YLabel.String = [];
end

%p = get(gcf,'position');
%set(gcf,'position',[p(1) p(2) 1000 500]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
figure(1); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/new_fig5a_L1Cspectraltrends_bigfont');
figure(2); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/new_fig5b_L1Cspectraltrends_unc_bigfont');
figure(3); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/new_fig5c_ERA5spectraltrends_bigfont');
figure(5); sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/new_fig5_montage_spectraltrends_bigfont');
%}
