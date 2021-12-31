addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/

disp('THIS IS 64x72 AVERAGES')

for ii=1:10; figure(ii); clf; end;

load llsmap5

%{
if ~exist('reconstruct_L3_spectra_geo2.mat')
  clear all
  driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2
end

if ~exist('reconstruct_era5_spectra_geo2.mat')
  clear all
  driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2
end

if ~exist('reconstruct_cmip6_spectra_geo2.mat')
  clear all
  driver_check_WV_T_RH_CMIP6_geo_and_spectral_rates2
end

L3 = load('reconstruct_L3_spectra_geo2.mat');
era5 = load('reconstruct_era5_spectra_geo2.mat');
cmip6 = load('reconstruct_cmip6_spectra_geo2.mat');
%}

obs = load('gather_tileCLRnight_rates_fits.mat');

check_WV_T_RH_ERA5_geo
check_WV_T_RH_AIRSL3_geo

error('kjglkjdg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
%%want a 2x2 tiled layout
ta = tiledlayout(1,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.25;

tafov(1) = nexttile;
pcolor(L3.rlatx,L3.plevsnwp,L3.rh2d_trendnwp); title('L3 24 levs')
shading interp;

tafov(2) = nexttile;
pcolor(L3.rlatx,L3.plevsx,L3.rh2d_trend(1:97,:)); title('rtp timeseries') 
shading interp;

tafov(3) = nexttile;
pcolor(L3.zonalrlat,L3.zonalplays,L3.zonalRHL3rate);  title('dT/dt,dWV/dt->dRH/dt')
shading interp;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;

set(tafov,'ydir','reverse');
set(tafov,'yscale','log');
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
set(tafov,'Ylim',[100 1000]);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

% Remove all ytick labels except for 1st column
for ii = [2 3]
   tafov(ii).YTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;

%title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(3),'ERA', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);

%% assign colorbar to one tile
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf;
%%want a 2x2 tiled layout
ta = tiledlayout(1,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.25;

tafov(1) = nexttile;
pcolor(era5.rlatx,era5.plevsnwp,era5.rh2d_trendnwp); title('ERA5 37 levs')
shading interp;

tafov(2) = nexttile;
pcolor(era5.rlatx,era5.plevsx,era5.rh2d_trend(1:97,:)); title('rtp timeseries') 
shading interp;

tafov(3) = nexttile;
pcolor(era5.zonalrlat,era5.zonalplays,era5.zonalRHERA5rate);  title('dT/dt,dWV/dt->dRH/dt')
shading interp;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;

set(tafov,'ydir','reverse');
set(tafov,'yscale','log');
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
set(tafov,'Ylim',[100 1000]);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

% Remove all ytick labels except for 1st column
for ii = [2 3]
   tafov(ii).YTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;

%title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(3),'ERA', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);

%% assign colorbar to one tile
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf;
%%want a 2x2 tiled layout
ta = tiledlayout(1,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.25;

tafov(1) = nexttile;
pcolor(cmip6.rlatx,cmip6.plevsnwp,cmip6.rh2d_trendnwp); title('CMIP6 19 levs')
shading interp;

tafov(2) = nexttile;
pcolor(cmip6.rlatx,cmip6.plevsx,cmip6.rh2d_trend(1:97,:)); title('rtp timeseries') 
shading interp;

tafov(3) = nexttile;
pcolor(cmip6.zonalrlat,cmip6.zonalplays,cmip6.zonalRHCMIP6rate);  title('dT/dt,dWV/dt->dRH/dt')
shading interp;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;

set(tafov,'ydir','reverse');
set(tafov,'yscale','log');
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
set(tafov,'Ylim',[100 1000]);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

% Remove all ytick labels except for 1st column
for ii = [2 3]
   tafov(ii).YTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;

%title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(3),'ERA', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);

%% assign colorbar to one tile
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf;
%%want a 2x2 tiled layout
ta = tiledlayout(1,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.15;

tafov(1) = nexttile;
pcolor(L3.rlatx,L3.plevsnwp,L3.t2d_trendnwp); title('L3 24 levs')
shading interp;

tafov(2) = nexttile;
pcolor(L3.rlatx,L3.plevsx,L3.t2d_trend(1:97,:)); title('rtp timeseries') 
shading interp;

tafov(3) = nexttile;
pcolor(L3.zonalrlat,L3.zonalplays,L3.zonalTL3rate);  title('dT/dt,dWV/dt->dT/dt')
shading interp;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;

set(tafov,'ydir','reverse');
set(tafov,'yscale','log');
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
set(tafov,'Ylim',[10 1000]);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

% Remove all ytick labels except for 1st column
for ii = [2 3]
   tafov(ii).YTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;

%title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(3),'ERA', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);

%% assign colorbar to one tile
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5); clf;
%%want a 2x2 tiled layout
ta = tiledlayout(1,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.15;

tafov(1) = nexttile;
pcolor(era5.rlatx,era5.plevsnwp,era5.t2d_trendnwp); title('ERA5 37 levs')
shading interp;

tafov(2) = nexttile;
pcolor(era5.rlatx,era5.plevsx,era5.t2d_trend(1:97,:)); title('rtp timeseries') 
shading interp;

tafov(3) = nexttile;
pcolor(era5.zonalrlat,era5.zonalplays,era5.zonalTERA5rate);  title('dT/dt,dWV/dt->dT/dt')
shading interp;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;

set(tafov,'ydir','reverse');
set(tafov,'yscale','log');
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
set(tafov,'Ylim',[10 1000]);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

% Remove all ytick labels except for 1st column
for ii = [2 3]
   tafov(ii).YTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;

%title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(3),'ERA', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);

%% assign colorbar to one tile
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6); clf;
%%want a 2x2 tiled layout
ta = tiledlayout(1,3);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

cmap = llsmap5;
cmax = 0.15;

tafov(1) = nexttile;
pcolor(cmip6.rlatx,cmip6.plevsnwp,cmip6.t2d_trendnwp); title('CMIP6 19 levs')
shading interp;

tafov(2) = nexttile;
pcolor(cmip6.rlatx,cmip6.plevsx,cmip6.t2d_trend(1:97,:)); title('rtp timeseries') 
shading interp;

tafov(3) = nexttile;
pcolor(cmip6.zonalrlat,cmip6.zonalplays,cmip6.zonalTCMIP6rate);  title('dT/dt,dWV/dt->dT/dt')
shading interp;

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;

set(tafov,'ydir','reverse');
set(tafov,'yscale','log');
set(tafov,'colormap',cmap,'CLim',[-cmax +cmax]);
set(tafov,'Ylim',[10 1000]);
%set(tafov,'shading','interp');

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';

% Remove all ytick labels except for 1st column
for ii = [2 3]
   tafov(ii).YTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;

%title(tafov(1),'UMBC', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(2),'AIRS L3', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);
%title(tafov(3),'ERA', 'Units', 'normalized', 'Position', [0.5, +0.9, 0]);

%% assign colorbar to one tile
cbh = colorbar(tafov(3));
cbh.Layout.Tile = 'south';

ColorOdrDef0 = get(gca,'ColorOrder'); 
ColorOdrDefNew = ColorOdrDef0(2:7,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7);
plot(L3.fchanx,L3.trend,cmip6.fchanx,cmip6.trend,era5.fchanx,era5.trend,L3.fchanx,nanmean(obs.rates,2),'linewidth',2);
%plot(L3.fchanx,L3.xtrend,cmip6.fchanx,cmip6.xtrend,era5.fchanx,era5.xtrend,L3.fchanx,nanmean(obs.rates,2),'linewidth',2);
  plotaxis2; hl = legend('L3','CMIP6','ERA5','AIRS L1 obs','location','best','fontsize',10); ylabel('dBT/dt');
axis([640 1640 -0.1 +0.05])
set(gca,'ColorOrder',ColorOdrDefNew);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8); clf
ta = tiledlayout(1,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
semilogy(L3.t_trend,L3.plevsx,cmip6.t_trend,cmip6.plevsx,era5.t_trend,era5.plevsx,'linewidth',2); 
  plotaxis2; ylim([10 1000]); hl = legend('L3','CMIP6','ERA5','location','best','fontsize',8);
xlim([-0.05 +0.05]);
title('from RTP 100 layers')
set(gca,'ColorOrder',ColorOdrDefNew);

tafov(2) = nexttile;
plot(L3.rh_trend,L3.plevsx,cmip6.rh_trend,cmip6.plevsx,era5.rh_trend,era5.plevsx,'linewidth',2); 
  plotaxis2; ylim([100 1000]); hl = legend('L3','CMIP6','ERA5','location','best','fontsize',8);
xlim([-0.25 +0.25]);
title('from RTP 100 layers')
set(gca,'ColorOrder',ColorOdrDefNew);

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;

set(tafov,'ydir','reverse');
tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;
tafov(1).XLabel.String = 'dT/dt (K/yr)'; tafov(1).YLabel.FontSize = 10;
tafov(2).XLabel.String = 'dRH/dt (percent/yr)'; tafov(1).YLabel.FontSize = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(9); clf
ta = tiledlayout(1,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
semilogy(nanmean(L3.t2d_trendnwp,2),L3.plevsnwp,nanmean(cmip6.t2d_trendnwp,2),cmip6.plevsnwp,nanmean(era5.t2d_trendnwp,2),era5.plevsnwp,'linewidth',2); 
  plotaxis2; ylim([10 1000]); hl = legend('L3','CMIP6','ERA5','location','best','fontsize',8);
xlim([-0.05 +0.05]);
title('from RAW Nlevs')
set(gca,'ColorOrder',ColorOdrDefNew);

tafov(2) = nexttile;
plot(nanmean(L3.rh2d_trendnwp,2),L3.plevsnwp,nanmean(cmip6.rh2d_trendnwp,2),cmip6.plevsnwp,nanmean(era5.rh2d_trendnwp,2),era5.plevsnwp,'linewidth',2); 
  plotaxis2; ylim([100 1000]); hl = legend('L3','CMIP6','ERA5','location','best','fontsize',8);
xlim([-0.25 +0.25]);
title('from RAW Nlevs')
set(gca,'ColorOrder',ColorOdrDefNew);

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;

set(tafov,'ydir','reverse');
tafov(1).YLabel.String = 'P(mb)'; tafov(1).YLabel.FontSize = 10;
tafov(1).XLabel.String = 'dT/dt (K/yr)'; tafov(1).YLabel.FontSize = 10;
tafov(2).XLabel.String = 'dRH/dt (percent/yr)'; tafov(1).YLabel.FontSize = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10); clf
plot(L3.rlatx,L3.st_trend,cmip6.rlatx,cmip6.st_trend,era5.rlatx,era5.st_trend,'linewidth',2);
set(gca,'ColorOrder',ColorOdrDefNew);
%hold on
%  plot(L3.rlatx,L3.bt1231_trend,'--',cmip6.rlatx,cmip6.bt1231_trend,'--',era5.rlatx,era5.bt1231_trend,'--','linewidth',2);
%hold off
  plotaxis2; hl = legend('L3','CMIP6','ERA5','location','best','fontsize',8); title('dST/dt');
