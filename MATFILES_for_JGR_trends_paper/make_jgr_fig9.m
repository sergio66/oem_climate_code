load fig9.mat

load llsmap5

figure(1); clf;
load llsmap5
pcolor(fig9_data.noMLS.rlat,fig9_data.noMLS.plays,fig9_data.noMLS.wvtrend); set(gca,'ydir','reverse'); colormap(llsmap5); shading interp; caxis([-1 +1]*0.01); colorbar
ylim([100 1000])
xlabel('Latitude'); ylabel('Pressure (mb)')
set(gca,'fontsize',14)

figure(2); clf;
load llsmap5
pcolor(fig9_data.yesMLS.rlat,fig9_data.yesMLS.plays,fig9_data.yesMLS.wvtrend); set(gca,'ydir','reverse'); colormap(llsmap5); shading interp; caxis([-1 +1]*0.01); colorbar
ylim([100 1000])
xlabel('Latitude'); ylabel('Pressure (mb)')
set(gca,'fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf
disp('Make sure you stretch Fig 3 to be as wide as Fig 1 + Fig 2')
disp('Make sure you stretch Fig 3 to be as wide as Fig 1 + Fig 2')
disp('Make sure you stretch Fig 3 to be as wide as Fig 1 + Fig 2')

ta = tiledlayout(1,2);

  tafov(1) = nexttile;
  pcolor(fig9_data.noMLS.rlat,fig9_data.noMLS.plays,fig9_data.noMLS.wvtrend); set(gca,'ydir','reverse'); colormap(llsmap5); shading interp; caxis([-1 +1]*0.01); 
  ylim([100 1000])
  xlabel('Latitude [deg]'); ylabel('Pressure (mb)')
  set(gca,'fontsize',14)


  tafov(2) = nexttile;
  pcolor(fig9_data.yesMLS.rlat,fig9_data.yesMLS.plays,fig9_data.yesMLS.wvtrend); set(gca,'ydir','reverse'); colormap(llsmap5); shading interp; caxis([-1 +1]*0.01); 
  ylim([100 1000])
  xlabel('Latitude [deg]'); %ylabel('Pressure (mb)')
  set(gca,'fontsize',14)
  
  tafov(2).YTickLabel = ' ';
  ta.Padding = 'none';
  ta.TileSpacing = 'compact';
  ta.Padding = 'compact';
  ta.TileSpacing = 'tight';

cb = colorbar('horizontal');;
cb.Layout.Tile = 'south';
text(-290,1200,'dWVfrac/dt \newline [yr^{-1}]','fontsize',12)
% sergioprintfig('/home/sergio/PAPERS/SUBMITPAPERS/trends_May2025/Figs_NoSmooth/newfig9');
