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
