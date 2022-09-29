addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load llsmap5

for ii = 1 : 64
  fin = ['SimulateTimeSeries/ERA5_ConstTracegas/reconstruct_era5_const_tracegas_spectra_geo_rlat' num2str(ii,'%02d') '_2014_09_2021_08.mat'];
  a = load(fin);
  zonalT1(:,ii) = nanmean(a.zonalTERA5rate,2);
  zonalT2(:,ii) = nanmean(a.thesave.t2d_xtrend,2);
  zonalWV2(:,ii) = nanmean(a.thesave.wv2d_xtrend,2);
  zonalOZ2(:,ii) = nanmean(a.thesave.oz2d_xtrend,2);
  semilogy(a.t_xconstr(1:97),a.plevsx); set(gca,'ydir','reverse'); ylim([0.1 1000]); 
end

figure(1); clf; colormap(llsmap5); pcolor(a.zonalrlat,a.plevsx,zonalT2(1:97,:)); set(gca,'ydir','reverse'); ylim([1 1000]);    set(gca,'yscale','log'); shading interp; colorbar; caxis([-1 +1]*0.15)
figure(2); clf; colormap(llsmap5); pcolor(a.zonalrlat,a.plevsx,zonalWV2(1:97,:)); set(gca,'ydir','reverse'); ylim([100 1000]); set(gca,'yscale','log'); shading interp; colorbar; caxis([-1 +1]*0.015)
figure(3); clf; colormap(llsmap5); pcolor(a.zonalrlat,a.plevsx,zonalOZ2(1:97,:)); set(gca,'ydir','reverse'); ylim([0.1 1000]); set(gca,'yscale','log'); shading interp; colorbar; caxis([-1 +1]*0.015)
