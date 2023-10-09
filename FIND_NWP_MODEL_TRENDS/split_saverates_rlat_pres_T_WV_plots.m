figure(41); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.T_z_lat.airsL3); shading interp 
  caxis([- 1 +1]*0.15); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(42); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.T_z_lat.umbc); shading interp 
  caxis([- 1 +1]*0.15); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(43); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.T_z_lat.climcapsL3); shading interp 
  caxis([- 1 +1]*0.15); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(44); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.T_z_lat.merra2); shading interp 
  caxis([- 1 +1]*0.15); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(45); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.T_z_lat.era5); shading interp 
  caxis([- 1 +1]*0.15); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')

figure(46); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.WVfrac_z_lat.airsL3); shading interp 
  caxis([- 1 +1]*0.015); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(47); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.WVfrac_z_lat.umbc); shading interp 
  caxis([- 1 +1]*0.015); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(48); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.WVfrac_z_lat.climcapsL3); shading interp 
  caxis([- 1 +1]*0.015); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(49); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.WVfrac_z_lat.merra2); shading interp 
  caxis([- 1 +1]*0.015); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(50); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.WVfrac_z_lat.era5); shading interp 
  caxis([- 1 +1]*0.015); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')

figure(51); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.RH_z_lat.airsL3); shading interp 
  caxis([- 1 +1]*0.5); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(52); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.RH_z_lat.umbc); shading interp 
  caxis([- 1 +1]*0.5); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(53); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.RH_z_lat.climcapsL3); shading interp 
  caxis([- 1 +1]*0.5); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(54); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.RH_z_lat.merra2); shading interp 
  caxis([- 1 +1]*0.5); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')
figure(55); pcolor(saverates_rlat_pres.trend_rlat64,saverates_rlat_pres.plays100,saverates_rlat_pres.RH_z_lat.era5); shading interp 
  caxis([- 1 +1]*0.5); colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000]); 
  xlabel('Latitude'); ylabel('Pressure (mb)')

iPrint = input('print these ? (-1 [default]/+1) : ');
if length(iPrint) == 0
  iPrint = -1;
end
if iPrint > 0
  dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/PAPER17_TRENDS/Figs/';
  dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
 
  figure(41); aslprint([dir0 '/separateplots_rlat_pres_Tz_airsL3.pdf']);
  figure(42); aslprint([dir0 '/separateplots_rlat_pres_Tz_umbc.pdf']);
  figure(43); aslprint([dir0 '/separateplots_rlat_pres_Tz_climcapsL3.pdf']);
  figure(44); aslprint([dir0 '/separateplots_rlat_pres_Tz_merra2.pdf']);
  figure(45); aslprint([dir0 '/separateplots_rlat_pres_Tz_era5.pdf']);
 
  figure(46); aslprint([dir0 '/separateplots_rlat_pres_WVfrac_airsL3.pdf']);
  figure(47); aslprint([dir0 '/separateplots_rlat_pres_WVfrac_umbc.pdf']);
  figure(48); aslprint([dir0 '/separateplots_rlat_pres_WVfrac_climcapsL3.pdf']);
  figure(49); aslprint([dir0 '/separateplots_rlat_pres_WVfrac_merra2.pdf']);
  figure(50); aslprint([dir0 '/separateplots_rlat_pres_WVfrac_era5.pdf']);
 
  figure(51); aslprint([dir0 '/separateplots_rlat_pres_RH_airsL3.pdf']);
  figure(52); aslprint([dir0 '/separateplots_rlat_pres_RH_umbc.pdf']);
  figure(53); aslprint([dir0 '/separateplots_rlat_pres_RH_climcapsL3.pdf']);
  figure(54); aslprint([dir0 '/separateplots_rlat_pres_RH_merra2.pdf']);
  figure(55); aslprint([dir0 '/separateplots_rlat_pres_RH_era5.pdf']);
end


