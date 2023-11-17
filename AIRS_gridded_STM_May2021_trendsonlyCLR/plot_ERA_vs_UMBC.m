figure(50); clf; aslmap(50,rlat65,rlon73,smoothn((reshape(era5.trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('dST/dt');     caxis([-1 +1]*0.15); colormap(llsmap5); title('ERA5 stemp trends')

%junk = reshape(permute(waterrate_ak0_era5,[3 1 2]),72,64,length(pavg));
%figure(42); pcolor(rlat,pjunk20,squeeze(nanmean(junk,1))'); shading interp;
%caxis([-1 +1]*0.15); set(gca,'ydir','reverse'); colorbar; set(gca,'yscale','linear'); ylim([100 1000]); title('ERA5 raw dWVfrac/dt 1/yr')
%caxis([-1 +1]*0.015); set(gca,'ydir','reverse'); colorbar; set(gca,'yscale','linear'); ylim([100 1000]); title('ERA5 dWVfrac/dt 1/yr')
%colormap(llsmap5)

figure(43); clf; pcolor(rlat,pjunk20,squeeze(nanmean(reshape(waterrate_ak0_era5,72,64,length(pavg)),1))'); shading interp;
caxis([-1 +1]*0.15); set(gca,'ydir','reverse'); colorbar; set(gca,'yscale','linear'); ylim([100 1000]); title('ERA5 raw dWVfrac/dt 1/yr')
caxis([-1 +1]*0.015); set(gca,'ydir','reverse'); colorbar; set(gca,'yscale','linear'); ylim([100 1000]); title('ERA5 dWVfrac/dt 1/yr')
colormap(llsmap5)

figure(44); clf; pcolor(rlat,pjunk20,squeeze(nanmean(reshape(waterrate_akF_era5,72,64,length(pavg)),1))'); shading interp;
caxis([-1 +1]*0.15); set(gca,'ydir','reverse'); colorbar; set(gca,'yscale','linear'); ylim([100 1000]); title('ERA5 * AK dWVfrac/dt 1/yr')
caxis([-1 +1]*0.015); set(gca,'ydir','reverse'); colorbar; set(gca,'yscale','linear'); ylim([100 1000]); title('ERA5 * AK dWVfrac/dt 1/yr')
colormap(llsmap5)

figure(45); clf; pcolor(rlat,pjunk20,squeeze(nanmean(reshape(mean_ak_wv,72,64,length(pavg)),1))'); shading interp;
set(gca,'ydir','reverse'); colorbar; set(gca,'yscale','linear'); ylim([100 1000]); title('AK WV')
colormap(jet); caxis([0 +1]*0.04); caxis([0 +1]*0.5); 

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(46); clf; pcolor(rlat,pjunk20,squeeze(nanmean(reshape(temprate_ak0_era5,72,64,length(pavg)),1))'); shading interp; 
caxis([-1 +1]*0.15); set(gca,'ydir','reverse'); colorbar; set(gca,'yscale','log'); ylim([1 1000]); title('ERA5 raw dT/dt K/yr')
colormap(llsmap5)

figure(47); clf; pcolor(rlat,pjunk20,squeeze(nanmean(reshape(temprate_akF_era5,72,64,length(pavg)),1))'); shading interp; 
caxis([-1 +1]*0.15); set(gca,'ydir','reverse'); colorbar; set(gca,'yscale','log'); ylim([1 1000]); title('ERA5 * AK dT/dt K/yr')
colormap(llsmap5)

figure(48); clf; pcolor(rlat,pjunk20,squeeze(nanmean(reshape(mean_ak_T,72,64,length(pavg)),1))'); shading interp;
set(gca,'ydir','reverse'); colorbar; set(gca,'yscale','log'); ylim([1 1000]); title('AK T')
colormap(jet); caxis([0 +1]*0.05); caxis([0 +1]*0.5); 

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(49); clf; pcolor(rlat,pjunk20,squeeze(nanmean(reshape(mean_ak_o3,72,64,length(pavg)),1))'); shading interp;
set(gca,'ydir','reverse'); colorbar; set(gca,'yscale','log'); ylim([0.1 1000]); title('AK O3')
colormap(jet); caxis([0 +1]*0.025); caxis([0 +1]*0.2); 

%%%%%%%%%%%%%%%%%%%%%%%%%

plot_T_WV_RH_era5_era5AK_umbc
plot_T_WV_RH_era5_xb_umbc

%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rlon')
  rlon73 = -180 : 5 : +180; rlon = meanvaluebin(rlon73);
end

%iCompare = input('Enter latbin over which to compare ERA5 vs UMBC trends (1:64, -1 to stop) : ');
iComparex = input('Enter rlat over which to compare ERA5 vs UMBC trends (-85 : +85, -9999 to stop) : ');
if iComparex > -91
  iCompare = find(rlat > iComparex,1);
else
  iCompare = -1;
end
fprintf(1,'you entered latitude %8.4f which corresponds to latbin %2i \n',iComparex,iCompare)
while iCompare > 0 & iCompare < 65
  do_compare_ERA5_UMBC_latbin
end
