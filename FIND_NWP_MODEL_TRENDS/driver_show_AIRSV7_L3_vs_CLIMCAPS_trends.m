addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load llsmap5

for ii = 1 : 6;
  figure(ii); clf
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iX = input('fastgrib (+1) or slow calc (-1) : ');
if iX > 0
  load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat
else
  load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug2022_20yr_desc.mat
end

figure(1); 
pcolor(thestats64x72.lats,Tlevs,squeeze(nanmean(thestats64x72.ptemprate,1))'); colormap(llsmap5); shading interp; colorbar;
set(gca,'ydir','reverse'); caxis([-1 +1]*0.15); title('AIRS L3 20 years dT/dt'); set(gca,'yscale','log'); ylim([10 1000])

figure(2)
pcolor(thestats64x72.lats,Qlevs,squeeze(nanmean(thestats64x72.waterrate,1))'); colormap(llsmap5); shading interp; colorbar; 
set(gca,'ydir','reverse'); caxis([-1 +1]*0.015); title('AIRS L3 20 years dWVfrac/dt'); ylim([100 1000])

figure(3)
pcolor(thestats64x72.lats,Qlevs,squeeze(nanmean(thestats64x72.RHrate,1))'); colormap(llsmap5); shading interp; colorbar; 
set(gca,'ydir','reverse'); caxis([-1 +1]*0.5); title('AIRS L3 20 years dRH/dt'); ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat

figure(4); 
pcolor(thestats64x72.lats,Tlevs/100,squeeze(nanmean(thestats64x72.ptemprate,1))'); colormap(llsmap5); shading interp; colorbar;
set(gca,'ydir','reverse'); caxis([-1 +1]*0.15); title('CLIMCAPS 20 years dT/dt'); set(gca,'yscale','log'); ylim([10 1000])

figure(5)
pcolor(thestats64x72.lats,Qlevs/100,squeeze(nanmean(thestats64x72.waterrate,1))'); colormap(llsmap5); shading interp; colorbar; 
set(gca,'ydir','reverse'); caxis([-1 +1]*0.015); title('CLIMCAPS 20 years dWVfrac/dt'); ylim([100 1000])

figure(6)
pcolor(thestats64x72.lats,Qlevs/100,100*squeeze(nanmean(thestats64x72.RHrate,1))'); colormap(llsmap5); shading interp; colorbar; 
set(gca,'ydir','reverse'); caxis([-1 +1]*0.5); title('CLIMCAPS 20 years dRH/dt'); ylim([100 1000])
