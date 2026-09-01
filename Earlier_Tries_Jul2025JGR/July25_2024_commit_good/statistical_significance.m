addpath /home/sergio/MATLABCODE/COLORMAP
load('/home/sergio/MATLABCODE/COLORMAP/LLS/llsmap5.mat')

load /asl/s1/sergio/JUNK//gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX100_50fatlayers_CLIMCAPS_MERRA2_AMIP6_feedback.mat
load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX3_50fatlayers_AIRSL3_ERA5_CMIP6_feedback.mat

woo = reshape(deltaT,101,72,64); whos woo; woo = squeeze(nanmean(woo,2)); whos woo
figure(1); pcolor(deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(2); pcolor(woo); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(3); pcolor(woo-deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); sum(sum(woo-deltaTlat'))

woo = reshape(deltaTunc,101,72,64); whos woo; woo = squeeze(nanmean(woo,2))/sqrt(72); whos woo
figure(2); pcolor(woo); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02)
woo = woo';

normcdf([-1.96 0 +1.96])
zscore = (abs(deltaTlat) - 0)./woo;
pvalue = 1-normcdf(zscore);
figure(3); pcolor(pvalue'); shading interp; colorbar; set(gca,'ydir','reverse'); 

figure(1); pcolor(deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dT/dt K/yr');
figure(2); pcolor(woo'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dT_{\sigma}/dt K/yr');
figure(3); pcolor(pvalue'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PValue');
