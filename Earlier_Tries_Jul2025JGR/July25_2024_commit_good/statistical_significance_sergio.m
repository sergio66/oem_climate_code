addpath /home/sergio/MATLABCODE/COLORMAP
load('/home/sergio/MATLABCODE/COLORMAP/LLS/llsmap5.mat')

%load /asl/s1/sergio/JUNK//gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX100_50fatlayers_CLIMCAPS_MERRA2_AMIP6_feedback.mat
load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX3_50fatlayers_AIRSL3_ERA5_CMIP6_feedback.mat

iX = 4;
tvalue = 0.025;
tvalue = 0.050;

%%%%%%%%%%%%%%%%%%%%%%%%%

iX = iX + 1;
xTtrend = reshape(deltaT,101,72,64); whos xTtrend; xTtrend = squeeze(nanmean(xTtrend,2)); whos xTtrend
figure(1); pcolor(deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(2); pcolor(xTtrend); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(3); pcolor(xTtrend-deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); sum(sum(xTtrend-deltaTlat'))

xTunc = reshape(deltaTunc,101,72,64); whos xTunc; xTunc = squeeze(nanmean(xTunc,2))/sqrt(72); whos xTunc
figure(2); pcolor(xTunc); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02)
xTunc = xTunc';

normcdf([-1.96 0 +1.96])
zscore = (abs(deltaTlat) - 0)./xTunc;
pvalueT = 1-normcdf(zscore);
figure(3); pcolor(pvalueT'); shading interp; colorbar; set(gca,'ydir','reverse'); 

figure(iX); pcolor(deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dT/dt K/yr');
figure(2); pcolor(xTunc'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dT_{\sigma}/dt K/yr');
figure(3); pcolor(pvalueT'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueT');

ppvalueT = pvalueT; ppvalueT(ppvalueT > tvalue) = NaN;
ppvalueT = zeros(size(pvalueT)); ppvalueT(pvalueT < tvalue) = 1;
figure(4); pcolor(ppvalueT'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueT');

figure(iX); hold on
for jj = 1 : 101
  for ii = 1 : 64
    if ppvalueT(ii,jj) == 1
      hold on; plot(ii,jj,'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iX); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iX = iX + 1;
xRHtrend = reshape(deltaRH,100,72,64); whos xRHtrend; xRHtrend = squeeze(nanmean(xRHtrend,2)); whos xRHtrend
figure(1); pcolor(deltaRHlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(2); pcolor(xRHtrend); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(3); pcolor(xRHtrend-deltaRHlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); nansum(nansum(xRHtrend-deltaRHlat'))

xRHunc = reshape(deltaRHunc,100,72,64); whos xRHunc; xRHunc = squeeze(nanmean(xRHunc,2))/sqrt(72); whos xRHunc
figure(2); pcolor(xRHunc); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02)
xRHunc = xRHunc';

normcdf([-1.96 0 +1.96])
zscore = (abs(deltaRHlat) - 0)./xRHunc;
pvalueRH = 1-normcdf(zscore);
figure(3); pcolor(pvalueRH'); shading interp; colorbar; set(gca,'ydir','reverse'); 

figure(iX); pcolor(deltaRHlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dRH/dt pc/yr');
figure(2); pcolor(xRHunc'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dRH_{\sigma}/dt pc/yr');
figure(3); pcolor(pvalueRH'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueRH');

ppvalueRH = pvalueRH; ppvalueRH(ppvalueRH > tvalue) = NaN;
ppvalueRH = zeros(size(pvalueRH)); ppvalueRH(pvalueRH < tvalue) = 1;
figure(4); pcolor(ppvalueRH'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueRH');

figure(iX); hold on
for jj = 1 : 100
  for ii = 1 : 64
    if ppvalueRH(ii,jj) == 1
      hold on; plot(ii,jj,'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iX); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iX = iX + 1;
xWVtrend = reshape(fracWV,101,72,64); whos xWVtrend; xWVtrend = squeeze(nanmean(xWVtrend,2)); whos xWVtrend
figure(1); pcolor(fracWVlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01)
figure(2); pcolor(xWVtrend); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01)
figure(3); pcolor(xWVtrend-fracWVlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01); nansum(nansum(xWVtrend-fracWVlat'))

xWVunc = reshape(fracWVunc,101,72,64); whos xWVunc; xWVunc = squeeze(nanmean(xWVunc,2))/sqrt(72); whos xWVunc
figure(2); pcolor(xWVunc); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02)
xWVunc = xWVunc';

normcdf([-1.96 0 +1.96])
zscore = (abs(fracWVlat) - 0)./xWVunc;
pvalueWV = 1-normcdf(zscore);
figure(3); pcolor(pvalueWV'); shading interp; colorbar; set(gca,'ydir','reverse'); 

figure(iX); pcolor(fracWVlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dWV/dt frac/yr');
figure(2); pcolor(xWVunc'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dWV_{\sigma}/dt frac/yr');
figure(3); pcolor(pvalueWV'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueWV');

ppvalueWV = pvalueWV; ppvalueWV(ppvalueWV > tvalue) = NaN;
ppvalueWV = zeros(size(pvalueWV)); ppvalueWV(pvalueWV < tvalue) = 1;
figure(4); pcolor(ppvalueWV'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueWV');

figure(iX); hold on
for jj = 1 : 101
  for ii = 1 : 64
    if ppvalueWV(ii,jj) == 1
      hold on; plot(ii,jj,'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iX); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
