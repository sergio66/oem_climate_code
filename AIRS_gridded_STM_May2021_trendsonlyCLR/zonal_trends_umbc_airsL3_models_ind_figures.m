junk = deltaT .* (ones(101,1) * xmaskLF);
for ii = 1 : length(rlat)
  boo = find(abs(p.rlat - rlat(ii)) < 0.5); findlat(ii) = length(boo);
  xdeltaTlat(ii,:) = nanmean(junk(:,boo),2);
end
figure(28); clf
pcolor(rlat,pavgLAY(1:97,3000),xdeltaTlat(:,1:97)'); 
pcolor(rlat,pavgLAY(1:97,3000),smoothn(xdeltaTlat(:,1:97)',1)); 
shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(xdeltaTlat)); junk = cos(rlat) * ones(1,101);
%area_wgtT = nansum(xdeltaTlat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgtT(1:97)*1000,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([T_Ylim 1000]); ccaxis([-0.15 +0.15],0.1/100); 
  colorbar('horizontal'); colormap(cmap); title(['Zonal d/dt T ' ocbstr ' Quantile' num2str(iQuantile,'%02d')]) %plotaxis2;
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_T_zonal_trends.pdf');

junk = deltaRH .* (ones(100,1) * xmaskLF);
for ii = 1 : length(rlat)
  boo = find(abs(p.rlat - rlat(ii)) < 0.5); findlat(ii) = length(boo);
  xdeltaRHlat(ii,:) = nanmean(junk(:,boo),2);
end
figure(29);
pcolor(rlat,pavgLAY(1:97,3000),xdeltaRHlat(:,1:97)'); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(xdeltaRHlat)); junk = cos(rlat) * ones(1,100);
area_wgtRH = nansum(xdeltaRHlat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgtRH(1:97)*100,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([WV_Ylim 1000]); ccaxis([-0.15 +0.15],0.1/100); colorbar('horizontal'); colormap(cmap); title(['Zonal d/dt RH ' ocbstr ' Quantile' num2str(iQuantile,'%02d')]) %plotaxis2;
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/umbc_rh_zonal_trends.pdf');

junk = fracWV .* (ones(101,1) * xmaskLF);
for ii = 1 : length(rlat)
  boo = find(abs(p.rlat - rlat(ii)) < 0.5); findlat(ii) = length(boo);
  xfracWVlat(ii,:) = nanmean(junk(:,boo),2);
end
figure(30); 
pcolor(rlat,pavgLAY(1:97,1000),xfracWVlat(:,1:97)'); shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log')
junk = zeros(size(xfracWVlat)); junk = cos(rlat) * ones(1,101);
area_wgt_xfracWV = nansum(xfracWVlat.*junk,1)./nansum(junk,1);
%hold on; plot(area_wgt_xfracWV(1:97)*10000,pavgLAY(1:97,1000),'color','r','linewidth',2); hold off
ylim([WV_Ylim 1000]); ccaxis([-2 +2]*1e-3,1e-5); colorbar('horizontal'); %plotaxis2;
title(['Zonal d/dt WVfrac ' ocbstr ' Quantile' num2str(iQuantile,'%02d')]) %plotaxis2;
colormap(cmap)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = 45;

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era.trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.15,0.1/100); ylim([T_Ylim 1000]); title([strNorD ' dT/dt ERA K/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era5.trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.15,0.1/100); ylim([T_Ylim 1000]); title([strNorD ' dT/dt ERA5 K/yr']);

boo = zeros(72,64,24); for ijunk = 1 : 24; boo(:,:,ijunk) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.* airsL3.thestats64x72.ptemprate; junk = squeeze(nanmean(junk,1))'; pcolor(rlat,airsL3.Tlevs,smoothn(junk(1:24,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.15,0.1/100); ylim([T_Ylim 1000]); title([strNorD ' dT/dt AIRS L3 K/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.* reshape(cmip6.trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.15,0.1/100); ylim([T_Ylim 1000]); title(['D/N dT/dt ' mip6str ' K/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era.trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.15,0.1/100); ylim([WV_Ylim 1000]); title([strNorD ' dRH/dt ERA percent/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era5.trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.15,0.1/100); ylim([WV_Ylim 1000]); title([strNorD ' dRH/dt ERA5 percent/yr']);

boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.* airsL3.thestats64x72.RHrate; junk = squeeze(nanmean(junk,1))'; pcolor(rlat,airsL3.Qlevs,smoothn(junk(1:12,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.15,0.1/100); ylim([WV_Ylim 1000]); title([strNorD ' dRH/dt AIRS L3 percent/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(cmip6.trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.15,0.1/100); ylim([WV_Ylim 1000]); title(['D/N dRH/dt ' mip6str ' percent/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era.trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.01,0.01/100); ylim([WV_Ylim 1000]); title([strNorD ' d(fracWV)/dt ERA 1/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(era5.trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.01,0.01/100); ylim([WV_Ylim 1000]); title([strNorD ' d(fracWV)/dt ERA5 1/yr']);

boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.* airsL3.thestats64x72.waterrate; junk = squeeze(nanmean(junk,1))'; pcolor(rlat,airsL3.Qlevs,smoothn(junk(1:12,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.01,0.01/100); ylim([WV_Ylim 1000]); title([strNorD ' d(fracWV)/dt AIRS L3 1/yr']);

boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end
iFig = iFig + 1; figure(iFig); clf;  junk = boo.*reshape(cmip6.trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); pcolor(rlat,pavgLAY(1:97,3000),smoothn(junk(1:97,:),1)); colorbar('southoutside'); colormap(llsmap5)
 shading interp; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ccaxis([-1 +1]*0.01,0.01/100); ylim([WV_Ylim 1000]); title(['D/N d(fracWV)/dt ' mip6str ' 1/yr']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
