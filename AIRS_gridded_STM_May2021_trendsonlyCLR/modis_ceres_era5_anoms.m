modis = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/MODIS_L3_AEROSOL_TRENDS/modis_L3_aerosol_trends_2002_0009_2024_0006.mat');

figure(11); clf; pcolor(2002 + modis.doy2002/365,rlat,squeeze(nanmean(modis.anom_deepblue,1))); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*3); title('MODIS Deep Blue anomaly')
figure(11); clf; pcolor(2002 + modis.doy2002/365,rlat,squeeze(nanmean(modis.anom_od,1))); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*3); title('MODIS AOD anomaly')

aod_anom = reshape(modis.anom_od,4608,262);
aod_anom_cosavg = sum(aod_anom.*(cos(reshape(YY,4608,1)*pi/180)*ones(1,262)))./sum(cos(reshape(YY,4608,1)*pi/180)*ones(1,262));
figure(11); clf; plot(2002 + modis.doy2002/365,aod_anom_cosavg); plotaxis2; xlim([2020 2025]); title('AOD')
pause(1)

deepblue_anom = reshape(modis.anom_deepblue,4608,262);
deepblue_anom_cosavg = sum(deepblue_anom.*(cos(reshape(YY,4608,1)*pi/180)*ones(1,262)))./sum(cos(reshape(YY,4608,1)*pi/180)*ones(1,262));
figure(11); clf; plot(2002 + modis.doy2002/365,deepblue_anom_cosavg); plotaxis2; xlim([2020 2025]); title('Deep Blue')
pause(1)

figure(11); clf; plot(2002 + modis.doy2002/365,aod_anom_cosavg,'k',2002 + modis.doy2002/365,deepblue_anom_cosavg,'b','linewidth',2); 
  plotaxis2; xlim([2020 2025]); hl = legend('AOD','Deep Blue','location','best'); title('MODIS AOD Cosine Avg')

disp('ret to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%

modis = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/MODIS_L3_AEROSOL_TRENDS/modis_L3_cloud_trends_2002_0009_2024_0006.mat');

figure(11); clf; pcolor(2002 + modis.doy2002/365,rlat,squeeze(nanmean(modis.anom_cldfrac,1))); shading interp; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('MODIS CldFrac anomaly')
figure(11); clf; pcolor(2002 + modis.doy2002/365,rlat,squeeze(nanmean(modis.anom_cldtop,1)));  shading interp; colorbar; colormap(usa2); caxis([-1 +1]*10);  title('MODIS CldTop anomaly')

aod_anom = reshape(modis.anom_cldfrac,4608,262);
aod_anom_cosavg = sum(aod_anom.*(cos(reshape(YY,4608,1)*pi/180)*ones(1,262)))./sum(cos(reshape(YY,4608,1)*pi/180)*ones(1,262));
figure(11); clf; plot(2002 + modis.doy2002/365,aod_anom_cosavg); plotaxis2; xlim([2020 2025]); title('Cld Frac')
pause(1)

deepblue_anom = reshape(modis.anom_cldtop,4608,262);
deepblue_anom_cosavg = sum(deepblue_anom.*(cos(reshape(YY,4608,1)*pi/180)*ones(1,262)))./sum(cos(reshape(YY,4608,1)*pi/180)*ones(1,262));
figure(11); clf; plot(2002 + modis.doy2002/365,deepblue_anom_cosavg); plotaxis2; xlim([2020 2025]); title('Cld Top')
pause(1)

figure(11); clf; plot(2002 + modis.doy2002/365,aod_anom_cosavg*1000,'k',2002 + modis.doy2002/365,deepblue_anom_cosavg,'b','linewidth',2); 
  plotaxis2; xlim([2020 2025]); hl = legend('Cld Frac','Cld Top','location','best'); title('MODIS AOD Cosine Avg')

figure(12); clf; 
  baboo = [5 6 7 9]; for ii = 1 : length(baboo); plot(2002+daysSince2002A/365,50*smooth(junk(iaChan(baboo(ii)),:),3),'linewidth',2); hold on; end; 
  plot(2002 + modis.doy2002/365,aod_anom_cosavg*1000,'k',2002 + modis.doy2002/365,deepblue_anom_cosavg,'b','linewidth',2); 
  hold off
  plotaxis2; hl = legend('800 cm-1','900 cm-1','960 cm-1','1231 cm-1','MODIS CldFrac * 1000','MODIS CldTop','location','best','fontsize',10); 
   title('Anomaly : MODIS CldFrac and CldTop \newline Window Channel * 50'); ylabel('no units, mb \newline  [50 * K]')
  xlim([2020 2025]);

disp('ret to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12); clf
ocean_heat_content

junk = btavg_cos;
figure(12); baboo = [5 6 7 9]; for ii = 1 : length(baboo); plot(2002+daysSince2002A/365,100*smooth(junk(iaChan(baboo(ii)),:),3),'linewidth',2); hold on; end; 
                               % plot(oceantime,ocean_0700.month_h22_WO*10-200,'rx-',oceantime,ocean_2000.month_h22_WO*10-200,'bx-'); hold off            
                               plot(oceantime,ocean_0700.month_h22_WO_anom*10-100,'rx-',oceantime,ocean_2000.month_h22_WO_anom*10-150,'bx-'); hold off            
  plotaxis2; hl = legend('800 cm-1','900 cm-1','960 cm-1','1231 cm-1','0-0700','0-2000','location','best','fontsize',10); 
   title('Anomaly : Ocean Heat Content \newline Window Channel *100 '); ylabel('[ZettaJoules (1 ZJ = 10^{21} J)] \newline  [100 * K]')
axis([2020 2025 -100 +100])
axis([2020 2025 -075 +075])
disp('ret to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12)
annual_solar_irradiance_data
disp('ret to continue'); pause

figure(12)
wind_speed_changes
disp('ret to continue'); pause

figure(12);
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/CERES_L3_TRENDS
quick_ceres_flux_cloud_anomaly_look

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coslat = cos(YY*pi/180)*ones(1,iNumAnomTimeSteps);

%freqind = h.vchan(1:2162);
%for ii = 1 : 4608
%  for jj = 1 : 496
%    junk = squeeze(btanomX(ii,:,jj));
%    fluxanomaly(ii,jj) = trapz(freqind,junk);
%  end
%end
if iDorA > 0
  fluxanomaly = fluxanomD;
else
  fluxanomaly = fluxanomA;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yymmSlope = 2020.50;
yymmSlope = 2002.75;
yymm = yyA + (mmA-1)/12;

iChanX = find(h.vchan(ind) >= 1231,1);
iChanX = find(h.vchan(ind) >= 0729,1);
iChanX = find(h.vchan(ind) >= 1419,1);
iChanX = find(h.vchan(ind) >= 1511,1);
for ii = 1 : 4608
  junk = squeeze(btanomX(ii,iChanX,:));
  boo = find(isfinite(junk));
  boo = find(isfinite(junk) & yymm' >= yymmSlope);
  P = polyfit(yymm(boo),junk(boo),1);
  slope_btanomX(ii) = P(1);
end

figure(13); clf
plot(rlat,nanmean(reshape(slope_btanomX,72,64),1)); plotaxis2; 
title(['BT' num2str(h.vchan(ind(iChanX))) ' anomaly --> trends \newline since ' num2str(yymmSlope)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 4608
  junk = fluxanomaly(ii,:);
  boo = find(isfinite(junk));
  boo = find(isfinite(junk) & yymm >= yymmSlope);
  P = polyfit(yymm(boo),junk(boo),1);
  slope_fluxanomaly(ii) = P(1);
end

figure(13); clf
plot(rlat,nanmean(reshape(slope_fluxanomaly,72,64),1),rlat,nanmean(reshape(slope_btanomX,72,64),1)); plotaxis2; 
hl = legend('Flux','Ind Channel','location','best');
title(['BT' num2str(h.vchan(ind(iChanX))) ' anomaly --> Flux Anomaly']); % \newline since ' num2str(yymmSlope)])

figure(14); clf
for ii = 1 : 496
  junk(ii) = sum(fluxanomaly(:,ii).*coslat(:,ii))/sum(coslat(:,ii));
end
junk = fluxanomaly.*coslat;
junk = sum(junk,1)./sum(coslat,1);

plot(yymm,smooth(junk,13)/2,'b',yymm_ceres,smooth(ceres_trend.anom_lw_clr/10,7),'g',yymm_ceres,smooth(ceres_trend.anom_lw/10,7),'r','linewidth',2); 
plotaxis2; title('CosAvg LW/MW IR Flux Anomaly W/m2')
legend('AIRS','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])
