disp('look at ../AIRS_gridded_STM_May2021_trendsonlyCLR/do_ceres_vs_umbc_anomalies_zonalavg.m')
disp('look at ../FIND_NWP_MODEL_TRENDS/driver_compare_read64x72flux72_save64zonalflux_ERA5_dRH_anom.m')

disp('first part from  ../FIND_NWP_MODEL_TRENDS/driver_read64x72flux72_save64zonalflux_ERA5_AIRSL3.m')

compare_OLR_from_monthly_2_OLR_from_sarta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%>> whos btanomX fluxanomaly
%  Name              Size                      Bytes  Class     Attributes%%
%  btanomX          64x2645x500            677120000  double
%  fluxanomaly      64x500x14                3584000  double   (:,:,1) === total flux

coslat = cos(rlat*pi/180)*ones(1,iNumAnomTimeSteps);

yymmSlope = 2020.50;
yymmSlope = 2002.75;

% iChanX = find(h.vchan(ind) >= 1231,1);
% iChanX = find(h.vchan(ind) >= 0729,1);
% iChanX = find(h.vchan(ind) >= 1419,1);
% iChanX = find(h.vchan(ind) >= 1511,1);
iChanX = find(f_ind >= 1231,1);
iChanX = find(f_ind >= 0729,1);
iChanX = find(f_ind >= 1419,1);
iChanX = find(f_ind >= 1511,1);
for ii = 1 : 64
  junk = squeeze(btanomX(ii,iChanX,:));
  boo = find(isfinite(junk));
  boo = find(isfinite(junk) & yymm' >= yymmSlope);
  P = polyfit(yymm(boo),junk(boo),1);
  slope_btanomX(ii) = P(1);
end

figure(13); clf
plot(rlat,slope_btanomX); plotaxis2; 
title(['BT' num2str(h.vchan(ind(iChanX))) ' anomaly --> trends \newline since ' num2str(yymmSlope)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk900 = squeeze(btanomX(:,i900,:)) * pi;

for ii = 1 : 64
  junk = squeeze(fluxanomaly(ii,:,1));
  boo = find(isfinite(junk));
  boo = find(isfinite(junk) & yymm >= yymmSlope);
  P = polyfit(yymm(boo),junk(boo),1);
  slope_fluxanomaly(ii) = P(1);

  junk = junk900(ii,:);
  boo = find(isfinite(junk));
  boo = find(isfinite(junk) & yymm >= yymmSlope);
  P = polyfit(yymm(boo),junk(boo),1);
  slope_bt900anomaly(ii) = P(1);

end

figure(13); clf
plot(rlat,slope_fluxanomaly,rlat,slope_btanomX); plotaxis2; 
hl = legend('Flux','Ind Channel','location','best');
title(['BT' num2str(h.vchan(ind(iChanX))) ' anomaly --> Flux Anomaly']); % \newline since ' num2str(yymmSlope)])

figure(14); clf
for ii = 1 : iNumAnomTimeSteps
  junk(ii) = sum(squeeze(fluxanomaly(:,ii,1)).*coslat(:,ii))/sum(coslat(:,ii));
end
junk = squeeze(fluxanomaly(:,:,1)).*coslat;
junk = sum(junk,1)./sum(coslat,1);
%% testing when I did not know conversion of 3.5788
%% plot(yymm,smooth(junk,13)/2,'b',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw_clr/10,7),'g',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw/10,7),'r','linewidth',2);
plot(yymm,smooth(junk,7),'b',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw_clr,7),'g',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw,7),'r','linewidth',2); 
plotaxis2; title('CosAvg LW/MW IR Flux Anomaly W/m2')
legend('AIRS','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

figure(15); 
plot(rlat,slope_bt900anomaly,'k',rlat,slope_fluxanomaly,'b',...
     rlat,nanmean(reshape(ceres_trend.trend_toa_lw_clr_t_4608,72,64)),'g',rlat,nanmean(reshape(ceres_trend.trend_toa_lw_all_4608,72,64)),'r','linewidth',2)
plotaxis2; title('Flux Trends W/m2/yr')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlabel('Latitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear junk900 junkairs
tropics = find(abs(rlat) <= 30);
tropicsYY = find(abs(YY) <= 30);
for ii = 1 : iNumAnomTimeSteps
  junk = squeeze(btanomX(:,i900,:)) * pi;
  junk900(ii) = nanmean(junk(tropics,ii));
  junkairs(ii) = nanmean(squeeze(fluxanomaly(tropics,ii,1)));
  if ii <= 259
    junckceres_clr(ii) = nanmean(ceres_trend.anom_toa_lw_clr_t_4608(tropicsYY,ii));
    junckceres_all(ii) = nanmean(ceres_trend.anom_toa_lw_all_4608(tropicsYY,ii));
  end
end
figure(16); clf
plot(yymm,smooth(junk900,7),'k',yymm,smooth(junkairs,7),'b',ceres_trend.yymm_ceres,smooth(junckceres_clr,7),'g',ceres_trend.yymm_ceres,smooth(junckceres_all,7),'r','linewidth',2); 
plotaxis2; title('Tropical LW/MW IR Flux Anomaly W/m2')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

tropics = find(abs(rlat) <= 30);
tropicsYY = find(abs(YY) <= 30);
for ii = 1 : iNumAnomTimeSteps
  junk = squeeze(btanomX(:,i900,:)) * pi;
  junk900(ii) = nanmean(junk(tropics,ii));
  junkairs(ii) = nanmean(squeeze(fluxanomaly(tropics,ii,1)));
  if ii <= 259
    junckceres_clr(ii) = nanmean(ceres_trend.anom_toa_lw_clr_t_4608(tropicsYY,ii));
    junckceres_all(ii) = nanmean(ceres_trend.anom_toa_lw_all_4608(tropicsYY,ii));
  end
end
figure(17); clf
plot(yymm,smooth(junk900,7),'k',yymm,smooth(junkairs,7),'b',ceres_trend.yymm_ceres,smooth(junckceres_clr,7),'g',ceres_trend.yymm_ceres,smooth(junckceres_all,7),'r','linewidth',2); 
plotaxis2; title('Tropical Ocean LW/MW IR Flux Anomaly W/m2')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

tropics = find(abs(rlat') <= 30);
tropicsYY = find(abs(YY) <= 30);
for ii = 1 : iNumAnomTimeSteps
  junk = squeeze(btanomX(:,i900,:)) * pi;
  junk900(ii) = nanmean(junk(tropics,ii));
  junkairs(ii) = nanmean(squeeze(fluxanomaly(tropics,ii,1)));
  if ii <= 259
    junckceres_clr(ii) = nanmean(ceres_trend.anom_toa_lw_clr_t_4608(tropicsYY,ii));
    junckceres_all(ii) = nanmean(ceres_trend.anom_toa_lw_all_4608(tropicsYY,ii));
  end
end
figure(18); clf
plot(yymm,smooth(junk900,7),'k',yymm,smooth(junkairs,7),'b',ceres_trend.yymm_ceres,smooth(junckceres_clr,7),'g',ceres_trend.yymm_ceres,smooth(junckceres_all,7),'r','linewidth',2); 
plotaxis2; title('Tropical Land LW/MW IR Flux Anomaly W/m2')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

midlat = find(abs(rlat) > 30 & abs(rlat) <= 60);
midlatYY = find(abs(YY) > 30 & abs(YY) <= 60);
for ii = 1 : iNumAnomTimeSteps
  junk = squeeze(btanomX(:,i900,:)) * pi;
  junk900(ii) = nanmean(junk(midlat,ii));
  junkairs(ii) = nanmean(squeeze(fluxanomaly(midlat,ii,1)));
  if ii <= 259
    junckceres_clr(ii) = nanmean(ceres_trend.anom_toa_lw_clr_t_4608(midlatYY,ii));
    junckceres_all(ii) = nanmean(ceres_trend.anom_toa_lw_all_4608(midlatYY,ii));
  end
end
figure(19); clf
plot(yymm,smooth(junk900,7),'k',yymm,smooth(junkairs,7),'b',ceres_trend.yymm_ceres,smooth(junckceres_clr,7),'g',ceres_trend.yymm_ceres,smooth(junckceres_all,7),'r','linewidth',2); 
plotaxis2; title('Midlat LW/MW IR Flux Anomaly W/m2')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

polar = find(abs(rlat) > 60);
polarYY = find(abs(YY) > 60);
for ii = 1 : iNumAnomTimeSteps
  junk = squeeze(btanomX(:,i900,:)) * pi;
  junk900(ii) = nanmean(junk(polar,ii));
  junkairs(ii) = nanmean(squeeze(fluxanomaly(polar,ii,1)));
  if ii <= 259
    junckceres_clr(ii) = nanmean(ceres_trend.anom_toa_lw_clr_t_4608(polarYY,ii));
    junckceres_all(ii) = nanmean(ceres_trend.anom_toa_lw_all_4608(polarYY,ii));
  end
end
figure(20); clf
plot(yymm,smooth(junk900,7),'k',yymm,smooth(junkairs,7),'b',ceres_trend.yymm_ceres,smooth(junckceres_clr,7),'g',ceres_trend.yymm_ceres,smooth(junckceres_all,7),'r','linewidth',2); 
plotaxis2; title('Polar LW/MW IR Flux Anomaly W/m2')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(22); clf
for ii = 1 : iNumAnomTimeSteps
  junk(ii) = sum(squeeze(fluxanomaly(:,ii,1)).*coslat(:,ii))/sum(coslat(:,ii));
end
junk = squeeze(fluxanomaly(:,:,1)).*coslat;
junk = sum(junk,1)./sum(coslat,1);
%% testing when I did not know conversion of 3.5788
%% plot(yymm,smooth(junk,13)/2,'b.-',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw_clr/10,7),'g.-',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw/10,7),'r.-','linewidth',2);
plot(yymm,smooth(junk,15),'b.-',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw_clr,7),'g.-',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw,7),'r.-','linewidth',2); 
hold on
plot(yymm_monthly(1:262),smooth(wahmean_airsL3,7)*3.5788,'c--',...
     yymm_monthly(1:264),smooth(wahmean_era5,7)*3.5788,'m--',...
     yymm_monthly(1:262),smooth(wahmean_airsL3_direct,7),'c',...
     yymm_monthly(1:264),smooth(wahmean_era5_direct,7),'m','linewidth',2); 
hold off;
plotaxis2; title('CosAvg LW/MW IR Flux Anomaly W/m2')
legend('AIRS','CERES LW CLR','CERES LW','AIRS L3 --> SARTA','ERA5 --> SARTA','AIRS L3 direct','ERA5 direct','location','best','fontsize',10);
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])
set(gca,'fontsize',10);

figure(23); clf
rr = 2:14;
  wah = squeeze(nansum(era5_sarta_flux.anomflux(:,rr,:),2));
  coswah = cos(rlat*pi/180)*ones(1,264);
  wahmean_era5 = nansum(wah.*coswah,1)./nansum(coswah,1);
  %%plot(smooth(squeeze(nanmean(era5_sarta_flux.anomflux(:,rr,:),1)),7)*3.5788); plotaxis2;
  %plot(1:264,smooth(wahmean_era5,7)*3.5788); plotaxis2;

  wah = squeeze(nansum(airsL3_sarta_flux.anomflux(:,rr,:),2));
  coswah = cos(rlat*pi/180)*ones(1,262);
  wahmean_airsL3 = nansum(wah.*coswah,1)./nansum(coswah,1);
  %%plot(smooth(squeeze(nanmean(airsL3_sarta_flux.anomflux(:,rr,:),1)),7)*3.5788); plotaxis2;
  %plot(1:262,smooth(wahmean_airsL3,7)*3.5788); plotaxis2;

  junk = squeeze(nansum(fluxanomaly(:,:,rr),3)).*coslat;
  junk = sum(junk,1)./sum(coslat,1);

  plot(yymm,smooth(junk,15),'b.-',yymm_monthly(1:262),smooth(wahmean_airsL3,7)*3.5788,'g--',yymm_monthly(1:264),smooth(wahmean_era5,7)*3.5788,'r--','linewidth',2);
  plotaxis2;
  legend('AIRS L1C','AIRS L3 --> SARTA','ERA5 --> SARTA','location','best','fontsize',10);
  title('Sum(All bands)')
set(gca,'fontsize',10);
  
figure(24);
for rr = 1 : 14
  wah = squeeze(era5_sarta_flux.anomflux(:,rr,:));
  coswah = cos(rlat*pi/180)*ones(1,264);
  wahmean_era5 = nansum(wah.*coswah,1)./nansum(coswah,1);
  %%plot(smooth(squeeze(nanmean(era5_sarta_flux.anomflux(:,rr,:),1)),7)*3.5788); plotaxis2;
  %plot(1:264,smooth(wahmean_era5,7)*3.5788); plotaxis2;

  wah = squeeze(airsL3_sarta_flux.anomflux(:,rr,:));
  coswah = cos(rlat*pi/180)*ones(1,262);
  wahmean_airsL3 = nansum(wah.*coswah,1)./nansum(coswah,1);
  %%plot(smooth(squeeze(nanmean(airsL3_sarta_flux.anomflux(:,rr,:),1)),7)*3.5788); plotaxis2;
  %plot(1:262,smooth(wahmean_airsL3,7)*3.5788); plotaxis2;

  junk = squeeze(fluxanomaly(:,:,rr)).*coslat;
  junk = sum(junk,1)./sum(coslat,1);

  plot(yymm,smooth(junk,15),'b.-',yymm_monthly(1:262),smooth(wahmean_airsL3,7)*3.5788,'g.-',yymm_monthly(1:264),smooth(wahmean_era5,7)*3.5788,'r.-','linewidth',2);
  plotaxis2;
  legend('AIRS L1C','AIRS L3 --> SARTA','ERA5 --> SARTA','location','best','fontsize',10);
  if rr == 1
    title('flux 605-2830 cm-1)')
  else
    title(['RRTM band ' num2str(rr-1) ' : ' num2str(RRTM_bands(rr-1)) '-' num2str(RRTM_bands(rr)) ' cm-1'])
  end
  set(gca,'fontsize',10);
  disp('ret to continue'); pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 14:24
  figure(ii); set(gca,'fontsize',10);
end

