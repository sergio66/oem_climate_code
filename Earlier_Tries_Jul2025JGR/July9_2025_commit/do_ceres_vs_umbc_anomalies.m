%>> whos btanomX fluxanomaly
%  Name              Size                      Bytes  Class     Attributes%%
%  btanomX          4608x444x500               677120000  double
%  fluxanomaly      4608x500x14                3584000  double   (:,:,1) === total flux

coslat = cos(YY*pi/180)*ones(1,iNumAnomTimeSteps);

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

junk900 = squeeze(btanomX(:,i900,:)) * pi;

for ii = 1 : 4608
  junk = fluxanomaly(ii,:);
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
plot(rlat,nanmean(reshape(slope_fluxanomaly,72,64),1),rlat,nanmean(reshape(slope_btanomX,72,64),1)); plotaxis2; 
hl = legend('Flux','Ind Channel','location','best');
title(['BT' num2str(h.vchan(ind(iChanX))) ' anomaly --> Flux Anomaly']); % \newline since ' num2str(yymmSlope)])

figure(14); clf
for ii = 1 : iNumAnomTimeSteps
  junk(ii) = sum(fluxanomaly(:,ii).*coslat(:,ii))/sum(coslat(:,ii));
end
junk = fluxanomaly.*coslat;
junk = sum(junk,1)./sum(coslat,1);
%% testing when I did not know conversion of 3.5788
%% plot(yymm,smooth(junk,13)/2,'b',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw_clr/10,7),'g',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw/10,7),'r','linewidth',2);
plot(yymm,smooth(junk,7),'b',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw_clr,7),'g',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_lw,7),'r','linewidth',2); 
plotaxis2; title('CosAvg LW/MW IR Flux Anomaly W/m2')
legend('AIRS','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

figure(15); 
plot(rlat,nanmean(reshape(slope_bt900anomaly,72,64)),'k',rlat,nanmean(reshape(slope_fluxanomaly,72,64)),'b',...
     rlat,nanmean(reshape(ceres_trend.trend_toa_lw_clr_t_4608,72,64)),'g',rlat,nanmean(reshape(ceres_trend.trend_toa_lw_all_4608,72,64)),'r','linewidth',2)
plotaxis2; title('Flux Trends W/m2/yr')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlabel('Latitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear junk900 junkairs
tropics = find(abs(YY) <= 30);
for ii = 1 : iNumAnomTimeSteps
  junk = squeeze(btanomX(:,i900,:)) * pi;
  junk900(ii) = nanmean(junk(tropics,ii));
  junkairs(ii) = nanmean(fluxanomaly(tropics,ii));
  if ii <= 259
    junckceres_clr(ii) = nanmean(ceres_trend.anom_toa_lw_clr_t_4608(tropics,ii));
    junckceres_all(ii) = nanmean(ceres_trend.anom_toa_lw_all_4608(tropics,ii));
  end
end
figure(16); clf
plot(yymm,smooth(junk900,7),'k',yymm,smooth(junkairs,7),'b',ceres_trend.yymm_ceres,smooth(junckceres_clr,7),'g',ceres_trend.yymm_ceres,smooth(junckceres_all,7),'r','linewidth',2); 
plotaxis2; title('Tropical LW/MW IR Flux Anomaly W/m2')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

tropics = find(abs(YY') <= 30 & p.landfrac <= 0.01);
tropics = find(abs(YY') <= 30);
for ii = 1 : iNumAnomTimeSteps
  junk = squeeze(btanomX(:,i900,:)) * pi;
  junk900(ii) = nanmean(junk(tropics,ii));
  junkairs(ii) = nanmean(fluxanomaly(tropics,ii));
  if ii <= 259
    junckceres_clr(ii) = nanmean(ceres_trend.anom_toa_lw_clr_t_4608(tropics,ii));
    junckceres_all(ii) = nanmean(ceres_trend.anom_toa_lw_all_4608(tropics,ii));
  end
end
figure(17); clf
plot(yymm,smooth(junk900,7),'k',yymm,smooth(junkairs,7),'b',ceres_trend.yymm_ceres,smooth(junckceres_clr,7),'g',ceres_trend.yymm_ceres,smooth(junckceres_all,7),'r','linewidth',2); 
plotaxis2; title('Tropical Ocean LW/MW IR Flux Anomaly W/m2')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

tropics = find(abs(YY') <= 30 & p.landfrac >= 0.99);
tropics = find(abs(YY') <= 30);
for ii = 1 : iNumAnomTimeSteps
  junk = squeeze(btanomX(:,i900,:)) * pi;
  junk900(ii) = nanmean(junk(tropics,ii));
  junkairs(ii) = nanmean(fluxanomaly(tropics,ii));
  if ii <= 259
    junckceres_clr(ii) = nanmean(ceres_trend.anom_toa_lw_clr_t_4608(tropics,ii));
    junckceres_all(ii) = nanmean(ceres_trend.anom_toa_lw_all_4608(tropics,ii));
  end
end
figure(18); clf
plot(yymm,smooth(junk900,7),'k',yymm,smooth(junkairs,7),'b',ceres_trend.yymm_ceres,smooth(junckceres_clr,7),'g',ceres_trend.yymm_ceres,smooth(junckceres_all,7),'r','linewidth',2); 
plotaxis2; title('Tropical Land LW/MW IR Flux Anomaly W/m2')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

midlat = find(abs(YY) > 30 & abs(YY) <= 60);
for ii = 1 : iNumAnomTimeSteps
  junk = squeeze(btanomX(:,i900,:)) * pi;
  junk900(ii) = nanmean(junk(midlat,ii));
  junkairs(ii) = nanmean(fluxanomaly(midlat,ii));
  if ii <= 259
    junckceres_clr(ii) = nanmean(ceres_trend.anom_toa_lw_clr_t_4608(midlat,ii));
    junckceres_all(ii) = nanmean(ceres_trend.anom_toa_lw_all_4608(midlat,ii));
  end
end
figure(19); clf
plot(yymm,smooth(junk900,7),'k',yymm,smooth(junkairs,7),'b',ceres_trend.yymm_ceres,smooth(junckceres_clr,7),'g',ceres_trend.yymm_ceres,smooth(junckceres_all,7),'r','linewidth',2); 
plotaxis2; title('Midlat LW/MW IR Flux Anomaly W/m2')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

polar = find(abs(YY) > 60);
for ii = 1 : iNumAnomTimeSteps
  junk = squeeze(btanomX(:,i900,:)) * pi;
  junk900(ii) = nanmean(junk(polar,ii));
  junkairs(ii) = nanmean(fluxanomaly(polar,ii));
  if ii <= 259
    junckceres_clr(ii) = nanmean(ceres_trend.anom_toa_lw_clr_t_4608(polar,ii));
    junckceres_all(ii) = nanmean(ceres_trend.anom_toa_lw_all_4608(polar,ii));
  end
end
figure(20); clf
plot(yymm,smooth(junk900,7),'k',yymm,smooth(junkairs,7),'b',ceres_trend.yymm_ceres,smooth(junckceres_clr,7),'g',ceres_trend.yymm_ceres,smooth(junckceres_all,7),'r','linewidth',2); 
plotaxis2; title('Polar LW/MW IR Flux Anomaly W/m2')
legend('\pi * AIRS 900cm-1','AIRS FLUX','CERES LW CLR','CERES LW','location','best');
xlim([2002 2025])
xticks([2002 2005:2:2022 2025])

for ii = 14:20
  figure(ii); set(gca,'fontsize',10);
end
