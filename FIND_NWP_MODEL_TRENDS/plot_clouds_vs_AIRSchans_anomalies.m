figure(1); clf
baboo = [5 6 7 9]; for ii = 1 : length(baboo); plot(2002+daysSince2002A/365,50*smooth(junk(iaChan(baboo(ii)),:),3),'linewidth',2); hold on; end;
  plot(2002 + modis_cloud.doy2002/365,aod_anom_cosavg*1000,'rx-',2002 + modis_cloud.doy2002/365,deepblue_anom_cosavg,'bx-','linewidth',4);
  hold off
  ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2);
  plotaxis2; hl = legend('800 cm-1','900 cm-1','960 cm-1','1231 cm-1','MODIS CldFrac * 1000','MODIS CldTop','location','best','fontsize',8);
  %title('Anomaly : MODIS CldFrac and CldTop \newline Window Channel * 50'); ylabel('no units, mb \newline  [50 * K]')
  xlim([2020 2025]);
set(gca,'fontsize',10)

figure(2); clf
baboo = [5 6 7 9]; for ii = 1 : length(baboo); plot(2002+daysSince2002A/365,100*smooth(junk(iaChan(baboo(ii)),:),3),'linewidth',2); hold on; end;
   % plot(ohc.oceantime,ohc.ocean_0700.month_h22_WO*10-200,'rx-',ohc.oceantime,ohc.ocean_2000.month_h22_WO*10-200,'bx-'); hold off
   plot(ohc.oceantime,ohc.ocean_0700.month_h22_WO_anom*10-100,'rx-',ohc.oceantime,ohc.ocean_2000.month_h22_WO_anom*10-150,'bx-','linewidth',4); hold off
  ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2);
  plotaxis2; hl = legend('800 cm-1','900 cm-1','960 cm-1','1231 cm-1','0-0700','0-2000','location','best','fontsize',8);
  %title('Anomaly : Ocean Heat Content \newline Window Channel *100 '); ylabel('[ZettaJoules (1 ZJ = 10^{21} J)] \newline  [100 * K]')
  axis([2020 2025 -100 +100])
  axis([2020 2025 -075 +075])
set(gca,'fontsize',10)

%% make_trends_anomalies_ceres_2022.m
figure(3); clf
plot(ceres_trend.yymm_ceres,smooth(ceres_trend.anom_cldtau,smn)*10,'r',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_cldarea,smn),'b',...
     ceres_trend.yymm_ceres,smooth(ceres_trend.anom_cldtemp,smn),'g',ceres_trend.yymm_ceres,smooth(ceres_trend.anom_cldpress,smn)/10,'k','linewidth',2); hold off
baboo = [5 6 7 9]; for ii = 1 : length(baboo); plot(2002+daysSince2002A/365,100*smooth(junk(iaChan(baboo(ii)),:),3),'linewidth',2); hold on; end;
plot(ceres_trend.yymm_ceres,50*smooth(ceres_trend.anom_cldtau,smn)*10,'rx-',ceres_trend.yymm_ceres,50*smooth(ceres_trend.anom_cldarea,smn),'bx-',...
     ceres_trend.yymm_ceres,50*smooth(ceres_trend.anom_cldtemp,smn),'gx-',ceres_trend.yymm_ceres,50*smooth(ceres_trend.anom_cldpress,smn)/10,'kx-','linewidth',4); hold off
  ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2);
  plotaxis2; hl = legend('800 cm-1','900 cm-1','960 cm-1','1231 cm-1','CERES \tau','CERES cldfrac','CERES cldTemp','CERES cldpress','location','best','fontsize',6);
  % title('Anomaly : CERES cld \newline Window Channel *100 '); ylabel('Mixed Units')
xlim([2020 2025])
set(gca,'fontsize',10)

figure(4); clf
plot(2002 + daysSince2002A/365,radanomD_fatbins.globalavg([i0667 i0723 i1231 i1419],:),'linewidth',2);
bonk = radanomD_fatbins.globalavg([i0667 i0723 i1231 i1419],:);
for ii = 1 : 4
  wonk(ii,:) = smooth(bonk(ii,:),2*smn+1);
end
plot(2002 + daysSince2002A/365,wonk,'linewidth',2);
bonk = reshape(era5_direct_flux.anom_stemp,72,64,264); bonk = squeeze(nanmean(bonk,1)); 
%bonk = nanmean(bonk,1);
bonk = sum(coslat.*bonk,1)./sum(coslat,1);
hold on; plot(yymm,smooth(bonk,smn),'b+-','linewidth',2);
  plotaxis2; 
  xlim([2002.75 2024.5])
   ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2);
  hl = legend('667 cm-1','723 cm-1','1231 cm-1','1419 cm-1','ERA5 SKT','location','best','fontsize',8);
set(gca,'fontsize',10)

%{
figure(1); sergioprint(['/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/anomalies_windowchan_modiscld'])
figure(2); sergioprint(['/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/anomalies_windowchan_ohc'])
figure(3); sergioprint(['/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/anomalies_windowchan_cerescld'])
figure(4); sergioprint(['/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/anomalies_windowchan_era5skt'])
%}
