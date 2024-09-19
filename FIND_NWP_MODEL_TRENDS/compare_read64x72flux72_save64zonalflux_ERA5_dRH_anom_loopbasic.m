for ii = 1 : 14
  %% easy
  % junk0 = squeeze(nanmean(a0.anomflux(:,ii,:),1)) * 3.5788; 
  % junk1 = squeeze(nanmean(a1.anomflux(:,ii,:),1)) * 3.5788;
  % junk2 = squeeze(nanmean(a2.anomflux(:,ii,:),1)) * 3.5788;
  % junk3 = squeeze(nanmean(a3.anomflux(:,ii,:),1)) * 3.5788;
  % junk4 = squeeze(nanmean(a4.anomflux(:,ii,:),1)) * 3.5788;
  % junk5 = squeeze(nanmean(a5.anomflux(:,ii,:),1)) * 3.5788;

  junk0 = squeeze(a0.anomflux(:,ii,:)) * 3.5788; 
  junk1 = squeeze(a1.anomflux(:,ii,:)) * 3.5788;
  junk2 = squeeze(a2.anomflux(:,ii,:)) * 3.5788;
  junk3 = squeeze(a3.anomflux(:,ii,:)) * 3.5788;
  junk4 = squeeze(a4.anomflux(:,ii,:)) * 3.5788;
  junk5 = squeeze(a5.anomflux(:,ii,:)) * 3.5788;
  junk6 = squeeze(a6.anomflux(:,ii,:)) * 3.5788;

  junk0trp = nanmean(junk0(tropics,:),1);
  junk1trp = nanmean(junk1(tropics,:),1);
  junk2trp = nanmean(junk2(tropics,:),1);
  junk3trp = nanmean(junk3(tropics,:),1);
  junk4trp = nanmean(junk4(tropics,:),1);
  junk5trp = nanmean(junk5(tropics,:),1);
  junk6trp = nanmean(junk6(tropics,:),1);

  junk0 = sum(junk0.*coslat)./sum(coslat);
  junk1 = sum(junk1.*coslat)./sum(coslat);
  junk2 = sum(junk2.*coslat)./sum(coslat);
  junk3 = sum(junk3.*coslat)./sum(coslat);
  junk4 = sum(junk4.*coslat)./sum(coslat);
  junk5 = sum(junk5.*coslat)./sum(coslat);
  junk6 = sum(junk6.*coslat)./sum(coslat);

  if ii > 1
    figure(1); plot(yymm,smooth(junk0,smn),'kx-',yymm,smooth(junk1,smn),'b',yymm,smooth(junk2,smn),'g',yymm,smooth(junk3,smn),'r',yymm,smooth(junk4,smn),'c',yymm,smooth(junk5,smn),'ms-','linewidth',2); 
      plotaxis2; ylabel('Flux [W/m2]'); title([num2str(ii) ' : ' num2str(RRTM_bands(ii-1)) ' to ' num2str(RRTM_bands(ii))]);
      ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('ALL','Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)
    figure(2); plot(yymm,smooth(junk0-junk1,smn),'b',yymm,smooth(junk0-junk2,smn),'g',yymm,smooth(junk0-junk3,smn),'r',yymm,smooth(junk0-junk4,smn),'c',yymm,smooth(junk0-junk5,smn),'ms-','linewidth',2); 
      plotaxis2; ylabel('\delta Flux : ClrSky-X [W/m2]'); title([num2str(ii) ' : ' num2str(RRTM_bands(ii-1)) ' to ' num2str(RRTM_bands(ii))]); 
      ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)

    figure(3); plot(yymm,smooth(junk0,smn),'kx-',yymm,smooth(junk1,smn),'b',yymm,smooth(junk2,smn),'g',yymm,smooth(junk3,smn),'r',yymm,smooth(junk4,smn),'c',yymm,smooth(junk5,smn),'ms-','linewidth',2); 
      plotaxis2; ylabel('Flux [W/m2]'); title([num2str(ii) ' : ' num2str(RRTM_bands(ii-1)) ' to ' num2str(RRTM_bands(ii))]);
      xlim([2020 2025]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('ALL','Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)
    figure(4); plot(yymm,smooth(junk0-junk1,smn),'b',yymm,smooth(junk0-junk2,smn),'g',yymm,smooth(junk0-junk3,smn),'r',yymm,smooth(junk0-junk4,smn),'c',yymm,smooth(junk0-junk5,smn),'ms-','linewidth',2); 
      plotaxis2; ylabel('\delta Flux : ClrSky-X [W/m2]'); title([num2str(ii) ' : ' num2str(RRTM_bands(ii-1)) ' to ' num2str(RRTM_bands(ii))]); 
      xlim([2020 2025]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)

    figure(5); plot(yymm,smooth(junk0trp,smn),'kx-',yymm,smooth(junk1trp,smn),'b',yymm,smooth(junk2trp,smn),'g',yymm,smooth(junk3trp,smn),'r',yymm,smooth(junk4trp,smn),'c',yymm,smooth(junk5trp,smn),'ms-','linewidth',2); 
      plotaxis2; ylabel('Flux [W/m2]'); title(['TRP ' num2str(ii) ' : ' num2str(RRTM_bands(ii-1)) ' to ' num2str(RRTM_bands(ii))]);
      xlim([2020 2025]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('ALL','Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)
    figure(6); plot(yymm,smooth(junk0trp-junk1trp,smn),'b',yymm,smooth(junk0trp-junk2trp,smn),'g',yymm,smooth(junk0trp-junk3trp,smn),'r',yymm,smooth(junk0trp-junk4trp,smn),'c',yymm,smooth(junk0trp-junk5trp,smn),'ms-','linewidth',2); 
      plotaxis2; ylabel('\delta Flux : ClrSky-X [W/m2]'); title(['TRP ' num2str(ii) ' : ' num2str(RRTM_bands(ii-1)) ' to ' num2str(RRTM_bands(ii))]); 
      xlim([2020 2025]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)

  else
    figure(1); plot(yymm,smooth(junk0,smn),'kx-',yymm,smooth(junk1,smn),'b',yymm,smooth(junk2,smn),'g',yymm,smooth(junk3,smn),'r',yymm,smooth(junk4,smn),'c',yymm,smooth(junk5,smn),'ms-','linewidth',2); 
      ylim([-1 +1]*2);
      plotaxis2; ylabel('Flux [W/m2]'); title('all flux')
      ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('ALL','Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)
    figure(2); plot(yymm,smooth(junk0-junk1,smn),'b',yymm,smooth(junk0-junk2,smn),'g',yymm,smooth(junk0-junk3,smn),'r',yymm,smooth(junk0-junk4,smn),'c',yymm,smooth(junk0-junk5,smn),'ms-','linewidth',2); 
      ylim([-1 +1]*2);
      plotaxis2; ylabel('\delta Flux : ClrSky-X [W/m2]'); title('all flux'); 
      ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)

    figure(3); plot(yymm,smooth(junk0,smn),'kx-',yymm,smooth(junk1,smn),'b',yymm,smooth(junk2,smn),'g',yymm,smooth(junk3,smn),'r',yymm,smooth(junk4,smn),'c',yymm,smooth(junk5,smn),'ms-','linewidth',2); 
      ylim([-1 +1]*2);
      plotaxis2; ylabel('Flux [W/m2]'); title('all flux')
      xlim([2020 2025]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('ALL','Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)
    figure(4); plot(yymm,smooth(junk0-junk1,smn),'b',yymm,smooth(junk0-junk2,smn),'g',yymm,smooth(junk0-junk3,smn),'r',yymm,smooth(junk0-junk4,smn),'c',yymm,smooth(junk0-junk5,smn),'ms-','linewidth',2); 
      ylim([-1 +1]*2);
      plotaxis2; ylabel('\delta Flux : ClrSky-X [W/m2]'); title('all flux'); 
      xlim([2020 2025]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)

    figure(5); plot(yymm,smooth(junk0trp,smn),'kx-',yymm,smooth(junk1trp,smn),'b',yymm,smooth(junk2trp,smn),'g',yymm,smooth(junk3trp,smn),'r',yymm,smooth(junk4trp,smn),'c',yymm,smooth(junk5trp,smn),'ms-','linewidth',2); 
      ylim([-1 +1]*2);
      plotaxis2; ylabel('Flux [W/m2]'); title('TRP all flux')
      xlim([2020 2025]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('ALL','Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)
    figure(6); plot(yymm,smooth(junk0trp-junk1trp,smn),'b',yymm,smooth(junk0trp-junk2trp,smn),'g',yymm,smooth(junk0trp-junk3trp,smn),'r',yymm,smooth(junk0trp-junk4trp,smn),'c',yymm,smooth(junk0trp-junk5trp,smn),'ms-','linewidth',2); 
      ylim([-1 +1]*2);
      plotaxis2; ylabel('\delta Flux : ClrSky-X [W/m2]'); title('TRP all flux'); 
      xlim([2020 2025]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('Const WV','Const CO2','Const T/ST','Const RH','Clouds','fontsize',8,'location','best');
      set(gca,'fontsize',10)
  end

  figure(1); title('All Latbins')
  figure(2); title('All Latbins')
  figure(3); title('All Latbins')
  figure(4); title('All Latbins')
  figure(5); title('TRP Latbins')
  figure(6); title('TRP Latbins')

  pause
end
