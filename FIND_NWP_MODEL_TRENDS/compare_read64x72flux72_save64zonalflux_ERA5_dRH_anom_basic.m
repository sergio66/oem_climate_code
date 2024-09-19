ii = 1;  junk0 = squeeze(a0.anomflux(:,ii,:)) * 3.5788; 
ii = 1;  junk5 = squeeze(a5.anomflux(:,ii,:)) * 3.5788;
  junk0trp = nanmean(junk0(tropics,:),1);
  junk5trp = nanmean(junk5(tropics,:),1);
  junk0 = sum(junk0.*coslat)./sum(coslat);
  junk5 = sum(junk5.*coslat)./sum(coslat);

figure(1)
  plot(yymm,smooth(era5.olr,smn),'r',yymm,smooth(junk5,smn),'m',yymm,smooth(era5.olr_clr,smn),'b',yymm,smooth(junk0,smn),'c','linewidth',2)
  plotaxis2; 
  ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
  legend('OLR CLD ERA5','OLR CLD SARTA','OLR CLR ERA5','OLR CLR SARTA','fontsize',8,'location','best');
  ylim([-1 +1]*1.25); ylabel('Flux W/m2')
  xlim([2002.75 2024.50]); set(gca,'fontsize',10)
figure(2)
  plot(yymm,smooth(era5.olr - junk5,smn),'r',yymm,smooth(era5.olr_clr - junk0,smn),'b','linewidth',2)
  plotaxis2;
  ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
  legend('OLR CLD: ERA5 - SARTA','OLR CLR : ERA5 - SARTA','fontsize',8,'location','best');
  ylim([-1 +1]*1.25); ylabel('\delta Flux : ClrSky-X W/m2')
  xlim([2002.75 2024.50]); set(gca,'fontsize',10)

figure(3)
  plot(yymm,smooth(era5.olr,smn),'r',yymm,smooth(junk5,smn),'m',yymm,smooth(era5.olr_clr,smn),'b',yymm,smooth(junk0,smn),'c','linewidth',2)
  xlim([2020 2025]); plotaxis2; ylabel('Flux [W/m2]');
  ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
  legend('OLR CLD ERA5','OLR CLD SARTA','OLR CLR ERA5','OLR CLR SARTA','fontsize',8,'location','best');
  ylim([-1 +1]*1.25); ylabel('Flux W/m2')
  set(gca,'fontsize',10)
figure(4)
  plot(yymm,smooth(era5.olr - junk5,smn),'r',yymm,smooth(era5.olr_clr - junk0,smn),'b','linewidth',2)
  xlim([2020 2025]); plotaxis2;
  ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
  legend('OLR CLD: ERA5 - SARTA','OLR CLR : ERA5 - SARTA','fontsize',8,'location','best');
  ylim([-1 +1]*1.25); ylabel('\delta Flux anomaly : ClrSky-X W/m2')
  set(gca,'fontsize',10)

figure(5)
  plot(yymm,smooth(era5.olr,smn),'r',yymm,smooth(junk5trp,smn),'m',yymm,smooth(era5.olr_clr,smn),'b',yymm,smooth(junk0trp,smn),'c','linewidth',2)
  xlim([2020 2025]); plotaxis2; ylabel('Flux [W/m2]');
  ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
  legend('OLR CLD ERA5','OLR CLD SARTA','OLR CLR ERA5','OLR CLR SARTA','fontsize',8,'location','best');
  ylim([-1 +1]*1.25); ylabel('TRP Flux W/m2')
  set(gca,'fontsize',10)
figure(6)
  plot(yymm,smooth(era5.olr - junk5trp,smn),'r',yymm,smooth(era5.olr_clr - junk0trp,smn),'b','linewidth',2)
  xlim([2020 2025]); plotaxis2;
  ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
  legend('OLR CLD: ERA5 - SARTA','OLR CLR : ERA5 - SARTA','fontsize',8,'location','best');
  ylim([-1 +1]*1.25); ylabel('TRP \delta Flux anomaly : ClrSky-X W/m2')
  set(gca,'fontsize',10)

figure(1); title('All Latbins')
figure(2); title('All Latbins')
figure(3); title('All Latbins')
figure(4); title('All Latbins')
figure(5); title('TRP Latbins')
figure(6); title('TRP Latbins')

