cosL1Clat = cos(rlat*pi/180) * ones(1,500);
cosL3lat  = cos(rlat*pi/180) * ones(1,262);
yymmL3 = yymm(1:262);

smn = 3;

rh_plot = ones(1,14);;

for ii = 1 : 14
  %% easy
  % junk0 = squeeze(nanmean(a0.anomflux(:,ii,:),1)) * 3.5788; 
  % junk1 = squeeze(nanmean(a1.anomflux(:,ii,:),1)) * 3.5788;
  % junk2 = squeeze(nanmean(a2.anomflux(:,ii,:),1)) * 3.5788;
  % junk3 = squeeze(nanmean(a3.anomflux(:,ii,:),1)) * 3.5788;
  % junk4 = squeeze(nanmean(a4.anomflux(:,ii,:),1)) * 3.5788;
  % junk5 = squeeze(nanmean(a5.anomflux(:,ii,:),1)) * 3.5788;
  % junkObsALL = squeeze(nanmean(fluxanomD(:,:,ii),1))' * 3.5788;
  % junkObsCLR = squeeze(nanmean(fluxanomA(:,:,ii),1))' * 3.5788;

  junk0 = squeeze(a0.anomflux(:,ii,:)) * 3.5788; 
  junk1 = squeeze(a1.anomflux(:,ii,:)) * 3.5788;
  junk2 = squeeze(a2.anomflux(:,ii,:)) * 3.5788;
  junk3 = squeeze(a3.anomflux(:,ii,:)) * 3.5788;
  junk4 = squeeze(a4.anomflux(:,ii,:)) * 3.5788;
  junk5 = squeeze(a5.anomflux(:,ii,:)) * 3.5788;
  junk6 = squeeze(a6.anomflux(:,ii,:)) * 3.5788;
  junkObsALL = squeeze(fluxanomD(:,:,ii)) * 3.5788;
  junkObsCLR = squeeze(fluxanomA(:,:,ii)) * 3.5788;
  junkAIRSL3 = squeeze(airsL3_sarta_flux.anomflux(:,ii,:)) * 3.5788;

  junk0trp = nanmean(junk0(tropics,:),1);
  junk1trp = nanmean(junk1(tropics,:),1);
  junk2trp = nanmean(junk2(tropics,:),1);
  junk3trp = nanmean(junk3(tropics,:),1);
  junk4trp = nanmean(junk4(tropics,:),1);
  junk5trp = nanmean(junk5(tropics,:),1);
  junk6trp = nanmean(junk6(tropics,:),1);
  junkObsALLtrp = nanmean(junkObsALL(tropics,:),1);
  junkObsCLRtrp = nanmean(junkObsCLR(tropics,:),1);
  junkAIRSL3trp = nanmean(junkAIRSL3(tropics,:),1);

  junk0 = sum(junk0.*coslat)./sum(coslat);
  junk1 = sum(junk1.*coslat)./sum(coslat);
  junk2 = sum(junk2.*coslat)./sum(coslat);
  junk3 = sum(junk3.*coslat)./sum(coslat);
  junk4 = sum(junk4.*coslat)./sum(coslat);
  junk5 = sum(junk5.*coslat)./sum(coslat);
  junk6 = sum(junk6.*coslat)./sum(coslat);
  junkObsALL = sum(junkObsALL.*cosL1Clat)./sum(cosL1Clat);
  junkObsCLR = sum(junkObsCLR.*cosL1Clat)./sum(cosL1Clat);
  junkAIRSL3 = sum(junkAIRSL3.*cosL3lat)./sum(cosL3lat);

%junk4 = junk6;
%junk4trp = junk6trp;

%junk4    = junk4    + (junk6-junk4)/20;
%junk4trp = junk4trp + (junk6trp-junk4trp)/20;

  if ii > 1
    figure(1); plot(yymm,smooth(junk0,smn),'m',yymmL1C,smooth(junkObsALL,2*smn+1),'b.-',yymmL1C,smooth(junkObsCLR,2*smn+1),'r.-',yymmL3,rh_plot(ii)*smooth(junkAIRSL3,smn),'g',yymm,smooth(junk5,smn),'c','linewidth',2); 
      plotaxis2; ylabel('Flux [W/m2]'); title([num2str(ii) ' : ' num2str(RRTM_bands(ii-1)) ' to ' num2str(RRTM_bands(ii))]);
      xlim([2002.5 2024.5]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('CLR sarta','Obs ALL','Obs CLR','AIRSv7 L3','ALL sarta','fontsize',8,'location','best');
      set(gca,'fontsize',10)

    figure(2); plot(yymm,smooth(junk0,smn),'m',yymmL1C,smooth(junkObsALL,2*smn+1),'bd-',yymmL1C,smooth(junkObsCLR,2*smn+1),'ro-',yymmL3,rh_plot(ii)*smooth(junkAIRSL3,smn),'g',yymm,smooth(junk5,smn),'c','linewidth',2); 
      plotaxis2; ylabel('Flux [W/m2]'); title([num2str(ii) ' : ' num2str(RRTM_bands(ii-1)) ' to ' num2str(RRTM_bands(ii))]);
      xlim([2020.5 2024.5]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('CLR sarta','Obs ALL','Obs CLR','AIRSv7 L3','ALL sarta','fontsize',8,'location','best');
      set(gca,'fontsize',10)

    figure(3); plot(yymm,smooth(junk0trp,smn),'m',yymmL1C,smooth(junkObsALLtrp,2*smn+1),'bd-',yymmL1C,smooth(junkObsCLRtrp,2*smn+1),'ro-',yymmL3,rh_plot(ii)*smooth(junkAIRSL3trp,smn),'g',yymm,smooth(junk5trp,smn),'c','linewidth',2); 
      plotaxis2; ylabel('Flux [W/m2]'); title(['TRP ' num2str(ii) ' : ' num2str(RRTM_bands(ii-1)) ' to ' num2str(RRTM_bands(ii))]);
      xlim([2020.5 2024.5]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('CLR sarta','Obs ALL','Obs CLR','AIRSv7 L3','ALL sarta','fontsize',8,'location','best');
      set(gca,'fontsize',10)

  else
    figure(1); plot(yymm,smooth(junk0,smn),'m',yymmL1C,smooth(junkObsALL,2*smn+1),'b.-',yymmL1C,smooth(junkObsCLR,2*smn+1),'r.-',yymmL3,rh_plot(ii)*smooth(junkAIRSL3,smn),'g',yymm,smooth(junk5,smn),'c','linewidth',2); 
      ylim([-1 +1]*2);
      plotaxis2; ylabel('Flux [W/m2]'); title('all flux')
      xlim([2002.5 2024.5]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('CLR sarta','Obs ALL','Obs CLR','AIRSv7 L3','ALL sarta','fontsize',8,'location','best');
      set(gca,'fontsize',10)

    figure(2); plot(yymm,smooth(junk0,smn),'m',yymmL1C,smooth(junkObsALL,2*smn+1),'bd-',yymmL1C,smooth(junkObsCLR,2*smn+1),'ro-',yymmL3,rh_plot(ii)*smooth(junkAIRSL3,smn),'g',yymm,smooth(junk5,smn),'c','linewidth',2); 
      ylim([-1 +1]*2);
      plotaxis2; ylabel('Flux [W/m2]'); title('all flux')
      xlim([2020.5 2024.5]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('CLR sarta','Obs ALL','Obs CLR','AIRSv7 L3','ALL sarta','fontsize',8,'location','best');
      set(gca,'fontsize',10)

    figure(3); plot(yymm,smooth(junk0trp,smn),'m',yymmL1C,smooth(junkObsALLtrp,2*smn+1),'bd-',yymmL1C,smooth(junkObsCLRtrp,2*smn+1),'ro-',yymmL3,rh_plot(ii)*smooth(junkAIRSL3trp,smn),'g',yymm,smooth(junk5trp,smn),'c','linewidth',2); 
      ylim([-1 +1]*2);
      plotaxis2; ylabel('Flux [W/m2]'); title('TRP all flux')
      xlim([2020.5 2024.5]); ax = axis; line([ht ht],[ax(3) ax(4)],'color','k','linewidth',2)
      legend('CLR sarta','Obs ALL','Obs CLR','AIRSv7 L3','ALL sarta','fontsize',8,'location','best');
      set(gca,'fontsize',10)
  end

  iPrint = input('print (-1 [default]/+1) : ');
  if length(iPrint) == 0
    iPrint = -1;
  end
  if iPrint > 0
    %figure(1); sergioprint(['/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/fluxes_era5_airsL3_band' num2str(ii,'%02d') '_obs_sarta']);
    figure(2); sergioprint(['/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/fluxes_era5_airsL3_band' num2str(ii,'%02d') '_zoom_obs_sarta']);
    %figure(3); sergioprint(['/home/sergio/PAPERS/AIRS/AIRS_STM_Sept2024/Trends/Figs/fluxes_era5_airsL3_band' num2str(ii,'%02d') '_zoom_topics_obs_sarta']);
  end   

end
