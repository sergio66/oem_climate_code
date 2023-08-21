iOffSet = 40;
figure(iOffSet+1); clf;
plot(fchanx,nanmean(cosYY.*obsrates.rates - cosYY.*obsrates.rates,2),'k',fchanx,nanmean(cosYY.*obsrates.rates - cosYY.*umbcL3.umbcL3_spectral_rates,2),'g',...
     fchanx,nanmean(cosYY.*obsrates.rates - cosYY.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYY.*obsrates.rates - cosYY.*climcapsL3.climcapsL3_spectral_rates,2),'c',......
     fchanx,nanmean(cosYY.*obsrates.rates - cosYY.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYY.*obsrates.rates - cosYY.*merra2.merra2_spectral_rates,2),'m','linewidth',2);
 xlim([640 1640]); ylim([-0.04 +0.04])
 plotaxis2;
 title('All lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Bias dBT/dt K/yr')

cosYYSPolar = cos(YY*pi/180)';
cosYYSPolar = ones(size(YY))';
cosYYSPolar(YY > -60) = 0;
cosYYSPolar = ones(2645,1) * cosYYSPolar;
figure(iOffSet+2); clf;
plot(fchanx,nanmean(cosYYSPolar.*obsrates.rates - cosYYSPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYSPolar.*obsrates.rates - cosYYSPolar.*umbcL3.umbcL3_spectral_rates,2),'g',...
     fchanx,nanmean(cosYYSPolar.*obsrates.rates - cosYYSPolar.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYYSPolar.*obsrates.rates - cosYYSPolar.*climcapsL3.climcapsL3_spectral_rates,2),'c',...
     fchanx,nanmean(cosYYSPolar.*obsrates.rates - cosYYSPolar.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYYSPolar.*obsrates.rates - cosYYSPolar.*merra2.merra2_spectral_rates,2),'m','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01])
 plotaxis2;
 title('SPolar lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Bias dBT/dt K/yr')

cosYYNPolar = cos(YY*pi/180)';
cosYYNPolar = ones(size(YY))';
cosYYNPolar(YY < +60) = 0;
cosYYNPolar = ones(2645,1) * cosYYNPolar;
figure(iOffSet+3); clf;
plot(fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*umbcL3.umbcL3_spectral_rates,2),'g',...
     fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*climcapsL3.climcapsL3_spectral_rates,2),'c',...
     fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*merra2.merra2_spectral_rates,2),'m','linewidth',2);
 xlim([640 1640]); ylim([-0.015 +0.015])
 plotaxis2;
 title('NPolar lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Bias dBT/dt K/yr')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(iOffSet+4); clf;
plot(fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*obsrates.rates'),'k',fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*umbcL3.umbcL3_spectral_rates'),'g',...
     fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*airsL3.airsL3_spectral_rates'),'b',fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*climcapsL3.climcapsL3_spectral_rates'),'c',......
     fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*era5.era5_spectral_rates'),'r',fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*merra2.merra2_spectral_rates'),'m','linewidth',2);
 xlim([640 1640]); ylim([0 +0.04])
  plotaxis2;
 title('All lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Std Dev dBT/dt K/yr')

cosYYSPolar = cos(YY*pi/180)';
cosYYSPolar = ones(size(YY))';
cosYYSPolar(YY > -60) = 0;
cosYYSPolar = ones(2645,1) * cosYYSPolar;
figure(iOffSet+5); clf;
plot(fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*obsrates.rates'),'k',fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*umbcL3.umbcL3_spectral_rates'),'g',...
     fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*airsL3.airsL3_spectral_rates'),'b',fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*climcapsL3.climcapsL3_spectral_rates'),'c',...
     fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*era5.era5_spectral_rates'),'r',fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*merra2.merra2_spectral_rates'),'m','linewidth',2);
 xlim([640 1640]); ylim([0 +0.03])
 plotaxis2;
 title('SPolar lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Std Dev dBT/dt K/yr')

cosYYNPolar = cos(YY*pi/180)';
cosYYNPolar = ones(size(YY))';
cosYYNPolar(YY < +60) = 0;
cosYYNPolar = ones(2645,1) * cosYYNPolar;
figure(iOffSet+6); clf;
plot(fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*obsrates.rates'),'k',fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*umbcL3.umbcL3_spectral_rates'),'g',...
     fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*airsL3.airsL3_spectral_rates'),'b',fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*climcapsL3.climcapsL3_spectral_rates'),'c',...
     fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*era5.era5_spectral_rates'),'r',fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*merra2.merra2_spectral_rates'),'m','linewidth',2);
 xlim([640 1640]); ylim([0 +0.03])
 plotaxis2;
 title('NPolar lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Std Dev dBT/dt K/yr')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(iOffSet+7); clf;
plot(fchanx,nanmean(cosYY.*obsrates.rates,2),'k',fchanx,nanmean(cosYY.*umbcL3.umbcL3_spectral_rates,2),'g',...
     fchanx,nanmean(cosYY.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYY.*climcapsL3.climcapsL3_spectral_rates,2),'c',......
     fchanx,nanmean(cosYY.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYY.*merra2.merra2_spectral_rates,2),'m','linewidth',2);
 xlim([640 1640]); ylim([-1 +1]*0.02)
 plotaxis2;
 title('All lats : (Obs or Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYSPolar = cos(YY*pi/180)';
cosYYSPolar = ones(size(YY))';
cosYYSPolar(YY > -60) = 0;
cosYYSPolar = ones(2645,1) * cosYYSPolar;
figure(iOffSet+8); clf;
plot(fchanx,nanmean(cosYYSPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYSPolar.*umbcL3.umbcL3_spectral_rates,2),'g',...
     fchanx,nanmean(cosYYSPolar.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYYSPolar.*climcapsL3.climcapsL3_spectral_rates,2),'c',...
     fchanx,nanmean(cosYYSPolar.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYYSPolar.*merra2.merra2_spectral_rates,2),'m','linewidth',2);
 xlim([640 1640]); ylim([-1 +1]*0.02)
 plotaxis2;
 title('SPolar lats : (Obs or Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYNPolar = cos(YY*pi/180)';
cosYYNPolar = ones(size(YY))';
cosYYNPolar(YY < +60) = 0;
cosYYNPolar = ones(2645,1) * cosYYNPolar;
figure(iOffSet+9); clf;
plot(fchanx,nanmean(cosYYNPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYNPolar.*umbcL3.umbcL3_spectral_rates,2),'g',...
     fchanx,nanmean(cosYYNPolar.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYYNPolar.*climcapsL3.climcapsL3_spectral_rates,2),'c',...
     fchanx,nanmean(cosYYNPolar.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYYNPolar.*merra2.merra2_spectral_rates,2),'m','linewidth',2);
 xlim([640 1640]); ylim([-1 +1]*0.02)
 plotaxis2;
 title('NPolar lats : (Obs or Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('just showing all 4608 lats Figs 40-45')
figure(iOffSet+10); clf;
  pcolor(fchanx,YY,obsrates.rates');                       shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('Obs Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+11); clf;
  pcolor(fchanx,YY,umbcL3.umbcL3_spectral_rates');         shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('UMBC Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+12); clf;
  pcolor(fchanx,YY,era5.era5_spectral_rates');             shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('ERA5 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+13); clf;
  pcolor(fchanx,YY,merra2.merra2_spectral_rates');         shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('MERRA2 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+14); clf;
  pcolor(fchanx,YY,airsL3.airsL3_spectral_rates');         shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('AIRS L3 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+15); clf;
  pcolor(fchanx,YY,climcapsL3.climcapsL3_spectral_rates'); shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('CLIMCAPS L3 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('doing zonal average Figs 50-55')
figure(iOffSet+20); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(obsrates.rates',72,64,2645),1)));                       
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('Obs Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+21); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(umbcL3.umbcL3_spectral_rates',72,64,2645),1)));         
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('UMBC Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+22); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(era5.era5_spectral_rates',72,64,2645),1)));             
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('ERA5 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+23); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(merra2.merra2_spectral_rates',72,64,2645),1)));         
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('MERRA2 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+24); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(airsL3.airsL3_spectral_rates',72,64,2645),1)));         
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('AIRS L3 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
figure(iOffSet+25); clf;
  pcolor(fchanx,nanmean(reshape(YY,72,64),1),squeeze(nanmean(reshape(climcapsL3.climcapsL3_spectral_rates',72,64,2645),1))); 
  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('CLIMCAPS L3 Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure(iOffSet+30); clf;
%  wah = obsrates.rates; wah = reshape(wah,2645,72,64); wah = squeeze(nanmean(wah,2));
%  pcolor(fchanx,nanmean(reshape(YY,72,64),1),wah);
%  shading flat; colorbar; colormap(usa2); caxis([-1 +1]*0.1); title('Obs Rates'); xlabel('Wavenumber'); ylabel('Latitude'); xlim([640 1640])
%keyboard_nowindow
