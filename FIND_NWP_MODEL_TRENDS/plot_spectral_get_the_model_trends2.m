figure(5); plot(fchanx,nanmean(cosYY.*obsrates.rates - cosYY.*obsrates.rates,2),'k',fchanx,nanmean(cosYY.*obsrates.rates - cosYY.*umbcL3.umbcL3_spectral_rates,2),'g',...
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
figure(6); plot(fchanx,nanmean(cosYYSPolar.*obsrates.rates - cosYYSPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYSPolar.*obsrates.rates - cosYYSPolar.*umbcL3.umbcL3_spectral_rates,2),'g',...
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
figure(7); plot(fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*umbcL3.umbcL3_spectral_rates,2),'g',...
                fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*climcapsL3.climcapsL3_spectral_rates,2),'c',...
                fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYYNPolar.*obsrates.rates - cosYYNPolar.*merra2.merra2_spectral_rates,2),'m','linewidth',2);
 xlim([640 1640]); ylim([-0.015 +0.015])
  plotaxis2;
 title('NPolar lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Bias dBT/dt K/yr')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8); plot(fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*obsrates.rates'),'k',fchanx,nanstd(cosYY'.*obsrates.rates' - cosYY'.*umbcL3.umbcL3_spectral_rates'),'g',...
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
figure(9); plot(fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*obsrates.rates'),'k',fchanx,nanstd(cosYYSPolar'.*obsrates.rates' - cosYYSPolar'.*umbcL3.umbcL3_spectral_rates'),'g',...
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
figure(10); plot(fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*obsrates.rates'),'k',fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*umbcL3.umbcL3_spectral_rates'),'g',...
                fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*airsL3.airsL3_spectral_rates'),'b',fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*climcapsL3.climcapsL3_spectral_rates'),'c',...
                fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*era5.era5_spectral_rates'),'r',fchanx,nanstd(cosYYNPolar'.*obsrates.rates' - cosYYNPolar'.*merra2.merra2_spectral_rates'),'m','linewidth',2);
 xlim([640 1640]); ylim([0 +0.03])
  plotaxis2;
 title('NPolar lats : (Obs-Cal)'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('Std Dev dBT/dt K/yr')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

