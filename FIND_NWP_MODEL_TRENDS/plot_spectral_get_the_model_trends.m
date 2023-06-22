figure(5); plot(fchanx,nanmean(cosYY.*obsrates.rates,2),'k',fchanx,nanmean(cosYY.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYY.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYY.*cmip6.cmip6_spectral_rates,2),'g','linewidth',2);
 xlim([640 1640]); ylim([-0.08 +0.04])
  plotaxis2;
 title('All lats'); hl = legend('OBS',iJorCstrSPECTRA,iEorMstrSPECTRA,iAorCstrSPECTRA,'location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYSPolar = cos(YY*pi/180)';
cosYYSPolar = ones(size(YY))';
cosYYSPolar(YY > -60) = 0;
cosYYSPolar = ones(2645,1) * cosYYSPolar;
figure(6); plot(fchanx,nanmean(cosYYSPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYSPolar.*airsL3.airsL3_spectral_rates,2),'b',...
                fchanx,nanmean(cosYYSPolar.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYYSPolar.*cmip6.cmip6_spectral_rates,2),'g','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01])
  plotaxis2;
 title('SPolar lats'); hl = legend('OBS',iJorCstrSPECTRA,iEorMstrSPECTRA,iAorCstrSPECTRA,'location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYNPolar = cos(YY*pi/180)';
cosYYNPolar = ones(size(YY))';
cosYYNPolar(YY < +60) = 0;
cosYYNPolar = ones(2645,1) * cosYYNPolar;
figure(7); plot(fchanx,nanmean(cosYYNPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYNPolar.*airsL3.airsL3_spectral_rates,2),'b',...
                fchanx,nanmean(cosYYNPolar.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYYNPolar.*cmip6.cmip6_spectral_rates,2),'g','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01])
  plotaxis2;
 title('NPolar lats'); hl = legend('OBS',iJorCstrSPECTRA,iEorMstrSPECTRA,iAorCstrSPECTRA,'location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); plot(fchanx,nanmean(cosYY.*obsrates.rates,2),'k',fchanx,nanmean(cosYY.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYY.*era5.era5_spectral_rates,2),'r','linewidth',2);
 xlim([640 1640]); ylim([-0.08 +0.04])
  plotaxis2;
 title('All lats'); hl = legend('OBS',iJorCstrSPECTRA,iEorMstrSPECTRA,iAorCstrSPECTRA,'location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYSPolar = cos(YY*pi/180)';
cosYYSPolar = ones(size(YY))';
cosYYSPolar(YY > -60) = 0;
cosYYSPolar = ones(2645,1) * cosYYSPolar;
figure(6); plot(fchanx,nanmean(cosYYSPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYSPolar.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYYSPolar.*era5.era5_spectral_rates,2),'r','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01])
  plotaxis2;
 title('SPolar lats'); hl = legend('OBS',iJorCstrSPECTRA,iEorMstrSPECTRA,'location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYNPolar = cos(YY*pi/180)';
cosYYNPolar = ones(size(YY))';
cosYYNPolar(YY < +60) = 0;
cosYYNPolar = ones(2645,1) * cosYYNPolar;
figure(7); plot(fchanx,nanmean(cosYYNPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYNPolar.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYYNPolar.*era5.era5_spectral_rates,2),'r','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01])
  plotaxis2;
 title('NPolar lats'); hl = legend('OBS',iJorCstrSPECTRA,iEorMstrSPECTRA,'location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5); plot(fchanx,nanmean(cosYY.*obsrates.rates,2),'k',fchanx,nanmean(cosYY.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYY.*climcapsL3.climcapsL3_spectral_rates,2),'c',......
                fchanx,nanmean(cosYY.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYY.*merra2.merra2_spectral_rates,2),'g','linewidth',2);
 xlim([640 1640]); ylim([-0.08 +0.04])
  plotaxis2;
 title('All lats'); hl = legend('OBS',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYSPolar = cos(YY*pi/180)';
cosYYSPolar = ones(size(YY))';
cosYYSPolar(YY > -60) = 0;
cosYYSPolar = ones(2645,1) * cosYYSPolar;
figure(6); plot(fchanx,nanmean(cosYYSPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYSPolar.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYYSPolar.*climcapsL3.climcapsL3_spectral_rates,2),'c',...
                fchanx,nanmean(cosYYSPolar.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYYSPolar.*merra2.merra2_spectral_rates,2),'g','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01])
  plotaxis2;
 title('SPolar lats'); hl = legend('OBS',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYNPolar = cos(YY*pi/180)';
cosYYNPolar = ones(size(YY))';
cosYYNPolar(YY < +60) = 0;
cosYYNPolar = ones(2645,1) * cosYYNPolar;
figure(7); plot(fchanx,nanmean(cosYYNPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYNPolar.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYYNPolar.*climcapsL3.climcapsL3_spectral_rates,2),'c',...
                fchanx,nanmean(cosYYNPolar.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYYNPolar.*merra2.merra2_spectral_rates,2),'g','linewidth',2);
 xlim([640 1640]); ylim([-0.015 +0.015])
  plotaxis2;
 title('NPolar lats'); hl = legend('OBS',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5); plot(fchanx,nanmean(cosYY.*obsrates.rates,2),'k',fchanx,nanmean(cosYY.*umbcL3.umbcL3_spectral_rates,2),'g',...
                fchanx,nanmean(cosYY.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYY.*climcapsL3.climcapsL3_spectral_rates,2),'c',......
                fchanx,nanmean(cosYY.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYY.*merra2.merra2_spectral_rates,2),'m','linewidth',2);
 xlim([640 1640]); ylim([-0.08 +0.04])
  plotaxis2;
 title('All lats'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYSPolar = cos(YY*pi/180)';
cosYYSPolar = ones(size(YY))';
cosYYSPolar(YY > -60) = 0;
cosYYSPolar = ones(2645,1) * cosYYSPolar;
figure(6); plot(fchanx,nanmean(cosYYSPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYSPolar.*umbcL3.umbcL3_spectral_rates,2),'g',...
                fchanx,nanmean(cosYYSPolar.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYYSPolar.*climcapsL3.climcapsL3_spectral_rates,2),'c',...
                fchanx,nanmean(cosYYSPolar.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYYSPolar.*merra2.merra2_spectral_rates,2),'m','linewidth',2);
 xlim([640 1640]); ylim([-0.01 +0.01])
  plotaxis2;
 title('SPolar lats'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

cosYYNPolar = cos(YY*pi/180)';
cosYYNPolar = ones(size(YY))';
cosYYNPolar(YY < +60) = 0;
cosYYNPolar = ones(2645,1) * cosYYNPolar;
figure(7); plot(fchanx,nanmean(cosYYNPolar.*obsrates.rates,2),'k',fchanx,nanmean(cosYYNPolar.*umbcL3.umbcL3_spectral_rates,2),'g',...
                fchanx,nanmean(cosYYNPolar.*airsL3.airsL3_spectral_rates,2),'b',fchanx,nanmean(cosYYNPolar.*climcapsL3.climcapsL3_spectral_rates,2),'c',...
                fchanx,nanmean(cosYYNPolar.*era5.era5_spectral_rates,2),'r',fchanx,nanmean(cosYYNPolar.*merra2.merra2_spectral_rates,2),'m','linewidth',2);
 xlim([640 1640]); ylim([-0.015 +0.015])
  plotaxis2;
 title('NPolar lats'); hl = legend('OBS','UMBC',iJorCstrSPECTRA,'CLIMCAPS L3',iEorMstrSPECTRA,'MERRA2','location','best','fontsize',8);
 xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

