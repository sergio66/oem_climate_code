iFig = 50;
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(nanmean(fracWV(1:100,mask),2),plays,'linewidth',2);                   ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
  hold on;                           subplot(121); semilogy(nanmean(fracWVunc(1:100,mask),2),plays,'linewidth',2);                ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac + unc pert');
                                     subplot(122); semilogy(1+nanstd(fracWV(1:100,mask),[],2),plays,'linewidth',2);               ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');
  hold on;                           subplot(122); semilogy(1+nanstd(fracWVunc(1:100,mask),[],2),plays,'linewidth',2);            ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(nanmean(deltaRH(1:100,mask),2),plays,'linewidth',2);                  ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
  hold on;                           subplot(121); semilogy(nanmean(deltaRHunc(1:100,mask),2),plays,'linewidth',2);               ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert + unc (%)');
                                     subplot(122); semilogy(nanstd(deltaRH(1:100,mask),[],2),plays,'linewidth',2);                ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');
  hold on;                           subplot(122); semilogy(nanstd(deltaRHunc(1:100,mask),[],2),plays,'linewidth',2);             ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(nanmean(deltaT(1:100,mask),2),plays,'linewidth',2);                   ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
  hold on;                           subplot(121); semilogy(nanmean(deltaTunc(1:100,mask),2),plays,'linewidth',2);                ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert + unc');
                                     subplot(122); semilogy(nanstd(deltaT(1:100,mask),[],2),plays,'linewidth',2);                 ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');
  hold on;                           subplot(122); semilogy(nanstd(deltaTunc(1:100,mask),[],2),plays,'linewidth',2);              ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(nanmean(deltaO3(1:97,mask),2),plays(1:97),'linewidth',2);             ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
  hold on;                           subplot(121); semilogy(nanmean(deltaO3unc(1:97,mask),2),plays(1:97),'linewidth',2);          ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert + unc');
                                     subplot(122); semilogy(nanstd(deltaO3(1:97,mask),[],2),plays(1:97),'linewidth',2);           ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');
  hold on;                           subplot(122); semilogy(nanstd(deltaO3unc(1:97,mask),[],2),plays(1:97),'linewidth',2);        ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');

if exist('airsL3')
  iFig = 50;
  Tlevs = airsL3.Tlevs;
  Qlevs = airsL3.Qlevs;

  boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = maskLFmatr'; end
  junk = airsL3.thestats64x72.waterrate.*boo; junk = reshape(junk,72*64,12);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,1)',Qlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
                                   subplot(122); hold on; semilogy(1+nanstd(junk,[],1)',Qlevs,'linewidth',2);     hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');

  boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = maskLFmatr'; end
  junk = airsL3.thestats64x72.RHrate.*boo; junk = reshape(junk,72*64,12);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,1)',Qlevs,'linewidth',2);        hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
                                   subplot(122); hold on; semilogy(nanstd(junk,[],1)',Qlevs,'linewidth',2);      hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');

  boo = zeros(72,64,12); for ijunk = 1 : 24; boo(:,:,ijunk) = maskLFmatr'; end
  junk = airsL3.thestats64x72.ptemprate.*boo; junk = reshape(junk,72*64,24);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,1)',Tlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
                                   subplot(122); hold on; semilogy(nanstd(junk,[],1)',Tlevs,'linewidth',2);       hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');

  boo = zeros(72,64,12); for ijunk = 1 : 24; boo(:,:,ijunk) = maskLFmatr'; end
  junk = airsL3.thestats64x72.ozonerate.*boo; junk = reshape(junk,72*64,24);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,1)',Tlevs,'linewidth',2);   hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
                                   subplot(122); hold on; semilogy(nanstd(junk,[],1)',Tlevs,'linewidth',2); hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');
end

if exist('cmip6')
  iFig = 50;
  Tlevs = plays;
  Qlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end

  junk = cmip6.trend_gas_1; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Qlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
                                   subplot(122); hold on; semilogy(1+nanstd(junk,[],2)',Qlevs,'linewidth',2);     hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');

  junk = cmip6.trend_RH; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Qlevs,'linewidth',2);        hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
                                   subplot(122); hold on; semilogy(nanstd(junk,[],2)',Qlevs,'linewidth',2);      hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');

  junk = cmip6.trend_ptemp; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
                                   subplot(122); hold on; semilogy(nanstd(junk,[],2)',Tlevs,'linewidth',2);       hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');

  junk = cmip6.trend_gas_3; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2);   hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
                                   subplot(122); hold on; semilogy(nanstd(junk,[],2)',Tlevs,'linewidth',2); hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');
end

if exist('era5')
  iFig = 50;
  Tlevs = plays;
  Qlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end

  junk = era5.trend_gas_1; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Qlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
                                   subplot(122); hold on; semilogy(1+nanstd(junk,[],2)',Qlevs,'linewidth',2);     hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');

  junk = era5.trend_RH; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Qlevs,'linewidth',2);        hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
                                   subplot(122); hold on; semilogy(nanstd(junk,[],2)',Qlevs,'linewidth',2);      hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');

  junk = era5.trend_ptemp; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
                                   subplot(122); hold on; semilogy(nanstd(junk,[],2)',Tlevs,'linewidth',2);       hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');

  junk = era5.trend_gas_3; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2);   hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
                                   subplot(122); hold on; semilogy(nanstd(junk,[],2)',Tlevs,'linewidth',2); hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');
end

iFig = 50;
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend('UMBC','UMBC+unc','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*2e-2)
                             subplot(122); hold on; plotaxis2; hl = legend('UMBC','UMBC+unc','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 0.01]/2+1)
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend('UMBC','UMBC+unc','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*0.25)
                             subplot(122); hold on; plotaxis2; hl = legend('UMBC','UMBC+unc','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 1]/2)
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend('UMBC','UMBC+unc','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*0.1) 
                             subplot(122); hold on; plotaxis2; hl = legend('UMBC','UMBC+unc','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 1]*0.05)
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend('UMBC','UMBC+unc','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*0.04)
                             subplot(122); hold on; plotaxis2; hl = legend('UMBC','UMBC+unc','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 1]*0.04)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('nwp_spectral_trends_umbc_only_unc_P')
  figure(9)

  %nwp_spectral_trends_umbc_only     = make_profile_spectral_trends_umbc_only(results,resultsWV,resultsT,resultsO3,pavg);
  %  sum(sum(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates-nwp_spectral_trends_umbc_only.umbc_spectral_rates))
  nwp_spectral_trends_umbc_only_unc_P = make_profile_spectral_trends_umbc_only(results+resultsunc,resultsWV+resultsWVunc,resultsT+resultsTunc,resultsO3+resultsO3unc,pavg);
  nwp_spectral_trends_umbc_only_unc_M = make_profile_spectral_trends_umbc_only(results-resultsunc,resultsWV-resultsWVunc,resultsT-resultsTunc,resultsO3-resultsO3unc,pavg);
  plot(h.vchan,nanmean(rates'),'k',h.vchan,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates'),'b',...
       h.vchan,nanmean(nwp_spectral_trends_umbc_only_unc_P.umbc_spectral_rates'),'c',h.vchan,nanmean(nwp_spectral_trends_umbc_only_unc_M.umbc_spectral_rates'),'m')
  axis([640 1640 -0.1 +0.1]); plotaxis2; hl = legend('Obs AIRS BT trends','OEM','OEM + unc','OEM - unc','location','best','fontsize',10);

  figure(9); clf
  plot(f,nanmean(rates(:,mask),2),'k',f,nanmean(fits(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),'linewidth',2)
  plot(f,nanmean(rates(:,mask),2),'k',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),'linewidth',2)
  plotaxis2; hl = legend('AIRS Obs','UMBC','AIRS L3','CMIP6','ERA5','location','best'); axis([640 1640 -0.1 0.05]); title('Spectral Rates');

  plot(f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates(:,mask),2),'x-',...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,mask),2),...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),f,nanmean(rates(:,mask),2),'linewidth',2)
  hold on
  plot(f,nanmean(nwp_spectral_trends_umbc_only_unc_P.umbc_spectral_rates,2),'color',[0.85 0.85 0.85],'linewidth',2);
  plot(f,nanmean(nwp_spectral_trends_umbc_only_unc_M.umbc_spectral_rates,2),'color',[0.95 0.95 0.95],'linewidth',2)
  hold off
  plotaxis2; hl = legend('UMBC','AIRS L3','CMIP6','ERA5','AIRS Obs','UMBC+unc','UMBC-unc','location','best','fontsize',8); axis([640 1640 -0.1 0.05]); title('Spectral Rates');

  disp('making d/dtRH = 0 calcs')
  %deltaTT = deltaT+deltaTunc;
  %deltaTT = deltaTunc;
  %for ii = 1 : 4608
  %  nlevs = p.nlevs(ii);
  %  deltaTT(nlevs:101,ii) = nan;
  %end
  compute_deltaRH.umbc_unc_M = const_rh_4608(results(:,6)'-resultsunc(:,6)',deltaT-(deltaTunc-deltaT),fracWV-(fracWVunc-fracWV),rates,mask); 
  compute_deltaRH.umbc_unc_P = const_rh_4608(results(:,6)'+resultsunc(:,6)',deltaT+(deltaTunc-deltaT),fracWV+(fracWVunc-fracWV),rates,mask); 
end
