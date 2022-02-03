%% see plot_profile_trends2.m

iFig = 34;
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(mncos100.*nanmean(coswgt100.*fracWV(1:100,xmask),2),plays,'linewidth',2);         ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
                                     subplot(122); semilogy(1+mncos100.*nanstd(coswgt100.*fracWV(1:100,xmask),[],2),plays,'linewidth',2);     ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(mncos100.*nanmean(coswgt100.*deltaRH(1:100,xmask),2),plays,'linewidth',2);        ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
                                     subplot(122); semilogy(mncos100.*nanstd(coswgt100.*deltaRH(1:100,xmask),[],2),plays,'linewidth',2);      ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(mncos100.*nanmean(coswgt100.*deltaT(1:100,xmask),2),plays,'linewidth',2);         ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
                                     subplot(122); semilogy(mncos100.*nanstd(coswgt100.*deltaT(1:100,xmask),[],2),plays,'linewidth',2);       ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(mncos097.*nanmean(coswgt097.*deltaO3(1:97,xmask),2),plays(1:97),'linewidth',2);   ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
                                     subplot(122); semilogy(mncos097.*nanstd(coswgt097.*deltaO3(1:97,xmask),[],2),plays(1:97),'linewidth',2); ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');

if exist('airsL3')
  iFig = 34;
  Tlevs = airsL3.Tlevs;
  Qlevs = airsL3.Qlevs;

  boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = xmaskLFmatr'; end
  junk = airsL3.thestats64x72.waterrate.*boo; junk = reshape(junk,72*64,12);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos012.*nanmean(coswgt012'.*junk(xmask,:),1)',Qlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
                                   subplot(122); hold on; semilogy(1+mncos012.*nanstd(coswgt012'.*junk(xmask,:),[],1)',Qlevs,'linewidth',2);     hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');

  boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = xmaskLFmatr'; end
  junk = airsL3.thestats64x72.RHrate.*boo; junk = reshape(junk,72*64,12);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos012.*nanmean(coswgt012'.*junk(xmask,:),1)',Qlevs,'linewidth',2);        hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
                                   subplot(122); hold on; semilogy(mncos012.*nanstd(coswgt012'.*junk(xmask,:),[],1)',Qlevs,'linewidth',2);      hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');

  boo = zeros(72,64,12); for ijunk = 1 : 24; boo(:,:,ijunk) = xmaskLFmatr'; end
  junk = airsL3.thestats64x72.ptemprate.*boo; junk = reshape(junk,72*64,24);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos024.*nanmean(coswgt024'.*junk(xmask,:),1)',Tlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
                                   subplot(122); hold on; semilogy(mncos024.*nanstd(coswgt024'.*junk(xmask,:),[],1)',Tlevs,'linewidth',2);       hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');

  boo = zeros(72,64,24); for ijunk = 1 : 24; boo(:,:,ijunk) = xmaskLFmatr'; end
  junk = airsL3.thestats64x72.ozonerate.*boo; junk = reshape(junk,72*64,24);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos024.*nanmean(coswgt024'.*junk(xmask,:),1)',Tlevs,'linewidth',2);   hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
                                   subplot(122); hold on; semilogy(mncos024.*nanstd(coswgt024'.*junk(xmask,:),[],1)',Tlevs,'linewidth',2); hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');
end

if exist('cmip6')
  iFig = 34;
  Tlevs = plays;
  Qlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end

  junk = cmip6.trend_gas_1; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Qlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
                                   subplot(122); hold on; semilogy(1+mncos100'.*nanstd(coswgt100.*junk(:,xmask),[],2)',Qlevs,'linewidth',2);     hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');

  junk = cmip6.trend_RH; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Qlevs,'linewidth',2);        hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
                                   subplot(122); hold on; semilogy(mncos100'.*nanstd(coswgt100.*junk(:,xmask),[],2)',Qlevs,'linewidth',2);      hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');

  junk = cmip6.trend_ptemp; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Tlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
                                   subplot(122); hold on; semilogy(mncos100'.*nanstd(coswgt100.*junk(:,xmask),[],2)',Tlevs,'linewidth',2);       hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');

  junk = cmip6.trend_gas_3; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Tlevs,'linewidth',2);   hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
                                   subplot(122); hold on; semilogy(mncos100'.*nanstd(coswgt100.*junk(:,xmask),[],2)',Tlevs,'linewidth',2); hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');
end

if exist('era5')
  iFig = 34;
  Tlevs = plays;
  Qlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = xmaskLFmatr'; end

  junk = era5.trend_gas_1; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Qlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
                                   subplot(122); hold on; semilogy(1+mncos100'.*nanstd(coswgt100.*junk(:,xmask),[],2)',Qlevs,'linewidth',2);     hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');

  junk = era5.trend_RH; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Qlevs,'linewidth',2);        hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
                                   subplot(122); hold on; semilogy(mncos100'.*nanstd(coswgt100.*junk(:,xmask),[],2)',Qlevs,'linewidth',2);      hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');

  junk = era5.trend_ptemp; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Tlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
                                   subplot(122); hold on; semilogy(mncos100'.*nanstd(coswgt100.*junk(:,xmask),[],2)',Tlevs,'linewidth',2);       hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');

  junk = era5.trend_gas_3; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(mncos100'.*nanmean(coswgt100.*junk(:,xmask),2)',Tlevs,'linewidth',2);   hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
                                   subplot(122); hold on; semilogy(mncos100'.*nanstd(coswgt100.*junk(:,xmask),[],2)',Tlevs,'linewidth',2); hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');
end

iFig = 34;
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend(ocbstr,'AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*1e-2)
                             subplot(122); hold on; plotaxis2; hl = legend(ocbstr,'AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 0.01]/2+1)
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend(ocbstr,'AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*0.25)
                             subplot(122); hold on; plotaxis2; hl = legend(ocbstr,'AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 1]/2)
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend(ocbstr,'AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*0.035) 
                             subplot(122); hold on; plotaxis2; hl = legend(ocbstr,'AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 1]*0.05)
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend(ocbstr,'AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*0.01/2)
                             subplot(122); hold on; plotaxis2; hl = legend(ocbstr,'AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 1]*0.01/2)
