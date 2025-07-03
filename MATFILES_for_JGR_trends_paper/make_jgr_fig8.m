load fig8.mat

  figure(1); clf
  plot(data8.rlat,nanmean(data8.umbc,1),'k',data8.rlat,nanmean(data8.airsv7,1),'b',data8.rlat,nanmean(data8.climcaps,1),'g',data8.rlat,nanmean(data8.era5,1),'r',data8.rlat,nanmean(data8.merra2),'m',data8.rlat,nanmean(data8.giss),'c','linewidth',2);
  xlim([-1 +1]*90); ylim([-0.06 +0.11]);
  plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','GISS','location','best','fontsize',8);

  %plot(data8.rlat,nanmean(data8.umbc,1),'k',data8.rlat,nanmean(data8.airsv7,1),'b',data8.rlat,nanmean(data8.climcaps,1),'g',data8.rlat,nanmean(data8.era5,1),'r',data8.rlat,nanmean(data8.merra2),'m','linewidth',2);
  %plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);

  xlabel('Latitude [deg]'); ylabel('dST/dt [K yr^{-1}]');

  %%%%%%%%%%%%%%%%%%%%%%%%%

  figure(2); clf
  plot(data8.rlat,nanmean(data8.ocean_umbc,1),'k',data8.rlat,nanmean(data8.ocean_airsv7,1),'b',data8.rlat,nanmean(data8.ocean_climcaps,1),'g',...
      data8.rlat,nanmean(data8.ocean_era5,1),'r',data8.rlat,nanmean(data8.ocean_merra2),'m',data8.rlat,nanmean(data8.ocean_giss),'c','linewidth',2);
  xlim([-1 +1]*90); ylim([-0.06 +0.11]);
  plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','GISS','location','best','fontsize',8);

  %plot(data8.rlat,nanmean(data8.ocean_umbc,1),'k',data8.rlat,nanmean(data8.ocean_airsv7,1),'b',data8.rlat,nanmean(data8.ocean_climcaps,1),'g',data8.rlat,nanmean(data8.ocean_era5,1),'r',data8.rlat,nanmean(data8.ocean_merra2),'m','linewidth',2);
  %plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);

  xlabel('Latitude [deg]'); ylabel('dST_{ocean}/dt [K yr^{-1}]')
