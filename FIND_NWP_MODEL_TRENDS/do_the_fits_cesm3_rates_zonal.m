warning off

zonkA = find(yy == StartY & mm == StartYM);
zonkB = find(yy == StopY  & mm == StopYM);
%zonk = 1:zonk;
zonk = zonkA : zonkB;
fprintf(1,'going from %4i/%2i to %4i/%2i .. found %3i time points \n',StartY,StartYM,StopY,StopYM,length(zonk))

thestats  = do_profilerate_fit(save_Q(:,:,zonk),save_T(:,:,zonk),save_stemp(:,zonk),...
                                   days(zonk),latbins);
xthestats = do_profilerate_fit(save_O3(:,:,zonk),save_RH(:,:,zonk),save_stemp(:,zonk),...
                                   days(zonk),latbins);

thestats.ozonerate    = xthestats.waterrate;
thestats.ozoneratestd = xthestats.waterratestd;

thestats.ozonerate        = xthestats.waterrate;
thestats.ozoneratestd     = xthestats.waterratestd;
thestats.ozonelag         = xthestats.waterlag;
thestats.ozoneratestd_lag = xthestats.waterratestd_lag;

thestats.RHrate        = xthestats.ptemprate;
thestats.RHratestd     = xthestats.ptempratestd;
thestats.RHlag         = xthestats.ptemplag;
thestats.RHratestd_lag = xthestats.ptempratestd_lag;

%clear xthestats
comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compute_cesm3_trends_desc_or_asc.m';

if iDorA > 0
  saver = ['save /asl/s1/sergio/CESM3/cesm3_zonal_rates_stats_' savestr_version '.mat thestats Tlevs Qlevs latbins save_lat zonk comment'];
end
eval(saver)

load llsmap5

figure(1); clf; plot(thestats.lats,thestats.stemprate); plotaxis2; title('CESM dST/dt')

junkP = squeeze(Tlevs); junkP = squeeze(junkP(1:58,90,140));
figure(2); clf; pcolor(thestats.lats,junkP/100,thestats.ptemprate'); shading interp; colormap(llsmap5); set(gca,'ydir','reverse');  set(gca,'yscale','log'); colorbar
title('CESM3 dT/dt K/yr')
caxis([-1 +1]*0.15)
ylim([10 1000])

junkP = squeeze(Tlevs); junkP = squeeze(junkP(1:58,90,140));
figure(3); clf; pcolor(thestats.lats,junkP/100,thestats.waterrate'); shading interp; colormap(llsmap5); set(gca,'ydir','reverse');  set(gca,'yscale','log'); colorbar
title('CESM3 dWVfrac/dt 1/yr')
caxis([-1 +1]*0.015)
ylim([100 1000])

junkP = squeeze(Tlevs); junkP = squeeze(junkP(1:58,90,140));
figure(4); clf; pcolor(thestats.lats,junkP/100,thestats.ozonerate'); shading interp; colormap(llsmap5); set(gca,'ydir','reverse');  set(gca,'yscale','log'); colorbar
title('CESM3 dO3frac/dt 1/yr')
caxis([-1 +1]*0.015)
ylim([10 1000])

junkP = squeeze(Tlevs); junkP = squeeze(junkP(1:58,90,140));
figure(5); clf; pcolor(thestats.lats,junkP/100,thestats.RHrate'); shading interp; colormap(llsmap5); set(gca,'ydir','reverse');  set(gca,'yscale','log'); colorbar
title('CESM3 dRH/dt percent/yr')
caxis([-1 +1]*0.15)
ylim([100 1000])
