figure(86); 
scatter_coast(p.rlon,p.rlat,100,nansum(Q0,1)); colormap jet
scatter_coast(p.rlon,p.rlat,100,nansum(fERA5_night.MSEtrend,1)); colormap(usa2); caxis([-1 +1]*5e6)
wah = nanmean(reshape(nansum(fERA5_night.MSEtrend,1),72,64),1);;
plot(rlat,wah); xlabel('rlat'); ylabel('ERA5 dMSE/dt')
junk = nansum(fERA5_night.MSEtrend,1)/1e6; 
aslmap(86,rlat65,rlon73,smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*2); 
title('ERA5 night dMSE/dt /1e6');

junk = nanmean(fERA5_night.ptemptrend(60:100,:),1); 
aslmap(87,rlat65,rlon73,smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.15); 
title('ERA5 night LowerTrop dT/dt');

junk = nanmean(fERA5_night.ptemptrend(60:100,:),1);
junk = junk.*fERA5_night.mmwtrend; 
aslmap(88,rlat65,rlon73,smoothn(reshape(junk,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.01); 
title('ERA5 night LowerTrop dT/dt .* dmmw/dt');

%%%%%%%%%%%%%%%%%%%%%%%%%

aslmap(89,rlat65,rlon73,smoothn(reshape(ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.5); title('CERES dILR/dt');

figure(89); clf
plot(fERA5_night.trend_stemp(:),ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608(:),'.')
[nn,nx,ny,nmeanILR,nstdILR] = myhist2d(fERA5_night.trend_stemp(:),ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608(:),-1:0.05:+1,-1:0.05:+1);
[nn,nx,ny,nmeanOLR,nstdOLR] = myhist2d(fERA5_night.trend_stemp(:),ceres_olr.ceres_trend.trend_toa_lw_clr_4608(:),-1:0.05:+1,-1:0.05:+1);
errorbar(-1:0.05:+1,nmeanOLR,nstdOLR,'color','b','linewidth',2);; hold on
errorbar(-1:0.05:+1,nmeanILR,nstdILR,'color','r','linewidth',2);; hold off
plotaxis2; xlabel('\delta SKT [K/yr]'); ylabel('\delta Flux [W/m2/yr]'); hl = legend('OLR','ILR','location','best','fontsize',10); title('LAND and OCEAN')

[nn,nx,ny,nmeanILR,nstdILR] = myhist2d(fERA5_night.trend_stemp(ocean),ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608(ocean),-1:0.01:+1,-1:0.01:+1);
[nn,nx,ny,nmeanOLR,nstdOLR] = myhist2d(fERA5_night.trend_stemp(ocean),ceres_olr.ceres_trend.trend_toa_lw_clr_4608(ocean),-1:0.01:+1,-1:0.01:+1);
errorbar(-1:0.01:+1,nmeanOLR,nstdOLR,'color','b','linewidth',2);; hold on
errorbar(-1:0.01:+1,nmeanILR,nstdILR,'color','r','linewidth',2);; hold off
plotaxis2; xlabel('\delta SKT [K/yr]'); ylabel('\delta Flux [W/m2/yr]'); hl = legend('OLR','ILR','location','best','fontsize',10); title('OCEAN')

%%%%%%%%%%%%%%%%%%%%%%%%%

wahooO = ceres_olr.ceres_trend.trend_toa_lw_clr_4608; wahooO(land) = nan;
wahooL = ceres_olr.ceres_trend.trend_toa_lw_clr_4608; wahooL(ocean) = nan;
wahoo = ceres_olr.ceres_trend.trend_toa_lw_clr_4608; 
plot(rlat,nanmean(reshape(wahoo,72,64),1),'k--',rlat,nanmean(reshape(wahooL,72,64),1),'m--',rlat,nanmean(reshape(wahooO,72,64),1),'c--','linewidth',2); 
hold on

wahooO = ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608; wahooO(land) = nan;
wahooL = ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608; wahooL(ocean) = nan;
wahoo = ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608; 
plot(rlat,nanmean(reshape(wahoo,72,64),1),'gx-',rlat,nanmean(reshape(wahooL,72,64),1),'r',rlat,nanmean(reshape(wahooO,72,64),1),'b','linewidth',2); 
hold off

plotaxis2;  hl = legend('OLR ALL','OLR Land','OLR Ocean','ILR ALL','ILR Land','ILR Ocean','location','best','fontsize',10); xlabel('Latitude'); ylabel('dILR/dt W/m2/yr'); xlim([-90 +90])
title('CERES OLR/ILR trends')

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(101); clf; 
wahoo = ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608; 
plot(rlat,nanmean(reshape(wahoo,72,64),1),'gx-',rlat,nanmean(reshape(fERA5_night.ilr_clr,72,64),1),'k',...
     rlat,nanmean(reshape(fERA5_night.ilr_adj,72,64),1),'b',rlat,nanmean(reshape(fERA5_night.ilr_Rld,72,64),1),'r','linewidth',2); 
plotaxis2;  hl = legend('CERES','ERA5 net?','ERA5 adj','Brutsaert','location','best','fontsize',10); xlabel('Latitude'); ylabel('dILR/dt W/m2/yr'); xlim([-90 +90])
title('ILR trends')

guess_delta_ILR
guess_delta_ILR_ERA5

if exist('fERA5_night')
  if isfield(fERA5_night,'ilr_adj')
    figure(101); clf
    plot(rlat,0.1 + 0.2*(sin(rlat*pi/180)).^2); 
    co2_ilr_adjust = 0.1 + 0.2*(sin(rlat*pi/180)).^2; %% this is ILR flux adjust due to CO2 changes .. 10% at humid tropics, 30% at dry poles
    brutsaert_era5 =  nanmean(reshape(fERA5_night.ilr_Rld,72,64),1);
    plot(rlat,brutsaert_era5,rlat,brutsaert_era5 .* (1 + co2_ilr_adjust'));

    wahoo = ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608;
    plot(rlat,nanmean(reshape(deltaILR,72,64),1),'k',...
         rlat,nanmean(reshape(wahoo,72,64),1),'gx-',rlat,nanmean(reshape(fERA5_night.ilr_clr,72,64),1),'c',rlat,nanmean(reshape(fERA5_night.ilr_adj,72,64),1),'b',...
         rlat,nanmean(reshape(fERA5_night.ilr_Rld,72,64),1),'r','linewidth',2);
    plotaxis2;  hl = legend('CHIRP\_A','CERES','ERA5 net?','ERA5 adj','Brutsaert','location','best','fontsize',10); xlabel('Latitude'); ylabel('dILR/dt W/m2/yr'); xlim([-90 +90])
    plot(rlat,nanmean(reshape(deltaILR,72,64),1),'k',...
         rlat,nanmean(reshape(wahoo,72,64),1),'gx-',rlat,nanmean(reshape(fERA5_night.ilr_clr,72,64),1),'c',rlat,nanmean(reshape(fERA5_night.ilr_adj,72,64),1),'b',...
         rlat,brutsaert_era5,'r','linewidth',2);
    plotaxis2;  hl = legend('CHIRP\_A','CERES','ERA5 net?','ERA5 adj','Brutsaert','location','best','fontsize',10); xlabel('Latitude'); ylabel('dILR/dt W/m2/yr'); xlim([-90 +90])

    plot(rlat,nanmean(reshape(deltaILR,72,64),1) .* (1 + co2_ilr_adjust') ,'k',...
         rlat,nanmean(reshape(wahoo,72,64),1),'gx-',rlat,nanmean(reshape(fERA5_night.ilr_adj,72,64),1),'b',...
         rlat,brutsaert_era5 .* (1 + co2_ilr_adjust'),'r','linewidth',2);
    plotaxis2;  hl = legend('CHIRP\_A with CO2 adj','CERES','ERA5 --> ERA5 + stemp adj','ERA5 Brutsaert with CO2 adj','location','best','fontsize',10); xlabel('Latitude'); 
    ylabel('dILR/dt W/m2/yr'); xlim([-90 +90])
    title('ILR trends')
    ylim([0 0.4])

    figure(102); clf
    aslmap(102,rlat65,rlon73,smoothn(reshape(deltaILR,72,64)',1), [-90 +90],[-180 +180]); title('CHIRP\_A d(ILR)/dt (W/m2/yr)'); caxis([-1 +1]*0.5); colormap(llsmap5)
    figure(102); clf
    co2_ilr_adjustx = 0.1 + 0.2*(sin(p.rlat*pi/180)).^2; %% this is ILR flux adjust due to CO2 changes .. 10% at humid tropics, 30% at dry poles
    aslmap(102,rlat65,rlon73,smoothn(reshape(deltaILR.*(1 + co2_ilr_adjustx),72,64)',1), [-90 +90],[-180 +180]); title('CHIRP\_A d(ILR)/dt (W/m2/yr)'); caxis([-1 +1]*0.5); colormap(llsmap5)

    figure(103); clf
    aslmap(103,rlat65,rlon73,smoothn(reshape(wahoo,72,64)',1), [-90 +90],[-180 +180]); title('CERES d(ILR)/dt (W/m2/yr)'); caxis([-1 +1]*0.5); colormap(llsmap5)

    figure(104); clf
    aslmap(104,rlat65,rlon73,smoothn(reshape(fERA5_night.ilr_Rld,72,64)',1), [-90 +90],[-180 +180]); title('ERA5 d(ILR)/dt (W/m2/yr)'); caxis([-1 +1]*0.5); colormap(llsmap5)

   figure(105); clf
   plot(rlat,nanmean(reshape(dTx2m_ERA5_orig,72,64),1),'m',rlat,nanmean(reshape(fERA5_night.t2m,72,64),1),'r',rlat,nanmean(reshape(dTx2m_orig,72,64),1),'k','linewidth',2);
   hold on
   plot(rlat,10*nanmean(reshape(dWVsurf_ERA5_orig,72,64),1),'c-.',rlat,10*nanmean(reshape(fERA5_night.frac_e2m,72,64),1),'b-.',rlat,10*nanmean(reshape(dWVsurf_orig,72,64),1),'g-.','linewidth',4);
   hold off

   ylim([0 0.075])
   plotaxis2;
   hl = legend('ERA5 interped T2m rate','ERA5 actual monthly T2m rate','CHIRP\_A interped T2m rate',...
               'ERA5 interped fracWV 2m rate  x10','ERA5 actual monthly fracWV 2m rate x10','CHIRP\_A interped fracWV 2m rate x10','location','best','fontsize',10);

  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
