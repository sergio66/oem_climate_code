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

if ~exist('iFig')
  iFig = 102;
end

guess_delta_ILR
guess_delta_ILR_ERA5
guess_delta_ILR_MERRA2
guess_delta_ILR_AIRSL3
guess_delta_ILR_CLIMCAPSL3

if exist('fERA5_night')
  if isfield(fERA5_night,'ilr_adj')
    disp(' << see ../AIRS_gridded_STM_May2021/driver_processs_aeri_uplook_40_4608.m  for polynomial coeffs for CO2 >> ')
    co2_ilr_adjustA = 0.10 + 0.20*(sin(rlat*pi/180)).^2; 
    co2_ilr_adjustB = 0.06 + 0.17*(sin(rlat*pi/180)).^2; 
    co2_ilr_adjustC = 0.0562 - 0.0595*(sin(rlat*pi/180)) + 0.1771*(sin(rlat*pi/180)).^2; 
    co2_ilr_adjustD = 0.0562 - 0.0486*(sin(rlat*pi/180)) + 0.1771*(sin(rlat*pi/180)).^2 - 0.0149*(sin(rlat*pi/180)).^3; 

    co2_ilr_adjust  = co2_ilr_adjustD;

    figure(101); clf
      brutsaert_era5 =  nanmean(reshape(fERA5_night.ilr_Rld,72,64),1);
      plot(rlat,nanmean(reshape(fERA5_night.ilr_adj,72,64),1),'b',rlat,brutsaert_era5,rlat,brutsaert_era5 .* (1 + co2_ilr_adjust'),'linewidth',2); 
      hl = legend('ERA5 ILR --> ERA5 ILR + stemp adj','ERA5 Brutsaert','ERA5 Brutsaert .* (1 + co2adjust)','location','best');

    figure(102); clf
    plot(rlat,co2_ilr_adjustA,rlat,co2_ilr_adjustB,rlat,co2_ilr_adjustC,rlat,co2_ilr_adjustD,'linewidth',2); 
    hl=legend('Orig40 0,1,0.2','Better40 0.05, 0.17','Quadratic 4608','Tertac 4608','location','best','fontsize',10);

    figure(103)
    deltaILR_chirp = deltaILR;    %% dT2m/dt
    deltaILR_chirp = deltaILR_ST; %% dSKT/dt
    wahoo = ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608;
    plot(rlat,nanmean(reshape(deltaILR_chirp,72,64),1),'k',...
         rlat,nanmean(reshape(wahoo,72,64),1),'gx-',rlat,nanmean(reshape(fERA5_night.ilr_clr,72,64),1),'c',rlat,nanmean(reshape(fERA5_night.ilr_adj,72,64),1),'b',...
         rlat,nanmean(reshape(fERA5_night.ilr_Rld,72,64),1),'r','linewidth',2);
    plotaxis2;  hl = legend('CHIRP\_A','CERES','ERA5 net?','ERA5 adj','Brutsaert','location','best','fontsize',10); xlabel('Latitude'); ylabel('dILR/dt W/m2/yr'); xlim([-90 +90])
    plot(rlat,nanmean(reshape(deltaILR_chirp,72,64),1),'k',...
         rlat,nanmean(reshape(wahoo,72,64),1),'gx-',rlat,nanmean(reshape(fERA5_night.ilr_clr,72,64),1),'c',rlat,nanmean(reshape(fERA5_night.ilr_adj,72,64),1),'b',...
         rlat,brutsaert_era5,'r','linewidth',2);
    plotaxis2;  hl = legend('CHIRP\_A','CERES','ERA5 net?','ERA5 adj','Brutsaert','location','best','fontsize',10); xlabel('Latitude'); ylabel('dILR/dt W/m2/yr'); xlim([-90 +90])

    plot(rlat,nanmean(reshape(deltaILR_chirp,72,64),1) .* (1 + co2_ilr_adjust') ,'k',...
         rlat,nanmean(reshape(wahoo,72,64),1),'gx-',rlat,nanmean(reshape(fERA5_night.ilr_adj,72,64),1),'b',...
         rlat,brutsaert_era5 .* (1 + co2_ilr_adjust'),'r','linewidth',2);
    plotaxis2;  hl = legend('CHIRP\_A with CO2 adj','CERES','ERA5 ILR --> ERA5 ILR + stemp adj','ERA5 Brutsaert with CO2 adj','location','best','fontsize',10); xlabel('Latitude'); 
    ylabel('dILR/dt W/m2/yr'); xlim([-90 +90])
    title('ILR trends')
    ylim([0 0.4])

    plot(rlat,nanmean(reshape(deltaILR_ST,72,64),1) .* (1 + co2_ilr_adjust') ,'b',rlat,nanmean(reshape(deltaILR,72,64),1) .* (1 + co2_ilr_adjust') ,'b--',...
         rlat,nanmean(reshape(wahoo,72,64),1),'gx-',...
         rlat,nanmean(reshape(fERA5_night.ilr_adj,72,64),1),'r',rlat,brutsaert_era5 .* (1 + co2_ilr_adjust'),'r--','linewidth',2);
    plotaxis2;  
    hl = legend('CHIRP\_A Brutsaert using dSKT/dt, CO2 adj','CHIRP\_A Brutsaert using dT2m/dt, CO2 adj',...
                 'CERES','ERA5 ILR --> ERA5 ILR + stemp adj','ERA5 Brutsaert with CO2 adj','location','best','fontsize',10); xlabel('Latitude'); 
    ylabel('dILR/dt W/m2/yr'); xlim([-90 +90])
    title('ILR trends')
    ylim([-0.2 0.4])

    figure(104); clf
    aslmap(104,rlat65,rlon73,smoothn(reshape(deltaILR_ST,72,64)',1), [-90 +90],[-180 +180]); title('CHIRP\_A d(ILR)/dt (W/m2/yr)'); caxis([-1 +1]*0.5); colormap(llsmap5)
    figure(104); clf
    co2_ilr_adjustx = 0.1 + 0.2*(sin(p.rlat*pi/180)).^2; %% this is ILR flux adjust due to CO2 changes .. 10% at humid tropics, 30% at dry poles
    aslmap(104,rlat65,rlon73,smoothn(reshape(deltaILR_ST.*(1 + co2_ilr_adjustx),72,64)',1), [-90 +90],[-180 +180]); title('CHIRP\_A d(ILR)/dt (W/m2/yr)'); caxis([-1 +1]*0.5); colormap(llsmap5)

    figure(105); clf
    aslmap(105,rlat65,rlon73,smoothn(reshape(wahoo,72,64)',1), [-90 +90],[-180 +180]); title('CERES d(ILR)/dt (W/m2/yr)'); caxis([-1 +1]*0.5); colormap(llsmap5)

    figure(106); clf
    aslmap(106,rlat65,rlon73,smoothn(reshape(fERA5_night.ilr_Rld,72,64)',1), [-90 +90],[-180 +180]); title('ERA5 d(ILR)/dt (W/m2/yr)'); caxis([-1 +1]*0.5); colormap(llsmap5)

   junkind = find(abs(p.rlat) < 60);
   %junkind = 1:4608;
   [mdT2m(1,1),sdT2m(1,1),mdT2m(1,2),sdT2m(1,2)] = weighted_mean_stddev(dTx2m_ERA5_orig(junkind),cos(p.rlat(junkind)*pi/180));
   [mdT2m(2,1),sdT2m(2,1),mdT2m(2,2),sdT2m(2,2)] = weighted_mean_stddev(fERA5_night.t2m(junkind),cos(p.rlat(junkind)*pi/180));
   [mdT2m(3,1),sdT2m(3,1),mdT2m(3,2),sdT2m(3,2)] = weighted_mean_stddev(dTx2m_MERRA2_orig(junkind),cos(p.rlat(junkind)*pi/180));
   [mdT2m(4,1),sdT2m(4,1),mdT2m(4,2),sdT2m(4,2)] = weighted_mean_stddev(dTx2m_AIRSL3_orig(junkind),cos(p.rlat(junkind)*pi/180));
   [mdT2m(5,1),sdT2m(5,1),mdT2m(5,2),sdT2m(5,2)] = weighted_mean_stddev(dTx2m_CLIMCAPSL3_orig(junkind),cos(p.rlat(junkind)*pi/180));
   [mdT2m(6,1),sdT2m(6,1),mdT2m(6,2),sdT2m(6,2)] = weighted_mean_stddev(dTx2m_orig(junkind),cos(p.rlat(junkind)*pi/180));
   [mdT2m(7,1),sdT2m(7,1),mdT2m(7,2),sdT2m(7,2)] = weighted_mean_stddev(results(junkind,6),cos(p.rlat(junkind)*pi/180));
   junkfrac = 0;
   junkfrac = 1;
   figure(107); clf
   hold off; plot(rlat,nanmean(reshape(dTx2m_ERA5_orig,72,64),1),'m',rlat,junkfrac*nanmean(reshape(fERA5_night.t2m,72,64),1),'r','linewidth',2); 
   hold on;  plot(rlat,nanmean(reshape(dTx2m_MERRA2_orig,72,64),1),'color',[1 1 1]*0.8,'linewidth',2);
   hold on;  plot(rlat,nanmean(reshape(dTx2m_AIRSL3_orig,72,64),1),'g',rlat,nanmean(reshape(dTx2m_CLIMCAPSL3_orig,72,64),1),'g--','linewidth',2);
   hold on;  plot(rlat,nanmean(reshape(dTx2m_orig,72,64),1),'b',rlat,junkfrac*nanmean(reshape(results(:,6),72,64),1),'cx-','linewidth',4);
   hold off
   ylim([-0.05 0.075])
   plotaxis2;
   hl = legend('ERA5 interped T2m rate','ERA5 actual monthly T2m rate','MERRA2 interped T2m rate','AIRSL3 interped T2m rate','CLIMCAPSL3 interped T2m rate',...
               'CHIRP\_A interped T2m rate','CHIRP\_A stemp rate','location','best','fontsize',8);

   junkind = find(abs(p.rlat) < 60);
   %junkind = 1:4608;
   [mdW2m(1,1),sdW2m(1,1),mdW2m(1,2),sdW2m(1,2)] = weighted_mean_stddev(dWVsurf_ERA5_orig(junkind),cos(p.rlat(junkind)*pi/180));
   [mdW2m(2,1),sdW2m(2,1),mdW2m(2,2),sdW2m(2,2)] = weighted_mean_stddev(fERA5_night.frac_e2m(junkind),cos(p.rlat(junkind)*pi/180));
   [mdW2m(3,1),sdW2m(3,1),mdW2m(3,2),sdW2m(3,2)] = weighted_mean_stddev(dWVsurf_MERRA2_orig(junkind),cos(p.rlat(junkind)*pi/180));
   [mdW2m(4,1),sdW2m(4,1),mdW2m(4,2),sdW2m(4,2)] = weighted_mean_stddev(dWVsurf_AIRSL3_orig(junkind),cos(p.rlat(junkind)*pi/180));
   [mdW2m(5,1),sdW2m(5,1),mdW2m(5,2),sdW2m(5,2)] = weighted_mean_stddev(dWVsurf_CLIMCAPSL3_orig(junkind),cos(p.rlat(junkind)*pi/180));
   [mdW2m(6,1),sdW2m(6,1),mdW2m(6,2),sdW2m(6,2)] = weighted_mean_stddev(dWVsurf_orig(junkind),cos(p.rlat(junkind)*pi/180));
   [mdW2m(7,1),sdW2m(7,1),mdW2m(7,2),sdW2m(7,2)] = weighted_mean_stddev(resultsWV(junkind,49),cos(p.rlat(junkind)*pi/180));
   figure(108); clf
   hold off; plot(rlat,10*nanmean(reshape(dWVsurf_ERA5_orig,72,64),1),'c-.',rlat,10*nanmean(reshape(fERA5_night.frac_e2m,72,64),1),'b-.','linewidth',4)
   hold on;  plot(rlat,10*nanmean(reshape(dWVsurf_orig,72,64),1),'g-.','linewidth',4);
   hold off
   ylim([0 0.075])
   plotaxis2;
   hl = legend('ERA5 interped fracWV 2m rate  x10','ERA5 actual monthly fracWV 2m rate x10','CHIRP\_A interped fracWV 2m rate x10','location','best','fontsize',8);

  end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
