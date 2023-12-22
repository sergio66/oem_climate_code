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

aslmap(89,rlat65,rlon73,smoothn(reshape(ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.5); title('CERES dILR/dt');

figure(89); clf
plot(fERA5_night.trend_stemp(:),ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608(:),'.')
[nn,nx,ny,nmeanILR,nstdILR] = myhist2d(fERA5_night.trend_stemp(:),ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608(:),-1:0.05:+1,-1:0.05:+1);
[nn,nx,ny,nmeanOLR,nstdOLR] = myhist2d(fERA5_night.trend_stemp(:),ceres_olr.ceres_trend.trend_toa_lw_clr_4608(:),-1:0.05:+1,-1:0.05:+1);
errorbar(-1:0.05:+1,nmeanOLR,nstdOLR,'color','b','linewidth',2);; hold on
errorbar(-1:0.05:+1,nmeanILR,nstdILR,'color','r','linewidth',2);; hold off
plotaxis2; xlabel('\delta SKT [K/yr]'); ylabel('\delta Flux [W/m2/yr]'); hl = legend('OLR','ILR','location','best','fontsize',10);


[nn,nx,ny,nmeanILR,nstdILR] = myhist2d(fERA5_night.trend_stemp(ocean),ceres_ilr.ceres_trend.trend_sfc_lw_clr_4608(ocean),-1:0.01:+1,-1:0.01:+1);
[nn,nx,ny,nmeanOLR,nstdOLR] = myhist2d(fERA5_night.trend_stemp(ocean),ceres_olr.ceres_trend.trend_toa_lw_clr_4608(ocean),-1:0.01:+1,-1:0.01:+1);
errorbar(-1:0.01:+1,nmeanOLR,nstdOLR,'color','b','linewidth',2);; hold on
errorbar(-1:0.01:+1,nmeanILR,nstdILR,'color','r','linewidth',2);; hold off
plotaxis2; xlabel('\delta SKT [K/yr]'); ylabel('\delta Flux [W/m2/yr]'); hl = legend('OLR','ILR','location','best','fontsize',10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
