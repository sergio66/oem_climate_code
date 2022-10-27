warning off

zonkA = find(yy == StartY & mm == StartYM);
zonkB = find(yy == StopY  & mm == StopYM);
%zonk = 1:zonk;
zonk = zonkA : zonkB;
fprintf(1,'going from %4i/%2i to %4i/%2i .. found %3i time points \n',StartY,StartYM,StopY,StopYM,length(zonk))

warning off
for ii = 1 : 72
  fprintf(1,'lonbin %2i of 72 \n',ii);
%  xthestats       = do_profilerate_fit_WV_T_O3_RH_ST(squeeze(save64x72_Q(:,ii,:,zonk)),squeeze(save64x72_T(:,ii,:,zonk)),squeeze(save64x72_RH(:,ii,:,zonk)),squeeze(save64x72_O3(:,ii,:,zonk)),...
%                                     squeeze(save64x72_stemp(:,ii,zonk)),squeeze(save64x72_RHSurf(:,ii,zonk)),squeeze(save64x72_TWetSurf(:,ii,zonk)),...
%                                     days(zonk),rlat);
  xthestats       = do_profilerate_fit_WV_T_O3_RH_ST(squeeze(save64x72_Q(:,ii,:,zonk)),squeeze(save64x72_T(:,ii,:,zonk)),squeeze(save64x72_RH(:,ii,:,zonk)),squeeze(save64x72_O3(:,ii,:,zonk)),...
                                     squeeze(save64x72_stemp(:,ii,zonk)),0*squeeze(save64x72_stemp(:,ii,zonk)),0*squeeze(save64x72_stemp(:,ii,zonk)),...
                                     days(zonk),rlat);

%  ythestats       = do_profilerate_fit_WV_T_O3_RH_ST(squeeze(save64x72_O3(:,ii,:,zonk)),squeeze(save64x72_T(:,ii,:,zonk)),squeeze(save64x72_RH(:,ii,:,zonk)),squeeze(save64x72_O3(:,ii,:,zonk)),...
%                                     squeeze(save64x72_stemp(:,ii,zonk)),squeeze(save64x72_RHSurf(:,ii,zonk)),squeeze(save64x72_TWetSurf(:,ii,zonk)),...
%                                     days(zonk),rlat);


  thestats64x72.lats = xthestats.lats;

  thestats64x72.waterrate(ii,:,:) = xthestats.waterrate;
  thestats64x72.waterratestd(ii,:,:) = xthestats.waterratestd;
  thestats64x72.waterlag(ii,:,:) = xthestats.waterlag;
  thestats64x72.waterratestd_lag(ii,:,:) = xthestats.waterratestd_lag;

  thestats64x72.ozonerate(ii,:,:) = xthestats.ozonerate;
  thestats64x72.ozoneratestd(ii,:,:) = xthestats.ozoneratestd;
  thestats64x72.ozonelag(ii,:,:) = xthestats.ozonelag;
  thestats64x72.ozoneratestd_lag(ii,:,:) = xthestats.ozoneratestd_lag;

  thestats64x72.RHrate(ii,:,:) = xthestats.RHrate;
  thestats64x72.RHratestd(ii,:,:) = xthestats.RHratestd;
  thestats64x72.RHlag(ii,:,:) = xthestats.RHlag;
  thestats64x72.RHratestd_lag(ii,:,:) = xthestats.RHratestd_lag;

  thestats64x72.ptemprate(ii,:,:) = xthestats.ptemprate;
  thestats64x72.ptempratestd(ii,:,:) = xthestats.ptempratestd;
  thestats64x72.ptemplag(ii,:,:) = xthestats.ptemplag;
  thestats64x72.ptempratestd_lag(ii,:,:) = xthestats.ptempratestd_lag;

  thestats64x72.stemprate(ii,:) = xthestats.stemprate;
  thestats64x72.stempratestd(ii,:) = xthestats.stempratestd;
  thestats64x72.stemplag(ii,:) = xthestats.stemplag;
  thestats64x72.stempratestd_lag(ii,:) = xthestats.stempratestd_lag;

  %thestats64x72.TWetSurfrate(ii,:) = xthestats.TWetSurfrate;
  %thestats64x72.TWetSurfratestd(ii,:) = xthestats.TWetSurfratestd;
  %thestats64x72.TWetSurflag(ii,:) = xthestats.TWetSurflag;
  %thestats64x72.TWetSurfratestd_lag(ii,:) = xthestats.TWetSurfratestd_lag;

  %thestats64x72.RHSurfrate(ii,:) = xthestats.RHSurfrate;
  %thestats64x72.RHSurfratestd(ii,:) = xthestats.RHSurfratestd;
  %thestats64x72.RHSurflag(ii,:) = xthestats.RHSurflag;
  %thestats64x72.RHSurfratestd_lag(ii,:) = xthestats.RHSurfratestd_lag;

end

addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load llsmap5
figure(1); pcolor(save_lat64x72,Tlevs/100,squeeze(nanmean(thestats64x72.ptemprate,1))'); shading interp; colorbar; caxis([-0.15 +0.15]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log')
  title('CESM2.X dT/dt K/yr');
figure(2); pcolor(save_lat64x72,Qlevs/100,squeeze(nanmean(thestats64x72.waterrate,1))'); shading interp; colorbar; caxis([-1 +1]*0.01); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','linear')
  title('CESM2.X d(WV)/<WV>/dt /yr');
figure(3); pcolor(save_lat64x72,Qlevs/100,squeeze(nanmean(thestats64x72.RHrate,1))'); shading interp; colorbar; caxis([-1 +1]*0.15); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','linear')
  title('CESM2.X d(RH)/dt /yr');

figure(1); ylim([1 1000]);
figure(2); ylim([100 1000]);
figure(3); ylim([100 1000]);

figure(4); pcolor(save_lon64x72,save_lat64x72,thestats64x72.stemprate'); shading interp; colorbar; caxis([-0.15 +0.15]); colormap(llsmap5);
figure(4); simplemap(thestats64x72.stemprate'); shading interp; colorbar; caxis([-0.15 +0.15]); colormap(llsmap5);
addpath /home/sergio/MATLABCODE/
addpath /asl/matlib/maps/
load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
%rlon = -180 : 5 : +180;  rlat = latB2;
%rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
%rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
figure(4); aslmap(4,rlat65,rlon73,smoothn(thestats64x72.stemprate',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('d/dt CESM3 K/yr'); colormap(llsmap5);
  title('CESM2.X d(stemp)/dt /yr');
caxis([-1 +1]*0.15)

figure(5); pcolor(save_lat64x72,Qlevs,squeeze(nanmean(save64x72_Q,[2 4]))'); shading interp; colorbar; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); title('CESM2.X <Q>')
%figure(5); pcolor(save_lat64x72,Qlevs,log10(squeeze(nanmean(save64x72_Q,[2 4]))')); shading interp; colorbar; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log')

comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compute_cesm3_trends.m';
saver = ['save /asl/s1/sergio/CESM3/cesm3_64x72_rates_stats_' savestr_version '_desc.mat thestats64x72 thestats64x72_other Tlevs Qlevs rlat rlon save_lon64x72 save_lat64x72 zonk comment'];
saver = ['save /asl/s1/sergio/CESM3/cesm3_64x72_rates_stats_' savestr_version '_desc.mat thestats64x72                     Tlevs Qlevs rlat rlon save_lon64x72 save_lat64x72 zonk comment'];
eval(saver)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_the_plots_64x72

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iAnom = input('do the anomaly (-1/+1) : ');
iAnom = -1;   %% code sucks
if iAnom > 0

  quick_junk_anom64x72
  anom_thestats64x72 = do_profilerate_fit_anom3(save64x72_Q(:,:,zonk),save64x72_O3(:,:,zonk),save64x72_T(:,:,zonk),save64x72_stemp(:,zonk),days(zonk),rlat);  
  anom_thestats64x72.days = days;
  if iL3orCLIMCAPS == +1
    saver = ['save /asl/s1/sergio/CESM3/fixedanomO3_airsL3_v7_64x72_rates_stats_' savestr_version '_all.mat thestats anom_thestats Tlevs Qlevs zonk'];  
  elseif iL3orCLIMCAPS == -1
    saver = ['save /asl/s1/sergio/CESM3/fixedanomO3_airsclimcaps_64x72_rates_stats_' savestr_version '_all.mat thestats anom_thestats Tlevs Qlevs zonk'];  
  end
  eval(saver)
  
  warning on

  addpath /home/sergio/MATLABCODE/COLORMAP
  figure(8); pcolor(squeeze(anom_thestats64x72.anomaly_temp(:,20,:)));
  title('equator T anom'); colorbar; colormap(usa2); caxis([-10 +10]); shading interp; %set(gca,'ydir','reverse')

  figure(9); pcolor(squeeze(anom_thestats64x72.anomaly_water(:,20,:)));
  title('equator WVfrac anom'); colorbar; colormap(usa2); caxis([-1 +1]); shading interp; %set(gca,'ydir','reverse')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
