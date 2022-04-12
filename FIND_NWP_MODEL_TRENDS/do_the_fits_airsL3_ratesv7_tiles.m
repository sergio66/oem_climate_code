warning off

zonkA = find(yy == StartY & mm == StartYM);
zonkB = find(yy == StopY  & mm == StopYM);
%zonk = 1:zonk;
zonk = zonkA : zonkB;
fprintf(1,'going from %4i/%2i to %4i/%2i .. found %3i time points \n',StartY,StartYM,StopY,StopYM,length(zonk))

%{
warning off
for ii = 1 : 72
  fprintf(1,'lonbin %2i of 72 \n',ii);
  junk = do_profilerate_fit_X(squeeze(save64x72_RH(:,ii,:,zonk)),days(zonk),rlat);
  thestats64x72.RHrate(ii,:,:) = junk.RHrate;
  thestats64x72.RHratestd(ii,:,:) = junk.RHratestd;
  thestats64x72.RHlag(ii,:,:) = junk.RHlag;
  thestats64x72.RHratestd_lag(ii,:,:) = junk.RHratestd_lag;
end
warning on
%}

warning off
for ii = 1 : 72
  fprintf(1,'lonbin %2i of 72 \n',ii);
  %xthestats       = do_profilerate_fit(squeeze(save64x72_Q(:,ii,:,zonk)),squeeze(save64x72_T(:,ii,:,zonk)),squeeze(save64x72_stemp(:,ii,zonk)),...
  %                                 days(zonk),rlat);
  if iL3orCLIMCAPS == +1
    xthestats       = do_profilerate_fit_WV_T_O3_RH_ST(squeeze(save64x72_Q(:,ii,:,zonk)),squeeze(save64x72_T(:,ii,:,zonk)),squeeze(save64x72_RH(:,ii,:,zonk)),squeeze(save64x72_O3(:,ii,:,zonk)),...
                                     squeeze(save64x72_stemp(:,ii,zonk)),squeeze(save64x72_RHSurf(:,ii,zonk)),squeeze(save64x72_TWetSurf(:,ii,zonk)),...
                                     days(zonk),rlat);

    xthestats_cld1  = do_profilerate_fit(squeeze(save64x72_cld_frac(:,ii,:,zonk)),squeeze(save64x72_cld_pres(:,ii,:,zonk)),squeeze(save64x72_stemp(:,ii,zonk)),...
                                   days(zonk),rlat);

  else
    %% wierdly, CLIMCAPS L3 does not have O3 so call RH twice, second time as a dummy
    xthestats       = do_profilerate_fit_WV_T_O3_RH_ST(squeeze(save64x72_Q(:,ii,:,zonk)),squeeze(save64x72_T(:,ii,:,zonk)),squeeze(save64x72_RH(:,ii,:,zonk)),squeeze(save64x72_RH(:,ii,:,zonk)),...
                                     squeeze(save64x72_stemp(:,ii,zonk)),squeeze(save64x72_RHSurf(:,ii,zonk)),squeeze(save64x72_TWetSurf(:,ii,zonk)),...
                                     days(zonk),rlat);
  end

  thestats64x72.lats = xthestats.lats;

  thestats64x72.waterrate(ii,:,:) = xthestats.waterrate;
  thestats64x72.waterratestd(ii,:,:) = xthestats.waterratestd;
  thestats64x72.waterlag(ii,:,:) = xthestats.waterlag;
  thestats64x72.waterratestd_lag(ii,:,:) = xthestats.waterratestd_lag;

  if iL3orCLIMCAPS == +1
    thestats64x72.ozonerate(ii,:,:) = xthestats.ozonerate;
    thestats64x72.ozoneratestd(ii,:,:) = xthestats.ozoneratestd;
    thestats64x72.ozonelag(ii,:,:) = xthestats.ozonelag;
    thestats64x72.ozoneratestd_lag(ii,:,:) = xthestats.ozoneratestd_lag;

    thestats64x72.cld_frac_rate(ii,:,:) = xthestats_cld1.waterrate;  thestats64x72.cld_frac_ratestd(ii,:,:) = xthestats_cld1.waterrate;
    thestats64x72.cld_pres_rate(ii,:,:) = xthestats_cld1.ptemprate;  thestats64x72.cld_pres_ratestd(ii,:,:) = xthestats_cld1.ptemprate;  %% fixed Sept 31, 2018 (or Oct -1, 2018)

  else
    thestats64x72.ozonerate(ii,:,:) = xthestats.ozonerate * 0;
    thestats64x72.ozoneratestd(ii,:,:) = xthestats.ozoneratestd * 0;
    thestats64x72.ozonelag(ii,:,:) = xthestats.ozonelag * 0;
    thestats64x72.ozoneratestd_lag(ii,:,:) = xthestats.ozoneratestd_lag * 0;
  end

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

  thestats64x72.TWetSurfrate(ii,:) = xthestats.TWetSurfrate;
  thestats64x72.TWetSurfratestd(ii,:) = xthestats.TWetSurfratestd;
  thestats64x72.TWetSurflag(ii,:) = xthestats.TWetSurflag;
  thestats64x72.TWetSurfratestd_lag(ii,:) = xthestats.TWetSurfratestd_lag;

  thestats64x72.RHSurfrate(ii,:) = xthestats.RHSurfrate;
  thestats64x72.RHSurfratestd(ii,:) = xthestats.RHSurfratestd;
  thestats64x72.RHSurflag(ii,:) = xthestats.RHSurflag;
  thestats64x72.RHSurfratestd_lag(ii,:) = xthestats.RHSurfratestd_lag;

end

if iL3orCLIMCAPS > 0
  %% remember we don't have cloud stuff here so this is sorta a waste
  for ii = 1 : 72
    %xthestats64x72_other = do_profilerate_fit_other_fraction(squeeze(save64x72_O3(:,ii,:,zonk)),squeeze(save64x72_olr(:,ii,zonk)),squeeze(save64x72_clrolr(:,ii,zonk)),days(zonk),rlat);
    xthestats64x72_other = do_profilerate_fit_other_fraction(squeeze(save64x72_CH4(:,ii,:,zonk)),squeeze(save64x72_olr(:,ii,zonk)),squeeze(save64x72_clrolr(:,ii,zonk)),days(zonk),rlat);
    thestats64x72.ch4rate(ii,:,:)        = xthestats.ozonerate;
    thestats64x72.ch4ratestd(ii,:,:)     = xthestats.ozoneratestd;
    thestats64x72.ch4lag(ii,:,:)         = xthestats.ozonelag;
    thestats64x72.ch4ratestd_lag(ii,:,:) = xthestats.ozoneratestd_lag;
    xthestats64x72_other = rmfield(xthestats64x72_other,'ozonerate');
    xthestats64x72_other = rmfield(xthestats64x72_other,'ozoneratestd');
    xthestats64x72_other = rmfield(xthestats64x72_other,'ozonelag');
    xthestats64x72_other = rmfield(xthestats64x72_other,'ozoneratestd_lag');
  
    xthestats64x72_other = do_profilerate_fit_other_fraction(squeeze(save64x72_CO(:,ii,:,zonk)),squeeze(save64x72_olr(:,ii,zonk)),squeeze(save64x72_clrolr(:,ii,zonk)),days(zonk),rlat);
    junk = save64x72_ice_od(:,zonk); junk = junk ./ (nanmean(junk')' * ones(1,length(zonk))+eps);
    thestats_cld2 = do_profilerate_fit_other_fraction(squeeze(save64x72_cld_frac(:,ii,:,zonk)),squeeze(save64x72_iceT(:,ii,zonk)),junk,days(zonk),rlat);
    thestats64x72.corate(ii,:,:)        = xthestats.ozonerate;
    thestats64x72.coratestd(ii,:,:)     = xthestats.ozoneratestd;
    thestats64x72.colag(ii,:,:)         = xthestats.ozonelag;
    thestats64x72.coratestd_lag(ii,:,:) = xthestats.ozoneratestd_lag;
    xthestats64x72_other = rmfield(xthestats64x72_other,'ozonerate');
    xthestats64x72_other = rmfield(xthestats64x72_other,'ozoneratestd');
    xthestats64x72_other = rmfield(xthestats64x72_other,'ozonelag');
    xthestats64x72_other = rmfield(xthestats64x72_other,'ozoneratestd_lag');
  
    thestats64x72_other.lats = xthestats64x72_other.lats;
    thestats64x72_other.olrrate(ii,:) = xthestats64x72_other.olrrate;
    thestats64x72_other.olrratestd(ii,:) = xthestats64x72_other.olrratestd;
    thestats64x72_other.clrolrrate(ii,:) = xthestats64x72_other.clrolrrate;
    thestats64x72_other.clrolrratestd(ii,:) = xthestats64x72_other.clrolrratestd;
   
    thestats64x72_other.iceT_rate(ii,:) = thestats_cld2.olrrate;       thestats64x72_other.iceT_ratestd(ii,:) = thestats_cld2.olrratestd;
    thestats64x72_other.ice_od_rate(ii,:) = thestats_cld2.clrolrrate;  thestats64x72_other.ice_od_ratestd(ii,:) = thestats_cld2.clrolrratestd;
  
    junk = squeeze(save64x72_icesze(:,ii,zonk));    junk = junk ./ (nanmean(junk')' * ones(1,length(zonk)));
    punk = squeeze(save64x72_liq_water(:,ii,zonk)); junk = junk ./ (nanmean(junk')' * ones(1,length(zonk)));    %% all NANS
    %thestats_cld3 = do_profilerate_fit_other_fraction(save64x72_cld_frac(:,:,zonk),junk,punk,days(zonk),rlat);
    %  thestats64x72_other.liq_water_rate = thestats_cld3.clrolrrate;  thestats64x72_other.liq_water_ratestd = thestats_cld3.clrolrratestd;
    thestats_cld3 = do_profilerate_fit_other_fraction(save64x72_cld_frac(:,:,zonk),junk,junk,days(zonk),rlat);
    thestats64x72_other.icesze_rate(ii,:)    = thestats_cld3.olrrate;      thestats64x72_other.icesze_ratestd(ii,:) = thestats_cld3.olrratestd;
    thestats64x72_other.liq_water_rate(ii,:) = nan*thestats_cld3.olrrate;  thestats64x72_other.liq_water_ratestd(ii,:) = nan*thestats_cld3.clrolrratestd;
  end
  warning on
else
  thestats64x72_other = 'nothing';
end

addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load llsmap5
figure(1); pcolor(save_lat64x72,Tlevs/100,squeeze(nanmean(thestats64x72.ptemprate,1))'); shading interp; colorbar; caxis([-0.15 +0.15]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log')
  title('AIRS L3 dT/dt K/yr');
figure(2); pcolor(save_lat64x72,Qlevs/100,squeeze(nanmean(thestats64x72.waterrate,1))'); shading interp; colorbar; caxis([-1 +1]*0.01); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log')
  title('AIRS L3 d(WV)/<WV>/dt /yr');
figure(3); pcolor(save_lat64x72,Qlevs/100,squeeze(nanmean(thestats64x72.RHrate*100,1))'); shading interp; colorbar; caxis([-1 +1]*0.15); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log')
  title('AIRS L3 d(RH)/dt /yr');

figure(1); ylim([10 1000]);
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
figure(4); aslmap(4,rlat65,rlon73,smoothn(thestats64x72.stemprate',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('d/dt AIRS L3 K/yr'); colormap(llsmap5);
  title('AIRS L3 d(stemp)/dt /yr');
caxis([-1 +1]*0.15)

figure(4); pcolor(save_lat64x72,Qlevs,squeeze(nanmean(save64x72_Q,[2 4]))'); shading interp; colorbar; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log')
figure(4); pcolor(save_lat64x72,Qlevs,log10(squeeze(nanmean(save64x72_Q,[2 4]))')); shading interp; colorbar; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log')

comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compute_AIRSL3_trends.m';
if iL3orCLIMCAPS == +1
  if iDorA > 0
    %saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_' savestr_version '_desc.mat thestats Tlevs Qlevs zonk'];
    saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_' savestr_version '_desc.mat thestats64x72 thestats64x72_other Tlevs Qlevs rlat rlon save_lon64x72 save_lat64x72 zonk comment'];
  else
    %saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_' savestr_version '_asc.mat thestats Tlevs Qlevs zonk'];
    saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_' savestr_version '_asc.mat thestats64x72 thestats64x72_other Tlevs Qlevs rlat rlon save_lon64x72 save_lat64x72 zonk comment'];
  end
elseif iL3orCLIMCAPS == -1
  if iDorA > 0
    %saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_' savestr_version '_desc.mat thestats Tlevs Qlevs zonk'];
    saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_' savestr_version '_desc.mat thestats64x72 thestats64x72_other Tlevs Qlevs rlat rlon save_lon64x72 save_lat64x72 zonk comment'];
  else
    %saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_' savestr_version '_asc.mat thestats Tlevs Qlevs zonk'];
    saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_stats_' savestr_version '_asc.mat thestats64x72 thestats64x72_other Tlevs Qlevs rlat rlon save_lon64x72 save_lat64x72 zonk comment'];
  end
end
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
    saver = ['save /asl/s1/sergio/AIRS_L3/fixedanomO3_airsL3_v7_64x72_rates_stats_' savestr_version '_all.mat thestats anom_thestats Tlevs Qlevs zonk'];  
  elseif iL3orCLIMCAPS == -1
    saver = ['save /asl/s1/sergio/AIRS_L3/fixedanomO3_airsclimcaps_64x72_rates_stats_' savestr_version '_all.mat thestats anom_thestats Tlevs Qlevs zonk'];  
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
