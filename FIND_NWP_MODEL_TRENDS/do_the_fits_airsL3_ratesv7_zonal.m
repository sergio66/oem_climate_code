warning off

zonkA = find(yy == StartY & mm == StartYM);
zonkB = find(yy == StopY  & mm == StopYM);
%zonk = 1:zonk;
zonk = zonkA : zonkB;
fprintf(1,'going from %4i/%2i to %4i/%2i .. found %3i time points \n',StartY,StartYM,StopY,StopYM,length(zonk))

thestatsOld    = do_profilerate_fit(save_Q(:,:,zonk),save_T(:,:,zonk),save_stemp(:,zonk),...
                                   days(zonk),latbins);
thestats       = do_profilerate_fit_WV_T_O3_RH_ST(save_Q(:,:,zonk),save_T(:,:,zonk),save_RH(:,:,zonk),...
                                   save_O3(:,:,zonk),save_stemp(:,zonk),save_RHSurf(:,zonk),save_TWetSurf(:,zonk),...
                                   days(zonk),latbins);
thestats_cld1  = do_profilerate_fit(save_cld_frac(:,:,zonk),save_cld_pres(:,:,zonk),save_stemp(:,zonk),...
                                   days(zonk),latbins);
  thestats.cld_frac_rate = thestats_cld1.waterrate;  thestats.cld_frac_ratestd = thestats_cld1.waterrate;
  thestats.cld_pres_rate = thestats_cld1.ptemprate;  thestats.cld_pres_ratestd = thestats_cld1.ptemprate;  %% fixed Sept 31, 2018 (or Oct -1, 2018)

figure(1);  caxis([-1 +1]*0.01)
figure(2);  caxis([-1 +1]*0.15)
figure(10); pcolor(thestats.lats,1:length(Airs_PQ),100*thestats.RHrate'); colormap(llsmap5); colorbar; shading interp; title('RH rate'); set(gca,'ydir','reverse');
  caxis([-1 +1]*0.15)

%{
%% upto Feb 2016
thestats_other = do_profilerate_fit_other_fraction(save_O3(:,:,zonk),save_olr(:,zonk),save_clrolr(:,zonk),days(zonk),latbins);
thestats_cld2 = do_profilerate_fit_other_fraction(save_cld_frac(:,:,zonk),save_iceT(:,zonk),save_ice_od(:,zonk),days(zonk),latbins);
  thestats_other.iceT_rate = thestats_cld2.olrrate;       thestats_other.iceT_ratestd = thestats_cld2.olrratestd;
  thestats_other.ice_od_rate = thestats_cld2.clrolrrate;  thestats_other.ice_od_ratestd = thestats_cld2.clrolrratestd;
thestats_cld3 = do_profilerate_fit_other_fraction(save_cld_frac(:,:,zonk),save_icesze(:,zonk),save_liq_water(:,zonk),days(zonk),latbins);
  thestats_other.icesze_rate = thestats_cld3.olrrate;       thestats_other.icesze_ratestd = thestats_cld3.olrratestd;
  thestats_other.liq_water_rate = thestats_cld3.clrolrrate;  thestats_other.liq_water_ratestd = thestats_cld3.clrolrratestd;
%}

if iL3orCLIMCAPS > 0
  %% after Mar 2016
  %thestats_other = do_profilerate_fit_other_fraction(save_O3(:,:,zonk),save_olr(:,zonk),save_clrolr(:,zonk),days(zonk),latbins);
  thestats_other = do_profilerate_fit_other_fraction(save_CH4(:,:,zonk),save_olr(:,zonk),save_clrolr(:,zonk),days(zonk),latbins);
    thestats.ch4rate        = thestats_other.ozonerate;
    thestats.ch4ratestd     = thestats_other.ozoneratestd;
    thestats.ch4lag         = thestats_other.ozonelag;
    thestats.ch4ratestd_lag = thestats_other.ozoneratestd_lag;
    thestats_other = rmfield(thestats_other,'ozonerate');
    thestats_other = rmfield(thestats_other,'ozoneratestd');
    thestats_other = rmfield(thestats_other,'ozonelag');
    thestats_other = rmfield(thestats_other,'ozoneratestd_lag');
  thestats_other = do_profilerate_fit_other_fraction(save_CO(:,:,zonk),save_olr(:,zonk),save_clrolr(:,zonk),days(zonk),latbins);
    thestats.corate        = thestats_other.ozonerate;
    thestats.coratestd     = thestats_other.ozoneratestd;
    thestats.colag         = thestats_other.ozonelag;
    thestats.coratestd_lag = thestats_other.ozoneratestd_lag;
    thestats_other = rmfield(thestats_other,'ozonerate');
    thestats_other = rmfield(thestats_other,'ozoneratestd');
    thestats_other = rmfield(thestats_other,'ozonelag');
    thestats_other = rmfield(thestats_other,'ozoneratestd_lag');
  
  junk = save_ice_od(:,zonk); junk = junk ./ (nanmean(junk')' * ones(1,length(zonk)));
  thestats_cld2 = do_profilerate_fit_other_fraction(save_cld_frac(:,:,zonk),save_iceT(:,zonk),junk,days(zonk),latbins);
    thestats_other.iceT_rate = thestats_cld2.olrrate;       thestats_other.iceT_ratestd = thestats_cld2.olrratestd;
    thestats_other.ice_od_rate = thestats_cld2.clrolrrate;  thestats_other.ice_od_ratestd = thestats_cld2.clrolrratestd;
  
  junk = save_icesze(:,zonk); junk = junk ./ (nanmean(junk')' * ones(1,length(zonk)));
  punk = save_liq_water(:,zonk); junk = junk ./ (nanmean(junk')' * ones(1,length(zonk)));    %% all NANS
  %thestats_cld3 = do_profilerate_fit_other_fraction(save_cld_frac(:,:,zonk),junk,punk,days(zonk),latbins);
  %  thestats_other.liq_water_rate = thestats_cld3.clrolrrate;  thestats_other.liq_water_ratestd = thestats_cld3.clrolrratestd;
  thestats_cld3 = do_profilerate_fit_other_fraction(save_cld_frac(:,:,zonk),junk,junk,days(zonk),latbins);
    thestats_other.icesze_rate = thestats_cld3.olrrate;         thestats_other.icesze_ratestd = thestats_cld3.olrratestd;
    thestats_other.liq_water_rate = nan*thestats_cld3.olrrate;  thestats_other.liq_water_ratestd = nan*thestats_cld3.clrolrratestd;

else
  thestats_other = 'this is CLIMCAPS so no cloud stuff';
end

comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compute_AIRSL3_trends.m';

if iL3orCLIMCAPS == +1
  if iDorA > 0
    %saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_' savestr_version '_desc.mat thestats Tlevs Qlevs zonk'];
    saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_' savestr_version '_desc.mat thestats thestats_other Tlevs Qlevs latbins save_lat zonk comment'];
  else
    %saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_' savestr_version '_asc.mat thestats Tlevs Qlevs zonk'];
    saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_' savestr_version '_asc.mat thestats thestats_other Tlevs Qlevs latbins save_lat zonk comment'];
  end
elseif iL3orCLIMCAPS == -1
  if iDorA > 0
    %saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_zonal_rates_stats_' savestr_version '_desc.mat thestats Tlevs Qlevs zonk'];
    saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_zonal_rates_stats_' savestr_version '_desc.mat thestats thestats_other Tlevs Qlevs latbins save_lat zonk comment'];
  else
    %saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_zonal_rates_stats_' savestr_version '_asc.mat thestats Tlevs Qlevs zonk'];
    saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_zonal_rates_stats_' savestr_version '_asc.mat thestats thestats_other Tlevs Qlevs latbins save_lat zonk comment'];
  end
end
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_the_plots_zonal

figure(1);  pcolor(thestats.lats,1:length(Airs_PQ),thestats.waterrate'); colormap(llsmap5); colorbar; shading interp; title('RH rate'); set(gca,'ydir','reverse'); caxis([-1 +1]*0.01)
figure(3);  pcolor(thestats.lats,1:length(Airs_PT),thestats.ptemprate'); colormap(llsmap5); colorbar; shading interp; title('RH rate'); set(gca,'ydir','reverse'); caxis([-1 +1]*0.15)
figure(10); pcolor(thestats.lats,1:length(Airs_PQ),100*thestats.RHrate'); colormap(llsmap5); colorbar; shading interp; title('RH rate'); set(gca,'ydir','reverse');  caxis([-1 +1]*1)

if iL3orCLIMCAPS == -1
  figure(1);  pcolor(thestats.lats,Airs_PQ/100,thestats.waterrate'); colormap(llsmap5); colorbar; shading interp; title('WV frac rate'); set(gca,'ydir','reverse'); caxis([-1 +1]*0.015); set(gca,'yscale','log'); ylim([100 1000])
  figure(3);  pcolor(thestats.lats,Airs_PT/100,thestats.ptemprate'); colormap(llsmap5); colorbar; shading interp; title('T rate'); set(gca,'ydir','reverse'); caxis([-1 +1]*0.15); set(gca,'yscale','log'); ylim([10 1000])
  figure(10); pcolor(thestats.lats,Airs_PQ/100,100*thestats.RHrate'); colormap(llsmap5); colorbar; shading interp; title('RH rate'); set(gca,'ydir','reverse');  caxis([-1 +1]*0.5); set(gca,'yscale','log'); ylim([100 1000])
else
  figure(1);  pcolor(thestats.lats,Airs_PQ,thestats.waterrate'); colormap(llsmap5); colorbar; shading interp; title('WV frac rate'); set(gca,'ydir','reverse'); caxis([-1 +1]*0.015); set(gca,'yscale','log'); ylim([100 1000])
  figure(3);  pcolor(thestats.lats,Airs_PT,thestats.ptemprate'); colormap(llsmap5); colorbar; shading interp; title('T rate'); set(gca,'ydir','reverse'); caxis([-1 +1]*0.15); set(gca,'yscale','log'); ylim([10 1000])
  figure(10); pcolor(thestats.lats,Airs_PQ,thestats.RHrate'); colormap(llsmap5); colorbar; shading interp; title('RH rate'); set(gca,'ydir','reverse');  caxis([-1 +1]*0.5); set(gca,'yscale','log'); ylim([100 1000])
end

if exist('quickSTrate')
  figure(2); pcolor(quickSTrate); shading interp; colorbar; caxis([-1 +1]*0.15); title('Surf Temp Rate'); colormap(usa2);
  addpath /umbc/xfs2/strow/asl/matlib/maps/aslmap.m
  addpath /home/sergio/MATLABCODE/COLORMAP/LLSMAPS
  load llsmap5
  figure(2); aslmap(2,-90:1:+90,-180:1:+180,quickSTrate,[-90 +90],[-180 +180]);  colormap(llsmap5); caxis([-0.15 +0.15]);
  title('stemp rates')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iAnom = input('do the anomaly (-1/+1) : ');
if iAnom > 0
  %anom_thestats = do_profilerate_fit_anom(save_Q(:,:,zonk),save_T(:,:,zonk),save_stemp(:,zonk),days(zonk),latbins);
  %anom_thestats.days = days;
  %saver = ['save /asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_stats_' savestr_version '_all.mat thestats anom_thestats Tlevs Qlevs zonk'];
  eval(saver)

  %% vers1
  %{
  anom_thestats = do_profilerate_fit_anom2(save_Q(:,:,zonk),save_T(:,:,zonk),save_stemp(:,zonk),days(zonk),latbins);  
  anom_thestats.days = days;
  %saver = ['save /asl/s1/sergio/AIRS_L3/anom_airsL3_v7_rates_stats_' savestr_version '_all.mat thestats anom_thestats Tlevs Qlevs zonk'];  
  saver = ['save /asl/s1/sergio/AIRS_L3/fixedanom_airsL3_v7_rates_stats_' savestr_version '_all.mat thestats anom_thestats Tlevs Qlevs zonk'];
  eval(saver)
  %}

  %% vers2
  %{
  anom_thestatsO3 = do_profilerate_fit_anom2(save_Q(:,:,zonk),save_O3(:,:,zonk),save_stemp(:,zonk),days(zonk),latbins);
  anom_thestats = do_profilerate_fit_anom2(save_Q(:,:,zonk),save_T(:,:,zonk),save_stemp(:,zonk),days(zonk),latbins);  
  anom_thestats.days = days;
  anom_thestats.anomaly_O3 = anom_thestatsO3.anomaly_temp;
  anom_thestats.timeseries_withouttrend_O3 = anom_thestatsO3.timeseries_withouttrend_temp;
  anom_thestats.O3rate = anom_thestatsO3.ptemprate;
  anom_thestats.O3ratestd = anom_thestatsO3.ptempratestd;
  saver = ['save /asl/s1/sergio/AIRS_L3/fixedanomO3_airsL3_v7_rates_stats_' savestr_version '_all.mat thestats anom_thestats Tlevs Qlevs zonk'];  
  eval(saver)
  %}

  %% vers3
  quick_junk_anom
  anom_thestats = do_profilerate_fit_anom3(save_Q(:,:,zonk),save_O3(:,:,zonk),save_T(:,:,zonk),save_stemp(:,zonk),days(zonk),latbins);  
  anom_thestats.days = days;
  if iL3orCLIMCAPS == +1
    saver = ['save /asl/s1/sergio/AIRS_L3/fixedanomO3_airsL3_v7_zonal_rates_stats_' savestr_version '_all.mat thestats anom_thestats Tlevs Qlevs zonk'];  
  elseif iL3orCLIMCAPS == -1
    saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/fixedanomO3_airsclimcaps_zonal_rates_stats_' savestr_version '_all.mat thestats anom_thestats Tlevs Qlevs zonk'];  
  end
  eval(saver)
  
  warning on

  addpath /home/sergio/MATLABCODE/COLORMAP
  figure(8); pcolor(squeeze(anom_thestats.anomaly_temp(:,20,:)));
  title('equator T anom'); colorbar; colormap(usa2); caxis([-10 +10]); shading interp; %set(gca,'ydir','reverse')

  figure(9); pcolor(squeeze(anom_thestats.anomaly_water(:,20,:)));
  title('equator WVfrac anom'); colorbar; colormap(usa2); caxis([-1 +1]); shading interp; %set(gca,'ydir','reverse')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
