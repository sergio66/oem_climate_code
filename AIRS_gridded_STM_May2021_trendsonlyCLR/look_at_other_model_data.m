%% now look at other model data

figure(28); ylim([100 1000]);
figure(29);
figure(30); ylim([100 1000]);
if dataset == 3
  figure(29); title('Zonal d/dt T  UMBC Extreme');
  figure(28); title('Zonal d/dt RH UMBC Extreme');      ylim([100 1000]);
  figure(30); title('Zonal d/dt WVfrac UMBC Extreme');  ylim([100 1000]);
end

redo_fig6_surfT_rate
redo_fig8_spectralrates_with_mask

addpath ../FIND_NWP_MODEL_TRENDS
driver_get_the_model_trends

%figure(8); plot_profile_trends          %% ERA_ERA5_AIRSL3,   with mask
%figure(8); plot_profile_trends2         %% CMIP6_ERA5_AIRSL3, with mask, but uses mean as mean(profile) and unc as std(profile)
figure(8); plot_profile_trends3_cmip6    %% CMIP6_ERA5_AIRSL3, with mask, but uses mean as mean(profile) and unc as mean(unc_profile),,,, I could call plot_profile_trends3_amip6 but I have cleared the variable

iX = input('DO you want to do atmospheric amplification (-1/+1) [-1 : default] : ');
if length(iX) == 0
  iX = -1;
end
if iX > 0
  atmospheric_amplification_plots2
end

iX = input('Do you want to also see profile/spectral trends with uncertainty??? (-1/+1) [-1 default] : ');
if length(iX) == 0
  iX = -1;
end
if iX > 0
  figure(8); plot_profile_trends2_with_unc %% redo UMBC with unc
end

iX = input('Do regressions of retrievals vs models? (-1/+1) [-1 default] ');
if length(iX) == 0
  iX = -1;
end
if iX > 0
  do_regressions
end

iX = input('compute feedbacks? (-1/+1) [default -1] ? ');
if length(iX) == 0
  iX = 0;
end
if iX > 0
  do_feedbacks
end

%plot_HadSurf_trends_36x72
iX = input('compare HadSurf to UMBC by converting from 36x72 to 4608(=64x72) (-1 no, default, +1 yes) ? ');
if length(iX) == 0
  iX = -1;
end
if iX > 0
  plot_HadSurf_trends_convert_to_64x72
end

iX = input('Plot wet bulb trends and PNAS2018 paper? (-/+1) [-1 default] : ');
if length(iX) == 0
  iX = -1;
end
if iX > 0
  %% wetbulb bad life
  addpath ../FIND_NWP_MODEL_TRENDS/
  find_wet_bulb_trends
  pnas2018_byrne_gorman
end

iX = input('Plot UMBC/AIRSL3/CMIP6/ERA5 zonal comparisons again (tiled plots) (-1/+1) [+1 default] ? ');
if length(iX) == 0
  iX = 1;
end
if iX > 0
  re_plot_zonal_trends_umbc_airsL3_models
end

iX = input('do you want to do spectral uncertainties? Takes a LOOOOONG time! (-1 no default, +1 yes) ?');
if iX > 0
  do_spectral_closure
end

pcld = p;
pcld = convert_rtp_to_cloudOD(h,pcld);
icecldOD = reshape(pcld.iceOD,72,64); icecldTOP = reshape(pcld.icetop,72,64); 
watercldOD = reshape(pcld.waterOD,72,64); watercldTOP = reshape(pcld.watertop,72,64); 
hold on; scatter(rlat,nanmean(icecldTOP,1),50,nanmean(icecldOD,1),'filled'); colorbar; set(gca,'ydir','reverse'); ylim([10 1000])
%hold on; scatter(rlat,nanmean(watercldTOP,1),50,nanmean(watercldOD,1),'filled'); colorbar; set(gca,'ydir','reverse'); ylim([10 1000])
hold on; scatter(rlat,nanmean(watercldTOP,1),50,nanmean(watercldOD,1)/10,'filled'); colorbar; set(gca,'ydir','reverse'); ylim([10 1000])
hold off
title('ice cldOD, water cldOD/10')
xlim([-90 +90]); grid; plotaxis2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
tons_of_aslprint_plots
%}
