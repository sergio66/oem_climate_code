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
figure(8); plot_ERA_ERA5_AIRSL3_CMIP6_trends
%figure(8); plot_profile_trends          %% ERA_ERA5_AIRSL3,   with mas
figure(8); plot_profile_trends2          %% CMIP6_ERA5_AIRSL3, with mask
figure(8); plot_profile_trends2_with_unc %% redo UMBC with unc

do_regressions

iFeedback = input('compute feedbacks? (-1 no, default, +1 yes) ? ');
if iFeedback > 0
  do_feedbacks
end

%plot_HadSurf_trends_36x72
iConvert = input('compare HadSurf to UMBC by converting from 36x72 to 4608(=64x72) (-1 no, default, +1 yes) ? ');
if iConvert > 0
  plot_HadSurf_trends_convert_to_64x72
end

iPlotWB = input('Plot wet bulb trends and PNAS2018 paper? (-1 no default, +1 yes) : ');
if iPlotWB > 0
  %% wetbulb bad life
  addpath ../FIND_NWP_MODEL_TRENDS/
  find_wet_bulb_trends
  pnas2018_byrne_gorman
end

iReplot = input('Plot UMBC/AIRSL3/CMIP6/ERA5 zonal comparisons again (tiled plots) (-1 no default, +1 yes) ? ');
if iReplot > 0
  re_plot_zonal_trends_umbc_airsL3_models
end

iUnc = input('do you want to do spectral uncertainties? Takes a LOOOOONG time! (-1 no default, +1 yes) ?');
if iUnc > 0
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
