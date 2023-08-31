%% see ../FIND_NWP_MODEL_TRENDS/driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m

%%% README >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%{

if you         have not saved off your ecRad runs, and      want to do fresh   regressions/save them .... do most of these steps, using [default] till (3) = do_feedbacks_wrt_globalSST.m, then DO NOT load the file you want to change
if you already have     saved off your ecRad runs, and just want to update the regressions/save them .... do all     these steps, using [default]


1) edit "driver_compute_feedbacks_from_smallFATfile.m"
      set eg a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 20
2) then run "driver_compute_feedbacks_from_smallFATfile.m"
3) when you get to "do_feedbacks_wrt_globalSST.m" just load in all the files eg
     /asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_20.mat
     /asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_20.mat
     /asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_20.mat
4) To recompute the regressions (polyfit/robutfit etc and any other new stuff) run 
       "show_olr_ecRad_feedback"  (which is a chip off the "compute_feedbacks_generic_ecRad.m" code)
     WARNING : if you do not have MERRA2/CLIMCAPSL3/AMIP6 it will get annoyed towards the end (when it is plotting) but else things should be fine
     This is eg for 05,10,15 years of trends
5) Save as needed using  "quick_save_olr_feedbacks_umbc_NWP_L3_XMIP6.m"!!!!

AND THEN WHEN HAPPY
AND THEN WHEN HAPPY
AND THEN WHEN HAPPY
clear all; for ii = 1 : 4; figure(ii); clf; end; iNumYears = 20; plot_show_olr_ecRad_feedback_globalsstfit
or
clear all; for ii = 1 : 5; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfit
%}
%%% README >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%%% README >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%{
if you do have [h,p] then can you run it through a generic ERA or ECM rtpmake??????

iNumYears = 20;
read_fileMean17years

%}
%%% README >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlib/maps
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/SHOWSTATS
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

a.topts.dataset = 10; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset10_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 05;  %% use CarbonTracker CO2 trends
a.topts.dataset = 11; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset11_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 10;  %% use CarbonTracker CO2 trends
a.topts.dataset = 12; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset12_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 15;  %% use CarbonTracker CO2 trends
a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 20;  %% use CarbonTracker CO2 trends ****

fprintf(1,'iNumYears = %2i .. reading in %s \n',iNumYears,strUMBC);
loader = ['load ' strUMBC];
eval(loader);

read_fileMean17years
h = hMean17years;
p = pMean17years;

% iSwap_ERA_2012_08_15 = +1;
iSwap_ERA_2012_08_15 = input('swap from 20 year AVG profile used for jacs .... to profiles on 2012/08/15 (-1 [default] / +1) : ');
if length(iSwap_ERA_2012_08_15) == 0
  iSwap_ERA_2012_08_15 = -1;
end
if iSwap_ERA_2012_08_15 > 0
  disp('swapping profile to ERA profile on 2012/08/15')
  %% 0130 = 1 hr 30 min = 01.50
  [h2,p2] = swap_mean_profile_TO_20XY_08_15(h,p,2012,01.50);
  [h3,p3] = swap_mean_profile_TO_20XY_Jan_to_Dec(h,p,2012,01.50);
  compare_profile_fields(p,p2)
  mmw0 = mmwater_rtp(h,p);
  mmw2 = mmwater_rtp(h2,p2);
  mmw3 = mmwater_rtp(h3,p3);

  %%%

  figure(1); scatter_coast(p.rlon,p.rlat,100,p.stemp-p2.stemp)
  figure(1); scatter_coast(p.rlon,p.rlat,100,mmw0-mmw2);

  semilogy(nanmean(p.ptemp'),nanmean(p.plevs'),nanmean(p2.ptemp'),nanmean(p.plevs'),nanmean(p3.ptemp'),nanmean(p.plevs')); set(gca,'ydir','reverse'); xlim([200 300]); ylim([1 1000]); plotaxis2;
  loglog(nanmean(p.gas_1'),nanmean(p.plevs'),nanmean(p2.gas_1'),nanmean(p.plevs'),nanmean(p3.gas_1'),nanmean(p.plevs')); set(gca,'ydir','reverse'); xlim([1e15 1e25]); ylim([1 1000]); plotaxis2;

  plot(nanmean(p.ptemp'-p2.ptemp'),nanmean(p.plevs'),'b',   nanstd(p.ptemp'-p2.ptemp'),nanmean(p.plevs'),'c',   nanmean(p.ptemp'-p3.ptemp'),nanmean(p.plevs'),'r', nanstd(p.ptemp'-p3.ptemp'),nanmean(p.plevs'),'m--'); 
    set(gca,'ydir','reverse'); xlim([-5 +10]); plotaxis2;
  plot(nanmean(p.gas_1'./p2.gas_1'),nanmean(p.plevs'),'b',1+nanstd(p.gas_1'./p2.gas_1'),nanmean(p.plevs'),'c--',nanmean(p.gas_1'./p3.gas_1'),nanmean(p.plevs'),'r',1+nanstd(p.gas_1'./p3.gas_1'),nanmean(p.plevs'),'m--'); 
    set(gca,'ydir','reverse'); xlim([-5 +10]); plotaxis2;

  %%%

  figure(1); scatter_coast(p.rlon,p.rlat,100,p3.stemp-p2.stemp)
  figure(1); scatter_coast(p.rlon,p.rlat,100,mmw3-mmw2);

  semilogy(nanmean(p3.ptemp'),nanmean(p3.plevs'),nanmean(p2.ptemp'),nanmean(p3.plevs')); set(gca,'ydir','reverse'); xlim([200 300]); ylim([1 1000]); plotaxis2;
  loglog(nanmean(p3.gas_1'),nanmean(p3.plevs'),nanmean(p2.gas_1'),nanmean(p3.plevs')); set(gca,'ydir','reverse'); xlim([1e15 1e25]); ylim([1 1000]); plotaxis2;

  plot(nanmean(p3.ptemp' - p2.ptemp'),nanmean(p3.plevs'),nanstd(p3.ptemp' - p2.ptemp'),nanmean(p3.plevs')); set(gca,'ydir','reverse'); xlim([-5 +10]); plotaxis2;
  plot(nanmean(p3.gas_1' ./ p2.gas_1'),nanmean(p3.plevs'),1 + nanstd(p3.gas_1' ./ p2.gas_1'),nanmean(p3.plevs')); set(gca,'ydir','reverse'); xlim([-5 +10]); plotaxis2;

  %%%

  [mean(p.stemp-p2.stemp) std(p.stemp-p2.stemp) mean(p.stemp-p3.stemp) std(p.stemp-p3.stemp)    mean(p3.stemp-p2.stemp) std(p3.stemp-p2.stemp)]
  [mean(mmw0-mmw2) std(mmw0-mmw2)               mean(mmw0-mmw3) std(mmw0-mmw3)                  mean(mmw3-mmw2) std(mmw3-mmw2)]

  h = h2;  p = p2;
  h = h3;  p = p3;
end

get_nan_bottom_layer
p = make_rtp_plays(p);
p.plays(p.plays <= eps) = NaN;
plays = nanmean(p.plays,2);
plays = plays(1:100);

if ~exist('deltaO3')
  disp('computing deltaO3')
  if ~exist('maskLF')
    disp('  warning : setting maskLF = ones(1,4608)')
    maskLF = ones(1,4608);
    maskLFmatr = reshape(maskLF,72,64)';
  end
  if ~exist('iNumLay')
    disp('  warning : setting iNumLay = 49')
    iNumLay = 49; 
  end
  if ~exist('xb')
    disp('  warning : setting xb = 0')
    xb = zeros(iNumLay,4608);
  end
  if ~exist('pavg')
    disp('  warning : setting pavg = pjunk20')
    pavg = pjunk20;
  end
  if ~exist('resultsunc')
    disp('  warning : setting resultsunc = results/100')
    resultsunc = results/100;
  end
  interp_resultsT_WV_O3_to_p
  get_deltaO3
end

if ~exist('xRH0')
  [xRH0,xRH1km0,xcolwater0] = layeramt2RH(h,p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = h.vchan;
mask = 1 : 4608;
mask = maskLF;
simple_get_the_model_trends_do_feedbacks
