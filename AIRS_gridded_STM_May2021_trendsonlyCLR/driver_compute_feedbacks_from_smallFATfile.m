%% see ../FIND_NWP_MODEL_TRENDS/driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m

%%% README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>>
%%% README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>>
%%% README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>>
%{
if you want to re-do the 05/10/15/20 year UMBC time series 
(because you have updated eg compute_feedbacks_regress_olr_ecRad_calcs.m)
then simply call 
  driver_quick_redo_regressions_umbc_feedbacks_timeseries.m

if you want to re-do the 20 year UMBC and ERA5/AIRSL3/CMIP6 and MERRA2/CLIMCAPSL3/AMIP6 
(because you have updated eg compute_feedbacks_regress_olr_ecRad_calcs.m)
then simply call 
  driver_quick_redo_regressions_20year_umbc_models.m
%}
%%% README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>>
%%% README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>>
%%% README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>> README QUICK REDO >>>>>>>>>>>>>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% README SLOW rerun ecRad etc >>>>>>>>>>>>>>>> README SLOW rerun ecRad etc >>>>>>>>>>>>>>>> README SLOW rerun ecRad etc >>>>>>>>>>>>>>>> README SLOW rerun ecRad etc >>>>>>>>>>>>>>>>
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
4) Depending on whether this is "raw" deltaX ecRad calcs, or "unc" deltaX + deltaXunc ecRad calcs, update_deltaT_WV_O3_SKT_using_uncertainty.m will be called <<<<<<<<<<<<<
5) To recompute the regressions (polyfit/robutfit etc and any other new stuff) run 
       "show_olr_ecRad_feedback"  (which is a chip off the "compute_feedbacks_generic_ecRad.m" code)
     WARNING : if you do not have MERRA2/CLIMCAPSL3/AMIP6 it will get annoyed towards the end (when it is plotting) but else things should be fine
     This is eg for 05,10,15 years of trends
6) Save as needed using  "quick_save_olr_feedbacks_umbc_NWP_L3_XMIP6.m"!!!!

%% Feldl, N., and G. H. Roe (2013), Four perspectives on climate feedbacks, Geophys. Res. Lett., 40, 4007–4011, doi:10.1002/grl.50711.
AND THEN WHEN HAPPY
AND THEN WHEN HAPPY
AND THEN WHEN HAPPY

FOR GLOBAL DEFN
clear all; for ii = 1 : 6; figure(ii); clf; end; iNumYears = 20; plot_show_olr_ecRad_feedback_globalsstfit_smooth2
or
clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfit_smooth    these two are symbolic links
clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfitsmooth     these two are symbolic links

clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfit_smooth2    these two are symbolic links
clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfitsm2         these two are symbolic links

FOR LOCAL DEFN
clear all; for ii = 1 : 6; figure(ii); clf; end; iNumYears = 20; plot_show_olr_ecRad_feedback_polyfit              %%% local feedbacks
clear all; for ii = 1 : 6; figure(ii); clf; end; iNumYears = 20; plot_show_olr_ecRad_feedback_robustfit            %%% local feedbacks
clear all; for ii = 1 : 6; figure(ii); clf; end; iNumYears = 20; plot_show_olr_ecRad_feedback_globalsstfit         %%% global feedbacks

clear all; for ii = 1 : 6; figure(ii); clf; end; iNumYears = 20; plot_show_olr_ecRad_feedback_polyfit_smooth       %%% local feedbacks
clear all; for ii = 1 : 6; figure(ii); clf; end; iNumYears = 20; plot_show_olr_ecRad_feedback_robustfit_smooth     %%% local feedbacks
clear all; for ii = 1 : 6; figure(ii); clf; end; iNumYears = 20; plot_show_olr_ecRad_feedback_globalsstfit_smooth  %%% global feedbacks

clear all; for ii = 1 : 6; figure(ii); clf; end; iNumYears = 20; plot_show_olr_ecRad_feedback_polyfit_smooth2      %%% local feedbacks
clear all; for ii = 1 : 6; figure(ii); clf; end; iNumYears = 20; plot_show_olr_ecRad_feedback_robustfit_smooth2    %%% local feedbacks   >>> used for paper
clear all; for ii = 1 : 6; figure(ii); clf; end; iNumYears = 20; plot_show_olr_ecRad_feedback_globalsstfit_smooth2 %%% global feedbacks  >>> used for paper

or 

clear all; for ii = 1 : 8; figure(ii); clf; end; iNumYears = 20; plot_unc_olr_ecRad_feedback_polyfit_smooth2      %%% local feedbacks   
clear all; for ii = 1 : 8; figure(ii); clf; end; iNumYears = 20; plot_unc_olr_ecRad_feedback_robustfit_smooth2    %%% local feedbacks   >>> used for paper, to show unc
clear all; for ii = 1 : 8; figure(ii); clf; end; iNumYears = 20; plot_unc_olr_ecRad_feedback_globalsstfit_smooth2 %%% global feedbacks  >>> used for paper, to show unc
 
or

clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_timeseries_robustfit           %%% local feedbacks
clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfit        %%% global feedbacks

clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_timeseries_robustfit_smooth    %%% local feedbacks
clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfitsmooth  %%% global feedbacks

clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_timeseries_robustfit_smooth2   %%% local feedbacks   >>> used for paper
clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfitsm2     %%% global feedbacks  >>> used for paper

or 

clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_seasonal_robustfit_smooth2   %%% local feedbacks   >>> used for paper
clear all; for ii = 1 : 6; figure(ii); clf; end; plot_show_olr_ecRad_feedback_umbc_seasonal_globalsstfitsm2     %%% global feedbacks  >>> used for paper

%}
%%% README SLOW rerun ecRad etc >>>>>>>>>>>>>>>> README SLOW rerun ecRad etc >>>>>>>>>>>>>>>> README SLOW rerun ecRad etc >>>>>>>>>>>>>>>> README SLOW rerun ecRad etc >>>>>>>>>>>>>>>>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlib/maps
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/SHOWSTATS
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

iRaw_or_Unc = -1; %% use raw profiles, and then perturbations+unc, in the ecRad or SARTA calcs
iRaw_or_Unc = +1; %% use raw profiles, and then perturbations,     in the ecRad or SARTA calcs

iGet_ERA5_AIRSL3_AMIP = +1;  %% to do UMBC and also load in eg MERRA2/ERA5 model trends
iGet_ERA5_AIRSL3_AMIP = -1;  %% default

iAllorSeasonal = +1;
a.topts.dataset = 10; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset10_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';      iNumYears = 05;  %% use CarbonTracker CO2 trends
a.topts.dataset = 12; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset12_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';      iNumYears = 15;  %% use CarbonTracker CO2 trends
a.topts.dataset = 11; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset11_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2v2.mat';    iNumYears = 10;  %% use CarbonTracker CO2 trends
a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';      iNumYears = 20;  %% use CarbonTracker CO2 trends ****, topts.iAdjLowerAtmWVfrac=0.25

%% use CarbonTracker CO2 trends ****, topts.iAdjLowerAtmWVfrac=0.25
a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_SON.mat';  iNumYears = 20;  iAllorSeasonal = -4;
a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_JJA.mat';  iNumYears = 20;  iAllorSeasonal = -3;
a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_MAM.mat';  iNumYears = 20;  iAllorSeasonal = -2;
a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_DJF.mat';  iNumYears = 20;  iAllorSeasonal = -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if ~exist('pavg')
  disp('  warning : setting pavg = pjunk20')
  pavg = pjunk20;
end

if ~exist('resultsTunc') 
  resultsTunc  = 0 * resultsT;
  if iRaw_or_Unc == -1
    error('resultsTunc DNE and you want to use it!!!')
  end
end
if ~exist('resultsWVunc') 
  resultsWVunc = 0 * resultsT;
  if iRaw_or_Unc == -1
    error('resultsWVunc DNE and you want to use it!!!')
  end
end
if ~exist('resultsO3unc') 
  resultsO3unc = 0 * resultsT;
  if iRaw_or_Unc == -1
    error('resultsO3unc DNE and you want to use it!!!')
  end
end
if ~exist('resultsunc') 
  resultsunc   = 0 * results;
  if iRaw_or_Unc == -1
    disp('WARNING resultsunc DNE and you want to use it!!!, set to 0.2 * stemp')
    resultsunc      = 0.005 * ones(4608,6);
    resultsunc(:,6) = 0.2 * results(:,6);
  end
end

if ~exist('maskLF')
  disp('  warning : setting maskLF = ones(1,4608)')
  maskLF = ones(1,4608);
end
maskLFmatr = reshape(maskLF,72,64)';

if ~exist('fracO3')
  disp('computing fracO3')
  if ~exist('iNumLay')
    disp('  warning : setting iNumLay = 49')
    iNumLay = 49; 
  end
  if ~exist('xb')
    disp('  warning : setting xb = 0')
    xb = zeros(iNumLay,4608);
  end
  if ~exist('resultsunc')
    disp('  warning : setting resultsunc = results/100')
    resultsunc = results/100;
  end

  interp_resultsT_WV_O3_to_p
  get_deltaO3
end

if ~exist('deltaO3')
  interp_resultsT_WV_O3_to_p
  [ppmvLAY1,ppmvAVG1,ppmvMAX1,pavgLAY1,tavgLAY1,ppmv500_1,ppmv75_1,ppmvSURF_1] = layers2ppmv(h,p,1:length(p.stemp),1);
  [ppmvLAY3,ppmvAVG3,ppmvMAX3,pavgLAY3,tavgLAY3,ppmv500_3,ppmv75_3,ppmvSURF_3] = layers2ppmv(h,p,1:length(p.stemp),3);

  [nlayO3,~] = size(ppmvLAY3);

  [ppmvLAYpert3,ppmvAVGpert3,ppmvMAXpert3,pavgLAYpert3,tavgLAYpert3,ppmv500pert3,ppmv75pert3,ppmvSURFpert3] = layers2ppmv(h,pert,1:length(p.stemp),3);
  [ppmvLAYpert_unc3,ppmvAVGpert_unc3,ppmvMAXpert_unc3,pavgLAYpert_unc3,tavgLAYpert_unc3,ppmv500pert_unc3,ppmv75pert_unc3,ppmvSURFpert_unc3] = layers2ppmv(h,pert_unc,1:length(p.stemp),3);
  deltaO3unc = ppmvLAYpert_unc3 - ppmvLAY3;
  deltaO3unc = deltaO3unc .* (ones(nlayO3,1) * maskLF);

  deltaO3 = ppmvLAYpert3 - ppmvLAY3;
  deltaO3 = deltaO3 .* (ones(nlayO3,1) * maskLF);
end

if iRaw_or_Unc == -1
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
