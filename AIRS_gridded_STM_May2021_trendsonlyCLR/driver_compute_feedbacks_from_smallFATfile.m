%% see ../FIND_NWP_MODEL_TRENDS/driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m

%%% README >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%{
if you already have saved off your ecRad runs, and just want to update the regressions/save them .... then 
1) set eg a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 20
2) just run "driver_compute_feedbacks_from_smallFATfile.m"
3) when you get to "do_feedbacks_wrt_globalSST.m" just load in all the files eg
     /asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_20.mat
     /asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_20.mat
     /asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_20.mat
4) Now run "show_olr_ecRad_feedback"  (which is a chip off the "compute_feedbacks_generic_ecRad.m" code)
5) Save as needed using  "quick_save_olr_feedbacks_umbc_NWP_L3_XMIP6.m"!!!!
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

a.topts.dataset = 12; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset12_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 15;  %% use CarbonTracker CO2 trends
a.topts.dataset = 11; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset11_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 10;  %% use CarbonTracker CO2 trends
a.topts.dataset = 10; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset10_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 05;  %% use CarbonTracker CO2 trends
a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 20;  %% use CarbonTracker CO2 trends

fprintf(1,'iNumYears = %2i .. reading in %s \n',iNumYears,strUMBC);
loader = ['load ' strUMBC];
eval(loader);

read_fileMean17years
h = hMean17years;
p = pMean17years;
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
