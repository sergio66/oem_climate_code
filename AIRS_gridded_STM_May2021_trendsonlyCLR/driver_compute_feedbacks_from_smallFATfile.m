%% see ../FIND_NWP_MODEL_TRENDS/driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m

addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

a.topts.dataset = 9; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 20;         %% use CarbonTracker CO2 trends

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

f = h.vchan;
mask = 1 : 4608;
mask = maskLF;
simple_get_the_model_trends_do_feedbacks
