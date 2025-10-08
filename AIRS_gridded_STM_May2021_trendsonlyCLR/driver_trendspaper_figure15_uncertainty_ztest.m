clear all
addpath /home/sergio/MATLABCODE_Git

%% see FIND_NWP_MODEL_TRENDS/driver_compare_trends_Day_vs_Night.m --> FIND_NWP_MODEL_TRENDS/get_umbc_day_night_name.m
iV3 = 6;   %% this is for the 400 mb AIRS_UMBC vs ERA5 Fig 13, and uncertainty z-scores Fig 14
iV3 = 9;   %% this is for the ERA5 trend and then the ERA5 trend retrieval, Fig 6

iSaveJGR = +1;
iSaveJGR = -1;

iSmooth = +1;    %% smoothn
iSmooth = +0.5;  %% smoothdata
iSmooth = 00;    %% no smooth


% iV3 = input('Enter iV3 = [9/default] ERA5 trend and retrieval comparison, makes Fig 6 of Wiley paper  [6]  makes Figs 13 (400 mb AIRS_UMBC vs ERA5) and Fig 14 (uncertainty z-scores) of Wiley paper : ');
% if iV3 == 9
%   disp('this is ERA5 trend simulations and retrievals')
%   umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_CAL_startwith0_50fatlayers_NoMODELS_MLS_fortrendspaper.mat';
%   umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_CAL_startwith0_50fatlayers_NoMODELS_MLS_fortrendspaper.mat';
% elseif iV3 == 6
%   %% test, should be same as iV3 == 2, uses MLS
%   umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss_MLS_day.mat';
%   umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss_MLS.mat';
% end
% 
% if iV3 == 9
%   iSmooth = 00;    %% no smooth  
% elseif iV3 == 6
%   iSmooth = +1;    %% smoothn
% end

iV3 = input('Enter iV3 = [9/default] ERA5 trend and retrieval comparison, makes Fig 6 of Wiley paper  [13] makes Figs 13 (400 mb AIRS_UMBC vs ERA5) and [14] Fig 14 (uncertainty z-scores) of Wiley paper : ');
if length(iV3) == 0
  iV3 = 9;
end

%% see ../FIND_NWP_MODEL_TRENDS/loop_driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m

if iV3 == 9 
  disp('this is Wiley fig06 : ERA5 trend simulations and retrievals')
  umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_CAL_startwith0_50fatlayers_NoMODELS_MLS_fortrendspaper.mat';
  umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_CAL_startwith0_50fatlayers_NoMODELS_MLS_fortrendspaper.mat';

  %umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_startwith0_50fatlayers_ERA5calcs_spectraltrends2.mat';
  %umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q16_newERA5_2021jacs_startwith0_50fatlayers_ERA5calcs_spectraltrends2.mat';  
   
elseif iV3 == 13
  disp('this is Wiley fig13 : 400 mb WV, AIRS UMBC vs ERA5 ')
  umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2D.mat';
  umbc_day_file = umbc_night_file;

elseif iV3 == 14
  %% test, should be same as iV3 == 2, uses MLS
  disp('this is Wiley fig14 : uncertainty and z score')
  umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss_MLS_day.mat';
  umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss_MLS.mat';
end

if iV3 == 9
  iSmooth = 00;    %% no smooth  
elseif iV3 == 13 | iV3 == 14
  iSmooth = +1;    %% smoothn
end

if iSmooth >= 1
  disp('using smoothn');
elseif iSmooth == 0.5
  disp('using smoothdata');
elseif iSmooth <= 0
  disp('no smooth');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aD = load(umbc_day_file);
aN = load(umbc_night_file);

results  = (aN.results + aD.results)*0.5;
resultsT  = (aN.resultsT + aD.resultsT)*0.5;
resultsWV = (aN.resultsWV + aD.resultsWV)*0.5;
resultsO3 = (aN.resultsO3 + aD.resultsO3)*0.5;

resultsunc  =  sqrt(aN.resultsunc.^2   + aD.resultsunc.^2)*0.5;
resultsTunc  = sqrt(aN.resultsTunc.^2  + aD.resultsTunc.^2)*0.5;
resultsWVunc = sqrt(aN.resultsWVunc.^2 + aD.resultsWVunc.^2)*0.5;
resultsO3unc = sqrt(aN.resultsO3unc.^2 + aD.resultsO3unc.^2)*0.5;

fracWV  = (aN.fracWV + aD.fracWV)*0.5;
deltaT  = (aN.deltaT + aD.deltaT)*0.5;
deltaRH = (aN.deltaRH + aD.deltaRH)*0.5;

fracWVunc  = sqrt(aN.fracWVunc.^2 + aD.fracWVunc.^2)*0.5;
deltaTunc  = sqrt(aN.deltaTunc.^2  + aD.deltaTunc.^2)*0.5;
deltaRHunc = sqrt(aN.deltaRHunc.^2 + aD.deltaRHunc.^2)*0.5;

iNumLay = aN.iNumLay;
pavg = aN.pavg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
quick_compare_era5_retrieval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
