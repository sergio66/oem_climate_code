clear showfeedbacks* strfeedbacks
clear all

figure(1); clf
figure(2); clf
figure(3); clf
figure(4); clf
figure(5); clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 20;  %% use CarbonTracker CO2 trends
fprintf(1,'iNumYears = %2i .. reading in %s \n',iNumYears,strUMBC);
%loader = ['load ' strUMBC];
%eval(loader);

read_fileMean17years
h = hMean17years;
p = pMean17years;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];
loader = ['load ' feedbacknameUMBC];
eval(loader)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
feedbacknameNWP_ERA5 = '/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_20.mat';
load(feedbacknameNWP_ERA5);
era5.trend_stemp = stemptrend.era5;   
aL3trend.stemp   = stemptrend.airsL3;
c6trend.stemp    = stemptrend.cmip6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

feedbacknameNWP_MERRA2 = '/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_20.mat';
load(feedbacknameNWP_MERRA2);
merra2.trend_stemp = stemptrend.era5;   
cL3trend.stemp     = stemptrend.airsL3;
a6trend.stemp      = stemptrend.cmip6;
amip6.stempjunk    = cL3trend.stemp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load latB64.mat
  %% I think before Sept 2023 I had this wrong
  %% mean weighted delta SST rate = 0.031962  0.003305  0.023971  0.019594 K/yr for 05/10/15/20 years   WRONG?? CORRECT?? 
  %% mean weighted delta SST rate = 0.069633  0.020002  0.028442  0.024870 K/yr for 05/10/15/20 years   CORRECT? WRONG??
  do_XX_YY_from_X_Y

disp(' ')
disp('PLEASE NOTE!!!!! calling show_olr_ecRad_feedback : do NOT do NOT do NOT    clear nwp_spectral_trends_cmip6_era5_airsL3_umbc or nwp_spectral_trends_amip6_merra2_climcapsL3_umbc');
disp('PLEASE NOTE!!!!! calling show_olr_ecRad_feedback : do NOT do NOT do NOT    clear nwp_spectral_trends_cmip6_era5_airsL3_umbc or nwp_spectral_trends_amip6_merra2_climcapsL3_umbc');
disp('PLEASE NOTE!!!!! calling show_olr_ecRad_feedback : do NOT do NOT do NOT    clear nwp_spectral_trends_cmip6_era5_airsL3_umbc or nwp_spectral_trends_amip6_merra2_climcapsL3_umbc');
disp('PLEASE NOTE!!!!! calling show_olr_ecRad_feedback : do NOT do NOT do NOT    clear nwp_spectral_trends_cmip6_era5_airsL3_umbc or nwp_spectral_trends_amip6_merra2_climcapsL3_umbc');
disp(' ')

iaComputeWhichFeedback = [0];     %% compute + plot feedbacks only
nwp_spectral_trends_cmip6_era5_airsL3_umbc = [];
nwp_spectral_trends_amip6_merra2_climcapsL3_umbc = [];

show_olr_ecRad_feedback
quick_save_olr_feedbacks_umbc_NWP_L3_XMIP6

disp('DONE, now showing results <RET> to continue'); pause
%plot_show_olr_ecRad_feedback_globalsstfit
%plot_show_olr_ecRad_feedback_globalsstfit_smooth
plot_show_olr_ecRad_feedback_globalsstfit_smooth2

