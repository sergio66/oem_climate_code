%% these are from driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m
%% these two are fig12 of JR paper, they have MLS
fumbc = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_MLS.mat';
fumbc = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss_MLS_day.mat';

%% these supposedly have no MLS
fumbc = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q05_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_feedback.mat';
fumbc = '/asl/s1/sergio/JUNK/test9_guessstartWV_Vers1_march22_2023.mat';
fumbc = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss.mat';
fumbc = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';
fumbc = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q05_newERA5_2021jacs_startwith_MLSL3_TOA_guessWV_dRH_zero_bot_50fatlayers.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% these are from get_umbc_day_night_name.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boo = load(fumbc,'rlat','resultsWV','pjunk20');
load llsmap5; booWV = squeeze(nanmean(reshape(boo.resultsWV,72,64,49),1)); figure(1); clf; pcolor(boo.rlat,boo.pjunk20,booWV'); 
set(gca,'ydir','reverse'); colormap(llsmap5); caxis([-1 +1]*0.01); ylim([100 1000]); shading interp; colorbar
