%%%% see loop_driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m
%%%% see loop_driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m
%%%% see loop_driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m

%% GULP
%umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2E_day.mat';
%umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2E.mat';

iV = input('Enter version of UMBC file (0) oringinally was in test.pdf (dRHdt = 0 ... trends paper) (1) for (Held/Jeevanjee) with d lnP/dSKT = 0.02 (constant) (2) for (Held/Jeevanjee) with zonally varying d lnP/dSKT from IPCC 2007 (3/default) YAYAYAYA clrjac in test.pdf (dRH/dt ~ 0.01) : ');
if length(iV) == 0
  iV = 0;
  iV = 3;
end

if iV == 0
  %% SEQN
  %% dRH = 0 (Sergio)
  umbc_day_file    = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2D_day.mat';  
  umbc_night_file  = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2D.mat';
  
elseif iV == 1
  %% dRH = f(dP/dt) (Held/Jeevanjee) with d lnP/dSKT = 0.02 (constant)
  umbc_day_file    = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_day.mat';
  umbc_night_file  = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F.mat';
  
elseif iV == 2
  %% dRH = f(dP/dt) (Held/Jeevanjee) with zonally varying d lnP/dSKT from IPCC 2007
  umbc_day_file    = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_vary_dlnPdST_day.mat';
  umbc_night_file  = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_vary_dlnPdST.mat';

elseif iV == 3
  iV3 = 1;
  iV3 = 2;  %% YAY, old default
  iV3 = 3;  %% YAY, new default???
  iV3 = 4;
  iV3 = 5;  %% testing is this same as iV3 = 2
  iV3 = 6;  %% testing is this same as iV3 = 2 BUT has MLS in UT/LS

  iV3 = 3;  %% YAY
  iV3 = 6;  %% testing is this same as iV3 = 2 BUT has MLS in UT/LS

  if iV3 == 1
    %% dRH = f(dP/dt) (Held/Jeevanjee) with zonally varying d lnP/dSKT from IPCC 2007, clrjac  -- BAD , have a bug
    umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_day.mat';
    umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac.mat';

  %%%%%%%%%%%%%%%%%%%%%%%%%
  elseif iV3 == 2
    %% ooops found bug in set_CO2_CH4_N2O_ESRL.m, calling arguments in jeevanjee_PBL_deltaRH_Ts  -- YAY YAY YAY
    umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_dayV2.mat';
    umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjacV2.mat';

  elseif iV3 == 3
    %% ooops found bug in set_CO2_CH4_N2O_ESRL.m, calling arguments in jeevanjee_PBL_deltaRH_Ts and now introducing emissivity trends -- YAY YAY const emiss, used in paper now?? in poster?
    umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_day_removeemisstrends.mat';
    umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_removeemisstrends.mat';

  %%%%%%%%%%%%%%%%%%%%%%%%%

  elseif iV3 == 4
    %% compare effects of chaging land emiss, mix of iV3 = 2,3
    umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjacV2.mat';
    umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_removeemisstrends.mat';

  %%%%%%%%%%%%%%%%%%%%%%%%%
  elseif iV3 == 5
    %% test, should be same as iV3 == 2
    umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss.mat';
    umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss.mat';

  elseif iV3 == 6
    %% test, should be same as iV3 == 2
    umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss_MLS_day.mat';
    umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss_MLS.mat';
  end

end  

fprintf(1,'umbc_day_file   = %s \n',umbc_day_file);
fprintf(1,'umbc_night_file = %s \n',umbc_night_file);
