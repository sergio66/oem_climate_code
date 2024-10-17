%%%% see loop_driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m
%%%% see loop_driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m
%%%% see loop_driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends.m

iV = input('Enter version of UMBC file (-1) GULP, old stuff, (0) oringinally was in test.pdf (dRHdt = 0 ... trends paper) (1) for (Held/Jeevanjee) with d lnP/dSKT = 0.02 (constant) (2) for (Held/Jeevanjee) with zonally varying d lnP/dSKT from IPCC 2007 (3/default) YAYAYAYA clrjac in test.pdf (dRH/dt ~ 0.01) : ');
if length(iV) == 0
  iV = 0;
  iV = 3;
end

iV3 = -1;  %% just for printing purposes at end

if iV == -1
  %% GULP
  umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2E_day.mat';
  umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2E.mat';

elseif iV == 0
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
  iV3 = 0;  %% Dec 26, 2023 commit .. but I had messed up args for jeevanjee_PLB_deltaRH_Ts(Ts,RHs,dlnPdT,trend_BT1231)  : used jeevanjee_PLB_deltaRH_Ts(Ts,RHs,trend_BT1231,dlnPdT)
  iV3 = 1;
  iV3 = 2;  %% YAY, old default
  iV3 = 3;  %% YAY, new default???, NOT FOR dcolWV/dt ... pretty bad!!!!!!!!!!!
  iV3 = 4;
  iV3 = 5;  %% testing is this same as iV3 = 2, except I have slightly switched  doing things in p.landfrac < 0.25 insteaad of eps, and also made sure args to jevanjee_PLB_deltaRH_Ts(a,b,c,d) are correct order
  iV3 = 6;  %% testing is this same as iV3 = 2 BUT has MLS in UT/LS, but is THE NEW DEFAUKT  not so great for dColWV/dt

  iV3 = 5;  %% YAY
  iV3 = 2;  %% YAY
  iV3 = 3;  %% YAY
  iV3 = 6;  %% testing is this same as iV3 = 2 BUT has MLS in UT/LS

  if iV3 == 0
    %% when I first found my jacobian mistake .. but then I made a mistake in args of jeevanjee_PLB_deltaRH_Ts(a,b,d,c) switched args
    %% thi is redone from  /asl/s1/sergio/JUNK/gitjunk5/oem_climate_code/AIRS_gridded_STM_May2021_trendsonlyCLR in March 2024 to test codes
    umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_gitjunk5_dec29_2023_day.mat';
    umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_gitjunk5_dec29_2023.mat';

  elseif iV3 == 1
    %% dRH = f(dP/dt) (Held/Jeevanjee) with zonally varying d lnP/dSKT from IPCC 2007, clrjac  -- BAD , have a bug, done Dec 27, 2023 .... but have the jevanjee_PLB_deltaRH_Ts(a,b,d,c) switched args
    umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_day.mat';
    umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac.mat';

  %%%%%%%%%%%%%%%%%%%%%%%%%
  elseif iV3 == 2
    %% ooops found bug in set_CO2_CH4_N2O_ESRL.m, calling arguments in jeevanjee_PBL_deltaRH_Ts  -- YAY YAY YAY, done Jan 1, 2024
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
    %% but see AIRS_gridded_STM_May2021_trendsonlyCLR/set_CO2_CH4_N2O_ESRL.m and search for <<<<<< till Dec 31, 2023 MAJOR DIFF >>>>>>
    %%   (a) ppp.landfrac(JOBJOBJOB) < eps     VS        ppp.landfrac(JOBJOBJOB) < 0.25
    %%   (b) jeevanjee_PBL_deltaRH_Ts(ppp.stemp(JOBJOBJOB),RHSurf(JOBJOBJOB)/100,dBT1231,dlnPdT) vs eevanjee_PBL_deltaRH_Ts(ppp.stemp(JOBJOBJOB),RHSurf(JOBJOBJOB)/100,dlnPdT,dBT1231);
    umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss_day.mat';
    umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss.mat';

  elseif iV3 == 6
    %% test, should be same as iV3 == 2, uses MLS
    umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss_MLS_day.mat';
    umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_NoMODELS_testjunk_constmiss_MLS.mat';
  end

end  

disp(' ')
fprintf(1,'iV,iV3 = %2i %2i umbc_day_file   = %s \n',iV,iV3,umbc_day_file);
fprintf(1,'iV,iV3 = %2i %2i umbc_night_file = %s \n',iV,iV3,umbc_night_file);
disp(' ')
