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
  %% dRH = f(dP/dt) (Held/Jeevanjee) with zonally varying d lnP/dSKT from IPCC 2007, clrjac
  umbc_day_file    = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_day.mat';
  umbc_night_file  = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac.mat';

end  
