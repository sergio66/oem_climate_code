iSave = -1; savename  = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithERA5_uncX3.mat';
iSave = +1; savename  = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX3_50fatlayers_AIRSL3_ERA5_CMIP6_feedback.mat';

iSave = -1;  savename = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithMLSL3_uncX3_50fatlayers_AIRSL3_ERA5_CMIP6_feedback.mat';
iSave = -1;  savename = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwithMLSL3_uncX100_50fatlayers_AIRSL3_ERA5_CMIP6_feedback.mat';
iSave = -1;  savename = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX100_50fatlayers_AIRSL3_ERA5_CMIP6_feedback.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'set_gather_savename_rates_AIRSL3_NWP_XMIP.m : MAIN FILE TO READ IN = %s \n',savename);
            savename2 = '/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX3_50fatlayers_CLIMCAPS_MERRA2_AMIP6_feedback.mat';

if iSave > 0
  junk = findstr(savename,'_2021jacs');
  savestr = savename(junk+9:end);
  savestr = savestr(1:end-4);
end