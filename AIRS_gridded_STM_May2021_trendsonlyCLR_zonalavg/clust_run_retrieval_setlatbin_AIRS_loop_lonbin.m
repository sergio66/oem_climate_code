%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS_loop_anomaly.m
%---------------------------------------------------------------------------

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil

%system_slurm_stats

t1x = tic;

%{
[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
scatter_coast(p.rlon,p.rlat,50,p.nemis)
mooO = find(p.nemis == 19);  i900_O = find(p.efreq(1:19,mooO(1)) >= 900,1);
mooL = find(p.nemis == 100); i900_L = find(p.efreq(1:100,mooL(1)) >= 900,1);
junk = [mooO mooL]; junkx = [p.emis(i900_O,mooO) p.emis(i900_L,mooL)];
scatter_coast(p.rlon(junk),p.rlat(junk),50,junkx)
%}

clear h ha p pa

%% this is the timestep : 1: 365 (coincidecne : there are 365 days/year and
%% I did 16 day averages .... so 365/16 steps per year ... and 2002-2018 is
%% 16 years so total of 365/16 * 16 = 365 steps

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 64 for the 64 latbins
%JOB = 33
%JOB = 47
%JOB = 22
%JOB = 23 %% bad

%%%%%%%%%% ANOM or RATES %%%%%%%%%%
% JOBTYPE = 1000;  %%% uncomment this when trying to fit for linear rates!!! fix change_important_topts_settings, and set <<< driver.i16daytimestep = -1 >>>;  iDoAnomalyOrRates = -1; below
%%%%%%%%%% ANOM or RATES %%%%%%%%%%

iDoAnomalyOrRates = +1;  %% do the anomalies
iDoAnomalyOrRates = -1;  %% do the trends/rates

%---------------------------------------------------------------------------
addpath /home/sergio/MATLABCODE/oem_pkg
addpath Plotutils
%---------------------------------------------------------------------------
% Doing debug?
 driver.debug = false;
 driver.debug_dir = '../Debug';

% Open debug file if desired
 if driver.debug
    writelog('open');
 end;
%---------------------------------------------------------------------------
% for this JOB timestep (1:388), loop over 64x72 grid points
%---------------------------------------------------------------------------
%JOB = 1000; iLon0 =  1; iLonE = 64*72;  %% trends
iLon0A = 1; iLonEA = 72;

if iDoAnomalyOrRates == -1
  JOBTYPE = 1000;
end

for iInd = 1

  disp(' ')
  fprintf(1,'latbin = %3i \n',JOB);

%------------------------------------------------------------------------
%% <<<<<<<    no real need to touch any of this  >>>>>>>>
  driver.iibin     = JOB;

  %%%%%%%%%% ANOM or RATES %%%%%%%%%%
  if iDoAnomalyOrRates == +1
    driver.i16daytimestep = JOB;  %% this is when doing anomaly
  elseif iDoAnomalyOrRates == -1  
    driver.i16daytimestep = -1;   %% for the rates, not anomalies, RUN BY HAND BY UN-COMMENTING THIS LINE and 
                                  %% on top JOB = 1000, in change_important_topts_settings.m also set topts.set_tracegas = -1;
  end
  %%%%%%%%%% ANOM or RATES %%%%%%%%%%
  driver.iDebugRatesUseNWP = 32; %% use AIRS L3 reconstructed trends
  driver.iDebugRatesUseNWP = 62; %% use CMIP6    reconstructed trends
  driver.iDebugRatesUseNWP = 52; %% use ERA      reconstructed trends
  driver.iDebugRatesUseNWP = -1; %% use AIRS observed trends

  iQuantile = 04;  %% 05% so very cloudy (hope SST jac can take care of that) and convection
  iQuantile = 14;  %% 
  iQuantile = 08;  %% 50% so has clouds (hope SST jac can take care of that) and convection -- this is 0.25 - 0.50
  iQuantile = 09;  %% 50% so has clouds (hope SST jac can take care of that) and convection -- this is 0.90 - 0.75
  iQuantile = 04;  %% quite cloudy (hopefully)
  iQuantile = 00;  %% mean     <<<<***** IF YOU SET THIS THEN topts.dataset is ignored, uses topts.dataset   = -3; *****>>>>
  iQuantile = 16;  %% hottest, for AIRS STM

  driver.NorD = -1; %% day, asc
  driver.NorD = +1; %% night, desc

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  driver.iQuantile = iQuantile;
  ix = iInd;
  %driver.iLon = iInd-iOffset;
  driver.iLat = JOB;

  driver.oem.dofit = true;
  driver.lls.dofit = false;
  driver.oem.nloop = 2;
  driver.oem.nloop = 3;
  driver.oem.nloop = 1;
  driver.oem.doplots = false;
%---------------------------------------------------------------------------
  change_important_topts_settings  % Override many settings and add covariance matrix
%------------------------------------------------------------------------

  topts.dataset   = -1;   %% (-1) AIRS 18 year quantile dataset, Sergio Aug 2021   2002/09-2020/08 
  topts.dataset   = +1;   %% (+1) AIRS 18 year quantile dataset, Strow  March 2021 2002/09-2020/08 
  topts.dataset   = +3;   %% (+3) AIRS 19 year extreme  dataset, Sergio Aug 2021   2002/09-2021/08
  topts.dataset   = -3;   %% (-3) AIRS 19 year mean     dataset, Sergio Aug 2021   2002/09-2020/08  AUTOMATIC USES Q00
  topts.dataset   = +2;   %% (+2) AIRS 19 year quantile dataset, Sergio Aug 2021   2002/09-2021/08

  %topts.iFixTG_NoFit = +1; %% dump out first scalar = CO2 boy certainly messes up CO2 even more!!!!!

  %% iNorD > 0 ==> night
  if topts.ocb_set == 0 & driver.i16daytimestep > 0 & driver.NorD > 0 & topts.dataset < 3
    driver.outfilename = ['OutputAnomaly_OBS/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(JOB,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif topts.ocb_set == 1 & driver.i16daytimestep > 0 & driver.NorD > 0 & topts.dataset < 3
    driver.outfilename = ['OutputAnomaly_CAL/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(JOB,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif driver.i16daytimestep < 0 & driver.NorD > 0 & topts.dataset < 3
    outdir = ['Output/Quantile' num2str(driver.iQuantile,'%02d')];
    driver.outfilename = [outdir  '/test' int2str(JOB) '.mat'];
  elseif driver.i16daytimestep < 0 & driver.NorD > 0 & topts.dataset == 3
    outdir = ['Output/Extreme/'];
    driver.outfilename = [outdir  '/test' int2str(JOB) '.mat'];
  %% iNorD < 0 ==> day
  elseif topts.ocb_set == 0 & driver.i16daytimestep > 0 & driver.NorD < 0 & topts.dataset < 3
    driver.outfilename = ['OutputAnomaly_OBS_Day/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(JOB,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif topts.ocb_set == 1 & driver.i16daytimestep > 0 & driver.NorD < 0 & topts.dataset < 3
    driver.outfilename = ['OutputAnomaly_CAL_Dat/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(JOB,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif driver.i16daytimestep < 0 & driver.NorD < 0 & topts.dataset < 3
    outdir = ['Output_Day/Quantile' num2str(driver.iQuantile,'%02d')];
    driver.outfilename = [outdir  '/test' int2str(JOB) '.mat'];
  elseif driver.i16daytimestep < 0 & driver.NorD < 0 & topts.dataset == 3
    outdir = ['Output_Day/Extreme/'];
    driver.outfilename = [outdir  '/test' int2str(JOB) '.mat'];
  end
  if ~exist(outdir)
    mker = ['!mkdir -p ' outdir];
    eval(mker);
  end

  oldstyle_removeCO2 = +1;
  if oldstyle_removeCO2 == -1
    error('jhave not coded up REMOVING CO2 from retrieval set')
  end

  if ~isfield(topts,'iFixTG_NoFit')
    oldstyle_removeCO2 = +1;
  elseif isfield(topts,'iFixTG_NoFit')
    if length(intersect(topts.iFixTG_NoFit,2)) >  0
      oldstyle_removeCO2 = -1;
    end
  end
  if oldstyle_removeCO2 > 0 
    if topts.iNlays_retrieve >= 97 & ~exist(driver.outfilename)
      [driver,aux] = strow_override_defaults_latbins_AIRS(driver,topts);
    elseif topts.iNlays_retrieve < 97 & ~exist(driver.outfilename)
      [driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,topts.iNlays_retrieve,topts);
    end
    %%[driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,topts.iNlays_retrieve,topts);
  else
    if topts.iNlays_retrieve >= 97 & ~exist(driver.outfilename)
      [driver,aux,co2rate] = strow_override_defaults_latbins_AIRS_removeCO2(driver,topts);
    elseif topts.iNlays_retrieve < 97 & ~exist(driver.outfilename)
      [driver,aux,co2rate] = strow_override_defaults_latbins_AIRS_fewlays_removeCO2(driver,topts.iNlays_retrieve,topts);
    end
 end

%---------------------------------------------------------------------------
  % Do the retrieval
  if ~exist(driver.outfilename)     
     driver = retrieval(driver,aux);
     if oldstyle_removeCO2 < 0
       %% now adjust everything!!!!!
       call_adjustCO2
     end
     figure(3); plot(1:2645,driver.rateset.rates,'b',1:2645,driver.oem.fit,'r'); title('AIRS_gridded_Oct2020_trendsonly')
  else
    fprintf(1,'%s already exists \n',driver.outfilename)
    driver.rateset.rates = zeros(2645,1);
  end
  
%error('the end IRS_gridded_Oct2020_trendsonly')

%---------------------------------------------------------------------------
  % Save retrieval output from this loop

  if isfield(topts,'iFixTz_NoFit')
    if topts.iFixTz_NoFit > 0
      junk = 1:length(driver.oem.finalrates)+length(aux.orig_temp_i);

      %% put ozone into correct (expected) spot
      driver.oem.finalrates(aux.orig_ozone_i) = driver.oem.finalrates(aux.orig_temp_i);
      driver.oem.finalsigs(aux.orig_ozone_i)  = driver.oem.finalsigs(aux.orig_temp_i);

      %% put fixed/unchanging T anomaly
      driver.oem.finalrates(aux.orig_temp_i) = aux.FixTz_NoFit;
      driver.oem.finalsigs(aux.orig_temp_i)  = 0;

      driver.jacobian.temp_i  = driver.jacobian.ozone_i;
      driver.jacobian.ozone_i = driver.jacobian.ozone_i + driver.jacobian.numlays;

    end
  end

  if isfield(topts,'iFixO3_NoFit')
    if topts.iFixO3_NoFit > 0

      %% put ozone into correct (expected) spot
      driver.oem.finalrates(aux.orig_ozone_i) = aux.FixO3_NoFit;
      driver.oem.finalsigs(aux.orig_ozone_i)  = 0;

      driver.jacobian.ozone_i = driver.jacobian.temp_i + driver.jacobian.numlays;
    end
  end

  if sum(abs(driver.rateset.rates)) > 0 & ~exist(driver.outfilename)
    fprintf(1,'saving to %s \n',driver.outfilename)
    save(driver.outfilename,'-struct','driver');
  elseif sum(abs(driver.rateset.rates)) < eps & ~exist(driver.outfilename)
    fprintf(1,'not saving %s since sum(abs(driver.rateset.rates)) = 0 \n',driver.outfilename);
  elseif sum(abs(driver.rateset.rates)) < eps & exist(driver.outfilename)
    fprintf(1,'not saving %s since it already exists\n',driver.outfilename);
   end

%   alld(JOB) = driver;

%---------------------------------------------------------------------------
% Some simple output
   if sum(abs(driver.rateset.rates)) > 0
     fprintf('Scalar Retrievals from OEM latbin %3i \n',JOB)
     if topts.co2lays == 1
       fprintf(1,'CO2   (ppm)   %5.3f  +- %5.3f \n',driver.oem.finalrates(1),driver.oem.finalsigs(1));
       fprintf(1,'N2O   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(2),driver.oem.finalsigs(2));
       fprintf(1,'CH4   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(3),driver.oem.finalsigs(3));
       fprintf(1,'CFC11 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(4),driver.oem.finalsigs(4));
       fprintf(1,'CFC12 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(5),driver.oem.finalsigs(5));
       fprintf(1,'SST   (K)    %5.3f  +- %5.3f \n',driver.oem.finalrates(6),driver.oem.finalsigs(6));
     elseif topts.co2lays == 3
       fprintf(1,'CO2 lower trop  (ppm)   %5.3f  +- %5.3f \n',driver.oem.finalrates(1),driver.oem.finalsigs(1));
       fprintf(1,'N2O   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(4),driver.oem.finalsigs(4));
       fprintf(1,'CH4   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(5),driver.oem.finalsigs(5));
       fprintf(1,'CFC11 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(6),driver.oem.finalsigs(6));
       fprintf(1,'CFC12 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(7),driver.oem.finalsigs(7));
       fprintf(1,'SST   (K)    %5.3f  +- %5.3f \n',driver.oem.finalrates(8),driver.oem.finalsigs(8));
     end

%% error('nyuk debug')

     %---------------------------------------------------------------------------
     % Pull interesting variable out for quick look
     if topts.co2lays == 1
       co2(JOB) = driver.oem.finalrates(1);
       co2_sigs(JOB) = driver.oem.finalsigs(1); 
       n2o(JOB) = driver.oem.finalrates(2); 
       n2o_sigs(JOB) = driver.oem.finalsigs(2); 
       ch4(JOB) = driver.oem.finalrates(3); 
       ch4_sigs(JOB) = driver.oem.finalsigs(3); 
       cfc11(JOB) = driver.oem.finalrates(4); 
       cfc11_sigs(JOB) = driver.oem.finalsigs(4); 
       cfc12(JOB) = driver.oem.finalrates(5); 
       cfc12_sigs(JOB) = driver.oem.finalsigs(5); 
       sst(JOB) = driver.oem.finalrates(6); 
       sst_sigs(JOB) = driver.oem.finalsigs(6); 
     elseif topts.co2lays == 3
       co2(JOB) = driver.oem.finalrates(1);
       co2_sigs(JOB) = driver.oem.finalsigs(1); 
       n2o(JOB) = driver.oem.finalrates(4); 
       n2o_sigs(JOB) = driver.oem.finalsigs(4); 
       ch4(JOB) = driver.oem.finalrates(5); 
       ch4_sigs(JOB) = driver.oem.finalsigs(5); 
       cfc11(JOB) = driver.oem.finalrates(6); 
       cfc11_sigs(JOB) = driver.oem.finalsigs(6); 
       cfc12(JOB) = driver.oem.finalrates(7); 
       cfc12_sigs(JOB) = driver.oem.finalsigs(7); 
       sst(JOB) = driver.oem.finalrates(8); 
       sst_sigs(JOB) = driver.oem.finalsigs(8); 
     end
   end

   % Plot Results
%{
  if topts.iNlays_retrieve >= 97
    plot_retrieval_latbins
  else
    plot_retrieval_latbins_fewlays
  end
   %disp('Hit return for next latitude'); pause
   pause(0.1)
%}

end % end of latbin loop  
%---------------------------------------------------------------------------
% Close debug file
if driver.debug
   writelog('close')
end
%--------------------------------------------------------------------------

%{
if topts.iNlays_retrieve >= 97
  plot_all_latbins_anom
else
  plot_all_latbins_fewlays_anom
end
%}

t2x = toc(t1x);
fprintf(1,'time taken (in seconds) for a %3i layer retrievals = %8.4f \n',topts.iNlays_retrieve,t2x)

if driver.i16daytimestep < 0
  plot_one_latlonbin_fewlays
end
