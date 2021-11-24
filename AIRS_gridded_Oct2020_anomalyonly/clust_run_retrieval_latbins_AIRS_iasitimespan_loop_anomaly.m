%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS_loop_anomaly.m
%---------------------------------------------------------------------------

addpath /home/sergio/MATLABCODE
system_slurm_stats

t1x = tic;

%% this is the timestep : 1: 365 (coincidecne : there are 365 days/year and
%% I did 16 day averages .... so 365/16 steps per year ... and 2002-2018 is
%% 16 years so total of 365/16 * 16 = 365 steps

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 64 for the 64 latbins
JOB = 32

%%%%%%%%%% ANOM or RATES %%%%%%%%%%
% JOBTYPE = 1000;  %%% uncomment this when trying to fit for linear rates!!! fix change_important_topts_settings, and set <<< driver.i16daytimestep = -1 >>>;  iDoAnomalyOrRates = -1; below
%%%%%%%%%% ANOM or RATES %%%%%%%%%%

iDoAnomalyOrRates = -1;  %% do the trends/rates
iDoAnomalyOrRates = +1;  %% do the anomalies

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

iLon0A = 1; iLonEA = 72;  %%%%%% DEFAULT
iLon0A = 1; iLonEA = 1;   %%%%%% DEBUG

if iDoAnomalyOrRates == -1
  JOBTYPE = 1000;
end

%iOffset = (JOB-1)*72;
%iLon0 = iLon0A + iOffset;  iLonE = iLonEA + iOffset; 

for iLon = iLon0A : iLonEA
  for iTime = 1 : 23*17     %% 23 timesteps/year x 17 years
    disp(' ')
    fprintf(1,'lat/lon = %4i %4i timestep = %3i \n',JOB,iLon,iTime);

    %------------------------------------------------------------------------
    %% <<<<<<<    no real need to touch any of this  >>>>>>>>
    driver.iibin     = iLon;

    %%%%%%%%%% ANOM or RATES %%%%%%%%%%
    if iDoAnomalyOrRates == +1
      driver.i16daytimestep = iTime;  %% this is when doing anomaly
    elseif iDoAnomalyOrRates == -1  
      driver.i16daytimestep = -1;   %% for the rates, not anomalies, RUN BY HAND BY UN-COMMENTING THIS LINE and 
                                  %% on top JOB = 1000, in change_important_topts_settings.m also set topts.set_tracegas = -1;
    end
    %%%%%%%%%% ANOM or RATES %%%%%%%%%%

    iLat = JOB;
    iy = iLat;
    ix = iLon;

    driver.latlon.latbin   = iLat;
    driver.latlon.lonbin   = iLon;
    driver.latlon.timestep = iTime;

    driver.oem.dofit = true;
    driver.lls.dofit = false;
    driver.oem.nloop = 2;
    driver.oem.nloop = 3;
    driver.oem.nloop = 1;
    driver.oem.doplots = false;
  %---------------------------------------------------------------------------
    change_important_topts_settings  % Override many settings and add covariance matrix
  
    if topts.ocb_set == 0 & driver.i16daytimestep > 0
      driver.outfilename = ['OutputAnomaly_OBS/' num2str(iLat,'%02d') '/' num2str(iLon,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
    elseif topts.ocb_set == 1 & driver.i16daytimestep > 0
      driver.outfilename = ['OutputAnomaly_CAL/' num2str(iLat,'%02d') '/' num2str(iLon,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
    elseif driver.i16daytimestep < 0
      driver.outfilename = ['Output/test' int2str(iLon) '.mat'];
    end
  
    if topts.iNlays_retrieve >= 97 & ~exist(driver.outfilename)
      [driver,aux] = strow_override_defaults_latbins_AIRS(driver,topts);
    elseif topts.iNlays_retrieve < 97 & ~exist(driver.outfilename)
      [driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,topts.iNlays_retrieve,topts);
    end
    %%[driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,topts.iNlays_retrieve,topts);
  
    %---------------------------------------------------------------------------
    % Do the retrieval
    if ~exist(driver.outfilename)
       driver = retrieval(driver,aux);
    else
      driver.rateset.rates = zeros(2645,1);
    end
    driver
    figure(3); plot(1:2645,driver.rateset.rates,'b',1:2645,driver.oem.fit,'r'); title('AIRS_gridded_Oct2020_trendsonly')
  
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
      save(driver.outfilename,'-struct','driver');
    elseif sum(abs(driver.rateset.rates)) < eps & ~exist(driver.outfilename)
      fprintf(1,'not saving %s since sum(abs(driver.rateset.rates)) = 0 \n',driver.outfilename);
    elseif sum(abs(driver.rateset.rates)) < eps & exist(driver.outfilename)
      fprintf(1,'not saving %s since it already exists\n',driver.outfilename);
     end
  
    %alld(JOB) = driver;

    %---------------------------------------------------------------------------
    % Some simple output
    if sum(abs(driver.rateset.rates)) > 0
      fprintf('Scalar Retrievals from OEM latbin %2i timestep %3i \n',iLon,JOB)
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

       %---------------------------------------------------------------------------
       % Pull interesting variable out for quick look
       if topts.co2lays == 1
         co2(JOB,iLon) = driver.oem.finalrates(1);
         co2_sigs(JOB,iLon) = driver.oem.finalsigs(1); 
         n2o(JOB,iLon) = driver.oem.finalrates(2); 
         n2o_sigs(JOB,iLon) = driver.oem.finalsigs(2); 
         ch4(JOB,iLon) = driver.oem.finalrates(3); 
         ch4_sigs(JOB,iLon) = driver.oem.finalsigs(3); 
         cfc11(JOB,iLon) = driver.oem.finalrates(4); 
         cfc11_sigs(JOB,iLon) = driver.oem.finalsigs(4); 
         cfc12(JOB,iLon) = driver.oem.finalrates(5); 
         cfc12_sigs(JOB,iLon) = driver.oem.finalsigs(5); 
         sst(JOB,iLon) = driver.oem.finalrates(6); 
         sst_sigs(JOB,iLon) = driver.oem.finalsigs(6); 
       elseif topts.co2lays == 3
         co2(JOB,iLon) = driver.oem.finalrates(1);
         co2_sigs(JOB,iLon) = driver.oem.finalsigs(1); 
         n2o(JOB,iLon) = driver.oem.finalrates(4); 
         n2o_sigs(JOB,iLon) = driver.oem.finalsigs(4); 
         ch4(JOB,iLon) = driver.oem.finalrates(5); 
         ch4_sigs(JOB,iLon) = driver.oem.finalsigs(5); 
         cfc11(JOB,iLon) = driver.oem.finalrates(6); 
         cfc11_sigs(JOB,iLon) = driver.oem.finalsigs(6); 
         cfc12(JOB,iLon) = driver.oem.finalrates(7); 
         cfc12_sigs(JOB,iLon) = driver.oem.finalsigs(7); 
         sst(JOB,iLon) = driver.oem.finalrates(8); 
         sst_sigs(JOB,iLon) = driver.oem.finalsigs(8); 
       end
     end
  end  %% end of lonbin loop
end    % end of latbin loop  
%---------------------------------------------------------------------------
% Close debug file
if driver.debug
   writelog('close')
end
%--------------------------------------------------------------------------

t2x = toc(t1x);
fprintf(1,'time taken (in seconds) for a %3i layer retrievals = %8.4f \n',topts.iNlays_retrieve,t2x)