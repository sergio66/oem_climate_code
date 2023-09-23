%% copied from clust_run_retrieval_setlatbin_AIRS_loop_lonbin.m

clear h ha p pa

%% this is the timestep : 1: 365 (coincidence : there are 365 days/year and
%% I did 16 day averages .... so 365/16 steps per year ... and 2002-2018 is
%% 16 years so total of 365/16 * 16 = 365 steps

if exist('JOBM')
  JOB = JOBM;
else
  JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 64 for the 64 latbins
end

%% JOB = 1 .. 10
%% JOB 1 == runs timesteps (1..454) and JOB 10 runs steps (4087 .. 4540)
JOB = input('Enter JOB (1..10) ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iDebug = -1;

ia_OorC_DataSet_Quantile = [+1 09 16  2]; %% ocb_set = 1 : MERRA2     cal fit, dataset = 9, iQuantile = 16  20 year rates
ia_OorC_DataSet_Quantile = [+1 09 16  5]; %% ocb_set = 1 : ERA5       cal fit, dataset = 9, iQuantile = 16  20 year rates
ia_OorC_DataSet_Quantile = [+1 09 16  3]; %% ocb_set = 1 : AIRSL3     cal fit, dataset = 9, iQuantile = 16  20 year rates, problems at eg latbin 1, lonbin 15-65
ia_OorC_DataSet_Quantile = [+1 09 16 -3]; %% ocb_set = 1 : CLIMCAPSL3 cal fit, dataset = 9, iQuantile = 16  20 year rates

ia_OorC_DataSet_Quantile = [+0 10 03 -9999]; %% ocb_set = 0 : obs fit, dataset = 10,iQuantile = 03    05 year rates, AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+0 12 03 -9999]; %% ocb_set = 0 : obs fit, dataset = 12,iQuantile = 03    15 year rates, AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+0 11 03 -9999]; %% ocb_set = 0 : obs fit, dataset = 11,iQuantile = 03    10 year rates, AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+0 09 03 -9999]; %% ocb_set = 0 : obs fit, dataset = 9, iQuantile = 03    20 year rates, AIRS obs Q(0.90-->1)

ia_OorC_DataSet_Quantile = [+2 09 03 -9999]; iNumAnomTimeSteps = 575; iNumAnomTiles = 20; iNumAnomJobsPerProc = 100; %% ocb_set = 2 : anomaly fit, dataset = 9, iQuantile = 03    25 year anomalies== > 20yrs* 23steps/yr = 575; AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+2 09 03 -9999]; iNumAnomTimeSteps = 454; iNumAnomTiles = 10; iNumAnomJobsPerProc =  72; %% ocb_set = 2 : anomaly fit, dataset = 9, iQuantile = 03    20 year anomalies== > 20yrs* 23steps/yr = 460; AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+2 09 03 -9999]; iNumAnomTimeSteps = 454; iNumAnomTiles = 19; iNumAnomJobsPerProc =  72; %% ocb_set = 2 : anomaly fit, dataset = 9, iQuantile = 03    20 year anomalies== > 20yrs* 23steps/yr = 460; AIRS obs Q(0.90-->1)

if ia_OorC_DataSet_Quantile(1) == 2
  disp('I suggest running test_run_retrieval_setlatbin_AIRS_loop_anomaly.m BEFOREHAND to see how many processors you need')
  disp('may be more than your usual 64!')
  clear mapperAnom2Processor
  iLocalTimeStep = 0;
  iPreviouslatnumber = 0;
  for input_spectrum_number = 1 : iNumAnomTimeSteps * iNumAnomTiles
    procnumber = floor((input_spectrum_number-1)/(iNumAnomJobsPerProc) + 1);
    latnumber = floor((input_spectrum_number-1)/(iNumAnomTimeSteps) + 1);
    if latnumber == iPreviouslatnumber
      iLocalTimeStep = iLocalTimeStep  + 1;
    else  
      iLocalTimeStep = 1;
    end
    iPreviouslatnumber = latnumber;
    fprintf(1,'input_spectrum_number : locallatbin localtimestep ---> procnumber = %4i : %4i %4i -----> %4i \n',input_spectrum_number,latnumber,iLocalTimeStep,procnumber);
    mapAnomData_to_processor(input_spectrum_number,1) = latnumber;
    mapAnomData_to_processor(input_spectrum_number,2) = iLocalTimeStep;
    mapAnomData_to_processor(input_spectrum_number,3) = procnumber;
  end
  clear latnumber iLocalTimeStep procnumber iPreviouslatnumber input_spectrum_number

  fprintf(1,' << ANOMALIES : with this configuration you need %3i processors! >>> \n',max(mapAnomData_to_processor(:,3)))
  fprintf(1,' << ANOMALIES : with this configuration you need %3i processors! >>> \n',max(mapAnomData_to_processor(:,3)))
  fprintf(1,' << ANOMALIES : with this configuration you need %3i processors! >>> \n',max(mapAnomData_to_processor(:,3)))
  disp('ret to continue'); pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% MAIN CODE %%%%%%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iDebug > 0
  if ia_OorC_DataSet_Quantile(1) < 2
    plot(1:4608,floor(((1:4608)-1)/72)+1)   %% this maps tile number and job (in sets of 72)
    JOB = floor((iDebug-1)/72+1);
  else
    plot(1:iNumAnomTimeSteps*iNumAnomTiles,floor(((1:iNumAnomTimeSteps*iNumAnomTiles)-1)/iNumAnomTimeSteps)+1)   %% this maps tile number and job (in sets of iNumTimeSteps)  note I only have iNumAnomTimeSteps*iNumAnomTiles timesteps
    JOB = floor((iDebug-1)/iNumAnomTimeSteps+1);
  end
end

%%%%%%%%%% ANOM or RATES %%%%%%%%%%
JOBTYPE = -1000;  %%% uncomment this when trying to fit for linear rates!!! fix change_important_topts_settings, and set <<< driver.i16daytimestep = -1 >>>;  iDoAnomalyOrRates = -1; below
%%%%%%%%%% ANOM or RATES %%%%%%%%%%

iDoAnomalyOrRates = +1;  %% do the anomalies
iDoAnomalyOrRates = -1;  %% do the trends/rates, default
if ia_OorC_DataSet_Quantile(1) == 2
  iDoAnomalyOrRates = +1;  %% do the anomalies
end

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
if iDoAnomalyOrRates == -1
  JOBTYPE = -1000;
  iLon0A = 1; iLonEA = 72;
  iOffset = (JOB-1)*72;
  iInd0 = iLon0A + iOffset;  iIndE = iLonEA + iOffset; 

  %%% this is new
  if iDebug > 0
    iInd0 = iDebug - (JOB-1)*72;
    iIndE = iInd0;
    fprintf(1,'iDebug = %4i   JOB (latbin) = %2i iInd0 (LonBin) = %2i \n',iDebug,JOB,iInd0)
    iInd0 = (JOB-1)*72 + iInd0;
    iIndE = iInd0;
  end

else
  JOBTYPE = +1000;
  mapperAnom2Processor = find(mapAnomData_to_processor(:,3) == JOB);
  iInd0 = mapperAnom2Processor(1);  iIndE = mapperAnom2Processor(end);
  fprintf(1,'JOB = %3i iInd0,iIndE = %4i %4i \n',JOB,iInd0,iIndE)
end

if ~exist('iForward')
  iForward = -1;   %% so start at 72 and go down to 01
  iForward = +1;   %% so start at 01 and go up   to 72
end
if iForward == +1
  iXX1 = iInd0; iXX2 = iIndE; idX = +1;
else
  iXX2 = iInd0; iXX1 = iIndE; idX = -1;
end

%for iInd = iIndE : -1 : iInd0
%for iInd = iInd0 + 21;
%for iInd = iInd0 : iIndE

for iInd = iXX1 : idX : iXX2

  %% so this is really iInd into 1:4608, 72 at a time
  %disp(' ')

%------------------------------------------------------------------------
%% <<<<<<<    no real need to touch any of this  >>>>>>>>
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called

  driver.ia_OorC_DataSet_Quantile = ia_OorC_DataSet_Quantile;
  driver.iibin     = iInd;

  driver.iAllorSeasonal = -1; %% DJF
  driver.iAllorSeasonal = -2; %% MAM
  driver.iAllorSeasonal = -3; %% JJA
  driver.iAllorSeasonal = -4; %% SON
  driver.iAllorSeasonal = +1; %% default
  if iDoAnomalyOrRates == +1
    driver.iAllorSeasonal = +1; %% default
  end

  %%%%%%%%%% ANOM or RATES %%%%%%%%%%
  if iDoAnomalyOrRates == +1
    driver.i16daytimestep = JOB;                               %% ORIG when doing only ONE anomaly time series
    driver.i16daytimestep = mapAnomData_to_processor(iInd,2);  %% NEW when doing about 10 anomaly time series
    driver.anomalylatbin  = mapAnomData_to_processor(iInd,1);  %% NEW when doing about 10 anomaly time series
  elseif iDoAnomalyOrRates == -1  
    driver.i16daytimestep = -1;   %% for the rates, not anomalies, RUN BY HAND BY UN-COMMENTING THIS LINE and 
                                  %% on top JOB = 1000, in change_important_topts_settings.m also set topts.set_tracegas = -1;
  end

  if iDoAnomalyOrRates < 0
    fprintf(1,'latbin = %3i lonbin = %3i   gridpoint = %4i i16daytimestep = %4i \n',JOB,(iInd-iInd0+1),iInd,driver.i16daytimestep);
  elseif iDoAnomalyOrRates > 0
    fprintf(1,'JOB = %5i gridpoint = %5i i16daytimestep = %5i i16dayLocalBin = %5i \n',JOB,iInd,driver.i16daytimestep,driver.anomalylatbin);
  end

  ix = iInd;
  if iDoAnomalyOrRates < 0
    driver.iLon = iInd-iOffset;
    driver.iLat = JOB;
  else
    driver.iLon = -1;
    driver.iLat = driver.anomalylatbin;
  end

end

