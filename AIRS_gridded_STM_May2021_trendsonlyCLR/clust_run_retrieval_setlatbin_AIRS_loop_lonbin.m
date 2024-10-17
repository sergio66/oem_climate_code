%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS_loop_anomaly.m
%---------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FOR CALCS (RETRIEVAL COMPARISONS TO ERA5 TRENDS)
%%  moo = load('Output_CAL/Quantile16/test4608.mat');
%%                   debug: 0
%%               debug_dir: '../Debug'
%%                   iibin: 4608
%%          i16daytimestep: -1
%%       iDebugRatesUseNWP: -1
%%                    NorD: 1
%%               iQuantile: 16   ---------->>>>>>
%%                    iLon: 72
%%                    iLat: 64
%%                     oem: [1x1 struct]
%%                     lls: [1x1 struct]
%%             outfilename: 'Output_CAL/Quantile16/test4608.mat'
%%                   topts: [1x1 struct]
%%                 rateset: [1x1 struct]
%%              jac_latbin: 64
%%      jac_indexINSIDEbin: 72
%%                jacobian: [1x1 struct]
%%                 qrenorm: [2.200000000000000e+00 1 5 1 1 1.000000000000000e-01 1.000000000000000e-02 1.000000000000000e-02 1.000000000000000e-02 1.000000000000000e-02 ... ]
%%  
%%  >> moo.topts
%%                 iSergioCO2: -1
%%      set_era5_cmip6_airsL3: -1
%%               UMBCvsERAjac: -1
%%              resetnorm2one: -1
%%                    dataset: 9    -------------->>>>>
%%                    co2lays: 1
%%                    numchan: 2645
%%                 chan_LW_SW: 0
%%                  descORasc: 1
%%               iFixTG_NoFit: -1
%%               iFixTz_NoFit: -1
%%               iFixO3_NoFit: -1
%%                offsetrates: -1
%%               set_tracegas: 1
%%                      iXJac: 0
%%                    ocb_set: 1   -------------->>>>>
%%          iDoStrowFiniteJac: 3
%%            iNlays_retrieve: 50
%%                     iChSet: 5
%%            obs_corr_matrix: -1
%%        tie_sst_lowestlayer: -1
%%                    invtype: 1
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTMISC/
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/LOADMIE

system_slurm_stats

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

%% this is the timestep : 1: 365 (coincidence : there are 365 days/year and
%% I did 16 day averages .... so 365/16 steps per year ... and 2002-2018 is
%% 16 years so total of 365/16 * 16 = 365 steps

if exist('JOBM')
  JOB = JOBM;
else
  JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 64 for the 64 latbins, 1:120 for anomalies
  if length(JOB) == 0
    JOB = 32;
    %JOB = 33;
    %JOB = 34;
    JOB = 29;
    JOB = 2;
    JOB = 45;
    JOB = 1;
  end
end

% JOB = 33
% JOB = 47
% JOB = 37
% JOB = 7
% JOB = 39
% JOB = 12
% JOB = 16
% JOB  = 15

% JOB = 70 %% for anomalies, > 64)

JOB0 = JOB;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iDebug = 2259; %% T
iDebug = 3439; %% NML
iDebug = 2268; %% T
iDebug = 3276; %% NML
iDebug = 0180;  %% SP terrible wiggles in lower atm
iDebug = 0108;  %% SP works terible

iDebug = 1233; %% SML
iDebug = 3233; %% NML
iDebug = 2233; %% T
iDebug = 2264; %% T
iDebug = 4483; %% NP

iDebug = -1;

%iDebug = 36;     %% check to make sure emiss trend comes in (over land)
%iDebug = 2000
%iDebug = 1598
%iDebug = 10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% see ../FIN_NWP_TRENDS/ : driver_compare_trends_Day_vs_Night.m --> prep_colWV_T_WV_trends_Day_vs_Night --> get_umbc_day_night_name.m
umbc_night_file  = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jac>> tartwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_removeemisstrends.mat';
oldresults = load(umbc_night_file);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
moo0 = load('SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin22/sarta_spectral_trends_const_tracegas_latbin22_2002_09_2022_08.mat');
moo1 = load('SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/all_4608_desc_2002_09_2022_08.mat');
dir0 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/';
fMERRA2 = load([dir0 'MERRA2_SARTA_SPECTRAL_RATES/all_4608_2002_09_2022_08.mat']);
umbc_night_file  = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjacV2.mat';

fUMBC_n = load(umbc_night_file,'rates','fits','resultsWV');
obs_n = fUMBC_n.rates;
fUMBC_n = fUMBC_n.fits;
plot(f,moo0.thesave.xtrend(:,20),'b',f,moo1.trend(:,1532),'r',f,fMERRA2.trend(:,1532),'g',f,obs_n(:,1532),'k'); plotaxis2; xlim([640 1630])  %% 20 is the lonbin, 22 is the latbin ... JOBJOBJOB = 1532 = (22-1)*72+20
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

driver.iTrendOrAnomaly = +1;   %% trends
driver.iTrendOrAnomaly = -1;   %% anomalies

if driver.iTrendOrAnomaly > 0
  %% JPL 2021 Science Team Meeting used dataset=4,quantile=16 and Princeton PCTS
  ia_OorC_DataSet_Quantile = [+0 04 16 -9999]; %% ocb_set = 0 : obs fit, dataset = 4, iQuantile = 16    19 year rates, AIRS obs Q(09.99), JPL Aprl 2022 meeting        04/23/2022 commit 30d2e554a97b34b0923ad58346d183a3c10d6bcb
  ia_OorC_DataSet_Quantile = [+0 05 50 -9999]; %% ocb_set = 0 : obs fit, dataset = 5, iQuantile = 50    12 year rates, AIRS obs Q(09.99), Princeton Aug 2022 meeting   09/04/2022 commit 0cb7d1fc6ca2485864b625b0590cbdbb7894e5ac
  %%
  ia_OorC_DataSet_Quantile = [+0 07 16 -9999]; %% ocb_set = 0 : obs fit, dataset = 7, iQuantile = 16    20 year rates, AIRS obs Q16
  ia_OorC_DataSet_Quantile = [+0 09 05 -9999]; %% ocb_set = 0 : obs fit, dataset = 9, iQuantile = 05    20 year rates, AIRS obs Q(0.97-->1)

  ia_OorC_DataSet_Quantile = [+1 09 16  2]; %% ocb_set = 1 : MERRA2     cal fit, dataset = 9,  iQuantile = 16  20 year rates
  ia_OorC_DataSet_Quantile = [+1 09 16  3]; %% ocb_set = 1 : AIRSL3     cal fit, dataset = 9,  iQuantile = 16  20 year rates, problems at eg latbin 1, lonbin 15-65
  ia_OorC_DataSet_Quantile = [+1 09 16 -3]; %% ocb_set = 1 : CLIMCAPSL3 cal fit, dataset = 9,  iQuantile = 16  20 year rates
  ia_OorC_DataSet_Quantile = [+1 09 16  5]; %% ocb_set = 1 : ERA5       cal fit, dataset = 9,  iQuantile = 16  20 year rates
  
  ia_OorC_DataSet_Quantile = [+0 14 03 -9999]; %% ocb_set = 0 : AIRSL1C obs fit, dataset = 14, iQuantile = 03    4 year rates, 2018/09-2022/08 AIRS obs Q(0.90-->1)
  
  %%%%% for trends paper START
  ia_OorC_DataSet_Quantile = [+1 09 16  5   ]; %% ocb_set = 1 : ERA5 cal fit,    dataset = 9,  iQuantile = 16   20 year rates
  
  ia_OorC_DataSet_Quantile = [+0 10 03 -9999]; %% ocb_set = 0 : AIRSL1C obs fit, dataset = 10, iQuantile = 03    05 year rates, 2002/09-2007/08 AIRS obs Q(0.90-->1)
  ia_OorC_DataSet_Quantile = [+0 12 03 -9999]; %% ocb_set = 0 : AIRSL1C obs fit, dataset = 12, iQuantile = 03    15 year rates, 2002/09-2012/08 AIRS obs Q(0.90-->1)
  ia_OorC_DataSet_Quantile = [+0 11 03 -9999]; %% ocb_set = 0 : AIRSL1C obs fit, dataset = 11, iQuantile = 03    10 year rates, 2002/09-2017/08 AIRS obs Q(0.90-->1)
  ia_OorC_DataSet_Quantile = [+0 09 03 -9999]; %% ocb_set = 0 : AIRSL1C obs fit, dataset = 09, iQuantile = 03    20 year rates, 2002/09-2022/08 AIRS obs Q(0.90-->1)
  ia_OorC_DataSet_Quantile = [+0 09 03 -9999]; %% ocb_set = 0 : AIRSL1C obs fit, dataset = 09, iQuantile = 03    20 year rates, 2002/09-2022/08 AIRS obs Q(0.90-->1)
  %%%%% for trends paper STOP 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  ia_OorC_DataSet_Quantile = [+2 30 01 -9999]; %% ocb_set = 0 : AMSU obs fit,    dataset = 09, iQuantile = 01    20 year anomalies, 2002/09-2022/08 AMSU obs Q(0.50-->1) -- technically this is "allsky average' but should be clear
  ia_OorC_DataSet_Quantile = [+0 30 01 -9999]; %% ocb_set = 0 : AMSU obs fit,    dataset = 09, iQuantile = 01    20 year rates,     2002/09-2022/08 AMSU obs Q(0.50-->1) -- technically this is "allsky average' but should be clear
  
  %%%%%%%%%%%%%%%%%%%%%%%%%

  ia_OorC_DataSet_Quantile = [+0 16 03 -9999]; %% ocb_set = 0 : AIRSL1C obs fit, dataset = 16,  iQuantile = 03   04 year rates, 2020/07-2024/06

elseif driver.iTrendOrAnomaly < 0
  set_anomaly_info
end

%%%%%%%%%%%%%%%%%%%%%%%%%

if ia_OorC_DataSet_Quantile(1) == 2
  %get_anomaly_processors  %% this sets mapAnomData_to_processor
  [jobjunk,mapAnomData_to_processor] = get_anomaly_processors(ia_OorC_DataSet_Quantile,iNumAnomTimeSteps,iNumAnomTiles,iNumAnomJobsPerProc,anomalydatafile,JOB0);
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

% <<< got rid of this here as it is reintroduced as  -(driver.iTrendOrAnomaly) !!!! >>>
% <<<<<<<     but I re-use it in change_important_topts_settings.m      >>>>>>>>>>
% iDoAnomalyOrRates = +2;  %% bootstrap anomaly?? see change_important_topts_settings.m
% iDoAnomalyOrRates = +1;  %% do the anomalies
% iDoAnomalyOrRates = -1;  %% do the trends/rates, default
% iDoAnomalyOrRates = -(driver.iTrendOrAnomaly);
% if ia_OorC_DataSet_Quantile(1) == 2
%   iDoAnomalyOrRates = +1;  %% do the anomalies
% end

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
if driver.iTrendOrAnomaly == +1
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
  mapperAnom2Processor = find(mapAnomData_to_processor(:,3) == JOB);    %% from get_anomaly_processors
  iInd0 = mapperAnom2Processor(1);  iIndE = mapperAnom2Processor(end);  %% from get_anomaly_processors
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
  disp(' ')

%------------------------------------------------------------------------
%% <<<<<<<    no real need to touch any of this  >>>>>>>>
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called

  driver.ia_OorC_DataSet_Quantile = ia_OorC_DataSet_Quantile;

  driver.iAllorSeasonal = -1; %% DJF
  driver.iAllorSeasonal = -2; %% MAM
  driver.iAllorSeasonal = -3; %% JJA
  driver.iAllorSeasonal = -4; %% SON
  driver.iAllorSeasonal = +1; %% default
  if driver.iTrendOrAnomaly == +1
    driver.iAllorSeasonal = +1; %% default
  end

  %%%%%%%%%% ANOM or RATES %%%%%%%%%%
  if driver.iTrendOrAnomaly == -1
    %driver.i16daytimestep = JOB;                              %% ORIG when doing only ONE anomaly time series
    driver.i16daytimestep = mapAnomData_to_processor(iInd,2);  %% NEW when doing about 10 anomaly time series, from get_anomaly_processors
    if (JOB <= iNumAnomTimeSteps)
      driver.anomalyinfo.global = +1;  %% this is global, over 64 latbins
    elseif (JOB > iNumAnomTimeSteps & JOB <= 2*iNumAnomTimeSteps)
      driver.anomalyinfo.global = +2;  %% this is tropical, over 64 latbins
    else
      driver.anomalyinfo.global = -1;  %% this is averaged over few latbins
    end 
    driver.anomalyinfo.timestep16day = driver.i16daytimestep;

    driver.anomalyinfo.fatprofilebin  = mapAnomData_to_processor(iInd,1);  %% NEW when doing about 10 anomaly time series, from get_anomaly_processors (1) GLOBAL (2) TROICAL (3-30) FAT LATBINS

    junkx = load(anomalydatafile,'usethese');
    driver.anomalyinfo.latbin = ceil(nanmean(junkx.usethese{mapAnomData_to_processor(iInd,1)}));
    driver.anomalyinfo.lonbin  = 36;

    driver.anomalyinfo.i4608eqv = (driver.anomalyinfo.latbin-1)*72 + driver.anomalyinfo.lonbin; 
    driver.anomalyinfo.timesteps = iNumAnomTimeSteps;
    driver.anomalyinfo.numtiles  = iNumAnomTiles;
    driver.anomalyinfo.location  = iInd;
    driver.iibin = (driver.anomalyinfo.latbin - 1)*72 + driver.anomalyinfo.lonbin; 

  elseif driver.iTrendOrAnomaly == +1  
    driver.i16daytimestep = -1;   %% for the rates, not anomalies, RUN BY HAND BY UN-COMMENTING THIS LINE and 
                                  %% on top JOB = 1000, in change_important_topts_settings.m also set topts.set_tracegas = -1;
    driver.iibin     = iInd;
  end

  if driver.iTrendOrAnomaly > 0
    fprintf(1,'latbin = %3i lonbin = %3i   saved trend lon/lat point = %4i of 4608 i16daytimestep = %4i \n',JOB,(iInd-iInd0+1),iInd,driver.i16daytimestep);
  elseif driver.iTrendOrAnomaly < 0
    fprintf(1,'JOB = %5i gridpoint = %5i i16daytimestep = %5i latbin/lonbin = %2i/%2i\n',JOB,iInd,driver.i16daytimestep,driver.anomalyinfo.latbin,driver.anomalyinfo.lonbin);
  end

  ix = iInd;
  if driver.iTrendOrAnomaly > 0
    driver.iLon = iInd-iOffset;
    driver.iLat = JOB;
  else
    %% see set_driver_jacfile.m, elseif driver.i16daytimestep > 0
    do_XX_YY_from_X_Y

    driver.anomalyinfo.datafile = anomalydatafile;

    if strfind(driver.anomalyinfo.datafile,'_tile_')
      junk = load(driver.anomalyinfo.datafile,'LatBin','LonBin');
      driver.iLon = junk.LonBin;
      driver.iLat = junk.LatBin;
    else  %% ~strfind(driver.anomalyinfo.datafile,'_tile_')

      %junk = load(driver.anomalyinfo.datafile,'usethese');
      %junk = junk.usethese{driver.anomalyinfo.latbin};
      %%YYmean = nanmean(YY(junk));
      %%junk = find(rlat65 >= YYmean,1) - 1;
      %junk = round(nanmean(junk));

      junk1 = jobjunk.latnumber(iInd-iXX1+1);
      junk = load(driver.anomalyinfo.datafile,'usethese');
      junk = ceil(nanmean(junk.usethese{junk1}));

      %{
      newLatGrid = [-90 -75 -60 [-55:5:+55] +60 +75 +90];;      from stand_alone_make_globalavg_and_N_average_anomalies.m
      so if eg get_anomaly_processors.m  says we are using input_spectrum_number : locallatbin localtimestep ---> procnumber = 8245 :   17  245 ----->   33 clustjob = 033
         locallatin = 17 ---> latitude = 10
         junk == 36 ==> rlat = meanvaluebin(rlat65)    rlat(36) = 9.625   YAYAYAYAYAYAYAYA
      %}

      driver.iLon = 36;
      driver.iLat = junk;
    end
  end

  driver.oem.dofit = true;
  driver.lls.dofit = false;
  driver.oem.nloop = 2;
  driver.oem.nloop = 3;
  driver.oem.nloop = 1;
  driver.oem.doplots = false;

  %%%%%%%%%% ANOM or RATES %%%%%%%%%%
  driver.iDebugRatesUseNWP = 32; %% use AIRS L3 constructed spectral trends from SARTA
  driver.iDebugRatesUseNWP = 62; %% use CMIP6   constructed spectral trends from SARTA
  driver.iDebugRatesUseNWP = 52; %% use ERA     constructed spectral trends from SARTA
  driver.iDebugRatesUseNWP = -1; %% use AIRS observed spectral trends >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%---------------------------------------------------------------------------
  %% need to run this before inputting topts
  change_important_topts_settings  % Override many settings and add covariance matrix
%------------------------------------------------------------------------

  % [150 = T(z)+ST 60 = WV(z) 100 = O3(z)]
  % [210 = surf temp, T(z)        and lowest WV(z)]
  % [214 = surf temp, lowest T(z) and lowest WV(z)]

  topts.iaSequential = [150 60 100 -1];            %% sequential, like SingleFootprint
  topts.iaSequential = [-1 150 60 100 -1];         %% sequential, like SingleFootprint
  topts.iaSequential = [-1 150 60 100 150 60];     %% sequential, like SingleFootprint
  topts.iaSequential = [210 150 60 100 150 60];    %% sequential, like SingleFootprint
  topts.iaSequential = [214 150 60 100 150 60];    %% sequential, like SingleFootprint
  topts.iaSequential = [150 60 100 150 60];        %% sequential, like SingleFootprint
  topts.iaSequential = [150 60];                   %% sequential, like SingleFootprint
  topts.iaSequential = [150 60 100 150 60];        %% sequential, like SingleFootprint
  topts.iaSequential = [214 150 60 100 ];          %% sequential, like SingleFootprint
  topts.iaSequential = [60 150 100];               %% sequential, like SingleFootprint

  topts.iaSequential = [214 150 60 100 150 60];    %% sequential, like SingleFootprint, decent column water results at tropics!, not so good polar WV results
  topts.iaSequential = -1;                         %% default one gulp, gives good results at poles but bad column water results at tropics!
  topts.iaSequential = [210 150 60 100 150 60];    %% sequential, like SingleFootprint, but now use 210 !!!! gives slightly better results than 214, decent column water results at tropics!, not so good polar WV results

  %% see driver_put_together_QuantileChoose_trends.m see driver_put_together_QuantileChoose_trends.m see driver_put_together_QuantileChoose_trends.m
  %% see driver_put_together_QuantileChoose_trends.m see driver_put_together_QuantileChoose_trends.m see driver_put_together_QuantileChoose_trends.m
  %% see driver_put_together_QuantileChoose_trends.m see driver_put_together_QuantileChoose_trends.m see driver_put_together_QuantileChoose_trends.m

  % quants = [0 0.01 0.02 0.03 0.04 0.05 0.10 0.25 0.50 0.75 0.9 0.95 0.96 0.97 0.98 0.99 1.00];
  % topts.dataset   = -1;   %% (-1) AIRS 18 year quantile dataset, Sergio Aug 2021   2002/09-2020/08 FULL 18 years
  % topts.dataset   = +1;   %% (+1) AIRS 18 year quantile dataset, Strow  March 2021 2002/09-2020/08 FULL 18 years
  % topts.dataset   = +2;   %% (+2) AIRS 19 year quantile dataset, Sergio Aug 2021   2002/09-2021/07 PARTIAL 19 years
  % topts.dataset   = +3;   %% (+3) AIRS 19 year extreme  dataset, Sergio Aug 2021   2002/09-2021/07 PARTIAL 19 years, EXTREME
  % topts.dataset   = -3;   %% (-3) AIRS 19 year mean     dataset, Sergio Aug 2021   2002/09-2020/08 AUTOMATIC USES Q00, MEAN
  % topts.dataset   = +4;   %% (+4) AIRS 19 year quantile dataset, Sergio Aug 2021   2002/09-2021/08 FULL 19 years ************************ JPL April 2021 Sounder Science Meeting
  % topts.dataset   = +5;   %% (+5) AIRS 12 year quantile dataset, Sergio Aug 2022   2002/09-2014/08 FULL 12 years
  % topts.dataset   = +6;   %% (+6) AIRS = CRIS NSR 07 year quantile dataset,        2012/05-2019/04 FULL 07 years
  % topts.dataset   = +7;   %% (+7) AIRS 20 year quantile dataset, Sergio Sep 2022   2002/09-2022/08 FULL 20 years ************************
  % topts.dataset   = +8;   %% (+8) AIRS = OCO2  07 year quantile dataset            2015/01-2021/12 OCO2 FULL 07 years

  % % quants = [0.50 0.80 0.90 0.95 0.97 1.00];
  % topts.dataset   = +09;   %% (+ 9) AIRS 20 year quantile dataset, Sergio Oct 2022   2002/09-2022/08 FULL 20 years, new way of doing quantile iQAX = 3  ************************
  % topts.dataset   = +10;   %% (+10) AIRS 05 year quantile dataset, Sergio May 2023   2002/09-2007/08 FULL 05 years, new way of doing quantile iQAX = 3  ************************
  % topts.dataset   = +11;   %% (+11) AIRS 10 year quantile dataset, Sergio May 2023   2002/09-2012/08 FULL 10 years, new way of doing quantile iQAX = 3  ************************
  % topts.dataset   = +12;   %% (+12) AIRS 15 year quantile dataset, Sergio May 2023   2002/09-2017/08 FULL 15 years, new way of doing quantile iQAX = 3  ************************

  % topts.dataset   = +13;   %% (+13) AIRS  8 year quantile dataset, Sergio May 2023    2002/09-2019/08 FULL 08 years, SW drifting, new way of doing quantile iQAX = 3  ************************
  % topts.dataset   = +14;   %% (+14) AIRS  4 year quantile dataset, Sergio May 2023    2018/09-2022/08 FULL 04 years, new way of doing quantile iQAX = 3  ************************
  % topts.dataset   = +15;   %% (+15) AIRS 14 year quantile dataset, Sergio May 2023    2008/01-2022/12 FULL 14 years, overlap with IASI for Sarah/Cathy, new way of doing quantile iQAX = 3  ************************

  % topts.dataset   = +16;   %% (+16) AIRS  4 year quantile dataset, Sergio May 2023    2020/07-2024/06 FULL 04 years, HOT HOT HOT new way of doing quantile iQAX = 3  ************************
  % topts.dataset   = +17;   %% (+17) AIRS 14 year quantile dataset, Sergio May 2023    2002/09-2024/06 FULL 22 years, new way of doing quantile iQAX = 3  ************************

  % topts.dataset   = +30;   %% (+30) AMSU 20 year average dataset, Sergio May 2024    2002/09-2022/08 FULL 20 years, one one quantile (avg), given by Stephen Leroy  ************************

  topts.numchan = 2645;                              %% 2645 AIRS channels
  topts.model   = ia_OorC_DataSet_Quantile(4);
  topts.dataset = ia_OorC_DataSet_Quantile(2);
  iQuantile     = ia_OorC_DataSet_Quantile(3);

  if topts.dataset == 30
    topts.iaSequential = -1;                         %% default one gulp, gives good results at poles but bad column water results at tropics!    
    topts.numchan = 13;                              %% 13 AMSU channels
    topts.iNlays_retrieve = 10;                      %% default, 10 AIRS lays thick since so few AMSU channels
  end

  if ia_OorC_DataSet_Quantile(1) == 2
    %% do the anomalies fast!
    topts.iaSequential = -1;                                          %% default one gulp, gives good results at poles but bad column water results at tropics!
    driver.actualanomalyimeStep = jobjunk.localtimestep(iInd-iXX1+1); %% this is used by oem_pkg/rodgers.m --> oem_pkg/common_rodgers_initializations1.m to remove CO2, N20, CH4 etc
    junk = [iXX1 idX  iXX2 iInd  driver.actualanomalyimeStep];
    fprintf(1,'     anomaly --- doing iInd = %5i : %2i : %5i; current iInd, iTimestep = %5i %5i \n',junk);
  end

  %%%%%%%%%%

  driver.removeEmisTrend = 0;  %% ignore      changing (LAND) emiss
  if ia_OorC_DataSet_Quantile(1) == 0 & ia_OorC_DataSet_Quantile(2) ~= 14
    %% ia_OorC_DataSet_Quantile(2) == 14 is 2018-2022 = 4 year dataset sp emiss can't have changed
    %% ia_OorC_DataSet_Quantile(2) == [everything else] = 14+ ears so changing land emiss
    driver.removeEmisTrend = +1; %% account for changing (LAND) emiss
  end
  driver.removeEmisTrend = 0;  %% ignore      changing (LAND) emiss default

  driver.NorD = -1; %% day, asc
  driver.NorD = +1; %% night, desc  DEFAULT
  driver.iQuantile = iQuantile;

  disp('see driver_put_together_QuantileChoose_trends.m  driver_put_together_QuantileChoose_trends.m  driver_put_together_QuantileChoose_trends.m')
  disp('see driver_put_together_QuantileChoose_trends.m  driver_put_together_QuantileChoose_trends.m  driver_put_together_QuantileChoose_trends.m')
  disp('see driver_put_together_QuantileChoose_trends.m  driver_put_together_QuantileChoose_trends.m  driver_put_together_QuantileChoose_trends.m')
  if abs(topts.dataset) == +1
    driver.iNumYears = 18;
  elseif abs(topts.dataset) >= 2 & abs(topts.dataset) <= 4
    driver.iNumYears = 18;
  elseif topts.dataset == 5
    driver.iNumYears = 12;
  elseif topts.dataset == 6
    driver.iNumYears = 07; %% oops 2012-2019
  elseif topts.dataset == 7
    driver.iNumYears = 20;
  elseif topts.dataset == 8
    driver.iNumYears = 07; %% oops 2015-2021
  elseif topts.dataset == 9
    driver.iNumYears = 20;
    disp('STANDARD 20 years, used for AIRS STM and trends paper')
  elseif topts.dataset == 10
    driver.iNumYears = 05;
    disp('STANDARD 05 years')
  elseif topts.dataset == 11
    driver.iNumYears = 10;
    disp('STANDARD 10 years')
  elseif topts.dataset == 12
    driver.iNumYears = 15;
    disp('STANDARD 15 years')
  elseif topts.dataset == 13
    driver.iNumYears = 20;
    disp('20 years, Q01-03 NEW DEFN ???')
    driver.iNumYears = 08;
    disp('08 years, Q01-03 NEW DEFN ???')
  elseif topts.dataset == 14
    driver.iNumYears = -04;  
    disp(' <<< 2018-2022 : warning ... the CO2 trends are from 2002 onwards, but this dataset is from 2018 >>>')
  elseif topts.dataset == 15
    driver.iNumYears = 14;  
    disp(' <<< 2008-2022 IASI overlap arning ... the CO2 trends are from 2002 onwards, but this dataset is from 2008 >>>')
  elseif topts.dataset == 16
    driver.iNumYears = 4;  
    disp(' <<< 2020-2024 AIRS hot 4 years, ... the CO2 trends are from 2002 onwards >>>')
  elseif topts.dataset == 17
    driver.iNumYears = 22;  
    disp(' <<< 2002-2024 AIRS for 22 years, ... the CO2 trends are from 2002 onwards >>>')
  elseif topts.dataset == 30
    driver.iNumYears = 20;  
    disp(' <<< AMSU 20 years allsky from Stephen Leroy >>> ')
  else
    fprintf(1,'topts.dataset = %2i \n',topts.dataset)
    error('unknown topts.dataset')
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% WARNING : set_CO2_CH4_N2O_ESRL.m has an internal setting iVers = 2 where both T and WV are adjusted in lower atmosphere, but that can be annulled by using topts.TfacAdjAtmosphericAmplification (default 0.5, set it to a smaller value or 0)
  iAdjLowerAtmWVfrac = 0;                             %% << WARNING <<< this does not  set WV in lower part of atmos >>> , depending on dBT1231/dt by using iAdjLoweAtmWVfrac >>
  iAdjLowerAtmWVfrac = 10;                            %% WARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!, this is for dRH/dt ~ 1 (see bk46)
  iAdjLowerAtmWVfrac = 05;                            %% WARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!, this is for dRH/dt ~ 0.5 (see bk46)
  iAdjLowerAtmWVfrac = 1;                             %% WARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!, this is for dRH/dt = 0 (see bk46), but I think tooooo much
  iAdjLowerAtmWVfrac = 0.25;                          %% WAARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!, this is for dRH/dt = 0 (see bk46), with an adjustment
  %% Nov 16, 2023 since the mmw trends in tropics are too large
  iAdjLowerAtmWVfrac = 0.25;                          %% WAARNING this also sets WV in lower part of atmos NEW
  iAdjLowerAtmWVfrac = 0.125;                         %% WAARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!, this is for dRH/dt = 0 (see bk46), with an adjustment, till Nov 2023
  iAdjLowerAtmWVfrac = 0.0625;                        %% WAARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!, this is for dRH/dt = 0 (see bk46), with an adjustment
  iAdjLowerAtmWVfrac = 1;                             %% WARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!, this is for dRH/dt = 0 (see bk46), but I think tooooo much
  topts.iAdjLowerAtmWVfrac = iAdjLowerAtmWVfrac;      %% WARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!

  topts.TfacAdjAtmosphericAmplification = 0.0;        %% this is an additional adjutment factor for a-priori WV believe it or not, TEST TO SEE WHAT HAPPENS
  topts.TfacAdjAtmosphericAmplification = 0.5;        %% this is an additional adjutment factor for a-priori WV believe it or not, default till Nov 2023
  topts.TfacAdjAtmosphericAmplification = 1.0;        %% new Dec 2023
  if length(topts.iaSequential) > 1
    topts.TfacAdjAtmosphericAmplification = 0.1;        %% this is an additional adjutment factor for a-priori WV believe it or not
    topts.TfacAdjAtmosphericAmplification = 0.25;       %% this is an additional adjutment factor for a-priori WV believe it or not
    topts.TfacAdjAtmosphericAmplification = 0.75;       %% this is an additional adjutment factor for a-priori WV believe it or not  
    topts.TfacAdjAtmosphericAmplification = 0.75;       %% this is an additional adjutment factor for a-priori WV believe it or not  
    topts.TfacAdjAtmosphericAmplification = 1.00;       %% this is an additional adjutment factor for a-priori WV believe it or not , pretty good BUT JUST LEAVE IT AS 1.0
  end

  %% do this if you want NO WV adj anywaehre, for LLS SARTA report for JPL
  %% else comment it out if you want to compare to ERA5 trends
  if driver.iNumYears < 0
    iAdjLowerAtmWVfrac = 0;
    topts.iAdjLowerAtmWVfrac = iAdjLowerAtmWVfrac;
    topts.TfacAdjAtmosphericAmplification = 0.00;
  end

  if driver.iTrendOrAnomaly > 0
    if driver.iLat <= 31
      %% we seem to need double or triple this factor in the Southern Hemisphere compared to Northern (too much ocean?????)
  
      %% this is /asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2.mat
      %% this is /asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2A.mat
      %% so let it be triple at S. Polar (JOB = 1) and then ramp down smoothly to unity at equator (JOB = 32)
      %% so the slope = (3-1)/(32-1) = -2/31
      intercept = 3;
  
      %% this is /asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2B.mat
      %% this is /asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2C.mat
      %% so let it be double at S. Polar (JOB = 1) and then ramp down smoothly to unity at equator (JOB = 32)
      %% so the slope = (2-1)/(32-1) = -1/31
      intercept = 2;
  
      %% this is /asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2D.mat
      %% so let it be double at S. Polar (JOB = 1) and then ramp down smoothly to unity at equator (JOB = 32)
      %% so the slope = (1-1)/(32-1) = 0/31
      intercept = 1;
  
      slope = (1-intercept)/31;
      yjunk = slope * (driver.iLat - 1) + intercept;
      topts.TfacAdjAtmosphericAmplification = topts.TfacAdjAtmosphericAmplification * yjunk;
      fprintf(1,'JOB = %4i latbin = %2i mult = %8.6f topts.TfacAdjAtmosphericAmplification = %8.6f \n',JOB,driver.iLat,yjunk,topts.TfacAdjAtmosphericAmplification);
    end
  end
  
  topts.tie_sst_lowestlayer = -1;  %% DEFAULT
  topts.tie_sst_lowestlayer = +1;  %% testing dataset=4,iQuantil=16,ocb_set=0 (the JPL SOunder meeting Apr 2022, 04/23/2022 commit 30d2e554a97b34b0923ad58346d183a3c10d6bcb

  %topts.iFixTG_NoFit = +1; %% dump out first scalar = CO2 boy certainly messes up CO2 even more!!!!!

  topts.iSergioCO2 = +1; %% fit for CO2/CH4/N2O
  topts.iSergioCO2 = -1; %% assume ESRL CO2/CH4
  if topts.dataset == 8 
    topts.iSergioCO2 = +1; %% fit for CO2/CH4/N2O
  end

  topts.ocb_set = ia_OorC_DataSet_Quantile(1);

  % << set_apriori_ERA5_MERRA2_or_AIRSL3_MLS_geophysical.m >>
  %   settings.set_era5_cmip6_airsL3_WV_T_O3 == -1   : set all
  %   settings.set_era5_cmip6_airsL3_WV_T_O3 == +1   : set WV
  %   settings.set_era5_cmip6_airsL3_WV_T_O3 == +2   : set T/ST
  %   settings.set_era5_cmip6_airsL3_WV_T_O3 == +3   : set O3
  %   settings.set_era5_cmip6_airsL3_WV_T_O3 == +4   : set T
  %   settings.set_era5_cmip6_airsL3_WV_T_O3 == +5   : set ST
  %   settings.set_era5_cmip6_airsL3_WV_T_O3 == +10  : set lower WV
  %   settings.set_era5_cmip6_airsL3_WV_T_O3 == +100 : set lower WV/upper WV with MLS

  topts.set_era5_cmip6_airsL3_WV_T_O3 = +2;  %% use T+ST
  topts.set_era5_cmip6_airsL3_WV_T_O3 = +40; %% use T only in the lower trop to "start things" in the polar region, based on BT1231
  topts.set_era5_cmip6_airsL3_WV_T_O3 = +4;  %% use T only
  topts.set_era5_cmip6_airsL3_WV_T_O3 = +1;  %% use WV only
  topts.set_era5_cmip6_airsL3_WV_T_O3 = +10; %% use WV only in the lower trop to "start things" in the polar region, based on BT1231
  topts.set_era5_cmip6_airsL3_WV_T_O3 = +100 %% set lower WV a priori, and then also upper WV with MLS   >>>>>
  topts.set_era5_cmip6_airsL3_WV_T_O3 = -1;  %% use WV/T/ST/O3, DEFAULT

  topts.set_era5_cmip6_airsL3 = 6;           %% use AMIP6    a priori
  topts.set_era5_cmip6_airsL3 = 3;           %% use AIRSL3   a priori
  topts.set_era5_cmip6_airsL3 = -3;          %% use CLIMCAPS a priori
  topts.set_era5_cmip6_airsL3 = 2;           %% use MERRA2   a priori
  topts.set_era5_cmip6_airsL3 = 5;           %% use ERA5     a priori
  topts.set_era5_cmip6_airsL3 = 8;           %% use MLS      a priori >>>>
  topts.set_era5_cmip6_airsL3 = 0;           %% use 0        a priori, DEFAULT

  if (length(intersect(ia_OorC_DataSet_Quantile(2),[10 11 12])) == 1) & topts.set_era5_cmip6_airsL3 == 8
    error('MLS UT/LS WV a-priori only for 20 years!!!!')
  end
  
  if ia_OorC_DataSet_Quantile(1) == 1 & topts.set_era5_cmip6_airsL3 ~= 0
    %% this is synthetic rates, and starting with big a-priori .... so do not mess with lower atm WV!!!
    disp('ia_OorC_DataSet_Quantile(1) == 1 & topts.set_era5_cmip6_airsL3 ~= 0 ==> synthetic rates + start with NWP trends as xb (a priori) means .... do not mess with PBL WV trends!!!!')
    disp('  making sure topts.iAdjLowerAtmWVfrac = 0 ')
    topts.iAdjLowerAtmWVfrac = 0;
  end

  if topts.set_era5_cmip6_airsL3 == 8   
    if topts.set_era5_cmip6_airsL3_WV_T_O3 == +10
      topts.set_era5_cmip6_airsL3_WV_T_O3 = +100;
      disp('you want MLS plus you want to set lower WV .... reset topts.set_era5_cmip6_airsL3_WV_T_O3 = +100')
    end
  end
  if topts.iAdjLowerAtmWVfrac > 0 & topts.set_era5_cmip6_airsL3 ~= 0
    disp('You have topts.iAdjLowerAtmWVfrac > 0 so you want to set PBL ... but you also want to set ERA5 or MERRA2 or AIRS L3 a-priori NOPE NOPE NOPE')
    topts.iAdjLowerAtmWVfrac = 0;             %% forget about PBL
    %topts.set_era5_cmip6_airsL3 = 0;         %% use 0 a priori for every layer ... but can still allow MLS input at UT/LS and above 300 mb
  end

  topts.iNlays_retrieve = 20; %% default, 5 AIRS lays thick
  topts.iNlays_retrieve = 50; %%          2 AIRS lays thick
  if topts.dataset == 30
    topts.iNlays_retrieve = 10;                      %% default, 10 AIRS lays thick since so few AMSU channels
  end

  %% iLatX is used in : build_cov_matrices.m to make the cov_setA, cov_srtB ... cov_set = wgtA cov_setA + wgtB cov_setB
  %%                  : get_jac_fast.m to set the iVersQRenorm and then qrenorm
  iLatX = 11; %% midlats
  iLatX = -1; %% do not any blending!!!     %%%%%%%%%%%%%%%%%%%%%%%%% TESTING TO GET OBS COV MATRICES SAME AS SYNTHETIC ERA5 COV MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%    

  topts.iLatX = iLatX;
  %% WARNING, when savesmallFATfile or savebigFATfile is called, topts.resetnorm2one will depend on which is the last file read in (could be anything, depending on the darn cluster)
  %% WARNING, when savesmallFATfile or savebigFATfile is called, topts.resetnorm2one will depend on which is the last file read in (could be anything, depending on the darn cluster)
  %% WARNING, when savesmallFATfile or savebigFATfile is called, topts.resetnorm2one will depend on which is the last file read in (could be anything, depending on the darn cluster)
  if driver.iTrendOrAnomaly > 0
    if driver.iLat <= iLatX | driver.iLat >= 64 - iLatX
      topts.resetnorm2one = -1;   %% DEFAULT, use eg 0.01 for T, 2.2/400 for CO2 etc
    else
      topts.resetnorm2one = +1;   %% use 1 for every param
    end
    %topts.resetnorm2one = +1;   %% use 1 for every param, this is new March 2023 and duplicates the Feb 9, 2023 commit
  end

  iChSet = topts.iChSet;
  iChSet = 2; %% new chans
  iChSet = 4; %% new chans + Tonga (high alt) + more LW
  iChSet = 1; %% old chans (default)
  iChSet = 3; %% new chans, but no CFC11   STROW PRISTINE SET, AMT 2019, also used for JPL April 2021 Sounder Science Meeting
  %iChSet = 5; %% new chans, + CO2 laser lines (window region, low altitude T)
  %iChSet = 6; %% SW T(z) chans + 800 - 1600 cm- 1 lines (window region, low altitude T)
  if topts.dataset == 4
    iChSet = 3; %% new chans, but no CFC11   STROW PRISTINE SET, AMT 2019, also used for JPL April 2021 Sounder Science Meeting
  end

  topts.iChSet = iChSet;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% NorD > 0 ==> night
  if topts.ocb_set == 0 & driver.i16daytimestep > 0 & driver.NorD > 0 & topts.dataset ~= 3
    driver.outfilename = ['OutputAnomaly/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(iInd,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif topts.ocb_set == 1 & driver.i16daytimestep > 0 & driver.NorD > 0 & topts.dataset ~= 3
    driver.outfilename = ['OutputAnomaly_CAL/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(iInd,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif topts.ocb_set == 2 & driver.i16daytimestep > 0 & driver.NorD > 0 %% NIGHT TIME ANOMALIES, NEW
    outdir = ['OutputAnomaly/Quantile' num2str(driver.iQuantile,'%02d')];
    driver.outfilename = [outdir  '/test' int2str(iInd) '.mat'];

  elseif topts.ocb_set == 0 & driver.i16daytimestep < 0 & driver.NorD > 0 & topts.dataset ~= 3
    outdir = ['Output/Quantile' num2str(driver.iQuantile,'%02d')];
    driver.outfilename = [outdir  '/test' int2str(iInd) '.mat'];
  elseif topts.ocb_set == 1 & driver.i16daytimestep < 0 & driver.NorD > 0 & topts.dataset ~= 3
    outdir = ['Output_CAL/Quantile' num2str(driver.iQuantile,'%02d')];
    driver.outfilename = [outdir  '/test' int2str(iInd) '.mat'];
  elseif topts.ocb_set == 0 & driver.i16daytimestep < 0 & driver.NorD > 0 & topts.dataset == 3 %% EXTREME
    outdir = ['Output/Extreme/'];
    driver.outfilename = [outdir  '/test' int2str(iInd) '.mat'];

  %% NorD < 0 ==> day
  elseif topts.ocb_set == 0 & driver.i16daytimestep > 0 & driver.NorD < 0 & topts.dataset ~= 3
    driver.outfilename = ['OutputAnomaly_OBS_Day/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(iInd,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif topts.ocb_set == 1 & driver.i16daytimestep > 0 & driver.NorD < 0 & topts.dataset ~= 3
    driver.outfilename = ['OutputAnomaly_CAL_Day/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(iInd,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];

  elseif topts.ocb_set == 0 & driver.i16daytimestep < 0 & driver.NorD < 0 & topts.dataset ~= 3
    outdir = ['Output_Day/Quantile' num2str(driver.iQuantile,'%02d')];
    driver.outfilename = [outdir  '/test' int2str(iInd) '.mat'];
  elseif topts.ocb_set == 1 & driver.i16daytimestep < 0 & driver.NorD < 0 & topts.dataset ~= 3
    outdir = ['Output_Day_CAL/Quantile' num2str(driver.iQuantile,'%02d')];
    driver.outfilename = [outdir  '/test' int2str(iInd) '.mat'];
  elseif topts.ocb_set == 0 & driver.i16daytimestep < 0 & driver.NorD < 0 & topts.dataset == 3 %% EXTREME
    outdir = ['Output_Day/Extreme/'];
    driver.outfilename = [outdir  '/test' int2str(iInd) '.mat'];

  end
  if ~exist(outdir)
    mker = ['!mkdir -p ' outdir];
    eval(mker);
  end

  if exist(driver.outfilename) & iDebug > 0
    disp(' ')
    fprintf(1,'WARNING iDebug = %5i and %s already exists, so deleting and continuing \n',iDebug,driver.outfilename);
    fprintf(1,'WARNING iDebug = %5i and %s already exists, so deleting and continuing \n',iDebug,driver.outfilename);
    fprintf(1,'WARNING iDebug = %5i and %s already exists, so deleting and continuing \n',iDebug,driver.outfilename);
    fprintf(1,'WARNING iDebug = %5i and %s already exists, so deleting and continuing \n',iDebug,driver.outfilename);
    fprintf(1,'WARNING iDebug = %5i and %s already exists, so deleting and continuing \n',iDebug,driver.outfilename);
    disp(' ')
    pause(2)
    rmer = ['!/bin/rm ' driver.outfilename];
    eval(rmer)
  end

  oldstyle_removeCO2 = +1;
  if ~isfield(topts,'iFixTG_NoFit')
    oldstyle_removeCO2 = +1;
  elseif isfield(topts,'iFixTG_NoFit')
    if length(intersect(topts.iFixTG_NoFit,2)) >  0
      oldstyle_removeCO2 = -1;
    end
  end

  if ~exist(driver.outfilename)
    if oldstyle_removeCO2 > 0 
      if topts.iNlays_retrieve >= 97 & ~exist(driver.outfilename)
        [driver,aux,topts] = strow_override_defaults_latbins_AIRS(driver,topts);
      elseif topts.iNlays_retrieve < 97 & ~exist(driver.outfilename)
        [driver,aux,topts] = strow_override_defaults_latbins_AIRS_fewlays(driver,topts.iNlays_retrieve,topts);
      end
      %%[driver,aux,topts] = strow_override_defaults_latbins_AIRS_fewlays(driver,topts.iNlays_retrieve,topts);
    else
      if topts.iNlays_retrieve >= 97 & ~exist(driver.outfilename)
        [driver,aux,co2rate] = strow_override_defaults_latbins_AIRS_removeCO2(driver,topts);
      elseif topts.iNlays_retrieve < 97 & ~exist(driver.outfilename)
        [driver,aux,co2rate] = strow_override_defaults_latbins_AIRS_fewlays_removeCO2(driver,topts.iNlays_retrieve,topts);
      end
    end
  else
    fprintf(1,'A outdir,outname = %s already exists \n',driver.outfilename)    
  end

%---------------------------------------------------------------------------
  fprintf(1,'A outdir,outname = %s \n',driver.outfilename)
  % Do the retrieval

  if ~exist(driver.outfilename)     

    wvmoo = driver.jacobian.water_i;
    plot(aux.f(driver.jacobian.chanset),sum(aux.m_ts_jac(driver.jacobian.chanset,wvmoo(end-2:end)),2),'b',...
       aux.f(driver.jacobian.chanset),aux.m_ts_jac(driver.jacobian.chanset,6)*10,'k',...
       aux.f(driver.jacobian.chanset),driver.rateset.rates(driver.jacobian.chanset)*50,'r')
    plotaxis2; hl = legend('lowest WVjac','ST jac','rateset','location','best');
    pause(0.1);

    driver = retrieval(driver,aux);
    if oldstyle_removeCO2 < 0
      %% now adjust everything!!!!!
      call_adjustCO2
    end

    junk  = load('h2645structure.mat');
    f2645 = junk.h.vchan;
    if topts.dataset == 30
      f2645 = aux.f;
    end
    figure(3); plot(1:length(f2645),driver.rateset.rates,'b',1:length(f2645),driver.oem.fit,'r'); title('AIRS_gridded_Oct2020_trendsonly')
    figure(3); plot(f2645,driver.rateset.rates,'b',f2645,driver.oem.fit,'r'); title('AIRS_gridded_Oct2020_trendsonly')

  else
    disp(' ')
    fprintf(1,'%s after all that strow_override_defaults, already exists \n',driver.outfilename)
    fprintf(1,'%s after all that, already exists \n',driver.outfilename)
    fprintf(1,'%s after all that, already exists \n',driver.outfilename)
    driver.rateset.rates = zeros(topts.numchan,1);
    disp(' ')
  end

% driver

  if ~exist('oem')
    oem = driver.oem;
  end
  if exist('oem')
    if isfield(oem,'cov_set')
      show_unc  
    end
  end

%---------------------------------------------------------------------------
  % Save retrieval output from this loop

  if isfield(topts,'iFixWV_NoFit')
    %% hardest, since WV is at the beginning
    if topts.iFixWV_NoFit > 0
      %% put ozone and temp into correct (expected) spots
      driver.oem.finalrates(aux.orig_ozone_i) = driver.oem.finalrates(aux.orig_temp_i);
      driver.oem.finalsigs(aux.orig_ozone_i)  = driver.oem.finalsigs(aux.orig_temp_i);

      driver.oem.finalrates(aux.orig_ozone_i) = driver.oem.finalrates(aux.orig_temp_i);
      driver.oem.finalsigs(aux.orig_ozone_i)  = driver.oem.finalsigs(aux.orig_temp_i);

      %% put fixed/unchanging T anomaly
      driver.oem.finalrates(aux.orig_temp_i) = aux.FixTz_NoFit;
      driver.oem.finalsigs(aux.orig_temp_i)  = 0;

      driver.oem.finalrates(aux.orig_temp_i) = aux.FixTz_NoFit;
      driver.oem.finalsigs(aux.orig_temp_i)  = 0;

      driver.jacobian.water_i = driver.jacobian.ozone_i;
      driver.jacobian.temp_i  = driver.jacobian.ozone_i;
      driver.jacobian.ozone_i = driver.jacobian.ozone_i + driver.jacobian.numlays;

    end
  end

  if isfield(topts,'iFixTz_NoFit')
    %% harder, since T is in the middle
    if topts.iFixTz_NoFit > 0
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
    %% easiest, since O3 is tail end
    if topts.iFixO3_NoFit > 0

      %% put ozone into correct (expected) spot
      driver.oem.finalrates(aux.orig_ozone_i) = aux.FixO3_NoFit;
      driver.oem.finalsigs(aux.orig_ozone_i)  = 0;

      driver.jacobian.ozone_i = driver.jacobian.temp_i + driver.jacobian.numlays;
    end
  end

%%%%   if sum(abs(driver.rateset.rates)) > 0 & ~exist(driver.outfilename)  till Aug 2023
  if nansum(abs(driver.rateset.rates)) >= 0 & ~exist(driver.outfilename)
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
     fprintf('Scalar Retrievals from OEM latbin %3i gridpoint %4i \n',JOB,iInd)
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

if (driver.iLat-1)*72 + driver.iLon == iDebug

  %disp('Hit return for next latitude'); pause
  pause(0.1)

  iNumYears = driver.iNumYears;
  read_fileMean17years

  print_cloud_params(hMean17years,pMean17years,driver.iibin); 

  if ~isfield(pMean17years,'plays')
    pMean17years = make_rtp_plays(pMean17years);
  end
  pMean17yearsx = pMean17years;

  junk = pMean17years.plays(1:101,driver.iibin);
  for iii = 1 : driver.jacobian.numlays
    junk2 = driver.jacobian.wvjaclays_used{iii}-6;
    playsRET(iii) = mean(junk(junk2));
  end

  pMean17yearsx.stemp(iDebug) = pMean17years.stemp(iDebug) + driver.oem.finalrates(6);
  junk = driver.oem.finalrates(driver.jacobian.ozone_i);
    njunk = pMean17years.nlevs(driver.iibin);
    junk2 = interp1(log10(playsRET(1:length(junk))),junk,log10(pMean17years.plays(1:njunk,driver.iibin)),[],'extrap');
    bad = find(isnan(junk2) | isinf(junk2)); junk2(bad) = 0;    
    pMean17yearsx.gas_3(1:njunk,driver.iibin) = pMean17years.gas_3(1:njunk,driver.iibin) .* (1+junk2);

  junk = driver.oem.finalrates(driver.jacobian.temp_i);
    njunk = pMean17years.nlevs(driver.iibin);
    junk2 = interp1(log10(playsRET(1:length(junk))),junk,log10(pMean17years.plays(1:njunk,driver.iibin)),[],'extrap');
    bad = find(isnan(junk2) | isinf(junk2)); junk2(bad) = 0;    
    pMean17yearsx.ptemp(1:njunk,driver.iibin) = pMean17years.ptemp(1:njunk,driver.iibin) + junk2;

  junk = driver.oem.finalrates(driver.jacobian.water_i);
    njunk = pMean17years.nlevs(driver.iibin);
    junk2 = interp1(log10(playsRET(1:length(junk))),junk,log10(pMean17years.plays(1:njunk,driver.iibin)),[],'extrap');
    bad = find(isnan(junk2) | isinf(junk2)); junk2(bad) = 0;    
    pMean17yearsx.gas_1(1:njunk,driver.iibin) = pMean17years.gas_1(1:njunk,driver.iibin) .* (1+junk2);
    semilogy(junk,playsRET(1:length(junk)),'b.-',junk2,pMean17years.plays(1:njunk,driver.iibin)); set(gca,'ydir','reverse');
    ylim([10 1050]); axx = axis; line([axx(1) axx(2)],[pMean17years.spres(driver.iibin) pMean17years.spres(driver.iibin)]);

  mmw0 = mmwater_rtp(hMean17years,pMean17years);
  mmwX = mmwater_rtp(hMean17years,pMean17yearsx);
  fprintf(1,'d(percent mmw)/d(ST) = %8.6f percent mm/K \n',100*(mmwX(driver.iibin)-mmw0(driver.iibin))/mmw0(driver.iibin)/driver.oem.finalrates(6))

  if topts.iNlays_retrieve >= 97
    %disp('ret to plot retrievals'); pause
    pause(1)
    plot_retrieval_latbins
  else
    %disp('ret to plot retrievals'); pause
    pause(1)
    plot_retrieval_latbins_fewlays
  end

  figure(5); clf
  plot(f,driver.oem.fitXcomponents); plotaxis2; 
  if topts.dataset < 30 
    xlim([640 1640]); 
    xlabel('AIRS Wavenumber cm-1')
  else
    xlabel('AMSU Frequency GHz')
  end
  title('Components of fit'); hl = legend('trace gases','ST','WV(z)','T(z)','O3(z)','location','best','fontsize',10);
  % printarray([aux.pavg; driver.oem.xb(driver.jacobian.water_i)']','[p(water) xb(water)]')

  figure(6); clf
  if topts.dataset < 30 
    plot(f,driver.rateset.rates); 
    xlim([640 1640]); 
    xlabel('AIRS Wavenumber cm-1')
  else
    plot(f,driver.rateset.rates,'.-'); 
    xlabel('AMSU Frequency GHz')
  end
  plotaxis2;     ylabel('dBT/dt [K/yr]')

  error('nnyuk iDebug')
end


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
