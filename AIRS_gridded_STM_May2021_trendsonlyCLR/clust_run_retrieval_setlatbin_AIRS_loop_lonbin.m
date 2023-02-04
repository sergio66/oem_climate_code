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
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil

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

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 64 for the 64 latbins
% JOB = 33
% JOB = 47
% JOB = 37
% JOB = 7
% JOB = 39


iDebug = 3233; %% NML
iDebug = 2233; %% T
iDebug = 2264; %% T
iDebug = 2268; %% T
iDebug = 1233; %% SML
%iDebug = 180
%iDebug = 754
%iDebug = 1;
%iDebug = 1665;
%iDebug = 2376;

%iDebug = 4483; %% NP
%iDebug = 4233; %% NP
%iDebug = 0249;  %% SP
%iDebug = 0233;  %% SP
%iDebug = 0108;  %% SP works nice
%iDebug = 0180;  %% SP terrible wiggles in lower atm

iDebug = 2259; %% T
iDebug = 3439; %% NML
iDebug = 3276; %% NML
iDebug = 2268; %% T
iDebug = -1;

%% JPL 2021 Science Team Meeting used dataset=4,quantile=16
ia_OorC_DataSet_Quantile = [+0 04 16]; %% ocb_set = 0 : obs fit, dataset = 4, iQuantile = 16    19 year rates, AIRS obs Q(09.99), JPL Aprl 2022 meeting        04/23/2022 commit 30d2e554a97b34b0923ad58346d183a3c10d6bcb
ia_OorC_DataSet_Quantile = [+0 05 50]; %% ocb_set = 0 : obs fit, dataset = 5, iQuantile = 50    12 year rates, AIRS obs Q(09.99), Princeton Aug 2022 meeting   09/04/2022 commit 0cb7d1fc6ca2485864b625b0590cbdbb7894e5ac
ia_OorC_DataSet_Quantile = [+0 09 05]; %% ocb_set = 0 : obs fit, dataset = 9, iQuantile = 05    20 year rates, AIRS obs Q(09.97-->1)
ia_OorC_DataSet_Quantile = [+1 09 16]; %% ocb_set = 1 : cal fit, dataset = 9, iQuantile = 16    20 year rates, ERA5 synthetic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% MAIN CODE %%%%%%% MAIN CODE %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iDebug > 0
  plot(1:4608,floor(((1:4608)-1)/72)+1)   %% this maps tile number and job (in sets of 72)
  JOB = floor((iDebug-1)/72+1);
end

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

%for iInd = iIndE : -1 : iInd0
%for iInd = iInd0 + 21;
for iInd = iInd0 : iIndE

  %% so this is really iInd into 1:4608, 72 at a time
  disp(' ')
  fprintf(1,'latbin = %3i lonbin = %3i gridpoint = %4i \n',JOB,(iInd-iInd0+1),iInd);

%------------------------------------------------------------------------
%% <<<<<<<    no real need to touch any of this  >>>>>>>>
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called
  %% DO NOT PUT TOPTS HERE, put TOPTS LATER after change_important_topts_settings is called

  driver.ia_OorC_DataSet_Quantile = ia_OorC_DataSet_Quantile;
  driver.iibin     = iInd;

  %%%%%%%%%% ANOM or RATES %%%%%%%%%%
  if iDoAnomalyOrRates == +1
    driver.i16daytimestep = JOB;  %% this is when doing anomaly
  elseif iDoAnomalyOrRates == -1  
    driver.i16daytimestep = -1;   %% for the rates, not anomalies, RUN BY HAND BY UN-COMMENTING THIS LINE and 
                                  %% on top JOB = 1000, in change_important_topts_settings.m also set topts.set_tracegas = -1;
  end
  %%%%%%%%%% ANOM or RATES %%%%%%%%%%
  driver.iDebugRatesUseNWP = 32; %% use AIRS L3 constructed spectral trends from SARTA
  driver.iDebugRatesUseNWP = 62; %% use CMIP6   constructed spectral trends from SARTA
  driver.iDebugRatesUseNWP = 52; %% use ERA     constructed spectral trends from SARTA
  driver.iDebugRatesUseNWP = -1; %% use AIRS observed spectral trends >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  ix = iInd;
  driver.iLon = iInd-iOffset;
  driver.iLat = JOB;

  driver.oem.dofit = true;
  driver.lls.dofit = false;
  driver.oem.nloop = 2;
  driver.oem.nloop = 3;
  driver.oem.nloop = 1;
  driver.oem.doplots = false;

%---------------------------------------------------------------------------
  %% need to run this before inputting topts
  change_important_topts_settings  % Override many settings and add covariance matrix
%------------------------------------------------------------------------

  topts.iaSequential = [150 60 100 -1];            %% sequential, like SingleFootprint
  topts.iaSequential = [-1 150 60 100 -1];         %% sequential, like SingleFootprint
  topts.iaSequential = [-1 150 60 100 150 60];     %% sequential, like SingleFootprint
  topts.iaSequential = [210 150 60 100 150 60];    %% sequential, like SingleFootprint
  topts.iaSequential = [214 150 60 100 150 60];    %% sequential, like SingleFootprint
  topts.iaSequential = [150 60 100 150 60];        %% sequential, like SingleFootprint
  topts.iaSequential = [150];                      %% sequential, like SingleFootprint
  topts.iaSequential = -1;                         %% default one gulp

  % quants = [0 0.01 0.02 0.03 0.04 0.05 0.10 0.25 0.50 0.75 0.9 0.95 0.96 0.97 0.98 0.99 1.00];
  topts.dataset   = -1;   %% (-1) AIRS 18 year quantile dataset, Sergio Aug 2021   2002/09-2020/08 FULL 18 years
  topts.dataset   = +1;   %% (+1) AIRS 18 year quantile dataset, Strow  March 2021 2002/09-2020/08 FULL 18 years
  topts.dataset   = +2;   %% (+2) AIRS 19 year quantile dataset, Sergio Aug 2021   2002/09-2021/07 PARTIAL 19 years
  topts.dataset   = +3;   %% (+3) AIRS 19 year extreme  dataset, Sergio Aug 2021   2002/09-2021/07 PARTIAL 19 years, EXTREME
  topts.dataset   = -3;   %% (-3) AIRS 19 year mean     dataset, Sergio Aug 2021   2002/09-2020/08 AUTOMATIC USES Q00, MEAN
  topts.dataset   = +4;   %% (+4) AIRS 19 year quantile dataset, Sergio Aug 2021   2002/09-2021/08 FULL 19 years ************************ JPL April 2021 Sounder Science Meeting
  topts.dataset   = +5;   %% (+5) AIRS 12 year quantile dataset, Sergio Aug 2022   2002/09-2014/08 FULL 12 years
  topts.dataset   = +6;   %% (+6) AIRS = CRIS NSR 07 year quantile dataset,        2012/05-2019/04 FULL 07 years
  topts.dataset   = +7;   %% (+7) AIRS 20 year quantile dataset, Sergio Sep 2022   2002/09-2022/08 FULL 20 years ************************
  topts.dataset   = +8;   %% (+8) AIRS = OCO2  07 year quantile dataset            2015/01-2021/12 OCO2 FULL 07 years
  % quants = [0.50 0.80 0.90 0.95 0.97 1.00];
  topts.dataset   = +9;   %% (+9) AIRS 20 year quantile dataset, Sergio Oct 2022   2002/09-2022/08 FULL 20 years, new way of douning quantile iQAX = 3  ************************

  topts.dataset   = ia_OorC_DataSet_Quantile(2);

  %%%%%%%%%%

  iQuantile = 04;  %% 05% so very cloudy (hope SST jac can take care of that) and convection
  iQuantile = 14;  %% 
  iQuantile = 08;  %% 50% so has clouds (hope SST jac can take care of that) and convection -- this is 0.25 - 0.50
  iQuantile = 09;  %% 50% so has clouds (hope SST jac can take care of that) and convection -- this is 0.90 - 0.75
  iQuantile = 04;  %% quite cloudy (hopefully)
  iQuantile = 00;  %% mean     <<<<***** IF YOU SET THIS THEN topts.dataset is ignored, uses topts.dataset   = -3; *****>>>>
  iQuantile = 50;  %% top 5 quantiles averaged (so some cloud and hottest)
  iQuantile = 04;  %% Q0.95, iQAX = 3, dataset = 9
  iQuantile = 01;  %% Q0.80, iQAX = 3, dataset = 9
  iQuantile = 03;  %% Q0.90, iQAX = 3, dataset = 9
  iQuantile = 16;  %% Q0.99 hottest, for AIRS STM, dataset = 7,9 (yeah the last is a fudge!) -- use this when fitting CAL
  iQuantile = 05;  %% Q0.97, iQAX = 3, dataset = 9

  iQuantile = ia_OorC_DataSet_Quantile(3);

  driver.NorD = -1; %% day, asc
  driver.NorD = +1; %% night, desc
  driver.iQuantile = iQuantile;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  iAdjLowerAtmWVfrac = 1;                             %% WARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!
  iAdjLowerAtmWVfrac = 0;                             %% WARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!
  topts.iAdjLowerAtmWVfrac = iAdjLowerAtmWVfrac;      %% WARNING this also sets WV in lower part of atmos, depending on dBT1231/dt by using iAdjLoweAtmWVfrac !!!!!

  topts.tie_sst_lowestlayer = -1;  %% DEFAULT
  topts.tie_sst_lowestlayer = +1;  %% testing dataset=4,iQuantil=16,ocb_set=0 (the JPL SOunder meeting Apr 2022, 04/23/2022 commit 30d2e554a97b34b0923ad58346d183a3c10d6bcb

  %topts.iFixTG_NoFit = +1; %% dump out first scalar = CO2 boy certainly messes up CO2 even more!!!!!

  topts.iSergioCO2 = +1; %% fit for CO2/CH4/N2O
  topts.iSergioCO2 = -1; %% assume ESRL CO2/CH4
  if topts.dataset == 8 
    topts.iSergioCO2 = +1; %% fit for CO2/CH4/N2O
  end

  topts.ocb_set = ia_OorC_DataSet_Quantile(1);

  topts.set_era5_cmip6_airsL3_WV_T_O3 = +2;  %% use T+ST
  topts.set_era5_cmip6_airsL3_WV_T_O3 = -1;  %% use WV/T+ST/O3
  topts.set_era5_cmip6_airsL3_WV_T_O3 = +40; %% use T only in the lower trop to "start things" in the polar region, based on BT1231
  topts.set_era5_cmip6_airsL3_WV_T_O3 = +4;  %% use T only

  topts.set_era5_cmip6_airsL3 = 8;           %% use MLS a priori
  topts.set_era5_cmip6_airsL3 = 5;           %% use ERA5 a priori
  topts.set_era5_cmip6_airsL3 = 0;           %% use 0 a priori

  topts.iNlays_retrieve = 20; %% default, 5 AIRS lays thick
  topts.iNlays_retrieve = 50; %%          2 AIRS lays thick

  topts.resetnorm2one = -1;   %% DEFAULT
  topts.resetnorm2one = +1;

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

  %% iNorD > 0 ==> night
  if topts.ocb_set == 0 & driver.i16daytimestep > 0 & driver.NorD > 0 & topts.dataset ~= 3
    driver.outfilename = ['OutputAnomaly_OBS/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(iInd,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif topts.ocb_set == 1 & driver.i16daytimestep > 0 & driver.NorD > 0 & topts.dataset ~= 3
    driver.outfilename = ['OutputAnomaly_CAL/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(iInd,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif topts.ocb_set == 0 & driver.i16daytimestep < 0 & driver.NorD > 0 & topts.dataset ~= 3
    outdir = ['Output/Quantile' num2str(driver.iQuantile,'%02d')];
    driver.outfilename = [outdir  '/test' int2str(iInd) '.mat'];
  elseif topts.ocb_set == 1 & driver.i16daytimestep < 0 & driver.NorD > 0 & topts.dataset ~= 3
    outdir = ['Output_CAL/Quantile' num2str(driver.iQuantile,'%02d')];
    driver.outfilename = [outdir  '/test' int2str(iInd) '.mat'];
  elseif topts.ocb_set == 0 & driver.i16daytimestep < 0 & driver.NorD > 0 & topts.dataset == 3 %% EXTREME
    outdir = ['Output/Extreme/'];
    driver.outfilename = [outdir  '/test' int2str(iInd) '.mat'];
  %% iNorD < 0 ==> day
  elseif topts.ocb_set == 0 & driver.i16daytimestep > 0 & driver.NorD < 0 & topts.dataset ~= 3
    driver.outfilename = ['OutputAnomaly_OBS_Day/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(iInd,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
  elseif topts.ocb_set == 1 & driver.i16daytimestep > 0 & driver.NorD < 0 & topts.dataset ~= 3
    driver.outfilename = ['OutputAnomaly_CAL_Dat/Quantile' num2str(driver.iQuantile,'%02d') '/' num2str(iInd,'%02d') '/anomtest_timestep' int2str(driver.i16daytimestep) '.mat'];
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
  fprintf(1,'A outdir,outname = %s \n',driver.outfilename)
  % Do the retrieval
  if ~exist(driver.outfilename)     
     driver = retrieval(driver,aux);
     if oldstyle_removeCO2 < 0
       %% now adjust everything!!!!!
       call_adjustCO2
     end
     junk  = load('h2645structure.mat');
     f2645 = junk.h.vchan;

     figure(3); plot(1:2645,driver.rateset.rates,'b',1:2645,driver.oem.fit,'r'); title('AIRS_gridded_Oct2020_trendsonly')
     figure(3); plot(f2645,driver.rateset.rates,'b',f2645,driver.oem.fit,'r'); title('AIRS_gridded_Oct2020_trendsonly')
  else
    disp(' ')
    fprintf(1,'%s after all that, already exists \n',driver.outfilename)
    fprintf(1,'%s after all that, already exists \n',driver.outfilename)
    fprintf(1,'%s after all that, already exists \n',driver.outfilename)
    driver.rateset.rates = zeros(2645,1);
    disp(' ')
  end
driver

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
  %[hMean17years,ha,pMean17years,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
  [hMean17years,ha,pMean17years,pa] = rtpread('/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/2012/FixedNAN/all4608_era5_full12months_Qcumulative09.rtp');
  print_cloud_params(hMean17years,pMean17years,driver.iibin); 

  if topts.iNlays_retrieve >= 97
    %disp('ret to plot retrievals'); pause
    pause(1)
    plot_retrieval_latbins
  else
    %disp('ret to plot retrievals'); pause
    pause(1)
    plot_retrieval_latbins_fewlays
  end
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
