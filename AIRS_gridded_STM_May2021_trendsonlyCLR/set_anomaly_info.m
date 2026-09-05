%% see     driver_stand_alone_make_globalavg_TWP_and_N_average_anomalies_zonalavg.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD ON TAKI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (ONLY global + tropics ONLY) x 500 timesteps = 01000 points to fit; if each processor needs to do 020 of these, we need 01000/020 = 50 procesors; or 750 processors to do everything
ia_OorC_DataSet_Quantile = [+2 17 03 -9999];
  iNumAnomTimeSteps = 500; iNumAnomTiles = 30; iNumAnomJobsPerProc =  050; 
  anomalydatafile = 'anomalyD_zonalavg_globalavg_and_TWPlat35lon66_and_28_averages_timeseries_Q03_numyears_22_iNumAnomTimeSteps_500.mat';  %% needs 500*30/250 = 60 processors   btavgAnomFinal = [2645x15000] = 500*28 anomaly time series + 1 global + 1 tropical
  anomalydatafile = 'anomalyD_zonalavg_globalavg_and_TWPlat35avg_and_28_averages_timeseries_Q03_numyears_22_iNumAnomTimeSteps_500.mat';    %% needs 500*30/250 = 60 processors   btavgAnomFinal = [2645x15000] = 500*28 anomaly time series + 1 global + 1 tropical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW ON CHIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ia_OorC_DataSet_Quantile = [+2 30 01 -9999]; %% ocb_set = 2 : AMSU obs fit,    dataset = 30, iQuantile = 01    20 year anomalies, 2002/09-2022/08 AMSU obs Q(0.50-->1) -- technically this is "allsky average' but should be clear

ia_OorC_DataSet_Quantile = [+2 20 01 -9999]; %% ocb_set = 2 : AIRS obs fit,    dataset = 20, iQuantile = 03    20 year anomalies, 2002/09-2022/08 AMSU obs Q(0.50-->1)
  iNumAnomTimeSteps = 523; iNumAnomTiles = 65; iNumAnomJobsPerProc =  550; 
  anomalydatafile = 'anomalyA_zonalavg_globalavg_and_64_averages_timeseries_Q03_numyears_23.00_iNumAnomTimeSteps_525.mat'; %% needs 523*65/550 = 62 processors   btavgAnomFinal = [2645x33995] = 523*64 anomaly time series + 523*1 glob   

