%% ocb_set = 2 : anomaly fit, dataset = 9, iQuantile = 03    20 year anomalies== > 20yrs* 23steps/yr = 460; AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+2 09 03 -9999]; iNumAnomTimeSteps = 454; iNumAnomTiles = 10; iNumAnomJobsPerProc =  72; 
  anomalydatafile = 'anomaly_globalavg_and_9_averages_timeseries_Q03.mat';   %% needs 454*10/72 = 64 processors       btavgAnomFinal: [2645x4540 double] ==> 454 timesteps x 10 anomaly time series

%% ocb_set = 2 : anomaly fit, dataset = 9, iQuantile = 03    20 year anomalies== > 20yrs* 23steps/yr = 460; AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+2 09 03 -9999]; iNumAnomTimeSteps = 454; iNumAnomTiles = 19; iNumAnomJobsPerProc =  72; 
  anomalydatafile = 'anomaly_globalavg_and_18_averages_timeseries_Q03.mat';  %% needs 454*19/72 = 120 processors      btavgAnomFinal: [2645x8626 double] ==> 454 timesteps x 19 anomaly time series

%% ocb_set = 2 : anomaly fit, dataset = 9, iQuantile = 03    20 year anomalies== > 20yrs* 23steps/yr = 460; AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+2 09 03 -9999]; iNumAnomTimeSteps = 454; iNumAnomTiles = 29; iNumAnomJobsPerProc =  110; 
  anomalydatafile = 'anomaly_globalavg_and_28_averages_timeseries_Q03.mat';  %% needs 454*29/110 = 120 processors     btavgAnomFinal: [2645x13166 double] ==> 454 timesteps x 99 anomaly time series

%% ocb_set = 2 : anomaly fit, dataset = 9, iQuantile = 03    20 year anomalies== > 20yrs* 23steps/yr = 460; AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+2 09 03 -9999]; iNumAnomTimeSteps = 454; iNumAnomTiles = 1; iNumAnomJobsPerProc =  20; 
  anomalydatafile = 'anomaly_tile_2515_timeseries_Q03.mat';  %% needs 454/20 = 23 processors                          btavgAnomFinal = [2645x454 double] ==> 454 x 1 

%% ocb_set = 2 : anomaly fit, dataset = 9, iQuantile = 03    20 year anomalies== > 20yrs* 23steps/yr = 460; AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+2 09 03 -9999]; iNumAnomTimeSteps = 454; iNumAnomTiles = 1; iNumAnomJobsPerProc =  20; 
  anomalydatafile = 'anomaly_tile_2515_timeseries_Q04.mat';  %% needs 454/20 = 23 processors                          btavgAnomFinal = [2645x454 double] ==> 454 x 1 

%% ocb_set = 2 : anomaly fit, dataset = 9, iQuantile = 03    20 year anomalies== > 20yrs* 23steps/yr = 460; AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+2 09 03 -9999]; iNumAnomTimeSteps = 454; iNumAnomTiles = 1; iNumAnomJobsPerProc =  20; 
  anomalydatafile = 'anomaly_tile_2515_timeseries_Q05.mat';  %% needs 454/20 = 23 processors                          btavgAnomFinal = [2645x454 double] ==> 454 x 1 

%% ocb_set = 2 : anomaly fit, dataset = 16, iQuantile = 03   22 year anomalies== > 22yrs* 23steps/yr = 500; AIRS obs Q(0.90-->1)
ia_OorC_DataSet_Quantile = [+2 17 03 -9999]; iNumAnomTimeSteps = 500; iNumAnomTiles = 64; iNumAnomJobsPerProc =  200; 
  anomalydatafile = '/asl/s1/sergio/JUNK/anomaly_zonalavg_ALL_Q03_numyears_22.00_iNumAnomTimeSteps_500_A.mat';          %% needs 500*64/200 =  160 processors     btavgAnomFinal = [64x2645x500] ==> 500 timesteps x 64 anomaly time series
%% (global + 28 lats) x 500 timesteps = 14500 points to fit; if each processor needs to do 250 of these, we need 14500/250 = 58 procesors
ia_OorC_DataSet_Quantile = [+2 17 03 -9999]; iNumAnomTimeSteps = 500; iNumAnomTiles = 29; iNumAnomJobsPerProc =  250; 
  anomalydatafile = 'anomalyD_zonalavg_globalavg_and_28_averages_timeseries_Q03_numyears_22_iNumAnomTimeSteps_500.mat';  %% needs 500*29/250 = 58 processors       btavgAnomFinal = [2645x14500] = 500*28 anomaly time series + 1 global

%% (global + tropics + 28 lats) x 500 timesteps = 15000 points to fit; if each processor needs to do 250 of these, we need 15000/250 = 60 procesors
ia_OorC_DataSet_Quantile = [+2 17 03 -9999]; iNumAnomTimeSteps = 500; iNumAnomTiles = 30; iNumAnomJobsPerProc =  250; 
  anomalydatafile = 'anomalyD_zonalavg_globalavg_and_tropics_and_28_averages_timeseries_Q03_numyears_22_iNumAnomTimeSteps_500.mat';  %% needs 500*30/250 = 60 processors   btavgAnomFinal = [2645x15000] = 500*28 anomaly time series + 1 global + 1 tropical
%% (ONLY global + tropics ONLY) x 500 timesteps = 01000 points to fit; if each processor needs to do 020 of these, we need 01000/020 = 50 procesors; or 750 processors to do everything
ia_OorC_DataSet_Quantile = [+2 17 03 -9999]; iNumAnomTimeSteps = 500; iNumAnomTiles = 30; iNumAnomJobsPerProc =  050; 
  anomalydatafile = 'anomalyD_zonalavg_globalavg_and_tropics_and_28_averages_timeseries_Q03_numyears_22_iNumAnomTimeSteps_500.mat';  %% needs 500*30/250 = 60 processors   btavgAnomFinal = [2645x15000] = 500*28 anomaly time series + 1 global + 1 tropical
