function [driver,settings] = set_driver_rateset_datafile(driver0,settings0);

%% iUseNWP = -1 for use AIRS obs/cal rates
%%         = +3 for AIRS L3, -3 for CLIMCAPS L3
%%         = +5 for ERA5, +2 for MERRA2
%%         = +6 for CMIP6, -6 for AMIP6
%%             1 for N, 2 for D

driver = driver0;
settings = settings0;

switch settings.dataset 
  case -1
    settings.iNumYears = 18;
  case +1
    settings.iNumYears = 18;
  case +2
    settings.iNumYears = 19;
  case +3
    settings.iNumYears = 19;
  case -3
    settings.iNumYears = 19;
  case +4
    settings.iNumYears = 19;
  case +5
    settings.iNumYears = 12;
  case +6
    settings.iNumYears = 07; %% though 2012/05 to 2019/04, CRIS
  case +7
    settings.iNumYears = 20;
  case +8
    settings.iNumYears = 07; %% though 2015/01 to 2021/12, OCO
  %%%%%%%%%%%%%%%%%%%%%%%%%
  case +9
    settings.iNumYears = 20;
  case +10
    settings.iNumYears = 05;
  case +11
    settings.iNumYears = 10;
  case +12
    settings.iNumYears = 15;
  %%%%%%%%%%%%%%%%%%%%%%%%%
  case +13
    settings.iNumYears = 20;
  case +14
    settings.iNumYears = -4.0;
  case +15
    settings.iNumYears = 15;
  %%%%%%%%%%%%%%%%%%%%%%%%%
  case +16
    settings.iNumYears = -4.1;
  case +17
    settings.iNumYears = 22;
  %%%%%%%%%%%%%%%%%%%%%%%%%
  case +30
    settings.iNumYears = 20;
end

%fprintf(1,' in set_driver_rateset_datafile.m : [settings.descORasc driver.i16daytimestep settings.dataset] = %3i %3i %3i \n',[settings.descORasc driver.i16daytimestep settings.dataset])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if settings.dataset == -2   
  disp('AIRS 16 year rates or anomalies, NO nu cal done')
  error('oops not done')

elseif abs(settings.dataset) >= 1 & abs(settings.dataset) <= 17
  disp('AIRS 07 or 18 or 19 or 20 year rates or anomalies, nu cal done in there')
  if settings.descORasc == +1 & driver.i16daytimestep < 0 & settings.dataset == 1    
    disp('doing Strow 18 yr gridded quantile rates')
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsonly_transpose.mat';  %% ah puts strides, as hoped/expected ie WRONG
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsonly.mat';           
      driver.rateset.datafile  = ['convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           

    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcs.mat';           
      driver.rateset.datafile  = 'AHAH';
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & settings.dataset == -1    
    disp('doing Sergio 18 yr gridded quantile rates')
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsonly_transpose.mat';  %% ah puts strides, as hoped/expected ie WRONG
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsonly.mat';           
      driver.rateset.datafile  = ['iType_-1_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           

    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcs.mat';           
      driver.rateset.datafile  = 'AHAH';
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & settings.dataset == 2    
    disp('doing Sergio 19 year gridded quantile rates')
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0
      driver.rateset.datafile  = ['iType_2_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           

    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcs.mat';           
      driver.rateset.datafile  = 'AHAH';
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & settings.dataset == 3
    disp('doing Sergio 19 year gridded extreme rates')
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0
      driver.rateset.datafile  = ['iType_3_extreme_convert_sergio_clearskygrid_obsonly.mat'];           

    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcs.mat';           
      driver.rateset.datafile  = 'AHAH';
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & settings.dataset == -3
    disp('doing Sergio 19 year gridded mean rates')
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0
      driver.rateset.datafile  = ['iType_-3_mean_convert_sergio_clearskygrid_obsonly.mat'];           

    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcs.mat';           
      driver.rateset.datafile  = 'AHAH';
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & settings.dataset == 4    
    disp('doing Sergio FULL 19 year gridded quantile rates')
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0
      driver.rateset.datafile  = ['iType_4_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcsAHAH.mat';           
      driver.rateset.datafile  = 'AHAH';
      strlatbin                = num2str(floor((driver.iibin-1)/72)+1,'%02d');
      driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_latbin' strlatbin '.mat'];                %% co2/n2o/ch4 change in time
      driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_const_tracegas_latbin' strlatbin '.mat']; %% co2/n2o/ch4 unchanging
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & settings.dataset == 5
    disp('doing Sergio FULL 13 year gridded quantile rates')
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0
      driver.rateset.datafile  = ['iType_5_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcsAHAH.mat';           
      driver.rateset.datafile  = 'AHAH';
      strlatbin                = num2str(floor((driver.iibin-1)/72)+1,'%02d');
      driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_latbin' strlatbin '.mat'];                %% co2/n2o/ch4 change in time
      driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_const_tracegas_latbin' strlatbin '.mat']; %% co2/n2o/ch4 unchanging
      driver.rateset.datafile  = 'AHAH';
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & settings.dataset == 6
    disp('doing Sergio FULL 07 year gridded quantile rates 2012/05-2019/04 CRIS NSR')
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0
      driver.rateset.datafile  = ['iType_6_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcsAHAH.mat';           
      driver.rateset.datafile  = 'AHAH';
      strlatbin                = num2str(floor((driver.iibin-1)/72)+1,'%02d');
      driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_latbin' strlatbin '.mat'];                %% co2/n2o/ch4 change in time
      driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_const_tracegas_latbin' strlatbin '.mat']; %% co2/n2o/ch4 unchanging
      driver.rateset.datafile  = 'AHAH';
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & settings.dataset == 7
    disp('doing Sergio FULL 20 year gridded quantile rates 2002/09-2022/08 ')
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0
      driver.rateset.datafile  = ['iType_7_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcsAHAH.mat';           
      driver.rateset.datafile  = 'AHAH';
      strlatbin                = num2str(floor((driver.iibin-1)/72)+1,'%02d');
      driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_latbin' strlatbin '.mat'];                %% co2/n2o/ch4 change in time
      driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_const_tracegas_latbin' strlatbin '.mat']; %% co2/n2o/ch4 unchanging
      driver.rateset.datafile  = 'AHAH';
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & settings.dataset == 8
    disp('doing Sergio FULL 07 year gridded quantile rates 2015/01-2021/12 OCO2')
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0
      driver.rateset.datafile  = ['iType_8_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcsAHAH.mat';           
      driver.rateset.datafile  = 'AHAH';
      driver.rateset.datafile  = ['iType_8_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep >= 0 & settings.dataset == 9
    disp('doing Sergio FULL N year anomalies')
    driver.rateset.datafile  = driver.anomalydatafile;   %%% ['anomaly_ALL_quantile_globalavg_and_18_averages_timeseries_Q' num2str(driver.iQuantile,'%02d') '.mat'];

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & (settings.dataset >= 9 & settings.dataset <= 12)
    disp('doing Sergio FULL 20/05/10/15 year gridded quantile rates 2002/09-2022/2007/2012/2017-08 , NEW WAY of doing quantile iQAX = 3')
    fprintf(1,'dataset = %2i where 09,10,11,12 are for 20,05,10,15 years of AIRS data .... \n',settings.dataset)
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0 & driver.iAllorSeasonal == +1
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & driver.iAllorSeasonal == -1
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '_DJF.mat'];           
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & driver.iAllorSeasonal == -2
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '_MAM.mat'];           
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & driver.iAllorSeasonal == -3
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '_JJA.mat'];           
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & driver.iAllorSeasonal == -4
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '_SON.mat'];           
    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcsAHAH.mat';           
      driver.rateset.datafile  = 'AHAH';
      disp(' settings.dataset == 9,10,11,12 but settings.ocb_set == 1 so set rateset = 2002/09-2021/08 clr (Q16)')
      strlatbin                = num2str(floor((driver.iibin-1)/72)+1,'%02d');

      if settings.model == 5
        driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_const_tracegas_latbin' strlatbin '.mat']; %% co2/n2o/ch4 unchanging
        driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_latbin' strlatbin '.mat'];                %% co2/n2o/ch4 change in time
        driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_const_tracegas_latbin' strlatbin '_2002_09_2022_08.mat']; %% co2/n2o/ch4 unchanging
      elseif settings.model == 2
        driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/MERRA2_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_latbin' strlatbin '_2002_09_2022_08.mat']; %% co2/n2o/ch4 changing
      elseif settings.model == 3
        driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/AIRSL3_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_latbin' strlatbin '_2002_09_2022_08.mat']; %% co2/n2o/ch4 changing
      elseif settings.model == -3
        driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/CLIMCAPSL3_SARTA_SPECTRAL_RATES/KCARTA_latbin' strlatbin '/sarta_spectral_trends_latbin' strlatbin '_2002_09_2022_08.mat']; %% co2/n2o/ch4 changing
      end
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & (settings.dataset >= 13 & settings.dataset <= 15)
    disp('doing Sergio FULL 20/05/15 year gridded quantile rates 2002/09-2022/08  2018/09-2022/08  2008/01-2022/12 , NEW WAY of doing quantile iQAX = 3')
    fprintf(1,'dataset = %2i where 13,14,15 are for 20,04,15 years of AIRS data .... \n',settings.dataset)
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 13
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_4_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 14
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 15
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & (settings.dataset >= 16 & settings.dataset <= 17)
    disp('doing Sergio FULL 22/04 year gridded quantile rates 2002/09-2024/08  2020/07-2024/06, NEW WAY of doing quantile iQAX = 3')
    fprintf(1,'dataset = %2i where 13,14,15 are for 22,04 years of AIRS data .... \n',settings.dataset)
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 16
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 17
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == -1 & driver.i16daytimestep < 0
    disp('doing ascending latbin rates')
    driver.rateset.datafile  = [];
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 2
    disp('doing ascending? descending? OBS ANOMALY')
    driver.rateset.datafile = driver.anomalydatafile;
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 1
    disp('doing descending CAL ANOMALY')
    driver.rateset.datafile = [];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & (settings.dataset == 30)
    disp('doing Sergio FULL 20 year AMSU AMSU AMSU gridded quantile rates 2002/09-2022/08 , just one case (all average) iQAX = 1')
    fprintf(1,'dataset = %2i for 20 years of AMSU data .... \n',settings.dataset)
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 30
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_AMSU_iQAX_' num2str(driver.iQuantile,'%02d') '.mat'];           
    end

end

