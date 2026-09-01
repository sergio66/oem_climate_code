function [driver,settings] = set_driver_rateset_datafile(driver0,settings0);

%% iUseNWP = -1 for use AIRS obs/cal rates
%%         = +3 for AIRS L3, -3 for CLIMCAPS L3
%%         = +5 for ERA5, +2 for MERRA2
%%         = +6 for CMIP6, -6 for AMIP6
%%             1 for N, 2 for D

driver = driver0;
settings = settings0;

switch settings.dataset
  %%% XXXXXXXXXXXXXXXXXXXXXXXXX ALL THIS LOST in Oct 2025 XXXXXXXXXXXXXXXXXXXXXXXXX %%%    
  %%% XXXXXXXXXXXXXXXXXXXXXXXXX ALL THIS LOST in Oct 2025 XXXXXXXXXXXXXXXXXXXXXXXXX %%%
  %%% XXXXXXXXXXXXXXXXXXXXXXXXX ALL THIS LOST in Oct 2025 XXXXXXXXXXXXXXXXXXXXXXXXX %%%
       
  % case -1
  %   settings.iNumYears = 18;
  % case +1
  %   settings.iNumYears = 18;
  % case +2
  %   settings.iNumYears = 19;
  % case +3
  %   settings.iNumYears = 19;
  % case -3
  %   settings.iNumYears = 19;
  % case +4
  %   settings.iNumYears = 19;
  % case +5
  %   settings.iNumYears = 12;
  % case +6
  %   settings.iNumYears = 07; %% though 2012/05 to 2019/04, CRIS
  % case +7
  %   settings.iNumYears = 20;
  % case +8
  %   settings.iNumYears = 07; %% though 2015/01 to 2021/12, OCO
  % %%%%%%%%%%%%%%%%%%%%%%%%%
  % case +9
  %   settings.iNumYears = 20;
  % case +10
  %   settings.iNumYears = 05;
  % case +11
  %   settings.iNumYears = 10;
  % case +12
  %   settings.iNumYears = 15;
  % %%%%%%%%%%%%%%%%%%%%%%%%%
  % case +13
  %   settings.iNumYears = 20;
  % case +14
  %   settings.iNumYears = -4.0;
  % case +15
  %   settings.iNumYears = 15;
  % %%%%%%%%%%%%%%%%%%%%%%%%%
  % case +16
  %   settings.iNumYears = -4.1;
  % case +17
  %   settings.iNumYears = 22;
  % case +18
  %   settings.iNumYears = 23;
  % %%%%%%%%%%%%%%%%%%%%%%%%%
  % case +30
  %   settings.iNumYears = 20;

  %%% XXXXXXXXXXXXXXXXXXXXXXXXX ALL THIS LOST in Oct 2025 XXXXXXXXXXXXXXXXXXXXXXXXX %%%    
  %%% XXXXXXXXXXXXXXXXXXXXXXXXX ALL THIS LOST in Oct 2025 XXXXXXXXXXXXXXXXXXXXXXXXX %%%
  %%% XXXXXXXXXXXXXXXXXXXXXXXXX ALL THIS LOST in Oct 2025 XXXXXXXXXXXXXXXXXXXXXXXXX %%%
    
  case +9
    settings.iNumYears = 20; %% 2002/09 to 2022/08
  case +19
    settings.iNumYears = 23; %% 2002/09 to 2025/08
  case +20
    settings.iNumYears = 20; %% 2003/01 to 2022/12
    
end

%fprintf(1,' in set_driver_rateset_datafile.m : [settings.descORasc driver.i16daytimestep settings.dataset] = %3i %3i %3i \n',[settings.descORasc driver.i16daytimestep settings.dataset])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if settings.dataset == -2   
  disp('AIRS 16 year rates or anomalies, NO nu cal done')
  error('oops not done')

elseif abs(settings.dataset) >= 1 & abs(settings.dataset) <= 29
  disp('AIRS 07 or 18 or 19 or 20 or 22 or 23 year rates or anomalies, nu cal done in there')
  fprintf(1,'  settings.dataset      = %2i \n',settings.dataset)
  fprintf(1,'  settings.ocb_set      = %2i \n',settings.ocb_set)
  fprintf(1,'  settings.descORasc    = %2i \n',settings.descORasc)
  fprintf(1,'  driver.i16daytimestep = %2i \n',driver.i16daytimestep)    
  % all_the_lost_rateset_datafile

  if settings.descORasc == +1 & driver.i16daytimestep < 0 & (settings.dataset >= 09 & settings.dataset <= 20)
    disp('doing Sergio FULL 22/04 year gridded quantile rates 2002/09-2024/08  2020/07-2024/06, NEW WAY of doing quantile iQAX = 3')
    fprintf(1,'dataset = %2i where 13,14,15 are for 22,04 years of AIRS data .... \n',settings.dataset)
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 16
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 17
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 18
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];
      
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 19
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 20
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];           
    elseif settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 9
      driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_obsonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];


    elseif settings.ocb_set == 1  & driver.i16daytimestep < 0 & settings.dataset == 9
      %driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_iQAX_3_convert_sergio_clearskygrid_calonly_Q' num2str(driver.iQuantile,'%02d') '.mat'];
      driver.rateset.datafile  = ['SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/STS/NIGHTorAVG/ERA5//20yrs/ERA5_spectraltrends_2002_09_2022_08.mat'];
      
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end

  elseif settings.descORasc == -1 & driver.i16daytimestep < 0
    disp('doing ascending latbin rates')
    driver.rateset.datafile  = [];
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 2
    disp('doing ascending? descending? OBS ANOMALY')
    driver.rateset.datafile = driver.anomalyinfo.datafile;
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 1
    disp('doing descending CAL ANOMALY')
    driver.rateset.datafile = [];
  end

%  elseif settings.descORasc == +1 & driver.i16daytimestep < 0 & (settings.dataset == 30)
%     disp('doing Sergio FULL 20 year AMSU AMSU AMSU gridded quantile rates 2002/09-2022/08 , just one case (all average) iQAX = 1')
%     fprintf(1,'dataset = %2i for 20 years of AMSU data .... \n',settings.dataset)
%     driver.rateset.datafile  = [];
%     if settings.ocb_set == 0  & driver.i16daytimestep < 0 & settings.dataset == 30
%       driver.rateset.datafile  = ['iType_' num2str(settings.dataset) '_AMSU_iQAX_' num2str(driver.iQuantile,'%02d') '.mat'];           
%     end
    
end

if ~isfield(driver.rateset,'datafile')
  fprintf(1,'OOOPS : set_driver_rateset_datafile.m did not set datafile \n')
  fprintf(1,'[settings.descORasc driver.i16daytimestep settings.dataset] = %3i %3i %3i \n',[settings.descORasc driver.i16daytimestep settings.dataset])
  error('please check')
end
if length(driver.rateset.datafile) == 0
  fprintf(1,'OOOPS : set_driver_rateset_datafile.m did not set datafile ... length(length(driver.rateset.datafile) == 0  \n')
  fprintf(1,'[settings.descORasc driver.i16daytimestep settings.dataset] = %3i %3i %3i \n',[settings.descORasc driver.i16daytimestep settings.dataset])
  error('please check')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
fprintf(1,' <<<< driver.rateset.datafile = %s >>>> \n',driver.rateset.datafile)
disp(' ')
