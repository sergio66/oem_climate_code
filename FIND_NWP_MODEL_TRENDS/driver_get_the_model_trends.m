%{
[sergio@taki-usr2 FIND_NWP_MODEL_TRENDS]$ ls -lt ../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries
total 120
drwxrwxr-x 2 sergio pi_strow 12288 Mar 31 23:38 AMIP6
drwxrwxr-x 2 sergio pi_strow 12288 Mar 31 21:04 CMIP6
drwxrwxr-x 2 sergio pi_strow 16384 Mar 31 17:46 MERRA2
drwxrwxr-x 2 sergio pi_strow 16384 Mar 31 11:48 CLIMCAPSL3
-rw-rw-r-- 1 sergio pi_strow  9761 Mar 31 07:27 Readme
drwxrwxr-x 2 sergio pi_strow 16384 Mar 31 07:23 AIRSL3
drwxrwxr-x 2 sergio pi_strow 12288 Mar 31 07:21 ERA5
%}

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% superceded driver_get_the_model_trends_orig
iJorC = input('(+1, default) J.Susskind AIRS L3   or (-1) C.Barnet CLIMCAPS  or (+3) CESM3 19 years : ');
if length(iJorC) == 0
  iJorC = +1;
end
if iJorC == -1
  iJorCstr = 'CLIMCAPS';
  iJorCstrSPECTRA = 'CLIMCAPSL3';
  iJorCFstr = '/reconstruct_climcapsL3_spectra_geo_rlat';
  disp('remember all your plots with title AIRSL3 should really be CLIMCAPSL3')
else
  iJorCstr = 'AIRSL3';
  iJorCstrSPECTRA = 'AIRSL3';
  iJorCFstr = '/reconstruct_airsL3_spectra_geo_rlat';
end

iEorM = input('(5,default) ERA5                   or (2) MERRA2 or (1) ERA-I  : ');
if length(iEorM) == 0
  iEorM = +5;
end
if iEorM == 1
  iEorMstr = 'ERA-I';
  disp('remember all your plots with title ERA5 should really be ERA-I')
elseif iEorM == 2
  iEorMstr = 'MERRA2';
  iEorMFstr = '/reconstruct_merra2_spectra_geo_rlat';
  disp('remember all your plots with title ERA5 should really be MERRA2')
elseif iEorM == 5
  iEorMstr = 'ERA5';
  iEorMFstr = '/reconstruct_era5_spectra_geo_rlat';
end
iEorMstrSPECTRA = iEorMstr;

iAorC = input('(-1,default) CMIP6 14 years        or (+1) AMIP6 14 years : ');
if length(iAorC) == 0
  iAorC = -1;
end
if iAorC == +3
  iAorCstr = 'CESM3';
  iAorCFstr = '/reconstruct_cesm3_spectra_geo_rlat';
  disp('remember all your plots with title CMIP6 should really be CESM3');
elseif iAorC == +1
  iAorCstr = 'AMIP6';
  iAorCFstr = '/reconstruct_amip6_spectra_geo_rlat';
  disp('remember all your plots with title CMIP6 should really be AMIP6');
elseif iAorC == -1
  iAorCstr = 'CMIP6';
  iAorCFstr = '/reconstruct_cmip6_spectra_geo_rlat';
end
iAorCstrSPECTRA = iAorCstr;

strMODELS = [iJorCstr '_' iEorMstr '_' iAorCstr];
if ~exist('iNumYears')
  disp('setting iNumYears == 20 in driver_get_the_model_trends')
  iNumYears = 20;
end
[airsL3,era5,cmip6] = driverchoose_AIRSvsNWPvsXMIP6(iJorC,iEorM,iAorC,iNumYears);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iSpectra = input('  now get model/L3 spectral rates???? (-1/+1 [default]) : ');
if length(iSpectra) == 0
  iSpectra = +1;
end

if iSpectra == 1
  %% see eg driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m
  dirout = '../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries';

  if ~exist('strUMBC')
    strUMBC = '/asl/s1/sergio/JUNK/test9_guessstartWV_Vers1_march22_2023.mat';
  end
  fprintf(1,'driver_get_the_model_trends.m is using obs rates from %s \n',strUMBC)

  obsrates = load(strUMBC,'rates');  
  obslat = load(strUMBC,'rlat');
  obslon = load(strUMBC,'rlon');

  [YY,XX] = meshgrid(obslat.rlat,obslon.rlon);
  YY = YY(:);
  cosYY = cos(YY*pi/180)';
  cosYY = ones(2645,1) * cosYY;

  get_spectral_get_the_model_trends

  airsL3.fchanx = fchanx;
  era5.fchanx = fchanx;
  cmip6.fchanx = fchanx;

  %plot_spectral_get_the_model_trends
  plot_spectral_get_the_model_trends2

end
  
