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

%% superceded driver_get_the_model_trends_orig
iJorC = input('(+1, default) J.Susskind AIRS L3  or (-1) C.Barnet CLIMCAPS  : ');
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

iEorM = input('(5,default) ERA5                  or (2) MERRA2 or (1) ERA-I  : ');
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

iAorC = input('(-1,default) CMIP6                or (+1) AMIP6 : ');
if length(iAorC) == 0
  iAorC = -1;
end
if iAorC == +1
  iAorCstr = 'AMIP6';
  iAorCFstr = '/reconstruct_amip6_spectra_geo_rlat';
  disp('remember all your plots with title CMIP6 should really be AMIP6');
elseif iAorC == -1
  iAorCstr = 'CMIP6';
  iAorCFstr = '/reconstruct_cmip6_spectra_geo_rlat';
end
iAorCstrSPECTRA = iAorCstr;

strMODELS = [iJorCstr '_' iEorMstr '_' iAorCstr];
[airsL3,era5,cmip6] = driverchoose_AIRSvsNWPvsXMIP6(iJorC,iEorM,iAorC,iNumYears);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iSpectra = input('  now get model/L3 spectral rates???? (-1/+1 default) : ');
if length(iSpectra) == 0
  iSpectra = +1;
end
if iSpectra == 1
  %% see eg driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m
  dirout = '../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries';

  %% see plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iJorCstrSPECTRA '/' iJorCFstr num2str(ii,'%02i') '.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    airsL3.airsL3_spectral_rates(:,ind) = thesave;
  end

  %% see plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2
  if iEorM == 5
    disp('loading in 2002/09 to 2022/08 spectral rates')
  end
  for ii = 1 : 64
    if iEorM ~= 5
      fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iEorMstrSPECTRA '/' iEorMFstr num2str(ii,'%02i') '.mat'];
    else
      fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iEorMstrSPECTRA '/' iEorMFstr num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    end
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    era5.era5_spectral_rates(:,ind) = thesave;
  end

  %% see plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iAorCstrSPECTRA '/' iAorCFstr num2str(ii,'%02i') '.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    cmip6.cmip6_spectral_rates(:,ind) = thesave;
  end

  airsL3.fchanx = fchanx;
  era5.fchanx = fchanx;
  cmip6.fchanx = fchanx;

  figure(5); plot(fchanx,nanmean(airsL3.airsL3_spectral_rates'),'b',fchanx,nanmean(era5.era5_spectral_rates'),'r',fchanx,nanmean(cmip6.cmip6_spectral_rates'),'k');
   xlim([640 1640]); ylim([-0.08 +0.04])
    plotaxis2;
   hl = legend(iJorCstrSPECTRA,iEorMstrSPECTRA,iAorCstrSPECTRA,'location','best','fontsize',8);
end
  
