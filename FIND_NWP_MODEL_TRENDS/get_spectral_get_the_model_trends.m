%% see plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2

iUMBC_SPECTRA = 'UMBC';
if ~exist('umbcXspectra')
  disp('get_spectral_get_the_model_trends.m : painfully reading in UMBC spectra')
  iUMBCstr = 'reconstruct_umbc_spectra_geo_rlat';
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iUMBC_SPECTRA '/' iUMBCstr num2str(ii,'%02i') '.mat'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iUMBC_SPECTRA '/' iUMBCstr num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    umbcL3.umbcL3_spectral_rates(:,ind) = thesave;
  end
else
  umbcXspectra = load(strUMBC,'fits','rates');
  umbcL3.umbcL3_spectral_rates = umbcXspectra.fits;
end

if iJorC == 1
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iJorCstrSPECTRA '/' iJorCFstr num2str(ii,'%02i') '.mat'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iJorCstrSPECTRA '/' iJorCFstr num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    airsL3.airsL3_spectral_rates(:,ind) = thesave;
  end
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/CLIMCAPSL3/' iJorCFstr num2str(ii,'%02i') '.mat'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/CLIMCAPSL3/reconstruct_climcapsL3_spectra_geo_rlat' num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    climcapsL3.climcapsL3_spectral_rates(:,ind) = thesave;
  end  

elseif iJorC == -1
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/AIRSL3/reconstruct_airsL3_spectra_geo_rlat' num2str(ii,'%02i') '.mat'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/AIRSL3/reconstruct_airsL3_spectra_geo_rlat' num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    airsL3.airsL3_spectral_rates(:,ind) = thesave;
  end
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/CLIMCAPSL3/' iJorCFstr num2str(ii,'%02i') '.mat'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/CLIMCAPSL3/reconstruct_climcapsL3_spectra_geo_rlat' num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    climcapsL3.climcapsL3_spectral_rates(:,ind) = thesave;
  end
end

%% see plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2
if iEorM == 5
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iEorMstrSPECTRA '/' iEorMFstr num2str(ii,'%02i') '.mat'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iEorMstrSPECTRA '/' iEorMFstr num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    era5.era5_spectral_rates(:,ind) = thesave;
  end
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/MERRA2/reconstruct_merra2_spectra_geo_rlat' iEorMFstr num2str(ii,'%02i') '.mat'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/MERRA2/reconstruct_merra2_spectra_geo_rlat' num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    merra2.merra2_spectral_rates(:,ind) = thesave;
  end
elseif iEorM == 2
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/ERA5/reconstruct_era5_spectra_geo_rlat' iEorMFstr num2str(ii,'%02i') '.mat'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/ERA5/reconstruct_era5_spectra_geo_rlat' num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    era5.era5_spectral_rates(:,ind) = thesave;
  end
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iEorMstrSPECTRA '/' iEorMFstr num2str(ii,'%02i') '.mat'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iEorMstrSPECTRA '/' iEorMFstr num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    merra2.merra2_spectral_rates(:,ind) = thesave;
  end
end

%% see plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2
for ii = 1 : 64
  fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iAorCstrSPECTRA '/' iAorCFstr num2str(ii,'%02i') '.mat'];
  %fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iAorCstrSPECTRA '/' iAorCFstr num2str(ii,'%02i') '_2002_09_2022_08.mat'];
  junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
  junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
  ind = (ii-1)*72 + (1:72);
  cmip6.cmip6_spectral_rates(:,ind) = thesave;
end
