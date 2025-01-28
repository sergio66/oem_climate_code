diroutX = dirout;
%diroutX = ' ';

ii = JOB;

disp('here I am doing this save')
if iType == 1
  saver = ['save ' diroutX 'reconstruct_umbc_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx *xconstr*'];
  %saver = [saver ' zonalrlat zonalplays zonalRHUMBCrate zonalTUMBCrate *xconstr*'];

%%%%%%%%%%

elseif iType == 2
  if ~exist('idRH')
    %% this is driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2.m
    saver = ['save ' diroutX 'reconstruct_merra2_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
    saver = [saver ' zonalrlat zonalplays zonalRHMERRA2rate zonalTMERRA2rate *xconstr*'];
  else
    %% this is driver_check_WV_T_RH_MERRA2_geo_and_spectral_rates2_deltaRH.m
    saver = ['save ' diroutX 'reconstruct_merra2_spectra_geo_idRH_' num2str(idRH) '_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
    saver = [saver ' zonalrlat zonalplays zonalRHMERRA2rate zonalTMERRA2rate *xconstr*'];
  end

%%%%%%%%%%

elseif iType == 5
  if ~exist('idRH')
    %% this is driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
    saver = ['save ' diroutX 'reconstruct_era5_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
    saver = [saver ' zonalrlat zonalplays zonalRHERA5rate zonalTERA5rate *xconstr*'];
  else
    %% this is driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_deltaRH.m
    saver = ['save ' diroutX 'reconstruct_era5_spectra_geo_idRH_' num2str(idRH) '_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
    saver = [saver ' zonalrlat zonalplays zonalRHERA5rate zonalTERA5rate *xconstr*'];
  end
elseif iType == 51
  %% this is driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_constracegas.m
  saver = ['save ' diroutX 'reconstruct_era5_const_tracegas_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHERA5rate zonalTERA5rate *xconstr*'];

%%%%%%%%%%

elseif iType == 6
  saver = ['save ' diroutX 'reconstruct_cmip6_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHCMIP6rate zonalTCMIP6rate *xconstr*'];
elseif iType == 61
  saver = ['save ' diroutX 'reconstruct_cmip6_const_tracegas_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHCMIP6rate zonalTCMIP6rate *xconstr*'];
elseif iType == 7
  saver = ['save ' diroutX 'reconstruct_amip6_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHAMIP6rate zonalTAMIP6rate *xconstr*'];

%%%%%%%%%%

elseif iType == 3
  saver = ['save ' diroutX 'reconstruct_airsL3_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHAIRSL3rate zonalTAIRSL3rate *xconstr*'];
elseif iType == 4
  saver = ['save ' diroutX 'reconstruct_climcapsL3_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHAIRSCLIMCAPSL3rate zonalTAIRSCLIMCAPSL3rate *xconstr*'];

end

foutname = findstr(saver,'.mat');
foutname = [saver(6:foutname) 'mat'];

if ~exist(foutname)
  saver = [saver ' plevsnwp plevsx dayOFtime'];
  fprintf(1,'%s \n',saver');
  eval(saver)
else
  fprintf(1,'%s already exists, not saving \n',foutname)
  eval(saver)
end

