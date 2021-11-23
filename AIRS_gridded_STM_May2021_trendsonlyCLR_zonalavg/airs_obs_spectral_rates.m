if ~exist('all64_airs_dbt_trends.mat')
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/DATAObsStats_StartSept2002_CORRECT_LatLon/'];
    fname = [fname '/LatBin' num2str(ii,'%02d') '/trends_zonalavg_fits_quantiles_LatBin' num2str(ii,'%02d') '_timesetps_001_429_V1.mat'];
    load(fname);
  
    airs_dbt_desc(:,ii)     = trends.dbt_desc(:,16);
    airs_dbt_desc_unc(:,ii) = trends.dbt_err_desc(:,16);
  
    airs_dbt_asc(:,ii)     = trends.dbt_asc(:,16);
    airs_dbt_asc_unc(:,ii) = trends.dbt_err_asc(:,16);
  end
  comment = 'see airs_obs_spectral_rates.m and get_rates.m');
  save all64_airs_dbt_trends.mat comment airs_dbt*
end
