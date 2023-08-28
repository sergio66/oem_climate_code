if ~exist('nwp_spectral_trends_cmip6_era5_airsL3_umbc')
  nwp_spectral_trends_cmip6_era5_airsL3_umbc = make_profile_spectral_trends(cmip6,era5,airsL3,results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,2,iVersJac,-1);
    aL3trend.stemp = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
    aL3trend.ptemp = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp;
    aL3trend.gas_1 = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1;
    aL3trend.gas_3 = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_3;

    c6trend.stemp = cmip6.trend_stemp;
    c6trend.ptemp = cmip6.trend_ptemp;
    c6trend.gas_1 = cmip6.trend_gas_1;
    c6trend.gas_3 = cmip6.trend_gas_3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('nwp_spectral_trends_amip6_merra2_climcapsL3_umbc')
  nwp_spectral_trends_amip6_merra2_climcapsL3_umbc = make_profile_spectral_trends(amip6,merra2,climcapsL3,results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,2,iVersJac,-1);
    cL3trend.stemp = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.stemp;
    cL3trend.ptemp = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.ptemp;
    cL3trend.gas_1 = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.gas_1;
    cL3trend.gas_3 = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.gas_3;

    a6trend.stemp = amip6.trend_stemp;
    a6trend.ptemp = amip6.trend_ptemp;
    a6trend.gas_1 = amip6.trend_gas_1;
    a6trend.gas_3 = amip6.trend_gas_3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iaComputeWhichFeedback = [0];     %% compute + plot feedbacks only
disp(' ')
disp('umbc')
umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1);
% dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
% figure(4); aslprint([dir0 'umbc_feedback_linearity_20yrs_UMBC.pdf'])

disp(' ')
disp('era5')
era5_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,era5.trend_ptemp,era5.trend_gas_1,era5.trend_gas_3,era5_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1);

disp(' ')
disp('airsL3')
airsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,aL3trend.stemp,aL3trend.ptemp,aL3trend.gas_1,aL3trend.gas_3,airsL3_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1);

disp(' ')
disp('cmip6')
cmip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,c6trend.stemp,c6trend.ptemp,c6trend.gas_1,c6trend.gas_3,cmip6_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1);

%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('merra2')
merra2_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,merra2.trend_stemp,merra2.trend_ptemp,merra2.trend_gas_1,merra2.trend_gas_3,merra2_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1);

disp(' ')
disp('climcapsL3')
climcapsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,cL3trend.stemp,cL3trend.ptemp,cL3trend.gas_1,cL3trend.gas_3,climcapsL3_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1);

disp(' ')
disp('amip6')
amip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,a6trend.stemp,a6trend.ptemp,a6trend.gas_1,a6trend.gas_3,amip6_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('POLYFIT results')
plot_show_olr_ecRad_feedback_polyfit

disp(' ')
disp('ROBUSTFIT results')
plot_show_olr_ecRad_feedback_robustfit

disp(' ')
disp('ROBUSTFIT results TROPICS')
plot_show_olr_ecRad_feedback_robustfit_tropics

