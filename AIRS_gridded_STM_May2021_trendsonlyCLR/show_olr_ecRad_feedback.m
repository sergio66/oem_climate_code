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

disp(' ')
disp('merra2')
merra2_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,merra2.trend_stemp,merra2.trend_ptemp,merra2.trend_gas_1,merra2.trend_gas_3,merra2_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1);

disp(' ')
disp('climcapsL3')
climcapsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,aL3trend.stemp,aL3trend.ptemp,aL3trend.gas_1,aL3trend.gas_3,climcapsL3_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1);

disp(' ')
disp('amip6')
amip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,c6trend.stemp,c6trend.ptemp,c6trend.gas_1,c6trend.gas_3,amip6_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_show_olr_ecRad_feedback

