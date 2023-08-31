junk = input('show_olr_ecRad_feedback.m : do you want to clear nwp_spectral_trends_cmip6_era5_airsL3_umbc and nwp_spectral_trends_amip6_merra2_climcapsL3_umbc (-1/+1 [default]) : ');
if length(junk) == 0
  junk = 1;
end
if junk > 0
  clear nwp_spectral_trends_*
end

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

if ~exist('nwp_spectral_trends_amip6_merra2_climcapsL3_umbc') & exist('amip6')
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

iaComputeWhichFeedback = [9999];  %% do the things for the direct comparisons
iaComputeWhichFeedback = [0];     %% compute + plot feedbacks only

disp(' ')
iaComputeWhichFeedback
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%

disp('umbc')
umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
% dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
% figure(4); aslprint([dir0 'umbc_feedback_linearity_20yrs_UMBC.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('era5')
era5_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,era5.trend_ptemp,era5.trend_gas_1,era5.trend_gas_3,era5_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'ERA5');

disp(' ')
disp('airsL3')
airsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,aL3trend.stemp,aL3trend.ptemp,aL3trend.gas_1,aL3trend.gas_3,airsL3_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'AIRS L3');

disp(' ')
disp('cmip6')
cmip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,c6trend.stemp,c6trend.ptemp,c6trend.gas_1,c6trend.gas_3,cmip6_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'CMIP6');

%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('amip6')
  disp(' ')
  disp('merra2')
  merra2_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,merra2.trend_stemp,merra2.trend_ptemp,merra2.trend_gas_1,merra2.trend_gas_3,merra2_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'MERRA2');
  
  disp(' ')
  disp('climcapsL3')
  climcapsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,cL3trend.stemp,cL3trend.ptemp,cL3trend.gas_1,cL3trend.gas_3,climcapsL3_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'CLIMCAPS L3');
  
  disp(' ')
  disp('amip6')
  amip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,a6trend.stemp,a6trend.ptemp,a6trend.gas_1,a6trend.gas_3,amip6_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'AMIP6');
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(' getting ready to show results of regressions : polyfit global, robustfit global, robustfit troipics')
disp('   warning : if you do not got MERRA2/CLIMCAPSL3/AMIPS this will eventually fail ... but can still run quick_save_olr_feedbacks_umbc_NWP_L3_XMIP6')

disp(' ')
disp('POLYFIT results GLOBAL')
plot_show_olr_ecRad_feedback_polyfit

disp(' ')
disp('GLOBAL(SST) FIT results GLOBAL')
plot_show_olr_ecRad_feedback_globalsstfit

disp(' ')
disp('ROBUSTFIT results GLOBAL')
plot_show_olr_ecRad_feedback_robustfit

disp(' ')
disp('ROBUSTFIT results TROPICS')
plot_show_olr_ecRad_feedback_robustfit_tropics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(umbc_spectral_olr,'perts9999')
  figure(10); clf
    plot(rlat,umbc_spectral_olr.feedback9999.atm_only_ecRad_polyfit_latbin + umbc_spectral_olr.feedback_ecRad.skt_ecRad_polyfit_latbin,'mx-',rlat,umbc_spectral_olr.feedback9999.atm_skt_ecRad_polyfit_latbin,'r',...
         rlat,era5_spectral_olr.feedback9999.atm_only_ecRad_polyfit_latbin + era5_spectral_olr.feedback_ecRad.skt_ecRad_polyfit_latbin,'cx-',rlat,era5_spectral_olr.feedback9999.atm_skt_ecRad_polyfit_latbin,'b','linewidth',2)
    xlim([-90 +90]); plotaxis2; hl = legend('UMBC atm+skt','UMBC','ERA5 atm+skt','ERA5','location','best','fontsize',10); xlabel('Latitude'); ylabel('Feedback \newline W/m2/K'); 

   %%%%%%%%%%%%%%%%%%%%%%%%%
   % see Fig 10
   % An Analytic Model for the Clear-Sky Longwave Feedback Daniel
   % D.B. Kolla , Nadir Jeevanjeeb , Nicholas J. Lutsko, JAS 2023
   % DOI 10.1175/JAS-D-22-0178.1

   figure(11); clf
    subplot(121)
    plot(umbc_spectral_olr.feedback9999.atm_skt_ecRad_polyfit_latbin,rlat,'b','linewidth',6); hold on;
    plot(umbc_spectral_olr.feedback_ecRad.skt_ecRad_polyfit_latbin,rlat,'k',umbc_spectral_olr.feedback9999.atm_only_ecRad_polyfit_latbin,rlat,'r','linewidth',2); hold off
    ylim([-90 +90]); plotaxis2; hl = legend('atm+skt','skt','atm','location','best','fontsize',8); ylabel('Latitude'); xlabel('Feedback \newline W/m2/K'); title('UMBC')

    subplot(122)
    plot(era5_spectral_olr.feedback9999.atm_skt_ecRad_polyfit_latbin,rlat,'b','linewidth',6); hold on;
    plot(era5_spectral_olr.feedback_ecRad.skt_ecRad_polyfit_latbin,rlat,'k',era5_spectral_olr.feedback9999.atm_only_ecRad_polyfit_latbin,rlat,'r','linewidth',2); hold off
    ylim([-90 +90]); plotaxis2; hl = legend('atm+skt','skt','atm','location','best','fontsize',8); ylabel('Latitude'); xlabel('Feedback \newline W/m2/K'); title('ERA5')

   %%%%%%%%%%%%%%%%%%%%%%%%%
   % see Fig 3
   % Schneider, D. P., Kay, J. E., & Hannay, C. (2022). Cloud and
   %  surface albedo feedbacks reshape 21st century warming in
   %  successive generations of an Earth System Model. Geophysical
   %  Research Letters, 49,
   %  e2022GL100653. https://doi. org/10.1029/2022GL100653

   % see Fig 1
   % Revisiting the Role of the Water Vapor and Lapse Rate Feedbacks
   %  in the Arctic Amplification of Climate Change EMMA BEER AND IAN
   %  EISENMAN, J. CLim 2022, v35, pg 2975-2988 DOI:
   %  10.1175/JCLI-D-21-0814.1

  figure(12); clf
    subplot(211)
    plot(rlat,umbc_spectral_olr.feedback_ecRad.planck_ecRad_polyfit_latbin,'r',rlat,umbc_spectral_olr.feedback_ecRad.lapse_ecRad_polyfit_latbin,'k',rlat,umbc_spectral_olr.feedback_ecRad.wv_ecRad_polyfit_latbin,'g','linewidth',2); hold on
    plot(rlat,umbc_spectral_olr.feedback_ecRad.planck_ecRad_polyfit_latbin + umbc_spectral_olr.feedback_ecRad.lapse_ecRad_polyfit_latbin + umbc_spectral_olr.feedback_ecRad.wv_ecRad_polyfit_latbin,'cx-',...
         rlat,umbc_spectral_olr.feedback9999.atm_skt_ecRad_polyfit_latbin,'b','linewidth',4); hold off
    xlim([-90 +90]); plotaxis2; hl = legend('planck','lapse','wv','planck+lapse+wv','atm+skt','location','best','fontsize',8); xlabel('Latitude'); ylabel('Feedback \newline W/m2/K'); title('UMBC')

    subplot(212)
    plot(rlat,era5_spectral_olr.feedback_ecRad.planck_ecRad_polyfit_latbin,'r',rlat,era5_spectral_olr.feedback_ecRad.lapse_ecRad_polyfit_latbin,'k',rlat,era5_spectral_olr.feedback_ecRad.wv_ecRad_polyfit_latbin,'g','linewidth',2); hold on
    plot(rlat,era5_spectral_olr.feedback_ecRad.planck_ecRad_polyfit_latbin + era5_spectral_olr.feedback_ecRad.lapse_ecRad_polyfit_latbin + era5_spectral_olr.feedback_ecRad.wv_ecRad_polyfit_latbin,'cx-',...
         rlat,era5_spectral_olr.feedback9999.atm_skt_ecRad_polyfit_latbin,'b','linewidth',4); hold off
    xlim([-90 +90]); plotaxis2; hl = legend('planck','lapse','wv','planck+lapse+wv','atm+skt','location','best','fontsize',8); xlabel('Latitude'); ylabel('Feedback \newline W/m2/K'); title('ERA5')
end

    
