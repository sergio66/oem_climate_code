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

iDefault = 0;
iaComputeWhichFeedback = [9999];  %% do the things for the direct comparisons
iaComputeWhichFeedback = [0];     %% compute + plot feedbacks only  DEFAULT

if iaComputeWhichFeedback ~= iDefault
  disp(' ')
  fprintf(1,'WARNING you have chosen iaComputeWhichFeedback = %5i instead of default %5i (== no computation. just plotting!!!!) \n',iaComputeWhichFeedback,iDefault)
  fprintf(1,'WARNING you have chosen iaComputeWhichFeedback = %5i instead of default %5i (== no computation. just plotting!!!!) \n',iaComputeWhichFeedback,iDefault)
  fprintf(1,'WARNING you have chosen iaComputeWhichFeedback = %5i instead of default %5i (== no computation. just plotting!!!!) \n',iaComputeWhichFeedback,iDefault)
  disp(' ')
end
disp(' ')
iaComputeWhichFeedback
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%

disp('umbc')
if iaComputeWhichFeedback ~= 0
  umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
else
  umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',[]    ,[]    ,[]    ,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
end
% dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
% figure(4); aslprint([dir0 'umbc_feedback_linearity_20yrs_UMBC.pdf'])
figure(4); figure(5); figure(6); figure(80); rett(2); %% disp('RET to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('era5')
if iaComputeWhichFeedback ~= 0
  era5_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,era5.trend_ptemp,era5.trend_gas_1,era5.trend_gas_3,era5_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'ERA5');
else
  era5_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,[]              ,[]              ,[]              ,era5_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'ERA5');
end
figure(4); figure(5); figure(6); figure(80); rett(2); %% disp('RET to continue'); pause

disp(' ')
disp('airsL3')
if iaComputeWhichFeedback ~= 0
  airsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,aL3trend.stemp,aL3trend.ptemp,aL3trend.gas_1,aL3trend.gas_3,airsL3_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'AIRS L3');
else
  airsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,aL3trend.stemp,[]            ,[]             ,[]           ,airsL3_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'AIRS L3');
end
figure(4); figure(5); figure(6); figure(80); rett(2); %% disp('RET to continue'); pause

disp(' ')
disp('cmip6')
if iaComputeWhichFeedback ~= 0
  cmip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,c6trend.stemp,c6trend.ptemp,c6trend.gas_1,c6trend.gas_3,cmip6_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'CMIP6');
else
  cmip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,c6trend.stemp,[]           ,[]           ,[]           ,cmip6_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'CMIP6');
end
figure(4); figure(5); figure(6); figure(80); rett(2); %% disp('RET to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('amip6')
  disp(' ')
  disp('merra2')
  if iaComputeWhichFeedback ~= 0
    merra2_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,merra2.trend_stemp,merra2.trend_ptemp,merra2.trend_gas_1,merra2.trend_gas_3,merra2_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'MERRA2');
  else
    merra2_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,merra2.trend_stemp,[]                ,[]                ,[]                ,merra2_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'MERRA2');
  end
  figure(4); figure(5); figure(6); figure(80); rett(2); %% disp('RET to continue'); pause  

  disp(' ')
  disp('climcapsL3')
  if iaComputeWhichFeedback ~= 0
    climcapsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,cL3trend.stemp,cL3trend.ptemp,cL3trend.gas_1,cL3trend.gas_3,climcapsL3_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'CLIMCAPS L3');
  else
    climcapsL3_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,cL3trend.stemp,[]            ,[]            ,[]            ,climcapsL3_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'CLIMCAPS L3');
  end
  figure(4); figure(5); figure(6); figure(80); rett(2); %% disp('RET to continue'); pause
  
  disp(' ')
  disp('amip6')
  if iaComputeWhichFeedback ~= 0
    amip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,a6trend.stemp,a6trend.ptemp,a6trend.gas_1,a6trend.gas_3,amip6_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'AMIP6');
  else
    amip6_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,a6trend.stemp,[]           ,[]           ,[]           ,amip6_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'AMIP6');
  end
  figure(4); figure(5); figure(6); figure(80); rett(2); %% disp('RET to continue'); pause
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(' getting ready to show results of regressions : polyfit global, robustfit global, robustfit troipics')
disp('   warning : if you do not got MERRA2/CLIMCAPSL3/AMIPS this will eventually fail ... but can still run quick_save_olr_feedbacks_umbc_NWP_L3_XMIP6')

iSmooth = 5;
iSmooth = 10;

disp(' ')
disp('POLYFIT results GLOBAL')
%plot_show_olr_ecRad_feedback_polyfit
plot_show_olr_ecRad_feedback_polyfit_smooth

disp(' ')
disp('GLOBAL(SST) FIT results GLOBAL')
%plot_show_olr_ecRad_feedback_globalsstfit
plot_show_olr_ecRad_feedback_globalsstfit_smooth

disp(' ')
disp('GLOBAL(SST) FIT results TROPICS')
plot_show_olr_ecRad_feedback_globalsstfit_tropics

disp(' ')
disp('ROBUSTFIT results GLOBAL')
%plot_show_olr_ecRad_feedback_robustfit
plot_show_olr_ecRad_feedback_robustfit_smooth

disp(' ')
disp('ROBUSTFIT results TROPICS')
plot_show_olr_ecRad_feedback_robustfit_tropics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(15); clf
subplot(211);
  plot(rlat,umbc_spectral_olr.feedback_ecRad.planck.polyfit_latbin + umbc_spectral_olr.feedback_ecRad.lapse.polyfit_latbin + umbc_spectral_olr.feedback_ecRad.wv.polyfit_latbin,'b',...
       rlat,umbc_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + umbc_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + umbc_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),'g',...
       rlat,umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,'r','linewidth',2)
       ylabel('UMBC \lambda'); plotaxis2; hl = legend('polyfit','robustfit','globalSKT','location','best','fontsize',8); axis([-90 +90 -5 2])

subplot(212);
  plot(rlat,era5_spectral_olr.feedback_ecRad.planck.polyfit_latbin + era5_spectral_olr.feedback_ecRad.lapse.polyfit_latbin + era5_spectral_olr.feedback_ecRad.wv.polyfit_latbin,'b',...
       rlat,era5_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + era5_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + era5_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),'g',...
       rlat,era5_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + era5_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + era5_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,'r','linewidth',2)
       ylabel('ERA5 \lambda'); plotaxis2; hl = legend('polyfit','robustfit','globalSKT','location','best','fontsize',8); axis([-90 +90 -5 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(umbc_spectral_olr,'perts9999')
  disp('Fig 10 checks that atm+skt == atm/skt')
  figure(10); clf
    plot(rlat,umbc_spectral_olr.perts9999.atm_only.polyfit_latbin + umbc_spectral_olr.feedback_ecRad.skt.polyfit_latbin,'mx-',rlat,umbc_spectral_olr.perts9999.atm_skt.polyfit_latbin,'r',...
         rlat,era5_spectral_olr.perts9999.atm_only.polyfit_latbin + era5_spectral_olr.feedback_ecRad.skt.polyfit_latbin,'cx-',rlat,era5_spectral_olr.perts9999.atm_skt.polyfit_latbin,'b','linewidth',2)
    xlim([-90 +90]); plotaxis2; hl = legend('UMBC atm+skt','UMBC','ERA5 atm+skt','ERA5','location','best','fontsize',10); xlabel('Latitude'); ylabel('Feedback \newline W/m2/K'); 

   %%%%%%%%%%%%%%%%%%%%%%%%%
   % see Fig 10
   % An Analytic Model for the Clear-Sky Longwave Feedback Daniel
   % D.B. Kolla , Nadir Jeevanjeeb , Nicholas J. Lutsko, JAS 2023
   % DOI 10.1175/JAS-D-22-0178.1

   figure(11); clf
    subplot(121)
    plot(umbc_spectral_olr.perts9999.atm_skt.polyfit_latbin,rlat,'b','linewidth',6); hold on;
    plot(umbc_spectral_olr.feedback_ecRad.skt.polyfit_latbin,rlat,'k',umbc_spectral_olr.perts9999.atm_only.polyfit_latbin,rlat,'r','linewidth',2); hold off
    ylim([-90 +90]); plotaxis2; hl = legend('atm+skt','skt','atm','location','best','fontsize',8); ylabel('Latitude'); xlabel('Feedback \newline W/m2/K'); title('UMBC')

    subplot(122)
    plot(era5_spectral_olr.perts9999.atm_skt.polyfit_latbin,rlat,'b','linewidth',6); hold on;
    plot(era5_spectral_olr.feedback_ecRad.skt.polyfit_latbin,rlat,'k',era5_spectral_olr.perts9999.atm_only.polyfit_latbin,rlat,'r','linewidth',2); hold off
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

  disp('Fig 12 checks that planck+lapse+wv === atm+skt')
  figure(12); clf
    subplot(211)
    plot(rlat,umbc_spectral_olr.feedback_ecRad.planck.polyfit_latbin,'r',rlat,umbc_spectral_olr.feedback_ecRad.lapse.polyfit_latbin,'k',rlat,umbc_spectral_olr.feedback_ecRad.wv.polyfit_latbin,'g','linewidth',2); hold on
    plot(rlat,umbc_spectral_olr.feedback_ecRad.planck.polyfit_latbin + umbc_spectral_olr.feedback_ecRad.lapse.polyfit_latbin + umbc_spectral_olr.feedback_ecRad.wv.polyfit_latbin,'cx-',...
         rlat,umbc_spectral_olr.perts9999.atm_skt.polyfit_latbin,'b','linewidth',4); hold off
    xlim([-90 +90]); plotaxis2; hl = legend('planck','lapse','wv','planck+lapse+wv','atm+skt','location','best','fontsize',8); xlabel('Latitude'); ylabel('Feedback \newline W/m2/K'); title('UMBC')

    subplot(212)
    plot(rlat,era5_spectral_olr.feedback_ecRad.planck.polyfit_latbin,'r',rlat,era5_spectral_olr.feedback_ecRad.lapse.polyfit_latbin,'k',rlat,era5_spectral_olr.feedback_ecRad.wv.polyfit_latbin,'g','linewidth',2); hold on
    plot(rlat,era5_spectral_olr.feedback_ecRad.planck.polyfit_latbin + era5_spectral_olr.feedback_ecRad.lapse.polyfit_latbin + era5_spectral_olr.feedback_ecRad.wv.polyfit_latbin,'cx-',...
         rlat,era5_spectral_olr.perts9999.atm_skt.polyfit_latbin,'b','linewidth',4); hold off
    xlim([-90 +90]); plotaxis2; hl = legend('planck','lapse','wv','planck+lapse+wv','atm+skt','location','best','fontsize',8); xlabel('Latitude'); ylabel('Feedback \newline W/m2/K'); title('ERA5')
end
