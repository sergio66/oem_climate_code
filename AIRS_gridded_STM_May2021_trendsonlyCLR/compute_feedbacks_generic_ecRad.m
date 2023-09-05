function x_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,trend_skt,trend_ptemp,trend_gas_1,trend_gas_3,x_spectral_olr0,iaComputeWhichFeedback,rlat65,rlon73,iPlotResults,caModelStr)

%% function x_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,trend_skt,trend_ptemp,trend_gas_1,trend_gas_3,x_spectral_olr0,iaComputeWhichFeedback,rlat65,rlon73,iPlotResults,caModelStr)
%%                                                           1 2   3        4           5           6           7             8               9               [10      11     12           13 ]
%% input
%%   h,p       = average profile used for jacobians
%%   results   = results(1:6,1:4608) which will have CO/N2O/CH4 trends in first three
%%   trend_skt   = 1x4608
%%   trend_ptemp = 100x4608
%%   trend_gas_1 = 100x4608
%%   trend_gas_3 = 100x4608
%%
%%   x_spectral_olr0        = "current"values of x_spectral_olr0, so you just overwrite what you want to compute
%%   iaComputeWhichFeedback = use compute_olr+superdriver_run_ecRad_rtp_loop_over_profiles to compute every feedback,   [[[ -1 for all ]]]    or
%%                                                               [[1      2     3   4] [5    6    ] [**9999**           0      ]]           to compute 
%%                                                               [[planck lapse o3 wv] [skt tz/co2] [simple all     LAMBDA ONLY]]
%%                                                                                                  [atm/skt/GHG                ]
%%   ** note 9999 only does ecRad calcs, then returns
%%
%%   [rlat65,rlon73,iPlotRresults] are optional args if you need to plot results
%% note if iaComputeWhichFeedback === 0 then you do not need to run ecRad which means do not need T(z),WV(z),O3(z) trends, only need SKT trends (and potential gas_2,gas_4,gas_6) trends ....
%%  so this is sufficient
%%  umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',[]    ,[]    ,[]    ,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% STARTING FROM SCRATCH  
     for complete list see  make_olr_ecrad_feedback.m

umbc_spectral_olr = struct;    %% so it has no fields
umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,-1,rlat65,rlon73,'UMBC');

era5_spectral_olr = struct;    %% so it has no fields
era5_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,era5.trend_ptemp,era5.trend_gas_1,era5.trend_gas_3,era5_spectral_olr,-1,'ERA5');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% RE-DOING EG planck,wv
     for complete list see  make_olr_ecrad_feedback.m
iaComputeWhichFeedback = [+1 +4]; %% compute planck (uniform pert) and wv ecRads
iaComputeWhichFeedback = [+1];    %% compute planck (uniform pert) ecRads
iaComputeWhichFeedback = [+1 +2]; %% compute planck (uniform pert) and lapse ecRads
iaComputeWhichFeedback = [+9999]; %% compute all (atm/skt/ghg) ecRads
umbc_spectral_olr    = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,'UMBC');
era5_spectral_olr    = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,era5.trend_ptemp,era5.trend_gas_1,era5.trend_gas_3,era5_spectral_olr,iaComputeWhichFeedback,'ERA5');
airsL3_spectral_olr  = compute_feedbacks_generic_ecRad(h,p,results,aL3trend.stemp,aL3trend.ptemp,aL3trend.gas_1,aL3trend.gas_3,airsL3_spectral_olr,iaComputeWhichFeedback,'AIRS L3');
cmip6_spectral_olr   = compute_feedbacks_generic_ecRad(h,p,results,c6trend.stemp,c6trend.ptemp,c6trend.gas_1,c6trend.gas_3,cmip6_spectral_olr,iaComputeWhichFeedback,'CMIP6');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% RE-DOING EG regressions and plots 
  for complete list see      show_olr_ecRad_feedback.m        
iaComputeWhichFeedback = [0];     %% compute + plot feedbacks only
disp(' ')
disp('umbc')
umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',[]    ,[]    ,[]    ,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');

disp(' ')
disp('era5')
era5_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,era5.trend_ptemp,era5.trend_gas_1,era5.trend_gas_3,era5_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'ERA5');
era5_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,era5.trend_stemp,[]              ,[]              ,[]               era5_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'ERA5');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 9
  error('need at least 9 args')
end

%% see “Simpson's Law” and the Spectral Cancellation of Climate Feedbacks
%% Nadir Jeevanjee1 , Daniel D. B. Koll2, and Nicholas Lutsko3
%% GRL eevanjee, N., Koll, D. D. B., & Lutsko, N. (2021). “Simpson's Law” and
%%  the spectral cancellation of climate feedbacks. Geophysical Research Letters, 48, e2021GL093699. https://doi. org/10.1029/2021GL093699

if length(fieldnames(x_spectral_olr0)) > 0
  x_spectral_olr = x_spectral_olr0;
  if length(intersect(iaComputeWhichFeedback,[0])) == 1
    disp('only redoing the final feedback calcs, no need for SARTA or ecRad, fast')
  else
    disp('overwriting input fields as needed, using SARTA and ecRad')
  end
else
  disp('starting feedback calcs afresh')
end

deltaST = trend_skt;     [mmbad,nnbad] = find(isinf(deltaST) | isnan(deltaST) > 1 | (isinf(deltaST) > 1));  [bad] = find(isinf(deltaST) | isnan(deltaST) > 1 | (isinf(deltaST) > 1));  deltaST(bad) = 0;
deltaT  = trend_ptemp;   [mmbad,nnbad] = find(isinf(deltaT) | isnan(deltaT) > 1 | (isinf(deltaT) > 1));     [bad] = find(isinf(deltaT) | isnan(deltaT) > 1 | (isinf(deltaT) > 1));     deltaT(bad) = 0;
fracWV  = trend_gas_1;   [mmbad,nnbad] = find(isinf(fracWV) | isnan(fracWV) > 1 | (isinf(fracWV) > 1));     [bad] = find(isinf(fracWV) | isnan(fracWV) > 1 | (isinf(fracWV) > 1));     fracWV(bad) = 0;
fracO3  = trend_gas_3;   [mmbad,nnbad] = find(isinf(fracO3) | isnan(fracO3) > 1 | (isinf(fracO3) > 1));     [bad] = find(isinf(fracO3) | isnan(fracO3) > 1 | (isinf(fracO3) > 1));     fracO3(bad) = 0;

cdRRTMback = ['cd ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iLambda_UseGlobalSST_regress = +1;   %% new way, use global avg SST, much safer (less likely to be 0)
iLambda_UseGlobalSST_regress = -1;   %% old way till March 2022, using computed SST per tile instead of glabal average; dangerous because if dSST = 0 oopsie when dividing

addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RobinHoganECMWF/ECRAD_ECMWF_version_of_flux/ecRad/create_ecrad_inputSergio/
addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/DRIVER_CODE_RRTM_Band17/

%% trying a faster method see     jaja = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);    driver_run_ecRad_rtp_loop_over_profiles.m iMethod = -1
%% scatter_coast(p.rlon,p.rlat,50,(x_spectral_olr.olr0_ecRad.clr-jaja.clr)./x_spectral_olr.olr0_ecRad.clr*100); colormap(usa2); caxis([-1 +1])
%% scatter_coast(p.rlon,p.rlat,50,(x_spectral_olr.olr0_rrtm-jaja.clr)./x_spectral_olr.olr0_rrtm*100); colormap(usa2); caxis([-1 +1])

%% iaComputeWhichFeedback == 0 ==> just in case there are updates to eg compute_feedbacks_ecRad_calcs.m, compute_feedbacks_regress_ecRad_calcs.m

if iaComputeWhichFeedback == -1
  %% if iaComputeWhichFeedback == =-1 then you are recomputing ALL SIX feedbacks, so may as well re-do base OLR calc
  disp('computing BASE OLR')
  px = p;
  x_spectral_olr.olr0 = compute_olr(h,px);
  %x_spectral_olr.olr0_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0); %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<
  x_spectral_olr.olr0_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
  eval(cdRRTMback);
end

indSST    = deltaST;
globalSST = nanmean(indSST);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iDoThis = -1;
if length(intersect(iaComputeWhichFeedback,9999)) == +1
  iDoThis = +1;
end
if iDoThis > 0
  %% all perts
  disp('perturbing Tsurf, T(z),WV(z),O3(z)) in one go, to compute feedback in one gulp!!!! and then exit code')

  px = p;
  
  %px.gas_2 = px.gas_2*(1+2.2/400);
  %if iLambda_UseGlobalSST == -1
  %  px.stemp = px.stemp + indSST;
  %else
  %  px.stemp = px.stemp + globalSST;
  %end
  
  px.gas_2 = px.gas_2 .* (ones(101,1)*(1+results(:,1)'/400));
  px.gas_4 = px.gas_4 .* (ones(101,1)*(1+results(:,2)'/300));
  px.gas_6 = px.gas_6 .* (ones(101,1)*(1+results(:,3)'/1800));
  px.stemp = px.stemp + indSST;
  px.ptemp(1:100,:) = px.ptemp(1:100,:) + trend_ptemp(1:100,:);
  fracJUNK = trend_gas_1(1:100,:); bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
    px.gas_1(1:100,:) = px.gas_1(1:100,:) .* (1 + fracJUNK);
  fracJUNK = trend_gas_3(1:100,:); bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
    px.gas_3(1:100,:) = px.gas_3(1:100,:) .* (1 + fracJUNK);
  x_spectral_olr.perts9999.atm_skt_ghg_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
  eval(cdRRTMback);
  
  boo0 = -(x_spectral_olr.perts9999.atm_skt_ghg_ecRad.clr - x_spectral_olr.olr0_ecRad.clr);
  scatter_coast(p.rlon,p.rlat,100,boo0./indSST); colormap jet
  junk0 = polyfit(indSST,boo0,1);
  [nn,nx,ny,nmean,nstd] = myhist2d(indSST,boo0,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
  str1 = ['atm/skt + GHG \newline d(OLR) = ' num2str(junk0(1)) ' d(SST) + ' num2str(junk0(2))];
  xlabel('dSST'); ylabel('d(OLR)'); 
  title(str1); fprintf(1,'%s \n',str1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% no tracegas perts
  
  px = p;
  
  %px.gas_2 = px.gas_2*(1+2.2/400);
  %if iLambda_UseGlobalSST == -1
  %  px.stemp = px.stemp + indSST;
  %else
  %  px.stemp = px.stemp + globalSST;
  %end
  
  %px.gas_2 = px.gas_2 .* (ones(101,1)*(1+results(:,1)'/400));
  %px.gas_4 = px.gas_4 .* (ones(101,1)*(1+results(:,2)'/300));
  %px.gas_6 = px.gas_6 .* (ones(101,1)*(1+results(:,3)'/1800));
  px.stemp = px.stemp + indSST;
  px.ptemp(1:100,:) = px.ptemp(1:100,:) + trend_ptemp(1:100,:);
  fracJUNK = trend_gas_1(1:100,:); bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
    px.gas_1(1:100,:) = px.gas_1(1:100,:) .* (1 + fracJUNK);
  fracJUNK = trend_gas_3(1:100,:); bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
    px.gas_3(1:100,:) = px.gas_3(1:100,:) .* (1 + fracJUNK);
  x_spectral_olr.perts9999.atm_skt_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
  eval(cdRRTMback);
  
  boo1 = -(x_spectral_olr.perts9999.atm_skt_ecRad.clr - x_spectral_olr.olr0_ecRad.clr);
  scatter_coast(p.rlon,p.rlat,100,boo1./indSST); colormap jet
  junk1 = polyfit(indSST,boo1,1);
  [nn,nx,ny,nmean,nstd] = myhist2d(indSST,boo1,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
  str2 = ['atm/skt only, no GHG \newline d(OLR) = ' num2str(junk1(1)) ' d(SST) + ' num2str(junk1(2))];
  xlabel('dSST'); ylabel('d(OLR)'); 
  title(str2); fprintf(1,'%s \n',str2);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% no tracegas perts, no SST perts ie just atmosphere
  
  px = p;
  
  %px.gas_2 = px.gas_2*(1+2.2/400);
  %if iLambda_UseGlobalSST == -1
  %  px.stemp = px.stemp + indSST;
  %else
  %  px.stemp = px.stemp + globalSST;
  %end
  
  %px.gas_2 = px.gas_2 .* (ones(101,1)*(1+results(:,1)'/400));
  %px.gas_4 = px.gas_4 .* (ones(101,1)*(1+results(:,2)'/300));
  %px.gas_6 = px.gas_6 .* (ones(101,1)*(1+results(:,3)'/1800));
  %px.stemp = px.stemp + indSST;
  px.ptemp(1:100,:) = px.ptemp(1:100,:) + trend_ptemp(1:100,:);
  fracJUNK = trend_gas_1(1:100,:); bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
    px.gas_1(1:100,:) = px.gas_1(1:100,:) .* (1 + fracJUNK);
  fracJUNK = trend_gas_3(1:100,:); bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
    px.gas_3(1:100,:) = px.gas_3(1:100,:) .* (1 + fracJUNK);
  x_spectral_olr.perts9999.atm_only_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
  eval(cdRRTMback);
  
  boo1 = -(x_spectral_olr.perts9999.atm_only_ecRad.clr - x_spectral_olr.olr0_ecRad.clr);
  scatter_coast(p.rlon,p.rlat,100,boo1./indSST); colormap jet
  junk1 = polyfit(indSST,boo1,1);
  [nn,nx,ny,nmean,nstd] = myhist2d(indSST,boo1,-0.15:0.01:+0.15,-0.4:0.05:+0.4);
  str2 = ['atm only, no SKT/GHG \newline d(OLR) = ' num2str(junk1(1)) ' d(SST) + ' num2str(junk1(2))];
  xlabel('dSST'); ylabel('d(OLR)'); 
  title(str2); fprintf(1,'%s \n',str2);

  return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Jeevanjee, N., Koll, D. D. B., & Lutsko, N. (2021). “Simpson's Law” and
%%% the spectral cancellation of climate  feedbacks. Geophysical Research
%%% Letters, 48, e2021GL093699. https://doi.org/10.1029/2021GL093699

%%% Assuming linearity in the finite differences, we then have the total feedback is summ of following three : 

% The Planck feedback is minus the (ΔTs-normalized) OLR response to a uniform change in surface and
% atmospheric temperatures, with qv held fixed at the initial profile.
% SO THIS IS UNIFORM PERTURBATION RESPONSE
if length(intersect(iaComputeWhichFeedback,[-1 1])) == 1
  disp('perturbing T(z),ST uniformly for planck OLR feedback')
  px = p;
  %px.gas_2 = px.gas_2*(1+2.2/400);
  %if iLambda_UseGlobalSST_regress == -1
    %%% REMEMBER ultimately EVERYTHING is regresse/deivided by indSST so use THIS as uniform perturbation!!!!
    px.stemp = px.stemp + indSST;
    px.ptemp = px.ptemp + ones(101,1)*indSST;
  %else
  %  px.stemp = px.stemp + globalSST;
  %  px.ptemp = px.ptemp + ones(101,4608)*globalSST;
  %end
  x_spectral_olr.planck = compute_olr(h,px);
  %x_spectral_olr.planck_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0);
  x_spectral_olr.planck_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
  eval(cdRRTMback);
end

% The lapse-rate feedback is minus the OLR response to the
% difference between the actual temperature response and the uniform
% Planck response, still holding qv fixed.
% SO HERE I COMPUTE ACTUAL RESPONSE TO T(z) and SKT, AND LATER SUBTRACT OUT THE UNIFORM PLANCK RESPONSE
if length(intersect(iaComputeWhichFeedback,[-1 2])) == 1
  disp('perturbing T(z),ST for lapse rate OLR feedback')
  px = p;
  %px.gas_2 = px.gas_2*(1+2.2/400);
  if iLambda_UseGlobalSST_regress == -1
    %% need the final state
    px.stemp = px.stemp + indSST;
    px.ptemp(1:100,:) = px.ptemp(1:100,:) + deltaT(1:100,:);
  else
    %% need the final state
    px.stemp = px.stemp + indSST;
    px.ptemp(1:100,:) = px.ptemp(1:100,:) + deltaT(1:100,:);
  end
  x_spectral_olr.lapse = compute_olr(h,px);   
  %x_spectral_olr.lapse_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0);
  x_spectral_olr.lapse_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
  eval(cdRRTMback);
end

% The water vapor feedback λwv is then minus the OLR response to the change in qv,
% holding temperatures fixed
if length(intersect(iaComputeWhichFeedback,[-1 4])) == 1
  disp('perturbing WV(z) for wv OLR feedback')
  px = p;
  %px.gas_2 = px.gas_2*(1+2.2/400);
  %% px.gas_1 = px.gas_1 .* (1 + fracWV);
  fracJUNK = fracWV(1:100,:); bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
  if iLambda_UseGlobalSST_regress == -1
    px.gas_1(1:100,:) = px.gas_1(1:100,:) .* (1 + fracJUNK);
  else
    px.gas_1(1:100,:) = px.gas_1(1:100,:) .* (1 + fracJUNK);
  end
  x_spectral_olr.wv = compute_olr(h,px);
  %x_spectral_olr.wv_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0);
  x_spectral_olr.wv_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
  eval(cdRRTMback);
end

%%%%%%%%%% 
if length(intersect(iaComputeWhichFeedback,[-1 3])) == 1
  disp('perturbing O3(z) for ozone OLR feedback')
  px = p;
  fracJUNK = fracO3(1:100,:); bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
  if iLambda_UseGlobalSST_regress == -1
    px.gas_3(1:100,:) = px.gas_3(1:100,:) .* (1 + fracJUNK);
  else
    px.gas_3(1:100,:) = px.gas_3(1:100,:) .* (1 + fracJUNK);
  end
  x_spectral_olr.o3 = compute_olr(h,px);
  %x_spectral_olr.o3_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0);
  x_spectral_olr.o3_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
  eval(cdRRTMback);
end

if length(intersect(iaComputeWhichFeedback,[-1 5])) == 1
  disp('perturbing SKT for OLR')
  px = p;
  %px.gas_2 = px.gas_2*(1+2.2/400);
  if iLambda_UseGlobalSST_regress == -1
    px.stemp = px.stemp + indSST;
  else
    px.stemp = px.stemp + globalSST;
  end
  x_spectral_olr.skt = compute_olr(h,px);
  %x_spectral_olr.skt_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0);
  x_spectral_olr.skt_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
  eval(cdRRTMback);
end
  
if length(intersect(iaComputeWhichFeedback,[-1 6])) == 1
  disp('perturbing CO2 and T(z) for OLR feedback')
  px = p;
  px.gas_2 = px.gas_2*(1+2.2/400);
  if iLambda_UseGlobalSST_regress == -1
    %% need the final state
    px.ptemp(1:100,:) = px.ptemp(1:100,:) + deltaT(1:100,:);
  else
    %% need the final state
    px.ptemp(1:100,:) = px.ptemp(1:100,:) + deltaT(1:100,:);
  end
  x_spectral_olr.ptemp_co2 = compute_olr(h,px);   
  %x_spectral_olr.ptemp_co2_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0);
  x_spectral_olr.ptemp_co2_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
  eval(cdRRTMback);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% COMPUTE INCOMPLETE FEEDBACKS FROM SARTA %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% remember this "spectral_olr" is only 645-2780 cm-1, at AIRS res!!!!!!
%% remember this "spectral_olr" is only 645-2780 cm-1, at AIRS res!!!!!!
%% remember this "spectral_olr" is only 645-2780 cm-1, at AIRS res!!!!!!

%% note Eq 8b of Jevanjee paper shows we need to use 
%%   x_spectral_olr.feedback_ecRad.lapse == x_spectral_olr.lapse - x_spectral_olr.planck
%% and not
%%   x_spectral_olr.feedback_ecRad.lapse == x_spectral_olr.lapse - x_spectral_olr.olr0

%% change radiance mW --> W and then multiply by pi for flux
ix1 = 1:2162; ix2 = 2163:2645;  %% basically have two bands of detectors!

x_spectral_olr = compute_feedbacks_regress_olr_sarta_calcs(x_spectral_olr,deltaST,iLambda_UseGlobalSST_regress,h);

if nargin <= 12  
  caModelStr = 'XYZ';
end
figure(71); scatter_coast(p.rlon,p.rlat,50,x_spectral_olr.feedback_sarta.planck); caxis([-4 0]*1);  colormap(jet);  title([caModelStr ' \lambda_{Planck} sarta'])
figure(72); scatter_coast(p.rlon,p.rlat,50,x_spectral_olr.feedback_sarta.lapse);  caxis([-5 +5]*1); colormap(usa2); title([caModelStr ' \lambda_{Lapse} sarta'])
figure(73); scatter_coast(p.rlon,p.rlat,50,x_spectral_olr.feedback_sarta.wv);     caxis([-2 +2]*1); colormap(usa2); title([caModelStr ' \lambda_{WV} sarta'])
figure(74); scatter_coast(p.rlon,p.rlat,50,x_spectral_olr.feedback_sarta.skt);    caxis([-2 0]*1);  colormap(jet);  title([caModelStr ' \lambda_{Skt} sarta'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% COMPUTE COMPLETE FEEDBACKS FROM ECRAD %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_spectral_olr = compute_feedbacks_regress_olr_ecRad_calcs(x_spectral_olr,deltaST,iLambda_UseGlobalSST_regress,caModelStr);

%% planck_ecRad ---> savenums(1,:) ./ indSST
%% lapse_ecRad ---> savenums(2,:)  ./ indSST
%% wv_ecRad ---> savenums(4,:)     ./ indSST
%% skt_ecRad ---> savenums(5,:)    ./ indSST

junkSKT = x_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
figure(75); scatter_coast(p.rlon,p.rlat,50,x_spectral_olr.feedback_ecRad.savenums(1,:)/junkSKT); caxis([-4 0]*1.5); colormap(jet);  title([caModelStr ' \lambda_{Planck} ecRad'])
figure(76); scatter_coast(p.rlon,p.rlat,50,x_spectral_olr.feedback_ecRad.savenums(2,:)/junkSKT); caxis([-5 +5]*2);  colormap(usa2); title([caModelStr ' \lambda_{Lapse} ecRad'])
figure(77); scatter_coast(p.rlon,p.rlat,50,x_spectral_olr.feedback_ecRad.savenums(4,:)/junkSKT); caxis([-2 +2]*2);  colormap(usa2); title([caModelStr ' \lambda_{WV} ecRad'])
figure(78); scatter_coast(p.rlon,p.rlat,50,x_spectral_olr.feedback_ecRad.savenums(5,:)/junkSKT); caxis([-2 0]*1);   colormap(jet);  title([caModelStr ' \lambda_{Skt} ecRad'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  DONE, EITHER EXIT OR MAKE PLOTS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 12
  iPlotResults = -1;
end

if iPlotResults > 0
  figsmap = colormap_soden_held_jclim2007;
  figsmap = usa2;

  %% rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = rlon73;  rlat = rlat65;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

  [Y,X] = meshgrid(rlat,rlon);
  X = X; Y = Y;
  YY = Y(:)';                     %% mean weighted delta SST rate = 0.031962  0.003305  0.023971  0.019594 K/yr for 05/10/15/20 years   WRONG
  YY = Y'; YY = YY(:); YY = YY';  %% mean weighted delta SST rate = 0.069633  0.020002  0.028442  0.024870 K/yr for 05/10/15/20 years   CORRECT

  figure(75); colormap(figsmap); caxis([-6 -2])   %% planck
  figure(77); colormap(figsmap); caxis([-3 +3])   %% wv

  maskLF = ones(1,4608);
  maskLFmatr = reshape(maskLF,72,64)';

  ns = 500;
  ns = 50;
  junkSKT = x_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;    
  wonk = x_spectral_olr.feedback_ecRad.savenums(1,:)/junkSKT; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
    aslmap(75,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);  colormap(jet);  caxis([-4 0]*1.5);  title([caModelStr '  \lambda_{Planck} GLOBAL SKT']);
  wonk = x_spectral_olr.feedback_ecRad.savenums(2,:)/junkSKT; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
    aslmap(76,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);  colormap(usa2); caxis([-5 5]*1);    title([caModelStr '  \lambda_{Lapse} GLOBAL SKT'])
  wonk = x_spectral_olr.feedback_ecRad.savenums(3,:)/junkSKT; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
    aslmap(77,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);  colormap(usa2); caxis([-2 2]*1);    title([caModelStr '  \lambda_{O3} GLOBAL SKT'])
  wonk = x_spectral_olr.feedback_ecRad.savenums(4,:)/junkSKT; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
    aslmap(78,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);  colormap(usa2); caxis([-2 2]*1);    title([caModelStr '  \lambda_{WV} GLOBAL SKT'])
  wonk = x_spectral_olr.feedback_ecRad.savenums(5,:)/junkSKT; wonk(wonk < -3) = NaN; wonk(wonk > +3) = NaN; 
    aslmap(79,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);  colormap(jet);  caxis([-2 0]*1.5);  title([caModelStr '  \lambda_{Skt} GLOBAL SKT'])
  wonk = x_spectral_olr.feedback_ecRad.savenums(1,:)/junkSKT + x_spectral_olr.feedback_ecRad.savenums(2,:)/junkSKT + x_spectral_olr.feedback_ecRad.savenums(4,:)/junkSKT; 
    wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
    aslmap(80,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);  colormap(jet);  caxis([-2 0]*1.5);  title([caModelStr '  \lambda_{ALL} GLOBAL SKT'])
  
  figure(75); colormap(figsmap); caxis([-10 0])    %% planck
  figure(76); colormap(figsmap); caxis([-3 +3])    %% lapse
  figure(77); colormap(figsmap); caxis([-1 +1])    %% o3
  figure(78); colormap(figsmap); caxis([-5 +5])    %% wv
  figure(79); colormap(figsmap); caxis([-3 0])     %% skt
  figure(80); colormap(figsmap); caxis([-4 +4])    %% all
  
  figure(75); colormap(figsmap); caxis([-6 -2])     %% planck
  figure(78); colormap(figsmap); caxis([-3 +3])     %% wv
  
  if iLambda_UseGlobalSST_regress == -1
    bad = find(abs(indSST) < 1e-4);
  else
    bad = [];
  end
  junklat = reshape(p.rlat,72,64); junklat = mean(junklat,1);
  figure(81); clf
    junk = x_spectral_olr.feedback_ecRad.savenums(1,:)/junkSKT; junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'r','linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(2,:)/junkSKT; junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'c','linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(4,:)/junkSKT; junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'b','linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(5,:)/junkSKT; junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'g','linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(1,:)/junkSKT + x_spectral_olr.feedback_ecRad.savenums(2,:)/junkSKT + x_spectral_olr.feedback_ecRad.savenums(3,:)/junkSKT + x_spectral_olr.feedback_ecRad.savenums(4,:)/junkSKT; 
      junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'kx-','linewidth',4); hold off
    ylim([-10 +10]/2); plotaxis2; title([caModelStr ' \lambda \newline (globalSKT/massaged for no NaNs etc)']); xlabel('Latitude'); hl = legend('Planck','Lapse','WV','SKT','ALL','location','best','fontsize',8);
  figure(81); clf  %% to match colors in Figs 4,5,6 of compute_feedbacks_regress_olr_ecRad_calcs
    junk = x_spectral_olr.feedback_ecRad.savenums(1,:)/junkSKT; junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(2,:)/junkSKT; junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(3,:)/junkSKT; junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(4,:)/junkSKT; junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(5,:)/junkSKT; junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(1,:)/junkSKT + x_spectral_olr.feedback_ecRad.savenums(2,:)/junkSKT + x_spectral_olr.feedback_ecRad.savenums(3,:)/junkSKT + x_spectral_olr.feedback_ecRad.savenums(4,:)/junkSKT; 
      junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'kx-','linewidth',4); hold off
    ylim([-10 +10]/2); plotaxis2; title([caModelStr ' \lambda \newline (globalSKT/massaged for no NaNs etc)']); xlabel('Latitude'); hl = legend('Planck','Lapse','O3','WV','SKT','ALL','location','best','fontsize',8);

  figure(82); clf  %% to match colors in Figs 4,5,6 of compute_feedbacks_regress_olr_ecRad_calcs
    junk = x_spectral_olr.feedback_ecRad.savenums(1,:)/junkSKT; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(2,:)/junkSKT; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(3,:)/junkSKT; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(4,:)/junkSKT; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(5,:)/junkSKT; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'linewidth',2); hold on
    junk = x_spectral_olr.feedback_ecRad.savenums(1,:)/junkSKT + x_spectral_olr.feedback_ecRad.savenums(2,:)/junkSKT + x_spectral_olr.feedback_ecRad.savenums(3,:)/junkSKT + x_spectral_olr.feedback_ecRad.savenums(4,:)/junkSKT; 
      junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'kx-','linewidth',4); hold off
    ylim([-10 +10]/2); plotaxis2; title([caModelStr ' \lambda (globalSKT RAW)']); xlabel('Latitude'); hl = legend('Planck','Lapse','O3','WV','SKT','ALL','location','best','fontsize',8);
  
  junk =  x_spectral_olr.feedback_ecRad.savenums(1,:)/junkSKT; junk(bad) = NaN;
  if ~exist('lps0')
    [mmw0,lps0] = mmwater_rtp_pstop_lapse(h,p);
  end
  figure(83); clf; plot(mean(lps0.lapse1(80:97,:),1),junk,'.');
  figure(83); clf; scatter(mean(lps0.lapse1(80:97,:),1),junk,10,p.rlat,'filled');      colorbar; colormap jet; axis([0 10 -5 -3])
  figure(83); clf; scatter(mean(lps0.lapse1(80:97,:),1),junk,10,abs(p.rlat),'filled'); colorbar; colormap jet; axis([0 10 -5 -3])
    xlabel('Lower Trop Lapse Rate K/km'); ylabel('\lambda_{Planck}');
  addpath /home/sergio/MATLABCODE/SHOWSTATS
  [n,nx,ny,nmean,nstd] = myhist2d(mean(lps0.lapse1(80:97,:),1),junk,0:0.25:10,-6:0.1:-3); errorbar(0:0.25:10,nmean,nstd); xlabel('Lower Trop Lapse Rate K/km'); ylabel('\lambda_{Planck}'); plotaxis2;

  figure(84); clf
  scatter(YY,trend_skt,10,YY,'filled'); colorbar; title([caModelStr]); xlabel('latitude'); ylabel('\delta SKT/yr')
  plotaxis2;

  pause(1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
