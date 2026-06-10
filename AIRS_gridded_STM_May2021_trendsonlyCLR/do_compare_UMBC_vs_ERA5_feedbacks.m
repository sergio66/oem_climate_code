function do_compare_UMBC_vs_ERA5_feedbacks(umbc_spectral_olr,era5_spectral_olr,umbc_skt_trends,era5_skt_trends);

%umbc_spectral_olr
%era5_spectral_olr

% do_compare_UMBC_vs_ERA5_feedbacks(umbc_spectral_olr,era5_spectral_olr,results(:,6),era5.trend_stemp);
%
% or
%
% do_compare_UMBC_vs_ERA5_feedbacks(umbc_spectral_delta_olr,era5_spectral_delta_olr,results(:,6),era5.trend_stemp);

do_XX_YY_from_X_Y
rlat = reshape(Y,72,64);
rlat = mean(rlat,1);

nsmooth = 5;
nsmooth = 11;
nsmooth = 15;

aslmap(101,rlat65,rlon73,smoothn((reshape(umbc_skt_trends',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(usa2); title('UMBC SKT trends')
aslmap(102,rlat65,rlon73,smoothn((reshape(era5_skt_trends',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(usa2); title('ERA5 SKT trends');
if isfield(umbc_spectral_olr,'feedback_ecRad')
  fprintf(1,'globally cosing avg SKT trends %.4f %.4f K/yr for UMBC and ERA5 \n',umbc_spectral_olr.feedback_ecRad.global_coslat_wgt_skt,era5_spectral_olr.feedback_ecRad.global_coslat_wgt_skt)
end

aslmap(103,rlat65,rlon73,smoothn((reshape(umbc_spectral_olr.olr0_ecRad.clr',72,64)') ,1), [-90 +90],[-180 +180]); colormap(jet); title('UMBC and ERA5 OLR 0')

%% note : see compute_feedbacks_generic_ecRad.m
%% if length(intersect(iaComputeWhichFeedback,9999)) == +1
%%   atm_skt_ghg_ecRad  is everything in one gulp : perturbing Tsurf,T(z),WV(z),O3(z),tracegas
%%   atm_skt_ecRad      is only perturbing SKT, T, WV, O3 so avoid it
%%   atm_only_ecRad     is only perturbing      T, WV, O3 so avoid it
%%   atmT_only_ecRad    is only perturbing      T
%%   ghg_only_ecRad     is only perturbing ghg
%% end
%% then   iaComputeWhichFeedback controls egs skt, t, wv, o3, lapse


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(umbc_spectral_olr,'perts8888') & isfield(era5_spectral_olr,'perts8888')
  wahU = umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.perts8888.atm_skt_ghg_ecRad.clr;  wahU = -wahU;
  wahE = era5_spectral_olr.olr0_ecRad.clr-era5_spectral_olr.perts8888.atm_skt_ghg_ecRad.clr;  wahE = -wahE;
  figure(104); clf; plot(1:4608,wahU,'b',1:4608,wahE,'r'); title('OLR trends (b) sarta (r) ecrad')
  
  wahU = reshape(umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.perts8888.atm_skt_ghg_ecRad.clr,72,64); wahU = -nanmean(wahU,1);
  wahE = reshape(era5_spectral_olr.olr0_ecRad.clr-era5_spectral_olr.perts8888.atm_skt_ghg_ecRad.clr,72,64); wahE = -nanmean(wahE,1);
  figure(105); clf; plot(rlat,wahU,'b',rlat,wahE,'r','linewidth',2); plotaxis2; legend('UMBC','ERA5','location','best') ; ylabel('dOLR/dt W/m2/yr'); 
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('plotting estimates of climate feedbacks ....')
dST = 0.0175; %% global surface warming K/yr
figure(104); clf; plot(1:4608,(nansum(umbc_spectral_olr.olr0,1)-nansum(umbc_spectral_olr.o3,1))/1e3*0.0025*pi*240/dST,'b',1:4608,(umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.o3_ecRad.clr)/dST,'r'); title('o3 (b) sarta (r) ecrad')
figure(105); clf; plot(1:4608,(nansum(umbc_spectral_olr.olr0,1)-nansum(umbc_spectral_olr.wv,1))/1e3*0.0025*pi*240/dST,'b',1:4608,(umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.wv_ecRad.clr)/dST,'r'); title('wv (b) sarta (r) ecrad')
figure(106); clf; plot(1:4608,(nansum(umbc_spectral_olr.olr0,1)-nansum(umbc_spectral_olr.lapse,1))/1e3*0.0025*pi*240/dST,'b',1:4608,(umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.lapse_ecRad.clr)/dST,'r'); title('lapse (b) sarta (r) ecrad')
figure(107); clf; plot(1:4608,(nansum(umbc_spectral_olr.olr0,1)-nansum(umbc_spectral_olr.ptemp_co2,1))/1e3*0.0025*pi*240/dST,'b',1:4608,(umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.ptemp_co2_ecRad.clr)/dST,'r'); title('tz (b) sarta (r) ecrad')
figure(108); clf; plot(1:4608,(nansum(umbc_spectral_olr.olr0,1)-nansum(umbc_spectral_olr.planck,1))/1e3*0.0025*pi*240/dST,'b',1:4608,(umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.planck_ecRad.clr)/dST,'r'); title('planck (b) sarta (r) ecrad')
figure(109); clf; plot(1:4608,(nansum(umbc_spectral_olr.olr0,1)-nansum(umbc_spectral_olr.skt,1))/1e3*0.0025*pi*240/dST,'b',1:4608,(umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.skt_ecRad.clr)/dST,'r'); title('skt (b) sarta (r) ecrad')

disp('ret to continue'); pause
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aslmap(104,rlat65,rlon73,smoothn((reshape(umbc_spectral_olr.feedback_ecRad.wv.individual',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*6); colormap(usa2); title('UMBC WV feedback')
aslmap(105,rlat65,rlon73,smoothn((reshape(era5_spectral_olr.feedback_ecRad.wv.individual',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 +1]*6); colormap(usa2); title('ERA5 WV feedback');

aslmap(106,rlat65,rlon73,smoothn((reshape(umbc_spectral_olr.feedback_ecRad.planck.individual',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 0]*5); colormap(usa2); title('UMBC planck feedback')
aslmap(107,rlat65,rlon73,smoothn((reshape(era5_spectral_olr.feedback_ecRad.planck.individual',72,64)') ,1), [-90 +90],[-180 +180]); caxis([-1 0]*5); colormap(usa2); title('ERA5 planck feedback');

figure(108)
blah1 = reshape(umbc_spectral_olr.feedback_ecRad.wv.individual',72,64); blah1 = nanmean(blah1,1);
blah2 = reshape(era5_spectral_olr.feedback_ecRad.wv.individual',72,64); blah2 = nanmean(blah2,1);
plot(rlat,smooth(blah1,nsmooth),rlat,smooth(blah2,nsmooth),'linewidth',2); plotaxis2; legend('umbc','era5'); title('WV feedback')
ylim([-1 5])

figure(109)
blah1 = reshape(umbc_spectral_olr.feedback_ecRad.planck.individual',72,64); blah1 = nanmean(blah1,1);
blah2 = reshape(era5_spectral_olr.feedback_ecRad.planck.individual',72,64); blah2 = nanmean(blah2,1);
plot(rlat,smooth(blah1,nsmooth),rlat,smooth(blah2,nsmooth),'linewidth',2); plotaxis2; legend('umbc','era5'); title('Planck feedback')
ylim([-20 10])

figure(110)
blah1 = reshape(umbc_spectral_olr.feedback_ecRad.lapse.individual',72,64); blah1 = nanmean(blah1,1);
blah2 = reshape(era5_spectral_olr.feedback_ecRad.lapse.individual',72,64); blah2 = nanmean(blah2,1);
plot(rlat,smooth(blah1,nsmooth),rlat,smooth(blah2,nsmooth),'linewidth',2); plotaxis2; legend('umbc','era5'); title('Lapse feedback')
ylim([-5 10])

figure(111)
blah1 = reshape(umbc_spectral_olr.feedback_ecRad.o3.individual',72,64); blah1 = nanmean(blah1,1);
blah2 = reshape(era5_spectral_olr.feedback_ecRad.o3.individual',72,64); blah2 = nanmean(blah2,1);
plot(rlat,smooth(blah1,nsmooth),rlat,smooth(blah2,nsmooth),'linewidth',2); plotaxis2; legend('umbc','era5'); title('O3 feedback')
ylim([-0.25 0.50])

figure(112)
blah1 = (umbc_spectral_olr.planck_ecRad.clr + umbc_spectral_olr.lapse_ecRad.clr + umbc_spectral_olr.wv_ecRad.clr + umbc_spectral_olr.o3_ecRad.clr + umbc_spectral_olr.skt_ecRad.clr)/5 - umbc_spectral_olr.olr0_ecRad.clr;
blah1 = reshape(blah1',72,64); blah1 = nanmean(blah1,1);
blah2 = (era5_spectral_olr.planck_ecRad.clr + era5_spectral_olr.lapse_ecRad.clr + era5_spectral_olr.wv_ecRad.clr + era5_spectral_olr.o3_ecRad.clr + era5_spectral_olr.skt_ecRad.clr)/5 - era5_spectral_olr.olr0_ecRad.clr;
blah2 = reshape(blah2',72,64); blah2 = nanmean(blah2,1);
plot(rlat,smooth(blah1,nsmooth),rlat,smooth(blah2,nsmooth),'linewidth',2); plotaxis2; legend('umbc','era5'); title('\Sigma feedbacks')

disp('ret to continue'); pause
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

blah1 = umbc_spectral_olr.wv_ecRad.clr - umbc_spectral_olr.olr0_ecRad.clr; blah1 = reshape(blah1',72,64)/umbc_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
blah2 = era5_spectral_olr.wv_ecRad.clr - umbc_spectral_olr.olr0_ecRad.clr; blah2 = reshape(blah2',72,64)/era5_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
aslmap(104,rlat65,rlon73,-smoothn(blah1',1), [-90 +90],[-180 +180]); caxis([-1 +1]*4); colormap(usa2); title('UMBC WV feedback')
aslmap(105,rlat65,rlon73,-smoothn(blah2',1), [-90 +90],[-180 +180]); caxis([-1 +1]*4); colormap(usa2); title('ERA5 WV feedback');

blah1 = umbc_spectral_olr.planck_ecRad.clr - umbc_spectral_olr.olr0_ecRad.clr; blah1 = reshape(blah1',72,64)/umbc_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
blah2 = era5_spectral_olr.planck_ecRad.clr - umbc_spectral_olr.olr0_ecRad.clr; blah2 = reshape(blah2',72,64)/era5_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
aslmap(106,rlat65,rlon73,-smoothn(blah1',1), [-90 +90],[-180 +180]); caxis([-1 +1]*4); colormap(usa2); title('UMBC PLANCK feedback')
aslmap(107,rlat65,rlon73,-smoothn(blah2',1), [-90 +90],[-180 +180]); caxis([-1 +1]*4); colormap(usa2); title('ERA5 PLANCK feedback');

figure(108)
blah1 = reshape(umbc_spectral_olr.wv_ecRad.clr-umbc_spectral_olr.olr0_ecRad.clr,72,64); blah1 = nanmean(blah1,1); blah1 = blah1/umbc_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
blah2 = reshape(era5_spectral_olr.wv_ecRad.clr-umbc_spectral_olr.olr0_ecRad.clr,72,64); blah2 = nanmean(blah2,1); blah2 = blah2/era5_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
plot(rlat,smooth(-blah1,nsmooth),rlat,smooth(-blah2,nsmooth),'linewidth',2); plotaxis2; legend('umbc','era5'); title('WV feedback')

figure(109)
blah1 = reshape(umbc_spectral_olr.planck_ecRad.clr-umbc_spectral_olr.olr0_ecRad.clr,72,64); blah1 = nanmean(blah1,1); blah1 = blah1/umbc_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
blah2 = reshape(era5_spectral_olr.planck_ecRad.clr-umbc_spectral_olr.olr0_ecRad.clr,72,64); blah2 = nanmean(blah2,1); blah2 = blah2/era5_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
plot(rlat,smooth(-blah1,nsmooth),rlat,smooth(-blah2,nsmooth),'linewidth',2); plotaxis2; legend('umbc','era5'); title('Planck feedback')

figure(110)
blah1 = reshape(umbc_spectral_olr.lapse_ecRad.clr-umbc_spectral_olr.olr0_ecRad.clr,72,64); blah1 = nanmean(blah1,1)/umbc_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
blah2 = reshape(era5_spectral_olr.lapse_ecRad.clr-umbc_spectral_olr.olr0_ecRad.clr,72,64); blah2 = nanmean(blah2,1)/era5_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
plot(rlat,smooth(-blah1,nsmooth),rlat,smooth(-blah2,nsmooth),'linewidth',2); plotaxis2; legend('umbc','era5'); title('Lapse feedback')

figure(111)
blah1 = reshape(umbc_spectral_olr.o3_ecRad.clr-umbc_spectral_olr.olr0_ecRad.clr,72,64); blah1 = nanmean(blah1,1)/umbc_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
blah2 = reshape(era5_spectral_olr.o3_ecRad.clr-umbc_spectral_olr.olr0_ecRad.clr,72,64); blah2 = nanmean(blah2,1)/era5_spectral_olr.feedback_ecRad.global_coslat_wgt_skt;
plot(rlat,smooth(-blah1,nsmooth),rlat,smooth(-blah2,nsmooth),'linewidth',2); plotaxis2; legend('umbc','era5'); title('O3 feedback')

figure(112)
blah1 = (umbc_spectral_olr.planck_ecRad.clr + umbc_spectral_olr.lapse_ecRad.clr + umbc_spectral_olr.wv_ecRad.clr + umbc_spectral_olr.o3_ecRad.clr + umbc_spectral_olr.skt_ecRad.clr)/5 - umbc_spectral_olr.olr0_ecRad.clr;
blah1 = reshape(blah1',72,64); blah1 = nanmean(blah1,1);
blah2 = (era5_spectral_olr.planck_ecRad.clr + era5_spectral_olr.lapse_ecRad.clr + era5_spectral_olr.wv_ecRad.clr + era5_spectral_olr.o3_ecRad.clr + era5_spectral_olr.skt_ecRad.clr)/5 - era5_spectral_olr.olr0_ecRad.clr;
blah2 = reshape(blah2',72,64); blah2 = nanmean(blah2,1);
plot(rlat,smooth(blah1,nsmooth),rlat,smooth(blah2,nsmooth),'linewidth',2); plotaxis2; legend('umbc','era5'); title('\Sigma feedbacks')

