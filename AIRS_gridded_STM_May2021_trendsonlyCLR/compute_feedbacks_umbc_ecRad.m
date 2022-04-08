%% see “Simpson's Law” and the Spectral Cancellation of Climate Feedbacks
%% Nadir Jeevanjee1 , Daniel D. B. Koll2, and Nicholas Lutsko3
%% GRL eevanjee, N., Koll, D. D. B., & Lutsko, N. (2021). “Simpson's Law” and
%%  the spectral cancellation of climate feedbacks. Geophysical Research Letters, 48, e2021GL093699. https://doi. org/10.1029/2021GL093699

cdRRTMback = ['cd ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RobinHoganECMWF/ECRAD_ECMWF_version_of_flux/ecRad/create_ecrad_inputSergio/
addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/DRIVER_CODE_RRTM_Band17/

px = p;
umbc_spectral_olr.olr0 = compute_olr(h,px);
%umbc_spectral_olr.olr0_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0); %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<
umbc_spectral_olr.olr0_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
eval(cdRRTMback);

%% trying a faster method see     jaja = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);    driver_run_ecRad_rtp_loop_over_profiles.m iMethod = -1
%% scatter_coast(p.rlon,p.rlat,50,(umbc_spectral_olr.olr0_ecRad.clr-jaja.clr)./umbc_spectral_olr.olr0_ecRad.clr*100); colormap(usa2); caxis([-1 +1])
%% scatter_coast(p.rlon,p.rlat,50,(umbc_spectral_olr.olr0_rrtm-jaja.clr)./umbc_spectral_olr.olr0_rrtm*100); colormap(usa2); caxis([-1 +1])
eval(cdRRTMback);

px = p;
%px.gas_2 = px.gas_2*(1+2.2/400);
px.stemp = px.stemp + results(:,6)';
umbc_spectral_olr.skt = compute_olr(h,px);
%umbc_spectral_olr.skt_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0);
umbc_spectral_olr.skt_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
eval(cdRRTMback);

px = p;
%px.gas_2 = px.gas_2*(1+2.2/400);
px.stemp = px.stemp + results(:,6)';
px.ptemp = px.ptemp + ones(101,1)*results(:,6)';
umbc_spectral_olr.planck = compute_olr(h,px);
%umbc_spectral_olr.planck_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0);
umbc_spectral_olr.planck_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
eval(cdRRTMback);

px = p;
%px.gas_2 = px.gas_2*(1+2.2/400);
px.stemp = px.stemp + results(:,6)';
px.ptemp = px.ptemp + deltaT;
umbc_spectral_olr.lapse = compute_olr(h,px);   
%umbc_spectral_olr.lapse_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0);
umbc_spectral_olr.lapse_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
eval(cdRRTMback);

px = p;
%px.gas_2 = px.gas_2*(1+2.2/400);
%% px.gas_1 = px.gas_1 .* (1 + fracWV);
fracJUNK = fracWV; bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
px.gas_1 = px.gas_1 .* (1 + fracJUNK);
umbc_spectral_olr.wv = compute_olr(h,px);
%umbc_spectral_olr.wv_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0);
umbc_spectral_olr.wv_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               
eval(cdRRTMback);

%%%%%%%%%%%%%%%%%%%%%%%%%

%% change radiance mW --> W and then multiply by pi for flux
ix1 = 1:2162; ix2 = 2163:2645;  %% basically have two bands of detectors!

%junk = pi/1000*sum(umbc_spectral_olr.planck - umbc_spectral_olr.olr0,1);
%junk = -junk./results(:,6)';
%umbc_spectral_olr.feedback.planck = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),umbc_spectral_olr.planck(ix1,:) - umbc_spectral_olr.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),umbc_spectral_olr.planck(ix2,:) - umbc_spectral_olr.olr0(ix2,:));
junk = junk1 + junk2;
junk = -junk./results(:,6)';
umbc_spectral_olr.feedback.planck = junk;

%% note Eq 8b of Jevanjee paper shows we need to use 
%%   umbc_spectral_olr.feedback.lapse == umbc_spectral_olr.lapse - umbc_spectral_olr.planck
%% and not
%%   umbc_spectral_olr.feedback.lapse == umbc_spectral_olr.lapse - umbc_spectral_olr.olr0
junk1 = pi/1000*trapz(h.vchan(ix1),umbc_spectral_olr.lapse(ix1,:) - umbc_spectral_olr.planck(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),umbc_spectral_olr.lapse(ix2,:) - umbc_spectral_olr.planck(ix2,:));
junk = junk1 + junk2;
junk = -junk./results(:,6)';
umbc_spectral_olr.feedback.lapse = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),umbc_spectral_olr.wv(ix1,:) - umbc_spectral_olr.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),umbc_spectral_olr.wv(ix2,:) - umbc_spectral_olr.olr0(ix2,:));
junk = junk1 + junk2;
junk = -junk./results(:,6)';
umbc_spectral_olr.feedback.wv = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),umbc_spectral_olr.skt(ix1,:) - umbc_spectral_olr.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),umbc_spectral_olr.skt(ix2,:) - umbc_spectral_olr.olr0(ix2,:));
junk = junk1 + junk2;
junk = -junk./results(:,6)';
umbc_spectral_olr.feedback.skt = junk;

figure(71); scatter_coast(p.rlon,p.rlat,50,umbc_spectral_olr.feedback.planck); caxis([-4 0]*1);  colormap(jet);  title('UMBC \lambda_{Planck}')
figure(72); scatter_coast(p.rlon,p.rlat,50,umbc_spectral_olr.feedback.lapse);  caxis([-5 +5]*1); colormap(usa2); title('UMBC \lambda_{Lapse}')
figure(73); scatter_coast(p.rlon,p.rlat,50,umbc_spectral_olr.feedback.wv);     caxis([-2 +2]*1); colormap(usa2); title('UMBC \lambda_{WV}')
figure(74); scatter_coast(p.rlon,p.rlat,50,umbc_spectral_olr.feedback.skt);    caxis([-2 0]*1);  colormap(jet);  title('UMBC \lambda_{Skt}')

umbc_spectral_olr.feedback.planck_ecRad = umbc_spectral_olr.planck_ecRad.clr-umbc_spectral_olr.olr0_ecRad.clr; 
  umbc_spectral_olr.feedback.planck_ecRad = -umbc_spectral_olr.feedback.planck_ecRad./results(:,6)';
umbc_spectral_olr.feedback.lapse_ecRad = umbc_spectral_olr.lapse_ecRad.clr-umbc_spectral_olr.planck_ecRad.clr;
  umbc_spectral_olr.feedback.lapse_ecRad = -umbc_spectral_olr.feedback.lapse_ecRad./results(:,6)';
umbc_spectral_olr.feedback.wv_ecRad = umbc_spectral_olr.wv_ecRad.clr-umbc_spectral_olr.olr0_ecRad.clr;
  umbc_spectral_olr.feedback.wv_ecRad = -umbc_spectral_olr.feedback.wv_ecRad./results(:,6)';
umbc_spectral_olr.feedback.skt_ecRad = umbc_spectral_olr.skt_ecRad.clr-umbc_spectral_olr.olr0_ecRad.clr;
  umbc_spectral_olr.feedback.skt_ecRad = -umbc_spectral_olr.feedback.skt_ecRad./results(:,6)';

figure(75); scatter_coast(p.rlon,p.rlat,50,umbc_spectral_olr.feedback.planck_ecRad); caxis([-4 0]*1.5);  colormap(jet);  title('UMBC \lambda_{Planck}')
figure(76); scatter_coast(p.rlon,p.rlat,50,umbc_spectral_olr.feedback.lapse_ecRad);  caxis([-5 +5]*2); colormap(usa2); title('UMBC \lambda_{Lapse}')
figure(77); scatter_coast(p.rlon,p.rlat,50,umbc_spectral_olr.feedback.wv_ecRad);     caxis([-2 +2]*2); colormap(usa2); title('UMBC \lambda_{WV}')
figure(78); scatter_coast(p.rlon,p.rlat,50,umbc_spectral_olr.feedback.skt_ecRad);    caxis([-2 0]*1);  colormap(jet);  title('UMBC \lambda_{Skt}')

wonk = umbc_spectral_olr.feedback.planck_ecRad; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = 500; aslmap(75,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('UMBC \lambda_{Planck}')
wonk = umbc_spectral_olr.feedback.lapse_ecRad; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = 500; aslmap(76,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);  colormap(usa2); caxis([-5 5]*1);    title('UMBC \lambda_{Lapse}')
wonk = umbc_spectral_olr.feedback.wv_ecRad; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = 500; aslmap(77,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('UMBC \lambda_{WV}')
wonk = umbc_spectral_olr.feedback.skt_ecRad; wonk(wonk < -3) = NaN; wonk(wonk > +3) = NaN; 
  ns = 500; aslmap(78,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);    colormap(jet);  caxis([-2 0]*1.5);  title('UMBC \lambda_{Skt}')

figure(75); colormap(colormap_soden_held_jclim2007); caxis([-6.5 -1])
figure(76); colormap(colormap_soden_held_jclim2007); caxis([-2 +2])
figure(77); colormap(colormap_soden_held_jclim2007); caxis([-0.5 +5])
figure(78); colormap(colormap_soden_held_jclim2007); caxis([-2 0])

figure(75); colormap(colormap_soden_held_jclim2007); caxis([-4 -3])
figure(77); colormap(colormap_soden_held_jclim2007); caxis([-1 +1])

bad = find(abs(results(:,6)) < 1e-4);
junklat = reshape(p.rlat,72,64); junklat = mean(junklat,1);
figure(79); clf
  junk = umbc_spectral_olr.feedback.planck_ecRad; junk(bad) = NaN; junk = reshape(junk,72,64);  junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'b','linewidth',2); hold on
  junk = umbc_spectral_olr.feedback.lapse_ecRad; junk(bad) = NaN; junk = reshape(junk,72,64);   junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'g','linewidth',2); hold on
  junk = umbc_spectral_olr.feedback.wv_ecRad; junk(bad) = NaN; junk = reshape(junk,72,64);      junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'r','linewidth',2); hold on
  junk = umbc_spectral_olr.feedback.skt_ecRad; junk(bad) = NaN; junk = reshape(junk,72,64);     junk = nanmean(junk,1); junk = smooth(junk,5); plot(junklat,junk,'k','linewidth',2); hold off
ylim([-10 +10]/2); plotaxis2; title('UMBC \lambda'); xlabel('Latitude'); hl = legend('Planck','Lapse','WV','SKT','location','best','fontsize',8);

junk = umbc_spectral_olr.feedback.planck_ecRad; junk(bad) = NaN;
if ~exist('lps0')
  [mmw0,lps0] = mmwater_rtp(h,p);
end
figure(79); clf; plot(mean(lps0.lapse1(80:97,:),1),junk,'.');
figure(79); clf; scatter(mean(lps0.lapse1(80:97,:),1),junk,10,p.rlat,'filled');      colorbar; colormap jet; axis([0 10 -5 -3])
figure(79); clf; scatter(mean(lps0.lapse1(80:97,:),1),junk,10,abs(p.rlat),'filled'); colorbar; colormap jet; axis([0 10 -5 -3])
  xlabel('Lower Trop Lapse Rate K/km'); ylabel('\lambda_{Planck}');
addpath /home/sergio/MATLABCODE/SHOWSTATS
[n,nx,ny,nmean,nstd] = myhist2d(mean(lps0.lapse1(80:97,:),1),junk,0:0.25:10,-6:0.1:-3); errorbar(0:0.25:10,nmean,nstd); xlabel('Lower Trop Lapse Rate K/km'); ylabel('\lambda_{Planck}'); plotaxis2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

