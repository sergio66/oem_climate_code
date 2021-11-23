%% see “Simpson's Law” and the Spectral Cancellation of Climate Feedbacks
%% Nadir Jeevanjee1 , Daniel D. B. Koll2, and Nicholas Lutsko3
%% GRL eevanjee, N., Koll, D. D. B., & Lutsko, N. (2021). “Simpson's Law” and
%%  the spectral cancellation of climate feedbacks. Geophysical Research Letters, 48, e2021GL093699. https://doi. org/10.1029/2021GL093699

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

px = p;
cmip6_spectral_olr.olr0 = compute_olr(h,px);

px = p;
px.stemp = px.stemp + cmip6.trend_stemp;
cmip6_spectral_olr.skt = compute_olr(h,px);

px = p;
px.stemp = px.stemp + cmip6.trend_stemp;
px.ptemp = px.ptemp + ones(101,1)*cmip6.trend_stemp;
cmip6_spectral_olr.planck = compute_olr(h,px);

px = p;
px.stemp = px.stemp + cmip6.trend_stemp;
px.ptemp(1:100,:) = px.ptemp(1:100,:) + cmip6.trend_ptemp(1:100,:);
cmip6_spectral_olr.lapse = compute_olr(h,px);   

px = p;
px.gas_1(1:100,:) = px.gas_1(1:100,:) .* (1 + cmip6.trend_gas_1(1:100,:));
cmip6_spectral_olr.wv = compute_olr(h,px);

%%%%%%%%%%%%%%%%%%%%%%%%%

%% change radiance mW --> W and then multiply by pi for flux
ix1 = 1:2162; ix2 = 2163:2645;  %% basically have two bands of detectors!

%junk = pi/1000*sum(cmip6_spectral_olr.planck - cmip6_spectral_olr.olr0,1);
%junk = -junk./cmip6.trend_stemp;
%cmip6_spectral_olr.feedback.planck = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),cmip6_spectral_olr.planck(ix1,:) - cmip6_spectral_olr.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),cmip6_spectral_olr.planck(ix2,:) - cmip6_spectral_olr.olr0(ix2,:));
junk = junk1 + junk2;
junk = -junk./cmip6.trend_stemp;
cmip6_spectral_olr.feedback.planck = junk;

%% note Eq 8b of Jevanjee paper shows we need to use 
%%   cmip6_spectral_olr.feedback.lapse == cmip6_spectral_olr.lapse - cmip6_spectral_olr.planck
%% and not
%%   cmip6_spectral_olr.feedback.lapse == cmip6_spectral_olr.lapse - cmip6_spectral_olr.olr0
junk1 = pi/1000*trapz(h.vchan(ix1),cmip6_spectral_olr.lapse(ix1,:) - cmip6_spectral_olr.planck(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),cmip6_spectral_olr.lapse(ix2,:) - cmip6_spectral_olr.planck(ix2,:));
junk = junk1 + junk2;
junk = -junk./cmip6.trend_stemp;
cmip6_spectral_olr.feedback.lapse = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),cmip6_spectral_olr.wv(ix1,:) - cmip6_spectral_olr.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),cmip6_spectral_olr.wv(ix2,:) - cmip6_spectral_olr.olr0(ix2,:));
junk = junk1 + junk2;
junk = -junk./cmip6.trend_stemp;
cmip6_spectral_olr.feedback.wv = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),cmip6_spectral_olr.skt(ix1,:) - cmip6_spectral_olr.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),cmip6_spectral_olr.skt(ix2,:) - cmip6_spectral_olr.olr0(ix2,:));
junk = junk1 + junk2;
junk = -junk./cmip6.trend_stemp;
cmip6_spectral_olr.feedback.skt = junk;

figure(71); scatter_coast(p.rlon,p.rlat,50,cmip6_spectral_olr.feedback.planck); caxis([-4 0]);  colormap(jet);  title('CMIP6 \lambda_{Planck}')
figure(72); scatter_coast(p.rlon,p.rlat,50,cmip6_spectral_olr.feedback.lapse);  caxis([-5 +5]); colormap(usa2); title('CMIP6 \lambda_{Lapse}')
figure(73); scatter_coast(p.rlon,p.rlat,50,cmip6_spectral_olr.feedback.wv);     caxis([-2 +2]); colormap(usa2); title('CMIP6 \lambda_{WV}')
figure(74); scatter_coast(p.rlon,p.rlat,50,cmip6_spectral_olr.feedback.skt);    caxis([-2 0]);  colormap(jet);  title('CMIP6 \lambda_{Skt}')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

