%% see “Simpson's Law” and the Spectral Cancellation of Climate Feedbacks
%% Nadir Jeevanjee1 , Daniel D. B. Koll2, and Nicholas Lutsko3
%% GRL Jeevanjee, N., Koll, D. D. B., & Lutsko, N. (2021). “Simpson's Law” and
%%  the spectral cancellation of climate feedbacks. Geophysical Research Letters, 48, e2021GL093699. https://doi. org/10.1029/2021GL093699

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

px = p;
airsL3_spectral_olr.olr0 = compute_olr(h,px);

px = p;
px.stemp = px.stemp + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
airsL3_spectral_olr.skt = compute_olr(h,px);

px = p;
px.stemp = px.stemp + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
px.ptemp(1:100,:) = px.ptemp(1:100,:) + ones(100,1)*nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
airsL3_spectral_olr.planck = compute_olr(h,px);

px = p;
px.stemp = px.stemp + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
px.ptemp(1:100,:) = px.ptemp(1:100,:) + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp(1:100,:);
airsL3_spectral_olr.lapse = compute_olr(h,px);   

px = p;
px.gas_1(1:100,:) = px.gas_1(1:100,:) .* (1 + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1(1:100,:));
airsL3_spectral_olr.wv = compute_olr(h,px);

%%%%%%%%%%%%%%%%%%%%%%%%%

%% change radiance mW --> W and then multiply by pi for flux
ix1 = 1:2162; ix2 = 2163:2645;  %% basically have two bands of detectors!

%junk = pi/1000*sum(airsL3_spectral_olr.planck - airsL3_spectral_olr.olr0,1);
%junk = -junk./nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
%airsL3_spectral_olr.feedback.planck = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),airsL3_spectral_olr.planck(ix1,:) - airsL3_spectral_olr.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),airsL3_spectral_olr.planck(ix2,:) - airsL3_spectral_olr.olr0(ix2,:));
junk = junk1 + junk2;
junk = -junk./nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
airsL3_spectral_olr.feedback.planck = junk;

%% note Eq 8b of Jevanjee paper shows we need to use 
%%   airsL3_spectral_olr.feedback.lapse == airsL3_spectral_olr.lapse - airsL3_spectral_olr.planck
%% and not
%%   airsL3_spectral_olr.feedback.lapse == airsL3_spectral_olr.lapse - airsL3_spectral_olr.olr0
junk1 = pi/1000*trapz(h.vchan(ix1),airsL3_spectral_olr.lapse(ix1,:) - airsL3_spectral_olr.planck(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),airsL3_spectral_olr.lapse(ix2,:) - airsL3_spectral_olr.planck(ix2,:));
junk = junk1 + junk2;
junk = -junk./nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
airsL3_spectral_olr.feedback.lapse = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),airsL3_spectral_olr.wv(ix1,:) - airsL3_spectral_olr.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),airsL3_spectral_olr.wv(ix2,:) - airsL3_spectral_olr.olr0(ix2,:));
junk = junk1 + junk2;
junk = -junk./nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
airsL3_spectral_olr.feedback.wv = junk;

junk1 = pi/1000*trapz(h.vchan(ix1),airsL3_spectral_olr.skt(ix1,:) - airsL3_spectral_olr.olr0(ix1,:));
junk2 = pi/1000*trapz(h.vchan(ix2),airsL3_spectral_olr.skt(ix2,:) - airsL3_spectral_olr.olr0(ix2,:));
junk = junk1 + junk2;
junk = -junk./nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
airsL3_spectral_olr.feedback.skt = junk;

figure(71); scatter_coast(p.rlon,p.rlat,50,airsL3_spectral_olr.feedback.planck); caxis([-4 0]);  colormap(jet);  title('AIRSL3 \lambda_{Planck}')
figure(72); scatter_coast(p.rlon,p.rlat,50,airsL3_spectral_olr.feedback.lapse);  caxis([-5 +5]); colormap(usa2); title('AIRSL3 \lambda_{Lapse}')
figure(73); scatter_coast(p.rlon,p.rlat,50,airsL3_spectral_olr.feedback.wv);     caxis([-2 +2]); colormap(usa2); title('AIRSL3 \lambda_{WV}')
figure(74); scatter_coast(p.rlon,p.rlat,50,airsL3_spectral_olr.feedback.skt);    caxis([-2 0]);  colormap(jet);  title('AIRSL3 \lambda_{Skt}')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

