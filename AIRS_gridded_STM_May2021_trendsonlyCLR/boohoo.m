compute_feedbacks_airsL3_ecRad ; pause(0.1)
compute_feedbacks_era5_ecRad   ; pause(0.1)
compute_feedbacks_cmip6_ecRad  ; pause(0.1)

ns0 = 500;
ns0 = 50;
ns0 = 1;
wonk = umbc_spectral_olr.feedback.planck_ecRad; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(75,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('UMBC \lambda_{Planck}')
wonk = umbc_spectral_olr.feedback.wv_ecRad; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = ns0; aslmap(76,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('UMBC \lambda_{WV}')

wonk = airsL3_spectral_olr.feedback.planck_ecRad; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(77,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('AIRSL3 \lambda_{Planck}')
wonk = airsL3_spectral_olr.feedback.wv_ecRad; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = ns0; aslmap(78,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('AIRSL3 \lambda_{WV}')

wonk = era5_spectral_olr.feedback.planck_ecRad; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(79,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('ERA5 \lambda_{Planck}')
wonk = era5_spectral_olr.feedback.wv_ecRad; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = ns0; aslmap(80,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('ERA5 \lambda_{WV}')

wonk = cmip6_spectral_olr.feedback.planck_ecRad; wonk(wonk < -10) = NaN; wonk(wonk > 0) = NaN; 
  ns = ns0; aslmap(81,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]); colormap(jet);  caxis([-4 0]*1.5);  title('CMIP6 \lambda_{Planck}')
wonk = cmip6_spectral_olr.feedback.wv_ecRad; wonk(wonk < -5) = NaN; wonk(wonk > +5) = NaN; 
  ns = ns0; aslmap(82,rlat65,rlon73,maskLFmatr.*smoothn((reshape(wonk,72,64)'),ns),[-90 +90],[-180 +180]);     colormap(usa2); caxis([-2 2]*1);    title('CMIP6 \lambda_{WV}')

figure(75); colormap(colormap_soden_held_jclim2007); caxis([-6 -3]); figure(76); colormap(colormap_soden_held_jclim2007); caxis([-2 +4])
figure(77); colormap(colormap_soden_held_jclim2007); caxis([-6 -3]); figure(78); colormap(colormap_soden_held_jclim2007); caxis([-2 +4])
figure(79); colormap(colormap_soden_held_jclim2007); caxis([-6 -3]); figure(80); colormap(colormap_soden_held_jclim2007); caxis([-2 +4])
figure(81); colormap(colormap_soden_held_jclim2007); caxis([-6 -3]); figure(82); colormap(colormap_soden_held_jclim2007); caxis([-2 +4])

figure(75); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(76); colormap(colormap_soden_held_jclim2007); caxis([-2 +4]*2)
figure(77); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(78); colormap(colormap_soden_held_jclim2007); caxis([-2 +4]*2)
figure(79); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(80); colormap(colormap_soden_held_jclim2007); caxis([-2 +4]*2)
figure(81); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(82); colormap(colormap_soden_held_jclim2007); caxis([-2 +4]*2)

figure(75); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(76); colormap(colormap_soden_held_jclim2007); caxis([-2 +2]*2)
figure(77); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(78); colormap(colormap_soden_held_jclim2007); caxis([-2 +2]*2)
figure(79); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(80); colormap(colormap_soden_held_jclim2007); caxis([-2 +2]*2)
figure(81); colormap(colormap_soden_held_jclim2007); caxis([-4.5 -3]); figure(82); colormap(colormap_soden_held_jclim2007); caxis([-2 +2]*2)

