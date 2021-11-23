clear wonk* mean_feedback
figure(2); clf

model = 1; 
  wonk11 = umbc_spectral_olr.feedback.planck_ecRad; wonk11(wonk11 < -10) = NaN; wonk11(wonk11 > +00) = NaN; wonk11 = reshape(wonk11,72,64); wonk11 = nanmean(wonk11,1); plot(rlat,wonk11); mean_feedback(model,1) = nanmean(wonk11); 
  wonk12 = umbc_spectral_olr.feedback.lapse_ecRad;  wonk12(wonk12 < -05) = NaN; wonk12(wonk12 > +05) = NaN; wonk12 = reshape(wonk12,72,64); wonk12 = nanmean(wonk12,1); plot(rlat,wonk12); mean_feedback(model,2) = nanmean(wonk12); 
  wonk13 = umbc_spectral_olr.feedback.wv_ecRad;     wonk13(wonk13 < -05) = NaN; wonk13(wonk13 > +05) = NaN; wonk13 = reshape(wonk13,72,64); wonk13 = nanmean(wonk13,1); plot(rlat,wonk13); mean_feedback(model,3) = nanmean(wonk13); 
  wonk14 = umbc_spectral_olr.feedback.skt_ecRad;    wonk14(wonk14 < -05) = NaN; wonk14(wonk14 > +00) = NaN; wonk14 = reshape(wonk14,72,64); wonk14 = nanmean(wonk14,1); plot(rlat,wonk14); mean_feedback(model,4) = nanmean(wonk14); 
model = 2; 
  wonk21 = airsL3_spectral_olr.feedback.planck_ecRad; wonk21(wonk21 < -10) = NaN; wonk21(wonk21 > +00) = NaN; wonk21 = reshape(wonk21,72,64); wonk21 = nanmean(wonk21,1); plot(rlat,wonk21); mean_feedback(model,1) = nanmean(wonk21); 
  wonk22 = airsL3_spectral_olr.feedback.lapse_ecRad;  wonk22(wonk22 < -05) = NaN; wonk22(wonk22 > +05) = NaN; wonk22 = reshape(wonk22,72,64); wonk22 = nanmean(wonk22,1); plot(rlat,wonk22); mean_feedback(model,2) = nanmean(wonk22); 
  wonk23 = airsL3_spectral_olr.feedback.wv_ecRad;     wonk23(wonk23 < -05) = NaN; wonk23(wonk23 > +05) = NaN; wonk23 = reshape(wonk23,72,64); wonk23 = nanmean(wonk23,1); plot(rlat,wonk23); mean_feedback(model,3) = nanmean(wonk23); 
  wonk24 = airsL3_spectral_olr.feedback.skt_ecRad;    wonk24(wonk24 < -05) = NaN; wonk24(wonk24 > +00) = NaN; wonk24 = reshape(wonk24,72,64); wonk24 = nanmean(wonk24,1); plot(rlat,wonk24); mean_feedback(model,4) = nanmean(wonk24); 
model = 3; 
  wonk31 = era5_spectral_olr.feedback.planck_ecRad; wonk31(wonk31 < -10) = NaN; wonk31(wonk31 > +00) = NaN; wonk31 = reshape(wonk31,72,64); wonk31 = nanmean(wonk31,1); plot(rlat,wonk31); mean_feedback(model,1) = nanmean(wonk31); 
  wonk32 = era5_spectral_olr.feedback.lapse_ecRad;  wonk32(wonk32 < -05) = NaN; wonk32(wonk32 > +05) = NaN; wonk32 = reshape(wonk32,72,64); wonk32 = nanmean(wonk32,1); plot(rlat,wonk32); mean_feedback(model,2) = nanmean(wonk32); 
  wonk33 = era5_spectral_olr.feedback.wv_ecRad;     wonk33(wonk33 < -05) = NaN; wonk33(wonk33 > +05) = NaN; wonk33 = reshape(wonk33,72,64); wonk33 = nanmean(wonk33,1); plot(rlat,wonk33); mean_feedback(model,3) = nanmean(wonk33); 
  wonk34 = era5_spectral_olr.feedback.skt_ecRad;    wonk34(wonk34 < -05) = NaN; wonk34(wonk34 > +00) = NaN; wonk34 = reshape(wonk34,72,64); wonk34 = nanmean(wonk34,1); plot(rlat,wonk34); mean_feedback(model,4) = nanmean(wonk34); 
model = 4; 
  wonk41 = cmip6_spectral_olr.feedback.planck_ecRad; wonk41(wonk41 < -10) = NaN; wonk41(wonk41 > +00) = NaN; wonk41 = reshape(wonk41,72,64); wonk41 = nanmean(wonk41,1); plot(rlat,wonk41); mean_feedback(model,1) = nanmean(wonk41); 
  wonk42 = cmip6_spectral_olr.feedback.lapse_ecRad;  wonk42(wonk42 < -05) = NaN; wonk42(wonk42 > +05) = NaN; wonk42 = reshape(wonk42,72,64); wonk42 = nanmean(wonk42,1); plot(rlat,wonk42); mean_feedback(model,2) = nanmean(wonk42); 
  wonk43 = cmip6_spectral_olr.feedback.wv_ecRad;     wonk43(wonk43 < -05) = NaN; wonk43(wonk43 > +05) = NaN; wonk43 = reshape(wonk43,72,64); wonk43 = nanmean(wonk43,1); plot(rlat,wonk43); mean_feedback(model,3) = nanmean(wonk43); 
  wonk44 = cmip6_spectral_olr.feedback.skt_ecRad;    wonk44(wonk44 < -05) = NaN; wonk44(wonk44 > +00) = NaN; wonk44 = reshape(wonk44,72,64); wonk44 = nanmean(wonk44,1); plot(rlat,wonk44); mean_feedback(model,4) = nanmean(wonk44); 

disp('           Planck          Lapse           WV            SKT')
fprintf(1,'UMBC     %9.6f      %9.6f     %9.6f    %9.6f \n',mean_feedback(1,:))
fprintf(1,'AIRS L3  %9.6f      %9.6f     %9.6f    %9.6f \n',mean_feedback(2,:))
fprintf(1,'ERA5     %9.6f      %9.6f     %9.6f    %9.6f \n',mean_feedback(3,:))
fprintf(1,'CMIP6    %9.6f      %9.6f     %9.6f    %9.6f \n',mean_feedback(4,:))

figure(2); clf
subplot(221); plot(rlat,[wonk11; wonk21; wonk31; wonk41],'linewidth',2); ylabel('\lambda Planck'); xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90]);
subplot(222); plot(rlat,[wonk12; wonk22; wonk32; wonk42],'linewidth',2); ylabel('\lambda Lapse');  xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90])
subplot(223); plot(rlat,[wonk13; wonk23; wonk33; wonk43],'linewidth',2); ylabel('\lambda WV');     xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90])
subplot(224); plot(rlat,[wonk14; wonk24; wonk34; wonk44],'linewidth',2); ylabel('\lambda SKT');    xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90])
