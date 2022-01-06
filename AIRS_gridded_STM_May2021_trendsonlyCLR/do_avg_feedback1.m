clear wonk* mean_feedback
figure(1); clf

model = 1; 
  wonk11 = umbc_spectral_olr.feedback.planck_ecRad; wonk11(wonk11 < -10) = NaN; wonk11(wonk11 > +00) = NaN; wonk11 = smooth(wonk11 .* cos(p.rlat*pi/180),72); mean_feedback(model,1) = nanmean(wonk11); plot(p.rlat,wonk11);
  wonk12 = umbc_spectral_olr.feedback.lapse_ecRad;  wonk12(wonk12 < -05) = NaN; wonk12(wonk12 > +05) = NaN; wonk12 = smooth(wonk12 .* cos(p.rlat*pi/180),72); mean_feedback(model,2) = nanmean(wonk12); plot(p.rlat,wonk12);
  wonk13 = umbc_spectral_olr.feedback.wv_ecRad;     wonk13(wonk13 < -05) = NaN; wonk13(wonk13 > +05) = NaN; wonk13 = smooth(wonk13 .* cos(p.rlat*pi/180),72); mean_feedback(model,3) = nanmean(wonk13); plot(p.rlat,wonk13);
  wonk14 = umbc_spectral_olr.feedback.skt_ecRad;    wonk14(wonk14 < -05) = NaN; wonk14(wonk14 > +00) = NaN; wonk14 = smooth(wonk14 .* cos(p.rlat*pi/180),72); mean_feedback(model,4) = nanmean(wonk14); plot(p.rlat,wonk14);
model = 2; 
  wonk21 = airsL3_spectral_olr.feedback.planck_ecRad; wonk21(wonk21 < -10) = NaN; wonk21(wonk21 > +00) = NaN; wonk21 = smooth(wonk21 .* cos(p.rlat*pi/180),72); mean_feedback(model,1) = nanmean(wonk21); plot(p.rlat,wonk21);
  wonk22 = airsL3_spectral_olr.feedback.lapse_ecRad;  wonk22(wonk22 < -05) = NaN; wonk22(wonk22 > +05) = NaN; wonk22 = smooth(wonk22 .* cos(p.rlat*pi/180),72); mean_feedback(model,2) = nanmean(wonk22); plot(p.rlat,wonk22);
  wonk23 = airsL3_spectral_olr.feedback.wv_ecRad;     wonk23(wonk23 < -05) = NaN; wonk23(wonk23 > +05) = NaN; wonk23 = smooth(wonk23 .* cos(p.rlat*pi/180),72); mean_feedback(model,3) = nanmean(wonk23); plot(p.rlat,wonk23);
  wonk24 = airsL3_spectral_olr.feedback.skt_ecRad;    wonk24(wonk24 < -05) = NaN; wonk24(wonk24 > +00) = NaN; wonk24 = smooth(wonk24 .* cos(p.rlat*pi/180),72); mean_feedback(model,4) = nanmean(wonk24); plot(p.rlat,wonk24);
model = 3; 
  wonk31 = era5_spectral_olr.feedback.planck_ecRad; wonk31(wonk31 < -10) = NaN; wonk31(wonk31 > +00) = NaN; wonk31 = smooth(wonk31 .* cos(p.rlat*pi/180),72); mean_feedback(model,1) = nanmean(wonk31); plot(p.rlat,wonk31);
  wonk32 = era5_spectral_olr.feedback.lapse_ecRad;  wonk32(wonk32 < -05) = NaN; wonk32(wonk32 > +05) = NaN; wonk32 = smooth(wonk32 .* cos(p.rlat*pi/180),72); mean_feedback(model,2) = nanmean(wonk32); plot(p.rlat,wonk32);
  wonk33 = era5_spectral_olr.feedback.wv_ecRad;     wonk33(wonk33 < -05) = NaN; wonk33(wonk33 > +05) = NaN; wonk33 = smooth(wonk33 .* cos(p.rlat*pi/180),72); mean_feedback(model,3) = nanmean(wonk33); plot(p.rlat,wonk33);
  wonk34 = era5_spectral_olr.feedback.skt_ecRad;    wonk34(wonk34 < -05) = NaN; wonk34(wonk34 > +00) = NaN; wonk34 = smooth(wonk34 .* cos(p.rlat*pi/180),72); mean_feedback(model,4) = nanmean(wonk34); plot(p.rlat,wonk34);
model = 4; 
  wonk41 = cmip6_spectral_olr.feedback.planck_ecRad; wonk41(wonk41 < -10) = NaN; wonk41(wonk41 > +00) = NaN; wonk41 = smooth(wonk41 .* cos(p.rlat*pi/180),72); mean_feedback(model,1) = nanmean(wonk41); plot(p.rlat,wonk41);
  wonk42 = cmip6_spectral_olr.feedback.lapse_ecRad;  wonk42(wonk42 < -05) = NaN; wonk42(wonk42 > +05) = NaN; wonk42 = smooth(wonk42 .* cos(p.rlat*pi/180),72); mean_feedback(model,2) = nanmean(wonk42); plot(p.rlat,wonk42);
  wonk43 = cmip6_spectral_olr.feedback.wv_ecRad;     wonk43(wonk43 < -05) = NaN; wonk43(wonk43 > +05) = NaN; wonk43 = smooth(wonk43 .* cos(p.rlat*pi/180),72); mean_feedback(model,3) = nanmean(wonk43); plot(p.rlat,wonk43);
  wonk44 = cmip6_spectral_olr.feedback.skt_ecRad;    wonk44(wonk44 < -05) = NaN; wonk44(wonk44 > +00) = NaN; wonk44 = smooth(wonk44 .* cos(p.rlat*pi/180),72); mean_feedback(model,4) = nanmean(wonk44); plot(p.rlat,wonk44);

for ii = 1 : 4
  for jj = 1 : 4
    str = ['junk = wonk' num2str(ii) num2str(jj) ';'];
    eval(str);
    junk = smooth(junk,4)';
    %junk = smooth(junk,8)';
    str = ['swonk' num2str(ii) num2str(jj) ' = junk;'];
    eval(str);
  end
end

disp('<<<<<<<<<<<<< ---- do_avg_feedback1.m ----- >>>>>>>>>>>>>>')
disp('           Planck          Lapse           WV            SKT')
fprintf(1,'UMBC     %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(1,:))
fprintf(1,'AIRS L3  %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(2,:))
fprintf(1,'ERA5     %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(3,:))
fprintf(1,'CMIP6    %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(4,:))

figure(1); clf
subplot(221); plot(p.rlat,[wonk11 wonk21 wonk31 wonk41],'linewidth',2); ylabel('\lambda Planck'); xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90])
subplot(222); plot(p.rlat,[wonk12 wonk22 wonk32 wonk42],'linewidth',2); ylabel('\lambda Lapse');  xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90])
subplot(223); plot(p.rlat,[wonk13 wonk23 wonk33 wonk43],'linewidth',2); ylabel('\lambda WV');     xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90])
subplot(224); plot(p.rlat,[wonk14 wonk24 wonk34 wonk44],'linewidth',2); ylabel('\lambda SKT');    xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90])
