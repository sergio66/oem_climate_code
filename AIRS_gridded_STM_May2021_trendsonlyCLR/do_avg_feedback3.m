addpath /home/sergio/MATLABCODE/PLOTTER/

clear wonk* mean_feedback

cosrlat = cos(rlat'*pi/180);
nmcos = 1/nanmean(cosrlat);

figure(2); clf

model = 1; 
  wonk11 = umbc_spectral_olr.feedback.planck_ecRad_polyfit_latbin;   plot(rlat,wonk11); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk11); 
  wonk12 = umbc_spectral_olr.feedback.lapse_ecRad_polyfit_latbin;    plot(rlat,wonk12); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk12); 
  wonk13 = umbc_spectral_olr.feedback.wv_ecRad_polyfit_latbin;       plot(rlat,wonk13); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk13); 
  wonk14 = umbc_spectral_olr.feedback.skt_ecRad_polyfit_latbin;      plot(rlat,wonk14); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk14); 
model = 2; 
  wonk21 = airsL3_spectral_olr.feedback.planck_ecRad_polyfit_latbin; plot(rlat,wonk21); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk21); 
  wonk22 = airsL3_spectral_olr.feedback.lapse_ecRad_polyfit_latbin;  plot(rlat,wonk22); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk22); 
  wonk23 = airsL3_spectral_olr.feedback.wv_ecRad_polyfit_latbin;     plot(rlat,wonk23); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk23); 
  wonk24 = airsL3_spectral_olr.feedback.skt_ecRad_polyfit_latbin;    plot(rlat,wonk24); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk24); 
model = 3; 
  wonk31 = era5_spectral_olr.feedback.planck_ecRad_polyfit_latbin;   plot(rlat,wonk31); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk31); 
  wonk32 = era5_spectral_olr.feedback.lapse_ecRad_polyfit_latbin;    plot(rlat,wonk32); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk32); 
  wonk33 = era5_spectral_olr.feedback.wv_ecRad_polyfit_latbin;       plot(rlat,wonk33); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk33); 
  wonk34 = era5_spectral_olr.feedback.skt_ecRad_polyfit_latbin;      plot(rlat,wonk34); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk34); 
model = 4; 
  wonk41 = cmip6_spectral_olr.feedback.planck_ecRad_polyfit_latbin;  plot(rlat,wonk41); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk41); 
  wonk42 = cmip6_spectral_olr.feedback.lapse_ecRad_polyfit_latbin;   plot(rlat,wonk42); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk42); 
  wonk43 = cmip6_spectral_olr.feedback.wv_ecRad_polyfit_latbin;      plot(rlat,wonk43); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk43); 
  wonk44 = cmip6_spectral_olr.feedback.skt_ecRad_polyfit_latbin;     plot(rlat,wonk44); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk44); 

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

disp('<<<<<<<<<<<<< ---- do_avg_feedback2cos.m ----- >>>>>>>>>>>>>>')
disp('           Planck          Lapse           WV            SKT')
fprintf(1,'UMBC     %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(1,:))
fprintf(1,'AIRS L3  %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(2,:))
fprintf(1,'ERA5     %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(3,:))
fprintf(1,'CMIP6    %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(4,:))

figure(2); clf
subplot(221); plot(rlat,[wonk11; wonk21; wonk31; wonk41],'linewidth',2); ylabel('\lambda Planck'); xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90]);
subplot(222); plot(rlat,[wonk12; wonk22; wonk32; wonk42],'linewidth',2); ylabel('\lambda Lapse');  xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90])
subplot(223); plot(rlat,[wonk13; wonk23; wonk33; wonk43],'linewidth',2); ylabel('\lambda WV');     xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90])
subplot(224); plot(rlat,[wonk14; wonk24; wonk34; wonk44],'linewidth',2); ylabel('\lambda SKT');    xlabel('latitude'); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); xlim([-90 +90])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf

ta = tiledlayout(2,2);
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile; plot(rlat,[wonk11; wonk21; wonk31; wonk41],'linewidth',2); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); 
tafov(2) = nexttile; plot(rlat,[wonk12; wonk22; wonk32; wonk42],'linewidth',2); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); 
tafov(3) = nexttile; plot(rlat,[wonk13; wonk23; wonk33; wonk43],'linewidth',2); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); 
tafov(4) = nexttile; plot(rlat,[wonk14; wonk24; wonk34; wonk44],'linewidth',2); plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); 

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;
tafov(3).FontSize = 10;
tafov(4).FontSize = 10;

set(tafov,'Xlim',[-90 +90]);

% Get rid of all extra space I can
ta.Padding = 'none';
ta.TileSpacing = 'none';
ta.TileSpacing = 'compact';

% Remove all ytick labels except for 1st column
%for ii = [2 4]
%   tafov(ii).YTickLabel = '';
%end
% Remove all xtick labels except for 2nd row
for ii = [1 2]
   tafov(ii).XTickLabel = '';
end

% Put in xlabel and ylable in the “middle”
tafov(1).YLabel.String = 'W/m2/K';   tafov(1).YLabel.FontSize = 10;
tafov(3).YLabel.String = '\lambda';  tafov(3).YLabel.FontSize = 10;
tafov(3).XLabel.String = 'Latitude'; tafov(3).XLabel.FontSize = 10;
tafov(4).XLabel.String = 'Latitude'; tafov(4).XLabel.FontSize = 10;

%% put titles
title(tafov(1),'\lambda Planck', 'Units', 'normalized', 'Position', [0.5, +1.025, 0]);
title(tafov(2),'\lambda Lapse', 'Units', 'normalized', 'Position', [0.5, +1.025, 0]);
title(tafov(3),'\lambda WV', 'Units', 'normalized', 'Position', [0.5, +1.025, 0]);
title(tafov(4),'\lambda SKT', 'Units', 'normalized', 'Position', [0.5, +1.025, 0]);

%% figure(2); figname = [printdir '/tiled_feedbackparams.pdf'];     aslprint(figname)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('do_avg_feedback3.m : remember feedback = lapse + planck + wv, do NOT INCLUDE skt as it is already in lapse')
disp('do_avg_feedback3.m : remember feedback = lapse + planck + wv, do NOT INCLUDE skt as it is already in lapse')
disp('do_avg_feedback3.m : remember feedback = lapse + planck + wv, do NOT INCLUDE skt as it is already in lapse')
disp(' ')

figure(4);
plot(rlat,[swonk11+swonk12+swonk13+swonk14*0; swonk21+swonk22+swonk23+swonk24*0; swonk31+swonk32+swonk33+swonk34*0; swonk41+swonk42+swonk43+swonk44*0],'linewidth',2); 
  plotaxis2; hl = legend('UMBC','AIRSL3','ERA5','CMIP6','location','best','fontsize',6); title('Sum \lambda')
junk = [swonk11+swonk12+swonk13+swonk14*0; swonk21+swonk22+swonk23+swonk24*0; swonk31+swonk32+swonk33+swonk34*0; swonk41+swonk42+swonk43+swonk44*0;];
globalavg_lambda = sum((ones(4,1)*cosrlat).*junk) ./ sum((ones(4,1)*cosrlat));
globalavg_lambda = sum((ones(4,1)*cosrlat).*junk,2) ./ sum((ones(4,1)*cosrlat),2);
fprintf(1,' do_avg_feedback3.m : global avg feedback = %8.6f %8.6f %8.6f %8.6f W/m2/K for UMBC/AIRSL3/ERA5/CMIP6 \n',globalavg_lambda)

figure(5);
model = 1; 
  wonk11 = umbc_spectral_olr.feedback.planck_ecRad;   
  wonk12 = umbc_spectral_olr.feedback.lapse_ecRad;    
  wonk13 = umbc_spectral_olr.feedback.wv_ecRad;       
  wonk14 = umbc_spectral_olr.feedback.skt_ecRad;      
  wonk1  = wonk11 + wonk12 + wonk13;
model = 2; 
  wonk21 = airsL3_spectral_olr.feedback.planck_ecRad; 
  wonk22 = airsL3_spectral_olr.feedback.lapse_ecRad;  
  wonk23 = airsL3_spectral_olr.feedback.wv_ecRad;     
  wonk24 = airsL3_spectral_olr.feedback.skt_ecRad;    
  wonk2  = wonk21 + wonk22 + wonk23;
model = 3; 
  wonk31 = era5_spectral_olr.feedback.planck_ecRad;   
  wonk32 = era5_spectral_olr.feedback.lapse_ecRad;    
  wonk33 = era5_spectral_olr.feedback.wv_ecRad;       
  wonk34 = era5_spectral_olr.feedback.skt_ecRad;      
  wonk3  = wonk31 + wonk32 + wonk33;
model = 4; 
  wonk41 = cmip6_spectral_olr.feedback.planck_ecRad;  
  wonk42 = cmip6_spectral_olr.feedback.lapse_ecRad;   
  wonk43 = cmip6_spectral_olr.feedback.wv_ecRad;      
  wonk44 = cmip6_spectral_olr.feedback.skt_ecRad;     
  wonk4  = wonk41 + wonk42 + wonk43;

%%%%%%%%%%%%%%%%%%%%%%%%%

model = 1; 
  wonk11 = umbc_spectral_olr.feedback.planck_ecRad;   junk(1)=mean(nanmean(wonk11(:))+nanmedian(wonk11(:))); junk(2)=nanstd(wonk11(:)); bad=find(wonk11 < (junk(1)-junk(2)) | wonk11 > (junk(1)+junk(2))); wonk11(bad)=NaN;
  wonk12 = umbc_spectral_olr.feedback.lapse_ecRad;    junk(1)=mean(nanmean(wonk12(:))+nanmedian(wonk12(:))); junk(2)=nanstd(wonk12(:)); bad=find(wonk12 < (junk(1)-junk(2)) | wonk12 > (junk(1)+junk(2))); wonk12(bad)=NaN;
  wonk13 = umbc_spectral_olr.feedback.wv_ecRad;       junk(1)=mean(nanmean(wonk13(:))+nanmedian(wonk13(:))); junk(2)=nanstd(wonk13(:)); bad=find(wonk13 < (junk(1)-junk(2)) | wonk13 > (junk(1)+junk(2))); wonk13(bad)=NaN;
  wonk14 = umbc_spectral_olr.feedback.skt_ecRad;      junk(1)=mean(nanmean(wonk14(:))+nanmedian(wonk14(:))); junk(2)=nanstd(wonk14(:)); bad=find(wonk14 < (junk(1)-junk(2)) | wonk14 > (junk(1)+junk(2))); wonk14(bad)=NaN;
  wonk1  = wonk11 + wonk12 + wonk13;
model = 2; 
  wonk21 = airsL3_spectral_olr.feedback.planck_ecRad; junk(1)=mean(nanmean(wonk21(:))+nanmedian(wonk21(:))); junk(2)=nanstd(wonk21(:)); bad=find(wonk21 < (junk(1)-junk(2)) | wonk21 > (junk(1)+junk(2))); wonk21(bad)=NaN;
  wonk22 = airsL3_spectral_olr.feedback.lapse_ecRad;  junk(1)=mean(nanmean(wonk22(:))+nanmedian(wonk22(:))); junk(2)=nanstd(wonk22(:)); bad=find(wonk22 < (junk(1)-junk(2)) | wonk22 > (junk(1)+junk(2))); wonk22(bad)=NaN;
  wonk23 = airsL3_spectral_olr.feedback.wv_ecRad;     junk(1)=mean(nanmean(wonk23(:))+nanmedian(wonk23(:))); junk(2)=nanstd(wonk23(:)); bad=find(wonk23 < (junk(1)-junk(2)) | wonk23 > (junk(1)+junk(2))); wonk23(bad)=NaN;
  wonk24 = airsL3_spectral_olr.feedback.skt_ecRad;    junk(1)=mean(nanmean(wonk24(:))+nanmedian(wonk24(:))); junk(2)=nanstd(wonk24(:)); bad=find(wonk24 < (junk(1)-junk(2)) | wonk24 > (junk(1)+junk(2))); wonk24(bad)=NaN;
  wonk2  = wonk21 + wonk22 + wonk23;
model = 3; 
  wonk31 = era5_spectral_olr.feedback.planck_ecRad;   junk(1)=mean(nanmean(wonk31(:))+nanmedian(wonk31(:))); junk(2)=nanstd(wonk31(:)); bad=find(wonk31 < (junk(1)-junk(2)) | wonk31 > (junk(1)+junk(2))); wonk31(bad)=NaN;
  wonk32 = era5_spectral_olr.feedback.lapse_ecRad;    junk(1)=mean(nanmean(wonk32(:))+nanmedian(wonk32(:))); junk(2)=nanstd(wonk32(:)); bad=find(wonk32 < (junk(1)-junk(2)) | wonk32 > (junk(1)+junk(2))); wonk32(bad)=NaN;
  wonk33 = era5_spectral_olr.feedback.wv_ecRad;       junk(1)=mean(nanmean(wonk33(:))+nanmedian(wonk33(:))); junk(2)=nanstd(wonk33(:)); bad=find(wonk33 < (junk(1)-junk(2)) | wonk33 > (junk(1)+junk(2))); wonk33(bad)=NaN;
  wonk34 = era5_spectral_olr.feedback.skt_ecRad;      junk(1)=mean(nanmean(wonk34(:))+nanmedian(wonk34(:))); junk(2)=nanstd(wonk34(:)); bad=find(wonk34 < (junk(1)-junk(2)) | wonk34 > (junk(1)+junk(2))); wonk34(bad)=NaN;
  wonk3  = wonk31 + wonk32 + wonk33;
model = 4; 
  wonk41 = cmip6_spectral_olr.feedback.planck_ecRad;  junk(1)=mean(nanmean(wonk41(:))+nanmedian(wonk41(:))); junk(2)=nanstd(wonk41(:)); bad=find(wonk41 < (junk(1)-junk(2)) | wonk41 > (junk(1)+junk(2))); wonk41(bad)=NaN;
  wonk42 = cmip6_spectral_olr.feedback.lapse_ecRad;   junk(1)=mean(nanmean(wonk42(:))+nanmedian(wonk42(:))); junk(2)=nanstd(wonk42(:)); bad=find(wonk42 < (junk(1)-junk(2)) | wonk42 > (junk(1)+junk(2))); wonk42(bad)=NaN;
  wonk43 = cmip6_spectral_olr.feedback.wv_ecRad;      junk(1)=mean(nanmean(wonk43(:))+nanmedian(wonk43(:))); junk(2)=nanstd(wonk43(:)); bad=find(wonk43 < (junk(1)-junk(2)) | wonk43 > (junk(1)+junk(2))); wonk43(bad)=NaN;
  wonk44 = cmip6_spectral_olr.feedback.skt_ecRad;     junk(1)=mean(nanmean(wonk44(:))+nanmedian(wonk44(:))); junk(2)=nanstd(wonk44(:)); bad=find(wonk44 < (junk(1)-junk(2)) | wonk44 > (junk(1)+junk(2))); wonk44(bad)=NaN;
  wonk4  = wonk41 + wonk42 + wonk43;


iS = 72; iS = 72/4; iS = 72*2; iS = 72/4;
plotoptions.str11 = 'UMBC'; plotoptions.str12 = 'AIRS L3'; plotoptions.str21 = 'ERA5'; plotoptions.str22 = 'CMIP6';
plotoptions.cx = [-6  -2]; plotoptions.maintitle = 'Planck'; aslmap_2x2tiledlayout(smooth(wonk11,iS),smooth(wonk21,iS),smooth(wonk31,iS),smooth(wonk41,iS),5,plotoptions); 
plotoptions.cx = [-20 +10]; plotoptions.maintitle = 'Total';  aslmap_2x2tiledlayout(smooth(wonk1,iS),smooth(wonk2,iS),smooth(wonk3,iS),smooth(wonk4,iS),5,plotoptions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6); clf
  kk = 1;
  for kk = 1 : 6
    errorbar(umbc_spectral_olr.ecRadresults.dbins,umbc_spectral_olr.ecRadresults.nmeanval(kk,:),umbc_spectral_olr.ecRadresults.nstdval(kk,:)); hold on
  end
  axis([-0.25 0.35 -1 +1]); plotaxis2; ylabel('OEM \newline ');
  hold off

  for kk = 1 : 6
    h2 = errorbar(era5_spectral_olr.ecRadresults.dbins,era5_spectral_olr.ecRadresults.nmeanval(kk,:),era5_spectral_olr.ecRadresults.nstdval(kk,:)); hold on
  end
  axis([-0.25 0.35 -1 +1]); plotaxis2; xlabel('\delta SKT (K)'); ylabel('ERA5 \newline '); 
  text(-0.3,0.25,'Forcing W/2','Rotation',90,'Fontsize',18)

  hl= legend('Planck','Lapse','O3','WV','SKT','T/CO2','location','northeast','fontsize',8); 
  hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
