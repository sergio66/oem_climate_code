addpath /home/sergio/MATLABCODE/PLOTTER/

clear wonk* mean_feedback

cosrlat = cos(rlat'*pi/180);
nmcos = 1/nanmean(cosrlat);

figure(2); clf

model = 1; 
  wonk11 = umbc_spectral_olr.feedback_ecRad.planck.polyfit_latbin;   plot(rlat,wonk11); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk11); 
  wonk12 = umbc_spectral_olr.feedback_ecRad.lapse.polyfit_latbin;    plot(rlat,wonk12); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk12); 
  wonk13 = umbc_spectral_olr.feedback_ecRad.wv.polyfit_latbin;       plot(rlat,wonk13); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk13); 
  wonk14 = umbc_spectral_olr.feedback_ecRad.skt.polyfit_latbin;      plot(rlat,wonk14); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk14); 
model = 2; 
  wonk21 = airsL3_spectral_olr.feedback_ecRad.planck.polyfit_latbin; plot(rlat,wonk21); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk21); 
  wonk22 = airsL3_spectral_olr.feedback_ecRad.lapse.polyfit_latbin;  plot(rlat,wonk22); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk22); 
  wonk23 = airsL3_spectral_olr.feedback_ecRad.wv.polyfit_latbin;     plot(rlat,wonk23); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk23); 
  wonk24 = airsL3_spectral_olr.feedback_ecRad.skt.polyfit_latbin;    plot(rlat,wonk24); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk24); 
model = 3; 
  wonk31 = era5_spectral_olr.feedback_ecRad.planck.polyfit_latbin;   plot(rlat,wonk31); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk31); 
  wonk32 = era5_spectral_olr.feedback_ecRad.lapse.polyfit_latbin;    plot(rlat,wonk32); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk32); 
  wonk33 = era5_spectral_olr.feedback_ecRad.wv.polyfit_latbin;       plot(rlat,wonk33); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk33); 
  wonk34 = era5_spectral_olr.feedback_ecRad.skt.polyfit_latbin;      plot(rlat,wonk34); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk34); 
model = 4; 
  wonk41 = cmip6_spectral_olr.feedback_ecRad.planck.polyfit_latbin;  plot(rlat,wonk41); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk41); 
  wonk42 = cmip6_spectral_olr.feedback_ecRad.lapse.polyfit_latbin;   plot(rlat,wonk42); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk42); 
  wonk43 = cmip6_spectral_olr.feedback_ecRad.wv.polyfit_latbin;      plot(rlat,wonk43); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk43); 
  wonk44 = cmip6_spectral_olr.feedback_ecRad.skt.polyfit_latbin;     plot(rlat,wonk44); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk44); 

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

disp('<<<<<<<<<<<<< ---- do_avg_feedback3.m ----- >>>>>>>>>>>>>>')
disp('     uses eg feedback_ecRad.planck.polyfit * cos(rlat)')
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

