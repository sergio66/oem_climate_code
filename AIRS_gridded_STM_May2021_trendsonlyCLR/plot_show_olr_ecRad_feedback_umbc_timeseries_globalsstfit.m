clear showfeedbacks* strfeedbacks

ix = 1; iNumYears = 05; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '05 years ';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted;
umbc05_spectral_olr = junk.umbc_spectral_olr;

ix = 2; iNumYears = 10; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '10 years ';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted;
umbc10_spectral_olr = junk.umbc_spectral_olr;

ix = 3; iNumYears = 15; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '15 years ';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted;
umbc15_spectral_olr = junk.umbc_spectral_olr;

ix = 4; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '20 years ';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted;
umbc20_spectral_olr = junk.umbc_spectral_olr;

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
ixx = ix;
showfeedbacks(1:ixx,7) = sum(showfeedbacks(1:ixx,[1 2 3 4]),2);

for ix = 1 : 4
  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},showfeedbacks(ix,[1 2 3 4 7]));
end
trends_paper_show = showfeedbacks(1:ixx,[1 2 3 4 7]);

figure(1); clf
bar(trends_paper_show');
ylabel('Feedback W/m2/K');
hl = legend('05','10','15','20','location','north','fontsize',10);
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rlat')
  load latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
  [Y,X] = meshgrid(rlat,rlon);
  X = X; Y = Y;
end

figure(2); clf;
subplot(221); 
plot(rlat,umbc05_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,rlat,umbc10_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,...
     rlat,umbc15_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,rlat,umbc20_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,...
     'linewidth',2);
  plotaxis2; %hl = legend('05','10','15','20','location','south','fontsize',8);
  title('Planck'); xlim([-90 +90])

subplot(222); 
plot(rlat,umbc05_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,rlat,umbc10_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,...
     rlat,umbc15_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,rlat,umbc20_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,...
     'linewidth',2);
  plotaxis2; %hl = legend('05','10','15','20','location','south','fontsize',8);
  title('Lapse'); xlim([-90 +90])

subplot(223); 
plot(rlat,umbc05_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,rlat,umbc10_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,...
     rlat,umbc15_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,rlat,umbc20_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,...
     'linewidth',2);
  plotaxis2; hl = legend('05','10','15','20','location','south','fontsize',8);
  title('Ozone'); xlim([-90 +90])

subplot(224); 
plot(rlat,umbc05_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,rlat,umbc10_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,...
     rlat,umbc15_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,rlat,umbc20_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,...
     'linewidth',2);
  plotaxis2; %hl = legend('05','10','15','20','location','south','fontsize',8);
  title('WV'); xlim([-90 +90])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3); clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf;
ta = tiledlayout(2,2,'TileSpacing','compact', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
plot(rlat,umbc05_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,rlat,umbc10_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,...
     rlat,umbc15_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,rlat,umbc20_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,...
     'linewidth',2);
  plotaxis2; box on; hl = legend('05','10','15','20','location','north','fontsize',8);
  xlim([-90 +90]); 

tafov(2) = nexttile;
plot(rlat,umbc05_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,rlat,umbc10_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,...
     rlat,umbc15_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,rlat,umbc20_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,...
     'linewidth',2);
  plotaxis2; box on;
  xlim([-90 +90])

tafov(3) = nexttile;
plot(rlat,umbc05_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,rlat,umbc10_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,...
     rlat,umbc15_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,rlat,umbc20_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,...
     'linewidth',2);
  plotaxis2; box on;
  xlim([-90 +90])

tafov(4) = nexttile;
plot(rlat,umbc05_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,rlat,umbc10_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,...
     rlat,umbc15_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,rlat,umbc20_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,...
     'linewidth',2);
  plotaxis2; box on;
  xlim([-90 +90])

% Remove all ytick labels except for 1st column
%for ii = [2 4]
%   tafov(ii).YTickLabel = '';
%   tafov(ii).YLabel.String = [];
%end

% Remove all xtick labels except for 3rd row
for ii = [1 2]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
end

xstr = 'Latitude';
ystr = 'W/m2/K';

tafov(1).YLabel.String = ystr; tafov(1).YLabel.FontSize = 18;
tafov(3).YLabel.String = ystr; tafov(3).YLabel.FontSize = 18;
tafov(3).XLabel.String = xstr; tafov(3).XLabel.FontSize = 18;
tafov(4).XLabel.String = xstr; tafov(4).XLabel.FontSize = 18;

plotoptions.str11 = 'Planck';
plotoptions.str12 = 'Lapse';
plotoptions.str21 = 'Ozone';
plotoptions.str22 = 'Water';

yposn = 0.8; %% inside, just below top line
yposn = 0.9; %% straddling top line
yposn = 1.0; %% outside, just above top line

title(tafov(1), plotoptions.str11, 'Units', 'normalized', 'Position', [0.5, yposn, 0]);
title(tafov(2), plotoptions.str12, 'Units', 'normalized', 'Position', [0.5, yposn, 0]);
title(tafov(3), plotoptions.str21, 'Units', 'normalized', 'Position', [0.5, yposn, 0]);
title(tafov(4), plotoptions.str22, 'Units', 'normalized', 'Position', [0.5, yposn, 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5); clf;
subplot(211)
plot(rlat,umbc05_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,rlat,umbc10_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,...
     rlat,umbc15_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,rlat,umbc20_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,...
     'linewidth',2);
  xlim([-90 +90]); ylim([-10 0])
  plotaxis2; %hl = legend('05','10','15','20','location','south','fontsize',8);
  title('Skt'); 

subplot(212)
plot(rlat,umbc05_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc05_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbc05_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbc05_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,...
     rlat,umbc10_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc10_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbc10_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbc10_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,...
     rlat,umbc15_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc15_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbc15_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbc15_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,...
     rlat,umbc20_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc20_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbc20_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbc20_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,...
     'linewidth',2);
  xlim([-90 +90]); ylim([-20 +20])
  plotaxis2; %hl = legend('05','10','15','20','location','south','fontsize',8);
  title('Sum Feedbacks'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
figure(1); aslprint([dir0 'feedbackparams_05_10_15_20yrs.pdf'])
%}

