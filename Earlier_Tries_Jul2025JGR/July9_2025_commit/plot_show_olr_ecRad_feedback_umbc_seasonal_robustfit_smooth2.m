clear showfeedbacks* strfeedbacks

ix = 1; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_DJF.mat']);
strfeedbacks{ix} = 'DJF';
%% the value
showfeedbacks_robustfit_all(ix,1,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all(ix,2,1) = junk.umbc_spectral_olr.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all(ix,3,1) = junk.umbc_spectral_olr.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all(ix,4,1) = junk.umbc_spectral_olr.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all(ix,5,1) = junk.umbc_spectral_olr.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all(ix,6,1) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all(ix,1,2) = junk.umbc_spectral_olr.feedback_ecRad.planck.robustfit_all(2);
showfeedbacks_robustfit_all(ix,2,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all(ix,3,2) = junk.umbc_spectral_olr.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all(ix,4,2) = junk.umbc_spectral_olr.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all(ix,5,2) = junk.umbc_spectral_olr.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all(ix,6,2) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.robustfit_all(2);
umbc05_spectral_olr = junk.umbc_spectral_olr;
umbc05_spectral_olr.deltaSKT = junk.results(:,6);

ix = 2; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_MAM.mat']);
strfeedbacks{ix} = 'MAM';
%% the value
showfeedbacks_robustfit_all(ix,1,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all(ix,2,1) = junk.umbc_spectral_olr.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all(ix,3,1) = junk.umbc_spectral_olr.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all(ix,4,1) = junk.umbc_spectral_olr.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all(ix,5,1) = junk.umbc_spectral_olr.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all(ix,6,1) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all(ix,1,2) = junk.umbc_spectral_olr.feedback_ecRad.planck.robustfit_all(2);
showfeedbacks_robustfit_all(ix,2,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all(ix,3,2) = junk.umbc_spectral_olr.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all(ix,4,2) = junk.umbc_spectral_olr.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all(ix,5,2) = junk.umbc_spectral_olr.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all(ix,6,2) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.robustfit_all(2);
umbc10_spectral_olr = junk.umbc_spectral_olr;
umbc10_spectral_olr.deltaSKT = junk.results(:,6);

ix = 3; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_JJA.mat']);
strfeedbacks{ix} = 'JJA';
%% the value
showfeedbacks_robustfit_all(ix,1,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all(ix,2,1) = junk.umbc_spectral_olr.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all(ix,3,1) = junk.umbc_spectral_olr.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all(ix,4,1) = junk.umbc_spectral_olr.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all(ix,5,1) = junk.umbc_spectral_olr.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all(ix,6,1) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all(ix,1,2) = junk.umbc_spectral_olr.feedback_ecRad.planck.robustfit_all(2);
showfeedbacks_robustfit_all(ix,2,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all(ix,3,2) = junk.umbc_spectral_olr.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all(ix,4,2) = junk.umbc_spectral_olr.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all(ix,5,2) = junk.umbc_spectral_olr.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all(ix,6,2) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.robustfit_all(2);
umbc15_spectral_olr = junk.umbc_spectral_olr;
umbc15_spectral_olr.deltaSKT = junk.results(:,6);

ix = 4; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_SON.mat']);
strfeedbacks{ix} = 'SON';
%% the value
showfeedbacks_robustfit_all(ix,1,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all(ix,2,1) = junk.umbc_spectral_olr.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all(ix,3,1) = junk.umbc_spectral_olr.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all(ix,4,1) = junk.umbc_spectral_olr.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all(ix,5,1) = junk.umbc_spectral_olr.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all(ix,6,1) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all(ix,1,2) = junk.umbc_spectral_olr.feedback_ecRad.planck.robustfit_all(2);
showfeedbacks_robustfit_all(ix,2,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all(ix,3,2) = junk.umbc_spectral_olr.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all(ix,4,2) = junk.umbc_spectral_olr.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all(ix,5,2) = junk.umbc_spectral_olr.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all(ix,6,2) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.robustfit_all(2);
umbc20_spectral_olr = junk.umbc_spectral_olr;
umbc20_spectral_olr.deltaSKT = junk.results(:,6);

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
ixx = ix;
showfeedbacks_robustfit_all(1:ixx,7,1) = sum(squeeze(showfeedbacks_robustfit_all(1:ixx,[1 2 3 4],1)),2);
junk = showfeedbacks_robustfit_all(1:ixx,[1 2 3 4],2);
junk = sqrt(sum(junk.*junk,2));
showfeedbacks_robustfit_all(1:ixx,7,2) = junk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rlat') | ~exist('Y')
  do_XX_YY_from_X_Y
end
%% oh oh may have this wrong before Sept 28 2023
%% mean weighted delta SST rate = 0.031962  0.003305  0.023971  0.019594 K/yr for DJF/MAM/JJA/SON seasons   WRONG
%% mean weighted delta SST rate = 0.069633  0.020002  0.028442  0.024870 K/yr for DJF/MAM/JJA/SON seasons   CORRECT

  coslat  = cos(YY*pi/180);
  indSST = umbc05_spectral_olr.deltaSKT'; boo(1) = sum(indSST .* coslat)/sum(coslat);
  indSST = umbc10_spectral_olr.deltaSKT'; boo(2) = sum(indSST .* coslat)/sum(coslat);
  indSST = umbc15_spectral_olr.deltaSKT'; boo(3) = sum(indSST .* coslat)/sum(coslat);
  indSST = umbc20_spectral_olr.deltaSKT'; boo(4) = sum(indSST .* coslat)/sum(coslat);
  fprintf(1,'mean weighted delta SST rate = %8.6f  %8.6f  %8.6f  %8.6f K/yr for DJF/MAM/JJA/SON seasons \n',boo)

disp('         Planck Lapse Ozone Water |  Total')
for ix = 1 : 4
  %fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},showfeedbacks_robustfit_all(ix,[1 2 3 4 7],1));
  junk = [showfeedbacks_robustfit_all(ix,[1 2 3 4 7],1) showfeedbacks_robustfit_all(ix,[1 2 3 4 7],2)];
  junk = junk([1 6 2 7 3 8 4 9 5 10]);
  fprintf(1,'%s %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f    %6.3f +/- %6.3f \n',strfeedbacks{ix},junk);
end
trends_paper_show = showfeedbacks_robustfit_all(1:ixx,[1 2 3 4 7]);

figure(1); clf
bar(trends_paper_show');
ylabel('Feedback W/m2/K');
hl = legend('05','10','15','20','location','south');
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('iSmooth')
  iSmooth = 5;
  iSmooth = 10;
end

if ~exist('rlat')
  load latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
  [Y,X] = meshgrid(rlat,rlon);
  X = X; Y = Y;
end

xsin = sin(rlat*pi/180);
  xtick = [-1 -sqrt(3)/2 -sqrt(2)/2 -1/2 -(0.25+0.01) 0 +(0.25+0.01) +1/2 +sqrt(2)/2 +sqrt(3)/2 +1]; %% -90 -60 -45 -30 -15 0 +15 +30 +45 +60 +90
  xtick = [-1            -sqrt(2)/2      -(0.25+0.01) 0 +(0.25+0.01)      +sqrt(2)/2            +1]; %% -90     -45     -15 0 +15     +45     +90
  xticklab = cellstr(num2str(round(180/pi*asin((xtick(:)))), '%d'));

figure(2); clf;
subplot(221); 
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     'linewidth',2);
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('Planck'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(222); 
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     'linewidth',2);
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('Lapse'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(223); 
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     'linewidth',2);
  plotaxis2; hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('Ozone'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(224); 
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     'linewidth',2);
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('WV'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf;
subplot(221); 
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,2),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,2),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,2),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,2),iSmooth),...
     'linewidth',2);
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('unc Planck'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(222); 
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,2),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,2),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,2),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,2),iSmooth),...
     'linewidth',2);
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('unc Lapse'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(223); 
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,2),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,2),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,2),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,2),iSmooth),...
     'linewidth',2);
  plotaxis2; hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('unc Ozone'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(224); 
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,2),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,2),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,2),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,2),iSmooth),...
     'linewidth',2);
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('unc WV'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf;
ta = tiledlayout(2,2,'TileSpacing','compact', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     'linewidth',2);
  plotaxis2; box on; hl = legend('DJF','MAM','JJA','SON','location','north','fontsize',8);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; 

tafov(2) = nexttile;
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(3) = nexttile;
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(4) = nexttile;
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

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
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.skt.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc10_spectral_olr.feedback_ecRad.skt.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.skt.robustfit_latbin(:,1),iSmooth),xsin,smooth(umbc20_spectral_olr.feedback_ecRad.skt.robustfit_latbin(:,1),iSmooth),...
     'linewidth',2);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; ylim([-3 0])
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('Skt'); 

subplot(212)
plot(xsin,smooth(umbc05_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + umbc05_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
                 umbc05_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1)     + umbc05_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc10_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + umbc10_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
                 umbc10_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1)     + umbc10_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc15_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + umbc15_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
                 umbc15_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1)     + umbc15_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc20_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + umbc20_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
                 umbc20_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1)     + umbc20_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     'linewidth',2);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; ylim([-3 +3])
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('Sum Feedbacks'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
figure(4); aslprint([dir0 'feedbackparams_seasonal_20yrs.pdf'])
%}

