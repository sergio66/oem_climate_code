addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath  /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS

load llsmap5

clear showfeedbacks* strfeedbacks
%clear all

for ii = 1 : 6; figure(ii); clf; end

ix = 1; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_DJF.mat']);
strfeedbacks{ix} = 'DJF';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted_all;
umbcDJF_spectral_olr = junk.umbc_spectral_olr;
umbcDJF_spectral_olr.deltaSKT = junk.results(:,6);

ix = 2; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_MAM.mat']);
strfeedbacks{ix} = 'MAM';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted_all;
umbcMAM_spectral_olr = junk.umbc_spectral_olr;
umbcMAM_spectral_olr.deltaSKT = junk.results(:,6);

ix = 3; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_JJA.mat']);
strfeedbacks{ix} = 'JJA';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted_all;
umbcJJA_spectral_olr = junk.umbc_spectral_olr;
umbcJJA_spectral_olr.deltaSKT = junk.results(:,6);

ix = 4; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_SON.mat']);
strfeedbacks{ix} = 'SON';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted_all;
umbcSON_spectral_olr = junk.umbc_spectral_olr;
umbcSON_spectral_olr.deltaSKT = junk.results(:,6);

ix = 5; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = 'ALL';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted_all;
umbcALL_spectral_olr = junk.umbc_spectral_olr;
umbcALL_spectral_olr.deltaSKT = junk.results(:,6);

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
ixx = ix;
showfeedbacks(1:ixx,7) = sum(showfeedbacks(1:ixx,[1 2 3 4]),2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('iSmooth')
  iSmooth = 5;
  iSmooth = 10;
end

if ~exist('rlat') | ~exist('Y')
  load latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
  [Y,X] = meshgrid(rlat,rlon);
  X = X; Y = Y;
end
%% oh oh may have this wrong before Sept 28 2023
%% mean weighted delta SST rate = 0.031962  0.003305  0.023971  0.019594 K/yr for DJF/MAM/JJA/SON seasons   WRONG
%% mean weighted delta SST rate = 0.069633  0.020002  0.028442  0.024870 K/yr for DJF/MAM/JJA/SON seasons   CORRECT
do_XX_YY_from_X_Y

  coslat  = cos(YY*pi/180);
  indSST = umbcDJF_spectral_olr.deltaSKT'; boo(1) = sum(indSST .* coslat)/sum(coslat);
  indSST = umbcMAM_spectral_olr.deltaSKT'; boo(2) = sum(indSST .* coslat)/sum(coslat);
  indSST = umbcJJA_spectral_olr.deltaSKT'; boo(3) = sum(indSST .* coslat)/sum(coslat);
  indSST = umbcSON_spectral_olr.deltaSKT'; boo(4) = sum(indSST .* coslat)/sum(coslat);
  indSST = umbcALL_spectral_olr.deltaSKT'; boo(5) = sum(indSST .* coslat)/sum(coslat);
fprintf(1,'mean weighted delta SST rate = %8.6f  %8.6f  %8.6f  %8.6f %8.6f K/yr for DJF/MAM/JJA/SON/ALL seasons \n',boo)

%%%%%%%%%%%%%%%%%%%%%%%%%

%% remember we do delta(OLR) so we need 3xolr0, and we actually do not need planck when summing
%% see compute_feedbacks_regress_olr_ecRad_calcs.m : 
%%  Planck = -F_planck + F_0
%%  Lapse  = -F_Lapse + F_planck
%%  ozone  = -F_ozone + F_0 
%%  water  = -F_water + F_0 
%% -------------------------
%%   SUM  = +3*F_0 + 0*F_planck -F_Lapse -F_ozone -F_water
%% -------------------------%% ------------------------- 
umbcDJF_spectral_olr.allsum = 3*umbcDJF_spectral_olr.olr0_ecRad.clr - (0*umbcDJF_spectral_olr.planck_ecRad.clr + umbcDJF_spectral_olr.lapse_ecRad.clr  + umbcDJF_spectral_olr.o3_ecRad.clr + umbcDJF_spectral_olr.wv_ecRad.clr);
umbcMAM_spectral_olr.allsum = 3*umbcMAM_spectral_olr.olr0_ecRad.clr - (0*umbcMAM_spectral_olr.planck_ecRad.clr + umbcMAM_spectral_olr.lapse_ecRad.clr  + umbcMAM_spectral_olr.o3_ecRad.clr + umbcMAM_spectral_olr.wv_ecRad.clr);
umbcJJA_spectral_olr.allsum = 3*umbcJJA_spectral_olr.olr0_ecRad.clr - (0*umbcJJA_spectral_olr.planck_ecRad.clr + umbcJJA_spectral_olr.lapse_ecRad.clr  + umbcJJA_spectral_olr.o3_ecRad.clr + umbcJJA_spectral_olr.wv_ecRad.clr);
umbcSON_spectral_olr.allsum = 3*umbcSON_spectral_olr.olr0_ecRad.clr - (0*umbcSON_spectral_olr.planck_ecRad.clr + umbcSON_spectral_olr.lapse_ecRad.clr  + umbcSON_spectral_olr.o3_ecRad.clr + umbcSON_spectral_olr.wv_ecRad.clr);
umbcALL_spectral_olr.allsum = 3*umbcALL_spectral_olr.olr0_ecRad.clr - (0*umbcALL_spectral_olr.planck_ecRad.clr + umbcALL_spectral_olr.lapse_ecRad.clr  + umbcALL_spectral_olr.o3_ecRad.clr + umbcALL_spectral_olr.wv_ecRad.clr);
z11 = umbcDJF_spectral_olr.allsum;
z12 = umbcMAM_spectral_olr.allsum;
zMID = umbcALL_spectral_olr.allsum;
z21 = umbcJJA_spectral_olr.allsum;
z22 = umbcSON_spectral_olr.allsum;
plotoptions.str11 = 'DJF'; plotoptions.str12 = 'MAM';
plotoptions.str21 = 'JJA'; plotoptions.str22 = 'SON';
plotoptions.ystr = 'Latitude'; plotoptions.xstr = 'Longitude';
plotoptions.cx = [-1 +1]*20; plotoptions.maintitle = '\lambda'; plotoptions.cmap = llsmap5; plotoptions.yReverseDir = -1; plotoptions.yLinearOrLog = +1;
aslmap_2x2tiledlayout(z11/boo(1),z12/boo(2),z21/boo(3),z22/boo(4),7,plotoptions);
plotoptions5.str11 = 'DJF'; plotoptions5.str12 = 'MAM';
plotoptions5.strzz = 'ALL'; 
plotoptions5.str21 = 'JJA'; plotoptions5.str22 = 'SON';
plotoptions5.ystr = 'Latitude'; plotoptions5.xstr = 'Longitude';
plotoption5.maintitle = '\lambda W/m2/K';
plotoptions5.cx = [-1 +1]*20; plotoptions5.maintitle = '\lambda'; plotoptions5.cmap = llsmap5; plotoptions5.yReverseDir = -1; plotoptions5.yLinearOrLog = +1;
aslmap_2x1x2tiledlayout(z11/boo(1),z12/boo(2),zMID/boo(5),z21/boo(3),z22/boo(4),7,plotoptions5);

%junk = [sum(umbcDJF_spectral_olr.allsum .* coslat)/sum(coslat) sum(umbcMAM_spectral_olr.allsum .* coslat)/sum(coslat) ...
%         sum(umbcJJA_spectral_olr.allsum .* coslat)/sum(coslat) sum(umbcSON_spectral_olr.allsum .* coslat)/sum(coslat)] ./ boo(1:4);
%fprintf(1,'feedbacks for DJF/MAM/JJA/SON = %8.5f %8.5f %8.5f %8.5f W/m2/K \n',junk);
junk = [sum(umbcDJF_spectral_olr.allsum .* coslat)/sum(coslat) sum(umbcMAM_spectral_olr.allsum .* coslat)/sum(coslat) sum(umbcJJA_spectral_olr.allsum .* coslat)/sum(coslat) ...
        sum(umbcSON_spectral_olr.allsum .* coslat)/sum(coslat) sum(umbcALL_spectral_olr.allsum .* coslat)/sum(coslat)] ./ boo(1:5);
fprintf(1,'feedbacks for DJF/MAM/JJA/SON/ALL = %8.5f %8.5f %8.5f %8.5f %8.5f W/m2/K \n',junk);

%%%%%%%%%%%%%%%%%%%%%%%%%

%disp('         Planck Lapse Ozone Water |  Total')
%for ix = 1 : 5
%  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f |  %5.2f \n',strfeedbacks{ix},showfeedbacks(ix,[1 2 3 4 7]));
%end
%trends_paper_show = showfeedbacks(1:ixx,[1 2 3 4 7]);

showfeedbacks(1:5,8) = junk; 
disp('         Planck Lapse Ozone Water |  Total')
for ix = 1 : 5
  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f |  %5.2f \n',strfeedbacks{ix},showfeedbacks(ix,[1 2 3 4 8]));
end
trends_paper_show = showfeedbacks(1:ixx,[1 2 3 4 8]);

figure(1); clf
bar(trends_paper_show');
ylabel('Feedback W/m2/K');
hl = legend('DJF','MAM','JJA','SON','ALL','location','north','fontsize',10);
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factor = 1;

xsin = sin(rlat*pi/180);
  xtick = [-1 -sqrt(3)/2 -sqrt(2)/2 -1/2 -(0.25+0.01) 0 +(0.25+0.01) +1/2 +sqrt(2)/2 +sqrt(3)/2 +1]; %% -90 -60 -45 -30 -15 0 +15 +30 +45 +60 +90
  xtick = [-1            -sqrt(2)/2      -(0.25+0.01) 0 +(0.25+0.01)      +sqrt(2)/2            +1]; %% -90     -45     -15 0 +15     +45     +90
  xticklab = cellstr(num2str(round(180/pi*asin((xtick(:)))), '%d'));

%% fixed on 9/2/2023
%% disp('since weighted global SST change for 10 years is about x10 smaller than rest, will divide by this factor')
%% disp('since weighted global SST change for 10 years is about x10 smaller than rest, will divide by this factor')
%% disp('since weighted global SST change for 10 years is about x10 smaller than rest, will divide by this factor')
%% factor = 0.1; 

figure(2); clf;
subplot(221); 
plot(xsin,smooth(umbcDJF_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcMAM_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcJJA_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcSON_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('Planck'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(222); 
plot(xsin,smooth(umbcDJF_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcMAM_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcJJA_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcSON_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('Lapse'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(223); 
plot(xsin,smooth(umbcDJF_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcMAM_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcJJA_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcSON_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('Ozone'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(224); 
plot(xsin,smooth(umbcDJF_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcMAM_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcJJA_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcSON_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('WV'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3); clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf;
ta = tiledlayout(2,2,'TileSpacing','compact', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
plot(xsin,smooth(umbcDJF_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcMAM_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcJJA_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcSON_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; 
  ylim([-40 +10])

tafov(2) = nexttile;
plot(xsin,smooth(umbcDJF_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcMAM_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcJJA_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcSON_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; box on; hl = legend('DJF','MAM','JJA','SON','location','north','fontsize',8);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');
  ylim([-05 +25])

tafov(3) = nexttile;
plot(xsin,smooth(umbcDJF_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcMAM_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcJJA_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcSON_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(4) = nexttile;
plot(xsin,smooth(umbcDJF_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcMAM_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcJJA_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcSON_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
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
figure(5); clf

xtick = [-1            -sqrt(2)/2      -(0.25+0.01) 0 +(0.25+0.01)      +sqrt(2)/2            +1]; %% -90     -45     -15 0 +15     +45     +90
  xtick = [-1 -sqrt(3)/2 -sqrt(2)/2 -1/2 -(0.25+0.01) 0 +(0.25+0.01) +1/2 +sqrt(2)/2 +sqrt(3)/2 +1]; %% -90 -60 -45 -30 -15 0 +15 +30 +45 +60 +90
  xticklab = cellstr(num2str(round(180/pi*asin((xtick(:)))), '%d'));

plot(xsin,smooth(umbcDJF_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbcDJF_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbcDJF_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbcDJF_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcMAM_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbcMAM_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbcMAM_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbcMAM_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcJJA_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbcJJA_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbcJJA_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbcJJA_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcSON_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbcSON_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbcSON_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbcSON_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  ylim([-1 +1]*7.5)
  plotaxis2; hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('Sum Feedbacks'); 
  set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex'); xlim([-1 +1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6); clf;
subplot(211)
plot(xsin,smooth(umbcDJF_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcMAM_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcJJA_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,iSmooth),xsin,smooth(umbcSON_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; ylim([-10 2])
  plotaxis2; %hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('Skt'); 

subplot(212)
plot(xsin,smooth(umbcDJF_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbcDJF_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbcDJF_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbcDJF_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcMAM_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbcMAM_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbcMAM_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbcMAM_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcJJA_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbcJJA_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbcJJA_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbcJJA_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbcSON_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbcSON_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbcSON_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbcSON_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
    xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; 
    ylim([-1 +1]*7.5)
  plotaxis2; hl = legend('DJF','MAM','JJA','SON','location','south','fontsize',8);
  title('Sum Feedbacks'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
figure(4); aslprint([dir0 'global_feedbackparams_seasonal_20yrs.pdf'])
figure(7); aslprint([dir0 'global_feedbackparams_seasonal_20yrs_latlon.pdf'])
%}
