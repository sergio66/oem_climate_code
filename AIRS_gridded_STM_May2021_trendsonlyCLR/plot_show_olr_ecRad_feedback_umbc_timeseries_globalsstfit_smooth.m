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
umbc05_spectral_olr.deltaSKT = junk.results(:,6);

ix = 2; iNumYears = 10; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '10 years ';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted;
umbc10_spectral_olr = junk.umbc_spectral_olr;
umbc10_spectral_olr.deltaSKT = junk.results(:,6);

ix = 3; iNumYears = 15; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '15 years ';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted;
umbc15_spectral_olr = junk.umbc_spectral_olr;
umbc15_spectral_olr.deltaSKT = junk.results(:,6);

ix = 4; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '20 years ';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt.globalSST_weighted;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2.globalSST_weighted;
umbc20_spectral_olr = junk.umbc_spectral_olr;
umbc20_spectral_olr.deltaSKT = junk.results(:,6);

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
YY = Y(:)';                     %% mean weighted delta SST rate = 0.031962  0.003305  0.023971  0.019594 K/yr for 05/10/15/20 years   WRONG
YY = Y'; YY = YY(:); YY = YY';  %% mean weighted delta SST rate = 0.069633  0.020002  0.028442  0.024870 K/yr for 05/10/15/20 years   CORRECT
  coslat  = cos(YY*pi/180);
  indSST = umbc05_spectral_olr.deltaSKT'; boo(1) = sum(indSST .* coslat)/sum(coslat);
  indSST = umbc10_spectral_olr.deltaSKT'; boo(2) = sum(indSST .* coslat)/sum(coslat);
  indSST = umbc15_spectral_olr.deltaSKT'; boo(3) = sum(indSST .* coslat)/sum(coslat);
  indSST = umbc20_spectral_olr.deltaSKT'; boo(4) = sum(indSST .* coslat)/sum(coslat);
  fprintf(1,'mean weighted delta SST rate = %8.6f  %8.6f  %8.6f  %8.6f K/yr for 05/10/15/20 years \n',boo)

disp('         Planck Lapse Ozone Water |  Total')
for ix = 1 : 4
  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f |  %5.2f \n',strfeedbacks{ix},showfeedbacks(ix,[1 2 3 4 7]));
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

factor = 1;

%% fixed on 9/2/2023
%% disp('since weighted global SST change for 10 years is about x10 smaller than rest, will divide by this factor')
%% disp('since weighted global SST change for 10 years is about x10 smaller than rest, will divide by this factor')
%% disp('since weighted global SST change for 10 years is about x10 smaller than rest, will divide by this factor')
%% factor = 0.1; 

figure(2); clf;
subplot(221); 
plot(rlat,smooth(umbc05_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc10_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),'x-',...
     rlat,smooth(umbc15_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc20_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; %hl = legend('05','10','15','20','location','south','fontsize',8);
  title('Planck'); xlim([-90 +90])

subplot(222); 
plot(rlat,smooth(umbc05_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc10_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),'x-',...
     rlat,smooth(umbc15_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc20_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; %hl = legend('05','10','15','20','location','south','fontsize',8);
  title('Lapse'); xlim([-90 +90])

subplot(223); 
plot(rlat,smooth(umbc05_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc10_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),'x-',...
     rlat,smooth(umbc15_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc20_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; hl = legend('05','10','15','20','location','south','fontsize',8);
  title('Ozone'); xlim([-90 +90])

subplot(224); 
plot(rlat,smooth(umbc05_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc10_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),'x-',...
     rlat,smooth(umbc15_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc20_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
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
plot(rlat,smooth(umbc05_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc10_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),'x-',...
     rlat,smooth(umbc15_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc20_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; box on;
  xlim([-90 +90]); 
  ylim([-40 +10])

tafov(2) = nexttile;
plot(rlat,smooth(umbc05_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc10_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),'x-',...
     rlat,smooth(umbc15_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc20_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; box on; hl = legend('05','10','15','20','location','north','fontsize',8);
  xlim([-90 +90])
  ylim([-05 +25])

tafov(3) = nexttile;
plot(rlat,smooth(umbc05_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc10_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),'x-',...
     rlat,smooth(umbc15_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc20_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  plotaxis2; box on;
  xlim([-90 +90])

tafov(4) = nexttile;
plot(rlat,smooth(umbc05_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc10_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),'x-',...
     rlat,smooth(umbc15_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc20_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
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
plot(rlat,smooth(umbc05_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc10_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,iSmooth),'x-',...
     rlat,smooth(umbc15_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,iSmooth),rlat,smooth(umbc20_spectral_olr.feedback_ecRad.skt.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  xlim([-90 +90]); ylim([-10 2])
  plotaxis2; %hl = legend('05','10','15','20','location','south','fontsize',8);
  title('Skt'); 

subplot(212)
plot(rlat,smooth(umbc05_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc05_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbc05_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbc05_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth)...
     rlat,smooth(umbc10_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc10_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbc10_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbc10_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),'x-',...
     rlat,smooth(umbc15_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc15_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbc15_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbc15_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     rlat,smooth(umbc20_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc20_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
          umbc20_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin     + umbc20_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     'linewidth',2);
  xlim([-90 +90]); ylim([-20 +20])
  plotaxis2; %hl = legend('05','10','15','20','location','south','fontsize',8);
  title('Sum Feedbacks'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iDebug10 = input('debug why 10 years is wierd? (-1/+1) : ');
if iDebug10 > 0
  figure(6); clf
  
  dbt = -1:0.01:+1; 
  dbt = -1:0.025:+1; 
  plot(dbt,histc(umbc05_spectral_olr.deltaSKT,dbt),dbt,histc(umbc10_spectral_olr.deltaSKT,dbt),dbt,histc(umbc15_spectral_olr.deltaSKT,dbt),dbt,histc(umbc20_spectral_olr.deltaSKT,dbt),'linewidth',2)
  semilogy(dbt,histc(umbc05_spectral_olr.deltaSKT,dbt),dbt,histc(umbc10_spectral_olr.deltaSKT,dbt),dbt,histc(umbc15_spectral_olr.deltaSKT,dbt),dbt,histc(umbc20_spectral_olr.deltaSKT,dbt),'linewidth',2)
  plotaxis2; hl = legend('05','10','15','20','location','best'); title('hist of d(SKT)'); xlabel('dSKT')
  
  plot(dbt,histc(coslat'.*umbc05_spectral_olr.deltaSKT,dbt),dbt,histc(coslat'.*umbc10_spectral_olr.deltaSKT,dbt),dbt,histc(coslat'.*umbc15_spectral_olr.deltaSKT,dbt),dbt,histc(coslat'.*umbc20_spectral_olr.deltaSKT,dbt),'linewidth',2)
  semilogy(dbt,histc(coslat'.*umbc05_spectral_olr.deltaSKT,dbt),dbt,histc(coslat'.*umbc10_spectral_olr.deltaSKT,dbt),dbt,histc(coslat'.*umbc15_spectral_olr.deltaSKT,dbt),dbt,histc(coslat'.*umbc20_spectral_olr.deltaSKT,dbt),'linewidth',2)
  plotaxis2; hl = legend('05','10','15','20','location','best'); title('hist of cos(rlat)*d(SKT)'); xlabel('cos(rlat)*dSKT')
  
  %plot(1:4608,umbc20_spectral_olr.deltaSKT,1:4608,umbc10_spectral_olr.deltaSKT)
  %plot(1:4608,umbc20_spectral_olr.planck_ecRad.clr-umbc20_spectral_olr.olr0_ecRad.clr,1:4608,umbc10_spectral_olr.planck_ecRad.clr-umbc10_spectral_olr.olr0_ecRad.clr)
  plot(YY,umbc05_spectral_olr.deltaSKT,YY,umbc10_spectral_olr.deltaSKT,YY,umbc15_spectral_olr.deltaSKT,YY,umbc20_spectral_olr.deltaSKT)
  plot(YY,umbc05_spectral_olr.deltaSKT,YY,umbc10_spectral_olr.deltaSKT)
  
  z11 = umbc05_spectral_olr.deltaSKT';
  z12 = umbc10_spectral_olr.deltaSKT';
  z21 = umbc15_spectral_olr.deltaSKT';
  z22 = umbc20_spectral_olr.deltaSKT';
  plotoptions.ystr = 'Latitude'; plotoptions.xstr = 'Latitude';
  plotoptions.cx = [-1 +1]*0.20; plotoptions.maintitle = 'dSKT/dt'; plotoptions.plotcolors = usa2; plotoptions.yReverseDir = -1; plotoptions.yLinearOrLog = +1; plotoptions.xLimits = [640  1640]; plotoptions.yLimits = [-90 +90];
  aslmap_2x2tiledlayout(z11,z12,z21,z22,7,plotoptions);
   
  %%%%%%%%%%%%%%%%%%%%%%%%%
  z11 = -(umbc05_spectral_olr.planck_ecRad.clr-umbc05_spectral_olr.olr0_ecRad.clr); %z11 = z11./umbc05_spectral_olr.deltaSKT';
  z12 = -(umbc10_spectral_olr.planck_ecRad.clr-umbc10_spectral_olr.olr0_ecRad.clr); %z12 = z11./umbc10_spectral_olr.deltaSKT';
  z21 = -(umbc15_spectral_olr.planck_ecRad.clr-umbc15_spectral_olr.olr0_ecRad.clr); %z21 = z21./umbc15_spectral_olr.deltaSKT';
  z22 = -(umbc20_spectral_olr.planck_ecRad.clr-umbc20_spectral_olr.olr0_ecRad.clr); %z22 = z22./umbc20_spectral_olr.deltaSKT';
  plotoptions.ystr = 'Latitude'; plotoptions.xstr = 'Latitude';
  plotoptions.cx = [-1 +1]*05; plotoptions.maintitle = '\delta(OLR) planck'; plotoptions.plotcolors = usa2; plotoptions.yReverseDir = -1; plotoptions.yLinearOrLog = +1; plotoptions.xLimits = [640  1640]; plotoptions.yLimits = [-90 +90];
  aslmap_2x2tiledlayout(z11,z12,z21,z22,8,plotoptions);
  
  z11 = -(umbc05_spectral_olr.planck_ecRad.clr-umbc05_spectral_olr.olr0_ecRad.clr); z11 = z11./umbc05_spectral_olr.deltaSKT';
  z12 = -(umbc10_spectral_olr.planck_ecRad.clr-umbc10_spectral_olr.olr0_ecRad.clr); z12 = z12./umbc10_spectral_olr.deltaSKT';
  z21 = -(umbc15_spectral_olr.planck_ecRad.clr-umbc15_spectral_olr.olr0_ecRad.clr); z21 = z21./umbc15_spectral_olr.deltaSKT';
  z22 = -(umbc20_spectral_olr.planck_ecRad.clr-umbc20_spectral_olr.olr0_ecRad.clr); z22 = z22./umbc20_spectral_olr.deltaSKT';
  plotoptions.ystr = 'Latitude'; plotoptions.xstr = 'Latitude';
  plotoptions.cx = [-1 0]*10; plotoptions.maintitle = '\lambda planck'; plotoptions.plotcolors = usa2; plotoptions.yReverseDir = -1; plotoptions.yLinearOrLog = +1; plotoptions.xLimits = [640  1640]; plotoptions.yLimits = [-90 +90];
  aslmap_2x2tiledlayout(z11,z12,z21,z22,9,plotoptions);
  
  zz11 = nanmean(reshape(z11,72,64),1); zz12 = nanmean(reshape(z12,72,64),1); zz21 = nanmean(reshape(z21,72,64),1); zz22 = nanmean(reshape(z22,72,64),1);
  figure(10); clf; plot(rlat,zz11,'b',rlat,zz12,'g',rlat,zz21,'r',rlat,zz22,'k','linewidth',2); hl = legend('05','10','15','20','location','best'); axis([-90 +90 -5 -2.5]); 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  z11 = -(umbc05_spectral_olr.wv_ecRad.clr-umbc05_spectral_olr.olr0_ecRad.clr); %z11 = z11./umbc05_spectral_olr.deltaSKT';
  z12 = -(umbc10_spectral_olr.wv_ecRad.clr-umbc10_spectral_olr.olr0_ecRad.clr); %z12 = z11./umbc10_spectral_olr.deltaSKT';
  z21 = -(umbc15_spectral_olr.wv_ecRad.clr-umbc15_spectral_olr.olr0_ecRad.clr); %z21 = z21./umbc15_spectral_olr.deltaSKT';
  z22 = -(umbc20_spectral_olr.wv_ecRad.clr-umbc20_spectral_olr.olr0_ecRad.clr); %z22 = z22./umbc20_spectral_olr.deltaSKT';
  plotoptions.ystr = 'Latitude'; plotoptions.xstr = 'Latitude';
  plotoptions.cx = [-1 +1]*05/2; plotoptions.maintitle = '\delta(OLR) wv'; plotoptions.plotcolors = usa2; plotoptions.yReverseDir = -1; plotoptions.yLinearOrLog = +1; plotoptions.xLimits = [640  1640]; plotoptions.yLimits = [-90 +90];
  aslmap_2x2tiledlayout(z11,z12,z21,z22,11,plotoptions);
  
  z11 = -(umbc05_spectral_olr.wv_ecRad.clr-umbc05_spectral_olr.olr0_ecRad.clr); z11 = z11./umbc05_spectral_olr.deltaSKT';
  z12 = -(umbc10_spectral_olr.wv_ecRad.clr-umbc10_spectral_olr.olr0_ecRad.clr); z12 = z12./umbc10_spectral_olr.deltaSKT';
  z21 = -(umbc15_spectral_olr.wv_ecRad.clr-umbc15_spectral_olr.olr0_ecRad.clr); z21 = z21./umbc15_spectral_olr.deltaSKT';
  z22 = -(umbc20_spectral_olr.wv_ecRad.clr-umbc20_spectral_olr.olr0_ecRad.clr); z22 = z22./umbc20_spectral_olr.deltaSKT';
  plotoptions.ystr = 'Latitude'; plotoptions.xstr = 'Latitude';
  plotoptions.cx = [-1 1]*10; plotoptions.maintitle = '\lambda wv'; plotoptions.plotcolors = usa2; plotoptions.yReverseDir = -1; plotoptions.yLinearOrLog = +1; plotoptions.xLimits = [640  1640]; plotoptions.yLimits = [-90 +90];
  aslmap_2x2tiledlayout(z11,z12,z21,z22,12,plotoptions);
  
  zz11 = nanmean(reshape(z11,72,64),1); zz12 = nanmean(reshape(z12,72,64),1); zz21 = nanmean(reshape(z21,72,64),1); zz22 = nanmean(reshape(z22,72,64),1);
  figure(13); clf; plot(rlat,smooth(zz11,10),'b',rlat,smooth(zz12,10),'g',rlat,smooth(zz21,10),'r',rlat,smooth(zz22,10),'k','linewidth',2); plotaxis2;
  hl = legend('05','10','15','20','location','best'); axis([-90 +90 -10 +10]); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
figure(1); aslprint([dir0 'feedbackparams_05_10_15_20yrs.pdf'])
%}

