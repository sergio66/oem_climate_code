clear showfeedbacks* strfeedbacks

if ~exist('iNumYears')
  disp('WARNING iNumYears DNE ... setting to 20')
  iNumYears = 20;
end

if ~exist('era5_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'era5_spectral_olr');
  era5_spectral_olr = junk.era5_spectral_olr;
end
ix = 1; junk = era5_spectral_olr;
strfeedbacks{ix} = 'ERA5       ';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck_ecRad_globalSST_weighted;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse_ecRad_globalSST_weighted;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3_ecRad_globalSST_weighted;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv_ecRad_globalSST_weighted;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt_ecRad_globalSST_weighted;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2_ecRad_globalSST_weighted;

if ~exist('merra2_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'merra2_spectral_olr');
  merra2_spectral_olr = junk.merra2_spectral_olr;
end
ix = 2; junk = merra2_spectral_olr;
strfeedbacks{ix} = 'MERRA2     ';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck_ecRad_globalSST_weighted;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse_ecRad_globalSST_weighted;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3_ecRad_globalSST_weighted;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv_ecRad_globalSST_weighted;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt_ecRad_globalSST_weighted;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2_ecRad_globalSST_weighted;

if ~exist('umbc_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'],'umbc_spectral_olr');
  umbc_spectral_olr = junk.umbc_spectral_olr;
end
ix = 3; junk = umbc_spectral_olr;
strfeedbacks{ix} = 'THIS WORK  ';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck_ecRad_globalSST_weighted;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse_ecRad_globalSST_weighted;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3_ecRad_globalSST_weighted;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv_ecRad_globalSST_weighted;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt_ecRad_globalSST_weighted;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2_ecRad_globalSST_weighted;

if ~exist('airsL3_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'airsL3_spectral_olr');
  airsL3_spectral_olr = junk.airsL3_spectral_olr;
end
ix = 4; junk = airsL3_spectral_olr;
strfeedbacks{ix} = 'AIRS L3    ';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck_ecRad_globalSST_weighted;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse_ecRad_globalSST_weighted;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3_ecRad_globalSST_weighted;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv_ecRad_globalSST_weighted;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt_ecRad_globalSST_weighted;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2_ecRad_globalSST_weighted;

if ~exist('climcapsL3_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'climcapsL3_spectral_olr');
  climcapsL3_spectral_olr = junk.climcapsL3_spectral_olr;
end
ix = 5; junk = climcapsL3_spectral_olr;
strfeedbacks{ix} = 'CLIMCAPS L3';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck_ecRad_globalSST_weighted;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse_ecRad_globalSST_weighted;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3_ecRad_globalSST_weighted;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv_ecRad_globalSST_weighted;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt_ecRad_globalSST_weighted;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2_ecRad_globalSST_weighted;

if ~exist('cmip6_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'cmip6_spectral_olr');
  cmip6_spectral_olr = junk.cmip6_spectral_olr;
end
ix = 6; junk = cmip6_spectral_olr;
strfeedbacks{ix} = 'CMIP6      ';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck_ecRad_globalSST_weighted;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse_ecRad_globalSST_weighted;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3_ecRad_globalSST_weighted;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv_ecRad_globalSST_weighted;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt_ecRad_globalSST_weighted;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2_ecRad_globalSST_weighted;

if ~exist('amip6_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'amip6_spectral_olr');
  amip6_spectral_olr = junk.amip6_spectral_olr;
end
ix = 7; junk = amip6_spectral_olr;
strfeedbacks{ix} = 'AMIP6      ';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck_ecRad_globalSST_weighted;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse_ecRad_globalSST_weighted;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3_ecRad_globalSST_weighted;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv_ecRad_globalSST_weighted;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt_ecRad_globalSST_weighted;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2_ecRad_globalSST_weighted;

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
ixx = ix;
showfeedbacks(1:ixx,7) = sum(showfeedbacks(1:ixx,[1 2 3 4]),2);

for ix = 1 : 5
  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},showfeedbacks(ix,[1 2 3 4 7]));
end
trends_paper_show = showfeedbacks(1:5,[1 2 3 4 7]);

figure(1); clf
bar(trends_paper_show')
ylabel('Feedback W/m2/K');
hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south');
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)

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
subplot(221); plot(rlat,era5_spectral_olr.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin,rlat,merra2_spectral_olr.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin,...
                     rlat,umbc_spectral_olr.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin,rlat,airsL3_spectral_olr.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin,...
                     rlat,climcapsL3_spectral_olr.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin);
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('Planck'); xlim([-90 +90])

subplot(222); plot(rlat,era5_spectral_olr.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin,rlat,merra2_spectral_olr.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin,...
                     rlat,umbc_spectral_olr.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin,rlat,airsL3_spectral_olr.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin,...
                     rlat,climcapsL3_spectral_olr.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin)
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('Lapse'); xlim([-90 +90])

subplot(223); plot(rlat,era5_spectral_olr.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin,rlat,merra2_spectral_olr.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin,...
                     rlat,umbc_spectral_olr.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin,rlat,airsL3_spectral_olr.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin,...
                     rlat,climcapsL3_spectral_olr.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin)
  plotaxis2; hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('Ozone'); xlim([-90 +90])

subplot(224); plot(rlat,era5_spectral_olr.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin,rlat,merra2_spectral_olr.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin,...
                     rlat,umbc_spectral_olr.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin,rlat,airsL3_spectral_olr.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin,...
                     rlat,climcapsL3_spectral_olr.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin)
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('WV'); xlim([-90 +90])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf;
ta = tiledlayout(2,2,'TileSpacing','compact', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
plot(rlat,era5_spectral_olr.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin,rlat,merra2_spectral_olr.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin,...
     rlat,umbc_spectral_olr.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin,rlat,airsL3_spectral_olr.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin,...
     rlat,climcapsL3_spectral_olr.feedback_ecRad.planck_ecRad_globalSST_weighted_latbin,'linewidth',2);
  plotaxis2; box on; hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','north','fontsize',8);
  xlim([-90 +90]); 

tafov(2) = nexttile;
plot(rlat,era5_spectral_olr.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin,rlat,merra2_spectral_olr.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin,...
     rlat,umbc_spectral_olr.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin,rlat,airsL3_spectral_olr.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin,...
     rlat,climcapsL3_spectral_olr.feedback_ecRad.lapse_ecRad_globalSST_weighted_latbin,'linewidth',2);
  plotaxis2; box on;
  xlim([-90 +90])

tafov(3) = nexttile;
plot(rlat,era5_spectral_olr.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin,rlat,merra2_spectral_olr.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin,...
     rlat,umbc_spectral_olr.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin,rlat,airsL3_spectral_olr.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin,...
     rlat,climcapsL3_spectral_olr.feedback_ecRad.o3_ecRad_globalSST_weighted_latbin,'linewidth',2);
  plotaxis2; box on;
  xlim([-90 +90])

tafov(4) = nexttile;
plot(rlat,era5_spectral_olr.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin,rlat,merra2_spectral_olr.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin,...
     rlat,umbc_spectral_olr.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin,rlat,airsL3_spectral_olr.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin,...
      rlat,climcapsL3_spectral_olr.feedback_ecRad.wv_ecRad_globalSST_weighted_latbin,'linewidth',2);
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

%{
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
figure(4); aslprint([dir0 'feedbackparams_20yrs_UMBC_ERA5_MERRA2_L3.pdf'])
%}

