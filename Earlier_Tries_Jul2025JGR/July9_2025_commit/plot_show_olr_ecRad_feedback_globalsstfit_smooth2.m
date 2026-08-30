addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath  /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS

load llsmap5

clear showfeedbacks* strfeedbacks
%clear all

if ~exist('iNumYears')
  disp('WARNING iNumYears DNE ... setting to 20')
  iNumYears = 20;
end

if ~exist('era5_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'era5_spectral_olr','stemptrend');
  era5_spectral_olr = junk.era5_spectral_olr;
  era5_spectral_olr.stemptrend = junk.stemptrend.era5;
end
ix = 1; junk = era5_spectral_olr;
strfeedbacks{ix} = 'ERA5       ';
thename{ix}      = 'era5';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_all;

if ~exist('merra2_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'merra2_spectral_olr','stemptrend');
  merra2_spectral_olr = junk.merra2_spectral_olr;
  merra2_spectral_olr.stemptrend = junk.stemptrend.era5;
end
ix = 2; junk = merra2_spectral_olr;
strfeedbacks{ix} = 'MERRA2     ';
thename{ix}      = 'merra2';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_all;

if ~exist('umbc_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'],'umbc_spectral_olr');
  umbc_spectral_olr = junk.umbc_spectral_olr;

  a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';      iNumYears = 20;  %% use CarbonTracker CO2 trends ****, topts.iAdjLowerAtmWVfrac=0.25
  junk = load(strUMBC,'results');
  umbc_spectral_olr.stemptrend = junk.results(:,6)';
end
ix = 3; junk = umbc_spectral_olr;
strfeedbacks{ix} = 'THIS WORK  ';
thename{ix}      = 'umbc';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_all;

if ~exist('airsL3_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'airsL3_spectral_olr','stemptrend');
  airsL3_spectral_olr = junk.airsL3_spectral_olr;
  airsL3_spectral_olr.stemptrend = junk.stemptrend.airsL3;
end
ix = 4; junk = airsL3_spectral_olr;
strfeedbacks{ix} = 'AIRS L3    ';
thename{ix}      = 'airsL3';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_all;

if ~exist('climcapsL3_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'climcapsL3_spectral_olr','stemptrend');
  climcapsL3_spectral_olr = junk.climcapsL3_spectral_olr;
  climcapsL3_spectral_olr.stemptrend = junk.stemptrend.airsL3;
end
ix = 5; junk = climcapsL3_spectral_olr;
strfeedbacks{ix} = 'CLIMCAPS L3';
thename{ix}      = 'climcapsL3';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_all;

if ~exist('cmip6_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'cmip6_spectral_olr','stemptrend');
  cmip6_spectral_olr = junk.cmip6_spectral_olr;
  cmip6_spectral_olr.stemptrend = junk.stemptrend.cmip6;
end
ix = 6; junk = cmip6_spectral_olr;
strfeedbacks{ix} = 'CMIP6      ';
thename{ix}      = 'cmip6';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_all;

if ~exist('amip6_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'amip6_spectral_olr','stemptrend');
  amip6_spectral_olr = junk.amip6_spectral_olr;
  amip6_spectral_olr.stemptrend = junk.stemptrend.cmip6;
end
ix = 7; junk = amip6_spectral_olr;
strfeedbacks{ix} = 'AMIP6      ';
thename{ix}      = 'amip6';
showfeedbacks(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_all;
showfeedbacks(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_all;
showfeedbacks(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_all;
showfeedbacks(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_all;
showfeedbacks(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_all;
showfeedbacks(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_all;

clear junk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear rlat
if ~exist('rlat')
  do_XX_YY_from_X_Y
end

if ~exist('iSmooth')
  iSmooth = 5;
  iSmooth = 10;
end

clear boo
  coslat  = cos(YY*pi/180);
  indSST = era5_spectral_olr.stemptrend; boo(1) = sum(indSST .* coslat)/sum(coslat);
  indSST = merra2_spectral_olr.stemptrend; boo(2) = sum(indSST .* coslat)/sum(coslat);
  indSST = umbc_spectral_olr.stemptrend; boo(3) = sum(indSST .* coslat)/sum(coslat);
  indSST = airsL3_spectral_olr.stemptrend; boo(4) = sum(indSST .* coslat)/sum(coslat);
  indSST = climcapsL3_spectral_olr.stemptrend; boo(5) = sum(indSST .* coslat)/sum(coslat);
  fprintf(1,'mean weighted delta SST rate = %8.6f  %8.6f  %8.6f  %8.6f  %8.6f K/yr for ERA5/MERRA2/UMBC/AIRSL3/CLIMCAPSL3 \n',boo)

%%%%%%%%%%%%%%%%%%%%%%%%%

boox = input('For feedback plots, Divide sum(all tile) by 1/meanSST or divide by 1 (-1[default]/+1) : ');
if length(boox) == 0
  boox = -1;
end
if boox == -1
  boox = boo;
  dascale = 1;
else
  boox = ones(size(boo));
  dascale = 1/20;
end

%% remember we do delta(OLR) so we need 3xolr0, and we actually do not need planck when summing
%% see compute_feedbacks_regress_olr_ecRad_calcs.m : 
%%  Planck = -F_planck + F_0
%%  Lapse  = -F_Lapse + F_planck
%%  ozone  = -F_ozone + F_0 
%%  water  = -F_water + F_0 
%% -------------------------
%%   SUM  = +3*F_0 + 0*F_planck -F_Lapse -F_ozone -F_water
%% -------------------------%% ------------------------- 
era5_spectral_olr.allsum       = 3*era5_spectral_olr.olr0_ecRad.clr - (0*era5_spectral_olr.planck_ecRad.clr + era5_spectral_olr.lapse_ecRad.clr  + era5_spectral_olr.o3_ecRad.clr + era5_spectral_olr.wv_ecRad.clr);
merra2_spectral_olr.allsum     = 3*merra2_spectral_olr.olr0_ecRad.clr - (0*merra2_spectral_olr.planck_ecRad.clr + merra2_spectral_olr.lapse_ecRad.clr  + merra2_spectral_olr.o3_ecRad.clr + merra2_spectral_olr.wv_ecRad.clr);
umbc_spectral_olr.allsum       = 3*umbc_spectral_olr.olr0_ecRad.clr - (0*umbc_spectral_olr.planck_ecRad.clr + umbc_spectral_olr.lapse_ecRad.clr  + umbc_spectral_olr.o3_ecRad.clr + umbc_spectral_olr.wv_ecRad.clr);
airsL3_spectral_olr.allsum     = 3*airsL3_spectral_olr.olr0_ecRad.clr - (0*airsL3_spectral_olr.planck_ecRad.clr + airsL3_spectral_olr.lapse_ecRad.clr  + airsL3_spectral_olr.o3_ecRad.clr + airsL3_spectral_olr.wv_ecRad.clr);
climcapsL3_spectral_olr.allsum = 3*climcapsL3_spectral_olr.olr0_ecRad.clr - (0*climcapsL3_spectral_olr.planck_ecRad.clr + climcapsL3_spectral_olr.lapse_ecRad.clr  + climcapsL3_spectral_olr.o3_ecRad.clr + climcapsL3_spectral_olr.wv_ecRad.clr);
junk = [sum(era5_spectral_olr.allsum .* coslat)/sum(coslat) sum(merra2_spectral_olr.allsum .* coslat)/sum(coslat) ...
        sum(umbc_spectral_olr.allsum .* coslat)/sum(coslat) sum(airsL3_spectral_olr.allsum .* coslat)/sum(coslat) sum(climcapsL3_spectral_olr.allsum .* coslat)/sum(coslat)]  ./ boo(1:5);

z11  = era5_spectral_olr.allsum;
z12  = merra2_spectral_olr.allsum;
zMID = umbc_spectral_olr.allsum;
z21  = airsL3_spectral_olr.allsum;
z22  = climcapsL3_spectral_olr.allsum;
  plotoptions5.str11 = 'ERA5'; plotoptions5.str12 = 'MERRA2';
  plotoptions5.strzz = 'UMBC'; 
  plotoptions5.str21 = 'AIRS L3'; plotoptions5.str22 = 'CLIMCAPS L3';
  plotoptions5.ystr = 'Latitude'; plotoptions5.xstr = 'Longitude';
  plotoptions5.maintitle = '\lambda W/m2/K';
  plotoptions5.cx = [-1 +1]*10*dascale; plotoptions5.cmap = llsmap5; plotoptions5.yReverseDir = -1; plotoptions5.yLinearOrLog = +1;
  aslmap_2x1x2tiledlayout(z11/boox(1),z12/boox(2),zMID/boox(3),z21/boox(4),z22/boox(5),7,plotoptions5);

z11  = -era5_spectral_olr.wv_ecRad.clr       + era5_spectral_olr.olr0_ecRad.clr;
z12  = -merra2_spectral_olr.wv_ecRad.clr     + merra2_spectral_olr.olr0_ecRad.clr;
zMID = -umbc_spectral_olr.wv_ecRad.clr       + umbc_spectral_olr.olr0_ecRad.clr;
z21  = -airsL3_spectral_olr.wv_ecRad.clr     + airsL3_spectral_olr.olr0_ecRad.clr;
z22  = -climcapsL3_spectral_olr.wv_ecRad.clr + climcapsL3_spectral_olr.olr0_ecRad.clr;
  plotoptions5.str11 = 'ERA5'; plotoptions5.str12 = 'MERRA2';
  plotoptions5.strzz = 'UMBC'; 
  plotoptions5.str21 = 'AIRS L3'; plotoptions5.str22 = 'CLIMCAPS L3';
  plotoptions5.ystr = 'Latitude'; plotoptions5.xstr = 'Longitude';
  plotoptions5.maintitle = '\lambda_{wv} W/m2/K';
  plotoptions5.cx = [-1 +1]*10*dascale; plotoptions5.cmap = llsmap5; plotoptions5.yReverseDir = -1; plotoptions5.yLinearOrLog = +1;
  aslmap_2x1x2tiledlayout(z11/boox(1),z12/boox(2),zMID/boox(3),z21/boox(4),z22/boox(5),8,plotoptions5);

z11  = -era5_spectral_olr.planck_ecRad.clr       + era5_spectral_olr.olr0_ecRad.clr;
z12  = -merra2_spectral_olr.planck_ecRad.clr     + merra2_spectral_olr.olr0_ecRad.clr;
zMID = -umbc_spectral_olr.planck_ecRad.clr       + umbc_spectral_olr.olr0_ecRad.clr;
z21  = -airsL3_spectral_olr.planck_ecRad.clr     + airsL3_spectral_olr.olr0_ecRad.clr;
z22  = -climcapsL3_spectral_olr.planck_ecRad.clr + climcapsL3_spectral_olr.olr0_ecRad.clr;
  plotoptions5.str11 = 'ERA5'; plotoptions5.str12 = 'MERRA2';
  plotoptions5.strzz = 'UMBC'; 
  plotoptions5.str21 = 'AIRS L3'; plotoptions5.str22 = 'CLIMCAPS L3';
  plotoptions5.ystr = 'Latitude'; plotoptions5.xstr = 'Longitude';
  plotoptions5.maintitle = '\lambda_{planck} W/m2/K';
  plotoptions5.cx = [-1 +1]*10*dascale; plotoptions5.cmap = llsmap5; plotoptions5.yReverseDir = -1; plotoptions5.yLinearOrLog = +1;
  aslmap_2x1x2tiledlayout(z11/boox(1),z12/boox(2),zMID/boox(3),z21/boox(4),z22/boox(5),9,plotoptions5);

z11  = -era5_spectral_olr.planck_ecRad.clr       + era5_spectral_olr.lapse_ecRad.clr;
z12  = -merra2_spectral_olr.planck_ecRad.clr     + merra2_spectral_olr.lapse_ecRad.clr;
zMID = -umbc_spectral_olr.planck_ecRad.clr       + umbc_spectral_olr.lapse_ecRad.clr;
z21  = -airsL3_spectral_olr.planck_ecRad.clr     + airsL3_spectral_olr.lapse_ecRad.clr;
z22  = -climcapsL3_spectral_olr.planck_ecRad.clr + climcapsL3_spectral_olr.lapse_ecRad.clr;
  plotoptions5.str11 = 'ERA5'; plotoptions5.str12 = 'MERRA2';
  plotoptions5.strzz = 'UMBC'; 
  plotoptions5.str21 = 'AIRS L3'; plotoptions5.str22 = 'CLIMCAPS L3';
  plotoptions5.ystr = 'Latitude'; plotoptions5.xstr = 'Longitude';
  plotoptions5.maintitle = '\lambda_{lapse} W/m2/K';
  plotoptions5.cx = [-1 +1]*10*dascale; plotoptions5.cmap = llsmap5; plotoptions5.yReverseDir = -1; plotoptions5.yLinearOrLog = +1;
  aslmap_2x1x2tiledlayout(z11/boox(1),z12/boox(2),zMID/boox(3),z21/boox(4),z22/boox(5),10,plotoptions5);

%junk = [sum(era5_spectral_olr.allsum .* coslat)/sum(coslat) sum(merra2_spectral_olr.allsum .* coslat)/sum(coslat) sum(umbc_spectral_olr.allsum .* coslat)/sum(coslat) sum(airsL3_spectral_olr.allsum .* coslat)/sum(coslat)] ./ boox(1:4);
%fprintf(1,'feedbacks for ERA5/MERRA2/AIRS L3/CLIMCAPS L3 = %8.5f %8.5f %8.5f %8.5f W/m2/K \n',junk);
junk = [sum(era5_spectral_olr.allsum .* coslat)/sum(coslat) sum(merra2_spectral_olr.allsum .* coslat)/sum(coslat) sum(umbc_spectral_olr.allsum .* coslat)/sum(coslat) ...
        sum(airsL3_spectral_olr.allsum .* coslat)/sum(coslat) sum(climcapsL3_spectral_olr.allsum .* coslat)/sum(coslat)] ./ boo(1:5);
fprintf(1,'feedbacks for ERA5/MERRA2/UMBC/AIRSL3/CLIMCAPSL3 = %8.5f %8.5f %8.5f %8.5f %8.5f W/m2/K \n',junk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cosYY = cos(YY*pi/180);
tropics = find(abs(YY) < 30);
poles   = find(abs(YY) >= 60);
midlats = find(abs(YY) >= 30 & abs(YY) < 60);

disp(' ')
disp('---------------|---------------------|---------------------|---------------------|---------------------')
disp('               |       GLOBAL        |         TROPICS             MIDLATS       |       POLAR         ')
disp(' <SKT>  K/yr   |    raw      coswgt  |    raw      coswgt  |    raw      coswgt  |    raw      coswgt  ')
disp('---------------|---------------------|---------------------|---------------------|---------------------')
for ix = 1 : 7
  str = ['boo = ' thename{ix} '_spectral_olr.stemptrend;'];
  eval(str);
  bah0 = sum(cosYY.*boo)/sum(cosYY);
  bahT = sum(cosYY(tropics).*boo(tropics))/sum(cosYY(tropics));
  bahM = sum(cosYY(midlats).*boo(midlats))/sum(cosYY(midlats));
  bahP = sum(cosYY(poles).*boo(poles))/sum(cosYY(poles));
  junkSKT(ix,:) = [mean(boo) bah0 mean(boo(tropics)) bahT mean(boo(midlats)) bahM mean(boo(poles)) bahP];
  fprintf(1,'  %12s | %8.3f   %8.3f | %8.3f   %8.3f | %8.3f   %8.3f | %8.3f   %8.3f \n',strfeedbacks{ix},junkSKT(ix,:))
end
disp('---------------|---------------------|---------------------|---------------------|---------------------')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
ixx = ix;
showfeedbacks(1:ixx,7) = sum(showfeedbacks(1:ixx,[1 2 3 4]),2);

disp('showing means from fits, 64 latbins')
%disp('         Planck Lapse Ozone Water |  Total')
%for ix = 1 : 5
%  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},showfeedbacks(ix,[1 2 3 4 7]));
%end
%trends_paper_show = showfeedbacks(1:5,[1 2 3 4 7]);

showfeedbacks(1:5,8) = junk; 
disp('         Planck Lapse Ozone Water |  Total')
for ix = 1 : 5
  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f |  %5.2f \n',strfeedbacks{ix},showfeedbacks(ix,[1 2 3 4 8]));
end
trends_paper_show = showfeedbacks(1:5,[1 2 3 4 8]);

figure(1); clf
bar(trends_paper_show')
ylabel('Feedback W/m2/K');
hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south');
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)

%%%%%%%%%%%%%%%%%%%%%%%%%
plot(era5_spectral_olr.stemptrend,era5_spectral_olr.olr0_ecRad.clr-era5_spectral_olr.planck_ecRad.clr,'.',...
     merra2_spectral_olr.stemptrend,merra2_spectral_olr.olr0_ecRad.clr-merra2_spectral_olr.planck_ecRad.clr,'.',...
     umbc_spectral_olr.stemptrend,umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.planck_ecRad.clr,'.',...
     airsL3_spectral_olr.stemptrend,airsL3_spectral_olr.olr0_ecRad.clr-airsL3_spectral_olr.planck_ecRad.clr,'.',...
     climcapsL3_spectral_olr.stemptrend,climcapsL3_spectral_olr.olr0_ecRad.clr-climcapsL3_spectral_olr.planck_ecRad.clr,'.')
plotaxis2; axis([-0.2 +0.2 -0.4 +0.4]); title('Planck')
sum_all4608(1,1) = sum((era5_spectral_olr.olr0_ecRad.clr-era5_spectral_olr.planck_ecRad.clr).*cos(YY*pi/180));
sum_all4608(1,2) = sum((merra2_spectral_olr.olr0_ecRad.clr-merra2_spectral_olr.planck_ecRad.clr).*cos(YY*pi/180));
sum_all4608(1,3) = sum((umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.planck_ecRad.clr).*cos(YY*pi/180));
sum_all4608(1,4) = sum((airsL3_spectral_olr.olr0_ecRad.clr-airsL3_spectral_olr.planck_ecRad.clr).*cos(YY*pi/180));
sum_all4608(1,5) = sum((climcapsL3_spectral_olr.olr0_ecRad.clr-climcapsL3_spectral_olr.planck_ecRad.clr).*cos(YY*pi/180));

plot(era5_spectral_olr.stemptrend,era5_spectral_olr.lapse_ecRad.clr-era5_spectral_olr.planck_ecRad.clr,'.',...
     merra2_spectral_olr.stemptrend,merra2_spectral_olr.lapse_ecRad.clr-merra2_spectral_olr.planck_ecRad.clr,'.',...
     umbc_spectral_olr.stemptrend,umbc_spectral_olr.lapse_ecRad.clr-umbc_spectral_olr.planck_ecRad.clr,'.',...
     airsL3_spectral_olr.stemptrend,airsL3_spectral_olr.lapse_ecRad.clr-airsL3_spectral_olr.planck_ecRad.clr,'.',...
     climcapsL3_spectral_olr.stemptrend,climcapsL3_spectral_olr.lapse_ecRad.clr-climcapsL3_spectral_olr.planck_ecRad.clr,'.')
plotaxis2; axis([-0.2 +0.2 -0.4 +0.4]); title('Lapse')
sum_all4608(2,1) = -sum((era5_spectral_olr.lapse_ecRad.clr-era5_spectral_olr.planck_ecRad.clr).*cos(YY*pi/180));
sum_all4608(2,2) = -sum((merra2_spectral_olr.lapse_ecRad.clr-merra2_spectral_olr.planck_ecRad.clr).*cos(YY*pi/180));
sum_all4608(2,3) = -sum((umbc_spectral_olr.lapse_ecRad.clr-umbc_spectral_olr.planck_ecRad.clr).*cos(YY*pi/180));
sum_all4608(2,4) = -sum((airsL3_spectral_olr.lapse_ecRad.clr-airsL3_spectral_olr.planck_ecRad.clr).*cos(YY*pi/180));
sum_all4608(2,5) = -sum((climcapsL3_spectral_olr.lapse_ecRad.clr-climcapsL3_spectral_olr.planck_ecRad.clr).*cos(YY*pi/180));

plot(era5_spectral_olr.stemptrend,era5_spectral_olr.olr0_ecRad.clr-era5_spectral_olr.o3_ecRad.clr,'.',...
     merra2_spectral_olr.stemptrend,merra2_spectral_olr.olr0_ecRad.clr-merra2_spectral_olr.o3_ecRad.clr,'.',...
     umbc_spectral_olr.stemptrend,umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.o3_ecRad.clr,'.',...
     airsL3_spectral_olr.stemptrend,airsL3_spectral_olr.olr0_ecRad.clr-airsL3_spectral_olr.o3_ecRad.clr,'.',...
     climcapsL3_spectral_olr.stemptrend,climcapsL3_spectral_olr.olr0_ecRad.clr-climcapsL3_spectral_olr.o3_ecRad.clr,'.')
plotaxis2; axis([-0.2 +0.2 -0.4 +0.4]); title('O3')
sum_all4608(3,1) = sum((era5_spectral_olr.olr0_ecRad.clr-era5_spectral_olr.o3_ecRad.clr).*cos(YY*pi/180));
sum_all4608(3,2) = sum((merra2_spectral_olr.olr0_ecRad.clr-merra2_spectral_olr.o3_ecRad.clr).*cos(YY*pi/180));
sum_all4608(3,3) = sum((umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.o3_ecRad.clr).*cos(YY*pi/180));
sum_all4608(3,4) = sum((airsL3_spectral_olr.olr0_ecRad.clr-airsL3_spectral_olr.o3_ecRad.clr).*cos(YY*pi/180));
sum_all4608(3,5) = sum((climcapsL3_spectral_olr.olr0_ecRad.clr-climcapsL3_spectral_olr.o3_ecRad.clr).*cos(YY*pi/180));

plot(era5_spectral_olr.stemptrend,era5_spectral_olr.olr0_ecRad.clr-era5_spectral_olr.wv_ecRad.clr,'.',...
     merra2_spectral_olr.stemptrend,merra2_spectral_olr.olr0_ecRad.clr-merra2_spectral_olr.wv_ecRad.clr,'.',...
     umbc_spectral_olr.stemptrend,umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.wv_ecRad.clr,'.',...
     airsL3_spectral_olr.stemptrend,airsL3_spectral_olr.olr0_ecRad.clr-airsL3_spectral_olr.wv_ecRad.clr,'.',...
     climcapsL3_spectral_olr.stemptrend,climcapsL3_spectral_olr.olr0_ecRad.clr-climcapsL3_spectral_olr.wv_ecRad.clr,'.')
plotaxis2; axis([-0.2 +0.2 -0.4 +0.4]); title('Wv')
sum_all4608(4,1) = sum((era5_spectral_olr.olr0_ecRad.clr-era5_spectral_olr.wv_ecRad.clr).*cos(YY*pi/180));
sum_all4608(4,2) = sum((merra2_spectral_olr.olr0_ecRad.clr-merra2_spectral_olr.wv_ecRad.clr).*cos(YY*pi/180));
sum_all4608(4,3) = sum((umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.wv_ecRad.clr).*cos(YY*pi/180));
sum_all4608(4,4) = sum((airsL3_spectral_olr.olr0_ecRad.clr-airsL3_spectral_olr.wv_ecRad.clr).*cos(YY*pi/180));
sum_all4608(4,5) = sum((climcapsL3_spectral_olr.olr0_ecRad.clr-climcapsL3_spectral_olr.wv_ecRad.clr).*cos(YY*pi/180));

junk = (era5_spectral_olr.olr0_ecRad.clr-era5_spectral_olr.planck_ecRad.clr) - (era5_spectral_olr.lapse_ecRad.clr-era5_spectral_olr.planck_ecRad.clr) + ...
       (era5_spectral_olr.olr0_ecRad.clr-era5_spectral_olr.o3_ecRad.clr) + (era5_spectral_olr.olr0_ecRad.clr-era5_spectral_olr.wv_ecRad.clr);
  sum_all4608(5,1) = sum(junk.*cos(YY*pi/180));
junk = (merra2_spectral_olr.olr0_ecRad.clr-merra2_spectral_olr.planck_ecRad.clr) - (merra2_spectral_olr.lapse_ecRad.clr-merra2_spectral_olr.planck_ecRad.clr) + ...
       (merra2_spectral_olr.olr0_ecRad.clr-merra2_spectral_olr.o3_ecRad.clr) + (merra2_spectral_olr.olr0_ecRad.clr-merra2_spectral_olr.wv_ecRad.clr);
  sum_all4608(5,2) = sum(junk.*cos(YY*pi/180));
junk = (umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.planck_ecRad.clr) - (umbc_spectral_olr.lapse_ecRad.clr-umbc_spectral_olr.planck_ecRad.clr) + ...
       (umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.o3_ecRad.clr) + (umbc_spectral_olr.olr0_ecRad.clr-umbc_spectral_olr.wv_ecRad.clr);
  sum_all4608(5,3) = sum(junk.*cos(YY*pi/180));
junk = (airsL3_spectral_olr.olr0_ecRad.clr-airsL3_spectral_olr.planck_ecRad.clr) - (airsL3_spectral_olr.lapse_ecRad.clr-airsL3_spectral_olr.planck_ecRad.clr) + ...
       (airsL3_spectral_olr.olr0_ecRad.clr-airsL3_spectral_olr.o3_ecRad.clr) + (airsL3_spectral_olr.olr0_ecRad.clr-airsL3_spectral_olr.wv_ecRad.clr);
  sum_all4608(5,4) = sum(junk.*cos(YY*pi/180));
junk = (climcapsL3_spectral_olr.olr0_ecRad.clr-climcapsL3_spectral_olr.planck_ecRad.clr) - (climcapsL3_spectral_olr.lapse_ecRad.clr-climcapsL3_spectral_olr.planck_ecRad.clr) + ...
       (climcapsL3_spectral_olr.olr0_ecRad.clr-climcapsL3_spectral_olr.o3_ecRad.clr) + (climcapsL3_spectral_olr.olr0_ecRad.clr-climcapsL3_spectral_olr.wv_ecRad.clr);
  sum_all4608(5,5) = sum(junk.*cos(YY*pi/180));

disp(' ')
disp('showing means from weighted rlat fits, all 4608 elements')
sum_all4608 = sum_all4608';
sum_all4608_deltaOLR = sum_all4608;
sum_all4608 = sum_all4608./(junkSKT(1:5,2)*ones(1,5))/sum(cos(YY*pi/180));
for ix = 1 : 5
  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},sum_all4608(ix,:));
end

disp('showing {cos{SKT)} {delta(OLR)=numerator}, {means}, {stddev}, {normalized deviations} from weighted rlat fits, all 4608 elements')
fprintf(1,'if sum(delta(OLR) looks large, remember we divide by sum(cos(YY*pi/180) which is %8.3f \n',sum(cos(YY*pi/180)))
wawoo  = sum_all4608_deltaOLR;
mwawoo = ones(5,1) * mean(sum_all4608_deltaOLR,1);
swawoo = ones(5,1) * std(sum_all4608_deltaOLR,1);
nwawoo = (sum_all4608_deltaOLR - mwawoo)./(swawoo);
for ix = 1 : 5
  kaboo = [junkSKT(ix,2) wawoo(ix,:) mwawoo(ix,:) swawoo(ix,:) nwawoo(ix,:)];
  fprintf(1,'%s || %7.4f ||  %7.2f %7.2f %7.2f %7.2f    %7.2f | %7.2f %7.2f %7.2f %7.2f    %7.2f | %7.2f %7.2f %7.2f %7.2f    %7.2f | %7.2f %7.2f %7.2f %7.2f    %7.2f \n',strfeedbacks{ix},kaboo)
end

%%%%%%%%%%%%%%%%%%%%%%%%%
return
clear cosavg
cosrlat = cos(rlat'*pi/180);

cosxrlat = cosrlat;
cosxrlat = ones(size(cosrlat));

disp(' ')
disp('THIS IS WIERD since we have already done the cosine weighting, and doing more cosine weighting')
disp('BUT IF YOU MAKE cosrlat == 1 then DOING the MEAN and ANSWERS come out GOOD')
disp('<<<< showing means from weighted rlat fits, 64 latbins >>> ')
disp('BUT IF YOU MAKE cosrlat == 1 then DOING the MEAN and ANSWERS come out GOOD')
disp('THIS IS WIERD since we have already done the cosine weighting, and doing more cosine weighting')

cosavg(1,1) = sum(era5_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(1,2) = sum(merra2_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(1,3) = sum(umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(1,4) = sum(airsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(1,5) = sum(climcapsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);

cosavg(2,1) = sum(era5_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(2,2) = sum(merra2_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(2,3) = sum(umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(2,4) = sum(airsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(2,5) = sum(climcapsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);

cosavg(3,1) = sum(era5_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(3,2) = sum(merra2_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(3,3) = sum(umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(3,4) = sum(airsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(3,5) = sum(climcapsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);

cosavg(4,1) = sum(era5_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(4,2) = sum(merra2_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(4,3) = sum(umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(4,4) = sum(airsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);
cosavg(4,5) = sum(climcapsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin .* cosxrlat) / sum(cosxrlat);

cosavg(5,:) = sum(cosavg(1:4,:));
cosavg(5,1) = sum((era5_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + era5_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
                   era5_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + era5_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin).* cosxrlat) / sum(cosxrlat);
cosavg(5,2) = sum((merra2_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + merra2_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
                   merra2_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + merra2_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin).* cosxrlat) / sum(cosxrlat);
cosavg(5,3) = sum((umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
                   umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin).* cosxrlat) / sum(cosxrlat);
cosavg(5,4) = sum((airsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + airsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
                   airsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + airsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin).* cosxrlat) / sum(cosxrlat);
cosavg(5,5) = sum((climcapsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + climcapsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
                   climcapsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + climcapsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin).* cosxrlat) / sum(cosxrlat);
cosavg = cosavg';

for ix = 1 : 5
  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},cosavg(ix,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xsin = sin(rlat*pi/180);
  xtick = [-1 -sqrt(3)/2 -sqrt(2)/2 -1/2 -(0.25+0.01) 0 +(0.25+0.01) +1/2 +sqrt(2)/2 +sqrt(3)/2 +1]; %% -90 -60 -45 -30 -15 0 +15 +30 +45 +60 +90
  xtick = [-1            -sqrt(2)/2      -(0.25+0.01) 0 +(0.25+0.01)      +sqrt(2)/2            +1]; %% -90     -45     -15 0 +15     +45     +90
  xticklab = cellstr(num2str(round(180/pi*asin((xtick(:)))), '%d'));

figure(2); clf;
subplot(221); plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
                   xsin,smooth(umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
                   xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth))
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('Planck'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(222); plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
                   xsin,smooth(umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
                   xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth))
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('Lapse'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(223); plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
                   xsin,smooth(umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
                   xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth))
  plotaxis2; hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('Ozone'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(224); plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
                   xsin,smooth(umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
                   xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth))
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('WV'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf;
ta = tiledlayout(2,2,'TileSpacing','compact', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),'linewidth',2);
  plotaxis2; box on; 
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; 

tafov(2) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');
  hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','north','fontsize',8);

tafov(3) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(4) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
      xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),'linewidth',2);
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

%% see  ~sergio/MATLABCODE/PLOTTER/TILEDPLOTS/profile_plots_2x1x2tiledlayout_wide.m
figure(5); clf;
ta = tiledlayout(3,4);
ta = tiledlayout(3,4,'TileSpacing','None', 'Padding','None');
ta = tiledlayout(3,4,'TileSpacing','compact', 'Padding','compact');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile([1,2]);
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin,iSmooth),'linewidth',2);
  plotaxis2; box on; 
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; 

tafov(2) = nexttile([1,2]);
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin,iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

nexttile([1,1]);
axis off
tafov(3) = nexttile([1,2]);
boo1 = era5_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + era5_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
       era5_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + era5_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin;
boo2 = merra2_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + merra2_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
       merra2_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + merra2_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin;
boo3 = umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
       umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin;
boo4 = airsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + airsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
       airsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + airsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin;
boo5 = climcapsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + climcapsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
       climcapsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + climcapsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin;
plot(xsin,smooth(boo1,iSmooth),xsin,smooth(boo2,iSmooth),xsin,smooth(boo3,iSmooth),xsin,smooth(boo4,iSmooth),xsin,smooth(boo5,iSmooth),'linewidth',2);
  plotaxis2; box on;
     hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','eastoutside','fontsize',8);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');
nexttile([1,1]);
axis off

tafov(4) = nexttile([1,2]);
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin,iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(5) = nexttile([1,2]);
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),...
      xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin,iSmooth),'linewidth',2);
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

tafov(1).YLabel.String = ystr; tafov(1).YLabel.FontSize = 10;
tafov(3).YLabel.String = ystr; tafov(3).YLabel.FontSize = 10;
tafov(4).YLabel.String = ystr; tafov(4).YLabel.FontSize = 18;
tafov(4).XLabel.String = xstr; tafov(3).XLabel.FontSize = 10;
tafov(5).XLabel.String = xstr; tafov(4).XLabel.FontSize = 10;

plotoptions.str11 = 'Planck';
plotoptions.str12 = 'Lapse';
plotoptions.strXY = 'LongWave';
plotoptions.str21 = 'Ozone';
plotoptions.str22 = 'Water';

yposn = 0.8; %% inside, just below top line
yposn = 0.9; %% straddling top line
yposn = 1.0; %% outside, just above top line

title(tafov(1), plotoptions.str11, 'Units', 'normalized', 'Position', [0.5, yposn, 0]);
title(tafov(2), plotoptions.str12, 'Units', 'normalized', 'Position', [0.5, yposn, 0]);
title(tafov(3), plotoptions.strXY, 'Units', 'normalized', 'Position', [0.5, yposn, 0]);
title(tafov(4), plotoptions.str21, 'Units', 'normalized', 'Position', [0.5, yposn, 0]);
title(tafov(5), plotoptions.str22, 'Units', 'normalized', 'Position', [0.5, yposn, 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xsin2 = sin(rlat*pi/180);
  xtick2 = [-1            -sqrt(2)/2      -(0.25+0.01) 0 +(0.25+0.01)      +sqrt(2)/2            +1]; %% -90     -45     -15 0 +15     +45     +90
  xtick2 = [-1 -sqrt(3)/2 -sqrt(2)/2 -1/2 -(0.25+0.01) 0 +(0.25+0.01) +1/2 +sqrt(2)/2 +sqrt(3)/2 +1]; %% -90 -60 -45 -30 -15 0 +15 +30 +45 +60 +90
  xtick2lab = cellstr(num2str(round(180/pi*asin((xtick2(:)))), '%d'));

figure(6); clf;
boo1 = era5_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + era5_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
       era5_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + era5_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin;
boo2 = merra2_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + merra2_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
       merra2_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + merra2_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin;
boo3 = umbc_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + umbc_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
       umbc_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + umbc_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin;
boo4 = airsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + airsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
       airsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + airsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin;
boo5 = climcapsL3_spectral_olr.feedback_ecRad.planck.globalSST_weighted_latbin + climcapsL3_spectral_olr.feedback_ecRad.lapse.globalSST_weighted_latbin + ...
       climcapsL3_spectral_olr.feedback_ecRad.o3.globalSST_weighted_latbin + climcapsL3_spectral_olr.feedback_ecRad.wv.globalSST_weighted_latbin;
plot(xsin2,smooth(boo1,iSmooth),xsin2,smooth(boo2,iSmooth),xsin2,smooth(boo3,iSmooth),xsin2,smooth(boo4,iSmooth),xsin2,smooth(boo5,iSmooth),'linewidth',2);
  xlim([-1 +1]); ylim([-1 +1]*7.5)
  plotaxis2; 
     hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','best','fontsize',8);
  set(gca,'Xtick',xtick2,'XTickLabel',xtick2lab,'TickLabelInterpreter','tex'); 
  title('Longwave \lambda feedback')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
addpath /asl/matlib/plotutils

dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
if dascale == +1
  figure(7);  aslprint([dir0 'global_feedbackparams_20yrs_ERA5_MERRA2_UMBC_AIRSL3_CLIMCAPSL3_total_correct_deltaSKT.pdf'])
  figure(8);  aslprint([dir0 'global_feedbackparams_20yrs_ERA5_MERRA2_UMBC_AIRSL3_CLIMCAPSL3_water_correct_deltaSKT.pdf'])
  figure(9);  aslprint([dir0 'global_feedbackparams_20yrs_ERA5_MERRA2_UMBC_AIRSL3_CLIMCAPSL3_planck_correct_deltaSKT.pdf'])
  figure(10); aslprint([dir0 'global_feedbackparams_20yrs_ERA5_MERRA2_UMBC_AIRSL3_CLIMCAPSL3_lapse_correct_deltaSKT.pdf'])
else
  figure(7);  aslprint([dir0 'global_feedbackparams_20yrs_ERA5_MERRA2_UMBC_AIRSL3_CLIMCAPSL3_total_unity_deltaSKT.pdf'])
  figure(8);  aslprint([dir0 'global_feedbackparams_20yrs_ERA5_MERRA2_UMBC_AIRSL3_CLIMCAPSL3_water_unity_deltaSKT.pdf'])
  figure(9);  aslprint([dir0 'global_feedbackparams_20yrs_ERA5_MERRA2_UMBC_AIRSL3_CLIMCAPSL3_planck_unity_deltaSKT.pdf'])
  figure(10); aslprint([dir0 'global_feedbackparams_20yrs_ERA5_MERRA2_UMBC_AIRSL3_CLIMCAPSL3_lapse_unity_deltaSKT.pdf'])
end
%}

