clear showfeedbacks* strfeedbacks

if exist('umbc_spectral_olr')
  disp('umbc_spectral_olr already exists, you may want to do   "clear *spectral_olr"    before running this')
end

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
%% the value
showfeedbacks_robustfit_all(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('merra2_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'merra2_spectral_olr');
  merra2_spectral_olr = junk.merra2_spectral_olr;
end
ix = 2; junk = merra2_spectral_olr;
strfeedbacks{ix} = 'MERRA2     ';
%% the value
showfeedbacks_robustfit_all(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('umbc_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'],'umbc_spectral_olr');
  umbc_spectral_olr = junk.umbc_spectral_olr;
end
ix = 3; junk = umbc_spectral_olr;
strfeedbacks{ix} = 'THIS WORK  ';
%% the value
showfeedbacks_robustfit_all(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('airsL3_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'airsL3_spectral_olr');
  airsL3_spectral_olr = junk.airsL3_spectral_olr;
end
ix = 4; junk = airsL3_spectral_olr;
strfeedbacks{ix} = 'AIRS L3    ';
%% the value
showfeedbacks_robustfit_all(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('climcapsL3_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'climcapsL3_spectral_olr');
  climcapsL3_spectral_olr = junk.climcapsL3_spectral_olr;
end
ix = 5; junk = climcapsL3_spectral_olr;
strfeedbacks{ix} = 'CLIMCAPS L3';
%% the value
showfeedbacks_robustfit_all(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('cmip6_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'cmip6_spectral_olr');
  cmip6_spectral_olr = junk.cmip6_spectral_olr;
end
ix = 6; junk = cmip6_spectral_olr;
strfeedbacks{ix} = 'CMIP6      ';
%% the value
showfeedbacks_robustfit_all(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('amip6_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'amip6_spectral_olr');
  amip6_spectral_olr = junk.amip6_spectral_olr;
end
ix = 7; junk = amip6_spectral_olr;
strfeedbacks{ix} = 'AMIP6      ';
%% the value
showfeedbacks_robustfit_all(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('showing means from fits, 64 latbins')

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
ixx = ix;
showfeedbacks_robustfit_all(1:ixx,7,1) = sum(squeeze(showfeedbacks_robustfit_all(1:ixx,[1 2 3 4],1)),2);
junk = showfeedbacks_robustfit_all(1:ixx,[1 2 3 4],2);
junk = sqrt(sum(junk.*junk,2));
showfeedbacks_robustfit_all(1:ixx,7,2) = junk;

ixshow = 5;
for ix = 1 : ixshow
  %fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},showfeedbacks_robustfit_all(ix,[1 2 3 4 7],1));
  junk = [showfeedbacks_robustfit_all(ix,[1 2 3 4 7],1) showfeedbacks_robustfit_all(ix,[1 2 3 4 7],2)];
  junk = junk([1 6 2 7 3 8 4 9 5 10]);
  fprintf(1,'%s %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f    %6.3f +/- %6.3f \n',strfeedbacks{ix},junk);
end
trends_paper_show = showfeedbacks_robustfit_all(1:ixshow,[1 2 3 4 7]);

figure(1); clf
bar(trends_paper_show')
ylabel('Feedback W/m2/K');
hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south');
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)
title('Global feedbacks')

%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rlat')
  load latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
  [Y,X] = meshgrid(rlat,rlon);
  X = X; Y = Y;
end

if ~exist('iSmooth')
  iSmooth = 5;
  iSmooth = 10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
clear cosavg
cosrlat = cos(rlat*pi/180);

cosxrlat = cosrlat;
cosxrlat = ones(size(cosrlat));
%cosxrlat = cosrlat.^0.5;

disp(' ')
disp('THIS IS WIERD since we have already done the cosine weighting, and doing more cosine weighting')
disp('BUT IF YOU MAKE cosrlat == 1 then DOING the MEAN and ANSWERS still are BAD, compared to plot_show_olr_ecRad_feedback_globalsstfit_smooth2')
disp('<<<< showing means from weighted rlat fits, 64 latbins >>> ')
disp('BUT IF YOU MAKE cosrlat == 1 then DOING the MEAN and ANSWERS come out GOOD')
disp('THIS IS WIERD since we have already done the cosine weighting, and doing more cosine weighting')

cosavg(1,1) = sum(era5_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(1,2) = sum(merra2_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(1,3) = sum(umbc_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(1,4) = sum(airsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(1,5) = sum(climcapsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);

cosavg(2,1) = sum(era5_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(2,2) = sum(merra2_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(2,3) = sum(umbc_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(2,4) = sum(airsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(2,5) = sum(climcapsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);

cosavg(3,1) = sum(era5_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(3,2) = sum(merra2_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(3,3) = sum(umbc_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(3,4) = sum(airsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(3,5) = sum(climcapsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);

cosavg(4,1) = sum(era5_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(4,2) = sum(merra2_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(4,3) = sum(umbc_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(4,4) = sum(airsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);
cosavg(4,5) = sum(climcapsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1) .* cosxrlat) / sum(cosxrlat);

cosavg(5,:) = sum(cosavg(1:4,:));
cosavg(5,1) = sum((era5_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + era5_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
                   era5_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) + era5_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1)).* cosxrlat) / sum(cosxrlat);
cosavg(5,2) = sum((merra2_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + merra2_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
                   merra2_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) + merra2_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1)).* cosxrlat) / sum(cosxrlat);
cosavg(5,3) = sum((umbc_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + umbc_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
                   umbc_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) + umbc_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1)).* cosxrlat) / sum(cosxrlat);
cosavg(5,4) = sum((airsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + airsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
                   airsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) + airsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1)).* cosxrlat) / sum(cosxrlat);
cosavg(5,5) = sum((climcapsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + climcapsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
                   climcapsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) + climcapsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1)).* cosxrlat) / sum(cosxrlat);
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
subplot(221); 
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('Planck'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(222); 
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('Lapse'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(223); 
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('Ozone'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(224); 
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
      xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('WV'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf;
subplot(221); 
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,2),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,2),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,2),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,2),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,2),iSmooth),'linewidth',2);
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('unc Planck'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(222); 
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,2),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,2),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,2),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,2),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,2),iSmooth),'linewidth',2);
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('unc Lapse'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(223); 
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,2),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,2),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,2),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,2),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,2),iSmooth),'linewidth',2);
  plotaxis2; hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('unc Ozone'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

subplot(224); 
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,2),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,2),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,2),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,2),iSmooth),...
      xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,2),iSmooth),'linewidth',2);
  plotaxis2; %hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',8);
  title('unc WV'); xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf;
ta = tiledlayout(2,2,'TileSpacing','compact', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on; hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','north','fontsize',8);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; 

tafov(2) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(3) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(4) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
      xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),'linewidth',2);
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
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on; 
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; 

tafov(2) = nexttile([1,2]);
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

nexttile([1,1]);
axis off
tafov(3) = nexttile([1,2]);
boo1 = era5_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + era5_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
       era5_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) + era5_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1);
boo2 = merra2_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + merra2_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
       merra2_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) + merra2_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1);
boo3 = umbc_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + umbc_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
       umbc_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) + umbc_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1);
boo4 = airsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + airsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
       airsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) + airsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1);
boo5 = climcapsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1) + climcapsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
       climcapsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1) + climcapsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1);
plot(xsin,smooth(boo1,iSmooth),xsin,smooth(boo2,iSmooth),xsin,smooth(boo3,iSmooth),xsin,smooth(boo4,iSmooth),xsin,smooth(boo5,iSmooth),'linewidth',2);
  plotaxis2; box on;
     hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','eastoutside','fontsize',8);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');
nexttile([1,1]);
axis off

tafov(4) = nexttile([1,2]);
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(5) = nexttile([1,2]);
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),'linewidth',2);
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

%{
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
figure(4); aslprint([dir0 'feedbackparams_20yrs_UMBC_ERA5_MERRA2_L3.pdf'])
%}

