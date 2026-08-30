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
showfeedbacks_robustfit_midlats(ix,1,1) = junk.feedback_ecRad.planck.robustfit_midlats(1);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_midlats(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,3,1) = junk.feedback_ecRad.o3.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,4,1) = junk.feedback_ecRad.wv.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,5,1) = junk.feedback_ecRad.skt.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(1);
%% the unc
showfeedbacks_robustfit_midlats(ix,1,2) = junk.feedback_ecRad.planck.robustfit_midlats(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_midlats(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,3,2) = junk.feedback_ecRad.o3.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,4,2) = junk.feedback_ecRad.wv.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,5,2) = junk.feedback_ecRad.skt.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(2);

if ~exist('merra2_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'merra2_spectral_olr');
  merra2_spectral_olr = junk.merra2_spectral_olr;
end
ix = 2; junk = merra2_spectral_olr;
strfeedbacks{ix} = 'MERRA2     ';
%% the value
showfeedbacks_robustfit_midlats(ix,1,1) = junk.feedback_ecRad.planck.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,3,1) = junk.feedback_ecRad.o3.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,4,1) = junk.feedback_ecRad.wv.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,5,1) = junk.feedback_ecRad.skt.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(1);
%% the unc
showfeedbacks_robustfit_midlats(ix,1,2) = junk.feedback_ecRad.planck.robustfit_midlats(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_midlats(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,3,2) = junk.feedback_ecRad.o3.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,4,2) = junk.feedback_ecRad.wv.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,5,2) = junk.feedback_ecRad.skt.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(2);

if ~exist('umbc_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'],'umbc_spectral_olr');
  umbc_spectral_olr = junk.umbc_spectral_olr;
end
ix = 3; junk = umbc_spectral_olr;
strfeedbacks{ix} = 'THIS WORK  ';
%% the value
showfeedbacks_robustfit_midlats(ix,1,1) = junk.feedback_ecRad.planck.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,3,1) = junk.feedback_ecRad.o3.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,4,1) = junk.feedback_ecRad.wv.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,5,1) = junk.feedback_ecRad.skt.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(1);
%% the unc
showfeedbacks_robustfit_midlats(ix,1,2) = junk.feedback_ecRad.planck.robustfit_midlats(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_midlats(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,3,2) = junk.feedback_ecRad.o3.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,4,2) = junk.feedback_ecRad.wv.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,5,2) = junk.feedback_ecRad.skt.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(2);

if ~exist('airsL3_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'airsL3_spectral_olr');
  airsL3_spectral_olr = junk.airsL3_spectral_olr;
end
ix = 4; junk = airsL3_spectral_olr;
strfeedbacks{ix} = 'AIRS L3    ';
%% the value
showfeedbacks_robustfit_midlats(ix,1,1) = junk.feedback_ecRad.planck.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,3,1) = junk.feedback_ecRad.o3.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,4,1) = junk.feedback_ecRad.wv.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,5,1) = junk.feedback_ecRad.skt.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(1);
%% the unc
showfeedbacks_robustfit_midlats(ix,1,2) = junk.feedback_ecRad.planck.robustfit_midlats(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_midlats(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,3,2) = junk.feedback_ecRad.o3.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,4,2) = junk.feedback_ecRad.wv.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,5,2) = junk.feedback_ecRad.skt.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(2);

if ~exist('climcapsL3_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'climcapsL3_spectral_olr');
  climcapsL3_spectral_olr = junk.climcapsL3_spectral_olr;
end
ix = 5; junk = climcapsL3_spectral_olr;
strfeedbacks{ix} = 'CLIMCAPS L3';
%% the value
showfeedbacks_robustfit_midlats(ix,1,1) = junk.feedback_ecRad.planck.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,3,1) = junk.feedback_ecRad.o3.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,4,1) = junk.feedback_ecRad.wv.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,5,1) = junk.feedback_ecRad.skt.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(1);
%% the unc
showfeedbacks_robustfit_midlats(ix,1,2) = junk.feedback_ecRad.planck.robustfit_midlats(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_midlats(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,3,2) = junk.feedback_ecRad.o3.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,4,2) = junk.feedback_ecRad.wv.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,5,2) = junk.feedback_ecRad.skt.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(2);

if ~exist('cmip6_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'cmip6_spectral_olr');
  cmip6_spectral_olr = junk.cmip6_spectral_olr;
end
ix = 6; junk = cmip6_spectral_olr;
strfeedbacks{ix} = 'CMIP6      ';
%% the value
showfeedbacks_robustfit_midlats(ix,1,1) = junk.feedback_ecRad.planck.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,3,1) = junk.feedback_ecRad.o3.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,4,1) = junk.feedback_ecRad.wv.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,5,1) = junk.feedback_ecRad.skt.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(1);
%% the unc
showfeedbacks_robustfit_midlats(ix,1,2) = junk.feedback_ecRad.planck.robustfit_midlats(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_midlats(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,3,2) = junk.feedback_ecRad.o3.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,4,2) = junk.feedback_ecRad.wv.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,5,2) = junk.feedback_ecRad.skt.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(2);

if ~exist('amip6_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'amip6_spectral_olr');
  amip6_spectral_olr = junk.amip6_spectral_olr;
end
ix = 7; junk = amip6_spectral_olr;
strfeedbacks{ix} = 'AMIP6      ';
%% the value
showfeedbacks_robustfit_midlats(ix,1,1) = junk.feedback_ecRad.planck.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,3,1) = junk.feedback_ecRad.o3.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,4,1) = junk.feedback_ecRad.wv.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,5,1) = junk.feedback_ecRad.skt.robustfit_midlats(1);
showfeedbacks_robustfit_midlats(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(1);
%% the unc
showfeedbacks_robustfit_midlats(ix,1,2) = junk.feedback_ecRad.planck.robustfit_midlats(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_midlats(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,3,2) = junk.feedback_ecRad.o3.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,4,2) = junk.feedback_ecRad.wv.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,5,2) = junk.feedback_ecRad.skt.robustfit_midlats(2);
showfeedbacks_robustfit_midlats(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_midlats(2);

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
ixx = ix;
ixshow = 5;
showfeedbacks_robustfit_midlats(1:ixx,7,1) = sum(squeeze(showfeedbacks_robustfit_midlats(1:ixx,[1 2 3 4],1)),2);
junk = showfeedbacks_robustfit_midlats(1:ixx,[1 2 3 4],2);
junk = sqrt(sum(junk.*junk,2));
showfeedbacks_robustfit_midlats(1:ixx,7,2) = junk;

for ix = 1 : ixshow
  %fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},showfeedbacks_robustfit_midlats(ix,[1 2 3 4 7],1));
  junk = [showfeedbacks_robustfit_midlats(ix,[1 2 3 4 7],1) showfeedbacks_robustfit_midlats(ix,[1 2 3 4 7],2)];
  junk = junk([1 6 2 7 3 8 4 9 5 10]);
  fprintf(1,'%s %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f    %6.3f +/- %6.3f \n',strfeedbacks{ix},junk);
end
trends_paper_show = showfeedbacks_robustfit_midlats(1:ixshow,[1 2 3 4 7]);

figure(1); clf
bar(trends_paper_show')
ylabel('Feedback W/m2/K');
hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south');
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)
title('Midlat feedbacks')
