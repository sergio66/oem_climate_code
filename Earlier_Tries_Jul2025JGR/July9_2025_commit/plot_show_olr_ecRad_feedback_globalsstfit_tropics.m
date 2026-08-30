clear showfeedbacks_globalSST_tropics* strfeedbacks

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
showfeedbacks_globalSST_tropics(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_tropics;

if ~exist('merra2_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'merra2_spectral_olr');
  merra2_spectral_olr = junk.merra2_spectral_olr;
end
ix = 2; junk = merra2_spectral_olr;
strfeedbacks{ix} = 'MERRA2     ';
showfeedbacks_globalSST_tropics(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_tropics;

if ~exist('umbc_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'],'umbc_spectral_olr');
  umbc_spectral_olr = junk.umbc_spectral_olr;
end
ix = 3; junk = umbc_spectral_olr;
strfeedbacks{ix} = 'THIS WORK  ';
showfeedbacks_globalSST_tropics(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_tropics;

if ~exist('airsL3_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'airsL3_spectral_olr');
  airsL3_spectral_olr = junk.airsL3_spectral_olr;
end
ix = 4; junk = airsL3_spectral_olr;
strfeedbacks{ix} = 'AIRS L3    ';
showfeedbacks_globalSST_tropics(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_tropics;

if ~exist('climcapsL3_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'climcapsL3_spectral_olr');
  climcapsL3_spectral_olr = junk.climcapsL3_spectral_olr;
end
ix = 5; junk = climcapsL3_spectral_olr;
strfeedbacks{ix} = 'CLIMCAPS L3';
showfeedbacks_globalSST_tropics(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_tropics;

if ~exist('cmip6_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'cmip6_spectral_olr');
  cmip6_spectral_olr = junk.cmip6_spectral_olr;
end
ix = 6; junk = cmip6_spectral_olr;
strfeedbacks{ix} = 'CMIP6      ';
showfeedbacks_globalSST_tropics(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_tropics;

if ~exist('amip6_spectral_olr')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '.mat'],'amip6_spectral_olr');
  amip6_spectral_olr = junk.amip6_spectral_olr;
end
ix = 7; junk = amip6_spectral_olr;
strfeedbacks{ix} = 'AMIP6      ';
showfeedbacks_globalSST_tropics(ix,1) = junk.feedback_ecRad.planck.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,2) = junk.feedback_ecRad.lapse.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,3) = junk.feedback_ecRad.o3.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,4) = junk.feedback_ecRad.wv.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,5) = junk.feedback_ecRad.skt.globalSST_weighted_tropics;
showfeedbacks_globalSST_tropics(ix,6) = junk.feedback_ecRad.ptemp_co2.globalSST_weighted_tropics;

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
ixx = ix;
showfeedbacks_globalSST_tropics(1:ixx,7) = sum(showfeedbacks_globalSST_tropics(1:ixx,[1 2 3 4]),2);

for ix = 1 : 5
  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},showfeedbacks_globalSST_tropics(ix,[1 2 3 4 7]));
end
trends_paper_show = showfeedbacks_globalSST_tropics(1:5,[1 2 3 4 7]);

figure(1); clf
bar(trends_paper_show')
ylabel('TROPICS Feedback W/m2/K');
hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south');
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)

