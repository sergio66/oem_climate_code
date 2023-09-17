clear showfeedbacks* strfeedbacks
clear *spectral_olr*

choose_max_ratio_plot

for ii=1:8; figure(ii); clf; end

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER

if exist('umbc_spectral_olr')
  disp('umbc_spectral_olr already exists, you may want to do   "clear *spectral_olr"    before running this')
end

if ~exist('iNumYears')
  disp('WARNING iNumYears DNE ... setting to 20')
  iNumYears = 20;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if ~exist('era5_spectral_olr_withunc')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '_unc_factor' num2str(maxratio,'%0.2f') '.mat'],'era5_spectral_olr');
  era5_spectral_olr_withunc = junk.era5_spectral_olr;
end
ix = 1; junk = era5_spectral_olr_withunc;
strfeedbacks{ix} = 'ERA5       ';
%% the value
showfeedbacks_robustfit_all_withunc(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all_withunc(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all_withunc(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all_withunc(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('merra2_spectral_olr_withunc')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '_unc_factor' num2str(maxratio,'%0.2f') '.mat'],'merra2_spectral_olr');
  merra2_spectral_olr_withunc = junk.merra2_spectral_olr;
end
ix = 2; junk = merra2_spectral_olr_withunc;
strfeedbacks{ix} = 'MERRA2     ';
%% the value
showfeedbacks_robustfit_all_withunc(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all_withunc(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all_withunc(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('umbc_spectral_olr_withunc')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '_unc_factor' num2str(maxratio,'%0.2f') '.mat'],'umbc_spectral_olr');
  umbc_spectral_olr_withunc = junk.umbc_spectral_olr;
end
ix = 3; junk = umbc_spectral_olr_withunc;
strfeedbacks{ix} = 'THIS WORK  ';
%% the value
showfeedbacks_robustfit_all_withunc(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all_withunc(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all_withunc(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('airsL3_spectral_olr_withunc')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '_unc_factor' num2str(maxratio,'%0.2f') '.mat'],'airsL3_spectral_olr');
  airsL3_spectral_olr_withunc = junk.airsL3_spectral_olr;
end
ix = 4; junk = airsL3_spectral_olr_withunc;
strfeedbacks{ix} = 'AIRS L3    ';
%% the value
showfeedbacks_robustfit_all_withunc(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all_withunc(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all_withunc(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('climcapsL3_spectral_olr_withunc')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '_unc_factor' num2str(maxratio,'%0.2f') '.mat'],'climcapsL3_spectral_olr');
  climcapsL3_spectral_olr_withunc = junk.climcapsL3_spectral_olr;
end
ix = 5; junk = climcapsL3_spectral_olr_withunc;
strfeedbacks{ix} = 'CLIMCAPS L3';
%% the value
showfeedbacks_robustfit_all_withunc(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all_withunc(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all_withunc(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('cmip6_spectral_olr_withunc')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_AIRSL3_ERA5_CMIP6_numyears_' num2str(iNumYears,'%02d') '_unc_factor' num2str(maxratio,'%0.2f') '.mat'],'cmip6_spectral_olr');
  cmip6_spectral_olr_withunc = junk.cmip6_spectral_olr;
end
ix = 6; junk = cmip6_spectral_olr_withunc;
strfeedbacks{ix} = 'CMIP6      ';
%% the value
showfeedbacks_robustfit_all_withunc(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all_withunc(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all_withunc(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

if ~exist('amip6_spectral_olr_withunc')
  junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6_numyears_' num2str(iNumYears,'%02d') '_unc_factor' num2str(maxratio,'%0.2f') '.mat'],'amip6_spectral_olr');
  amip6_spectral_olr_withunc = junk.amip6_spectral_olr;
end
ix = 7; junk = amip6_spectral_olr_withunc;
strfeedbacks{ix} = 'AMIP6      ';
%% the value
showfeedbacks_robustfit_all_withunc(ix,1,1) = junk.feedback_ecRad.planck.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,2,1) = junk.feedback_ecRad.lapse.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,3,1) = junk.feedback_ecRad.o3.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,4,1) = junk.feedback_ecRad.wv.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,5,1) = junk.feedback_ecRad.skt.robustfit_all(1);
showfeedbacks_robustfit_all_withunc(ix,6,1) = junk.feedback_ecRad.ptemp_co2.robustfit_all(1);
%% the unc
showfeedbacks_robustfit_all_withunc(ix,1,2) = junk.feedback_ecRad.planck.robustfit_all(2);   %% (1 == value, 2 == unc)
showfeedbacks_robustfit_all_withunc(ix,2,2) = junk.feedback_ecRad.lapse.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,3,2) = junk.feedback_ecRad.o3.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,4,2) = junk.feedback_ecRad.wv.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,5,2) = junk.feedback_ecRad.skt.robustfit_all(2);
showfeedbacks_robustfit_all_withunc(ix,6,2) = junk.feedback_ecRad.ptemp_co2.robustfit_all(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
disp('raw calcs, just using trends')

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
hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',10);
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)
title('Global feedbacks')
ylim([-1 +1]*5);

%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
disp('unc calcs, using trends + unc trends')

ixx = ix;
showfeedbacks_robustfit_all_withunc(1:ixx,7,1) = sum(squeeze(showfeedbacks_robustfit_all_withunc(1:ixx,[1 2 3 4],1)),2);
junk = showfeedbacks_robustfit_all_withunc(1:ixx,[1 2 3 4],2);
junk = sqrt(sum(junk.*junk,2));
showfeedbacks_robustfit_all_withunc(1:ixx,7,2) = junk;

ixshow = 5;
for ix = 1 : ixshow
  %fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},showfeedbacks_robustfit_all_withunc(ix,[1 2 3 4 7],1));
  junk = [showfeedbacks_robustfit_all_withunc(ix,[1 2 3 4 7],1) showfeedbacks_robustfit_all_withunc(ix,[1 2 3 4 7],2)];
  junk = junk([1 6 2 7 3 8 4 9 5 10]);
  fprintf(1,'%s %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f    %6.3f +/- %6.3f \n',strfeedbacks{ix},junk);
end
trends_paper_show_withunc = showfeedbacks_robustfit_all_withunc(1:ixshow,[1 2 3 4 7]);

disp(' ')
disp('difference between the two')
for ix = 1 : 5
  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},abs(showfeedbacks_robustfit_all_withunc(ix,[1 2 3 4 7])-showfeedbacks_robustfit_all(ix,[1 2 3 4 7])));
end


figure(2); clf
bar(trends_paper_show_withunc')
ylabel('Feedback W/m2/K');
hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south','fontsize',10);
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)
title('Global feedbacks with unc added')
ylim([-1 +1]*5);

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

if ~exist('iSmooth')
  iSmooth = 5;
  iSmooth = 10;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xsin = sin(rlat*pi/180);
  xtick = [-1 -sqrt(3)/2 -sqrt(2)/2 -1/2 -(0.25+0.01) 0 +(0.25+0.01) +1/2 +sqrt(2)/2 +sqrt(3)/2 +1]; %% -90 -60 -45 -30 -15 0 +15 +30 +45 +60 +90
  xtick = [-1            -sqrt(2)/2      -(0.25+0.01) 0 +(0.25+0.01)      +sqrt(2)/2            +1]; %% -90     -45     -15 0 +15     +45     +90
  xticklab = cellstr(num2str(round(180/pi*asin((xtick(:)))), '%d'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf;
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
ystr1 = 'W/m2/K';
ystr3 = 'trends only';

tafov(1).YLabel.String = ystr1; tafov(1).YLabel.FontSize = 18;
tafov(3).YLabel.String = ystr3; tafov(3).YLabel.FontSize = 18;
tafov(3).XLabel.String = xstr;  tafov(3).XLabel.FontSize = 18;
tafov(4).XLabel.String = xstr;  tafov(4).XLabel.FontSize = 18;

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

figure(4); clf;
ta = tiledlayout(2,2,'TileSpacing','compact', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
plot(xsin,smooth(era5_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on; hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','north','fontsize',8);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; 

tafov(2) = nexttile;
plot(xsin,smooth(era5_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(3) = nexttile;
plot(xsin,smooth(era5_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(4) = nexttile;
plot(xsin,smooth(era5_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(merra2_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),xsin,smooth(airsL3_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
      xsin,smooth(climcapsL3_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),'linewidth',2);
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
ystr1 = 'W/m2/K';
ystr3 = 'trends+unc';

tafov(1).YLabel.String = ystr1; tafov(1).YLabel.FontSize = 18;
tafov(3).YLabel.String = ystr3; tafov(3).YLabel.FontSize = 18;
tafov(3).XLabel.String = xstr;  tafov(3).XLabel.FontSize = 18;
tafov(4).XLabel.String = xstr;  tafov(4).XLabel.FontSize = 18;

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
ta = tiledlayout(2,2,'TileSpacing','compact', 'Padding','None');
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

tafov(1) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1)-era5_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(merra2_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1)-merra2_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1)-umbc_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(airsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1)-airsL3_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.planck.robustfit_latbin(:,1)-climcapsL3_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on; hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','north','fontsize',8);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');; 

tafov(2) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1)-era5_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(merra2_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1)-merra2_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1)-umbc_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(airsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1)-airsL3_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.lapse.robustfit_latbin(:,1)-climcapsL3_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(3) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1)-era5_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(merra2_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1)-merra2_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1)-umbc_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(airsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1)-airsL3_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.o3.robustfit_latbin(:,1)-climcapsL3_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1),iSmooth),'linewidth',2);
  plotaxis2; box on;
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');

tafov(4) = nexttile;
plot(xsin,smooth(era5_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1)-era5_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(merra2_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1)-merra2_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(umbc_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1)-umbc_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(airsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1)-airsL3_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),...
     xsin,smooth(climcapsL3_spectral_olr.feedback_ecRad.wv.robustfit_latbin(:,1)-climcapsL3_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1),iSmooth),'linewidth',2);
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
ystr1 = 'W/m2/K';
ystr3 = '\delta';

tafov(1).YLabel.String = ystr1; tafov(1).YLabel.FontSize = 18;
tafov(3).YLabel.String = ystr3; tafov(3).YLabel.FontSize = 18;
tafov(3).XLabel.String = xstr;  tafov(3).XLabel.FontSize = 18;
tafov(4).XLabel.String = xstr;  tafov(4).XLabel.FontSize = 18;

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

boo1x = era5_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1) + era5_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
       era5_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1) + era5_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1);
boo2x = merra2_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1) + merra2_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
       merra2_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1) + merra2_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1);
boo3x = umbc_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1) + umbc_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
       umbc_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1) + umbc_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1);
boo4x = airsL3_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1) + airsL3_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
       airsL3_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1) + airsL3_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1);
boo5x = climcapsL3_spectral_olr_withunc.feedback_ecRad.planck.robustfit_latbin(:,1) + climcapsL3_spectral_olr_withunc.feedback_ecRad.lapse.robustfit_latbin(:,1) + ...
       climcapsL3_spectral_olr_withunc.feedback_ecRad.o3.robustfit_latbin(:,1) + climcapsL3_spectral_olr_withunc.feedback_ecRad.wv.robustfit_latbin(:,1);

figure(6)
plot(xsin,smooth(boo1,iSmooth),xsin,smooth(boo2,iSmooth),xsin,smooth(boo3,iSmooth),xsin,smooth(boo4,iSmooth),xsin,smooth(boo5,iSmooth),'linewidth',2);
  plotaxis2;
     hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','best','fontsize',8);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');
title('only retrieved trends')

figure(7)
plot(xsin,smooth(boo1x,iSmooth),xsin,smooth(boo2x,iSmooth),xsin,smooth(boo3x,iSmooth),xsin,smooth(boo4x,iSmooth),xsin,smooth(boo5x,iSmooth),'linewidth',2);
  plotaxis2;
     hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','best','fontsize',8);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');
title('retrieved trends + unc trends')

figure(8)
plot(xsin,smooth(boo1-boo1x,iSmooth),xsin,smooth(boo2-boo2x,iSmooth),xsin,smooth(boo3-boo3x,iSmooth),xsin,smooth(boo4-boo4x,iSmooth),xsin,smooth(boo5-boo5x,iSmooth),'linewidth',2);
  plotaxis2;
     hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','best','fontsize',8);
  xlim([-1 +1]); set(gca,'XTick',xtick,'XTickLabel',xticklab,'TickLabelInterpreter','tex');
title('difference raw-with unc')

figure(6); xlabel('Latitude'); ylabel('\lambda = \Sigma_{i} \lambda_{i}','interpreter','tex'); ylim([-1 +1]*2.5); ylim([-2.5 1.0])
figure(7); xlabel('Latitude'); ylabel('\lambda = \Sigma_{i} \lambda_{i}','interpreter','tex'); ylim([-1 +1]*2.5); ylim([-2.5 1.0])
figure(8); xlabel('Latitude'); ylabel('\lambda = \Sigma_{i} \lambda_{i}','interpreter','tex'); ylim([-1 +1]*1.0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

