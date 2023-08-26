clear showfeedbacks

ix = 1; junk = era5_spectral_olr;
showfeedbacks(ix,1) = junk.feedback.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.feedback.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.feedback.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.feedback.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.feedback.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.feedback.ptemp_co2_ecRad_polyfit;

ix = 2; junk = merra2_spectral_olr;
showfeedbacks(ix,1) = junk.feedback.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.feedback.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.feedback.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.feedback.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.feedback.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.feedback.ptemp_co2_ecRad_polyfit;

ix = 3; junk = umbc_spectral_olr;
showfeedbacks(ix,1) = junk.feedback.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.feedback.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.feedback.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.feedback.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.feedback.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.feedback.ptemp_co2_ecRad_polyfit;

ix = 4; junk = airsL3_spectral_olr;
showfeedbacks(ix,1) = junk.feedback.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.feedback.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.feedback.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.feedback.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.feedback.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.feedback.ptemp_co2_ecRad_polyfit;

ix = 5; junk = climcapsL3_spectral_olr;
showfeedbacks(ix,1) = junk.feedback.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.feedback.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.feedback.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.feedback.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.feedback.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.feedback.ptemp_co2_ecRad_polyfit;

ix = 6; junk = cmip6_spectral_olr;
showfeedbacks(ix,1) = junk.feedback.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.feedback.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.feedback.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.feedback.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.feedback.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.feedback.ptemp_co2_ecRad_polyfit;

ix = 7; junk = amip6_spectral_olr;
showfeedbacks(ix,1) = junk.feedback.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.feedback.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.feedback.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.feedback.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.feedback.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.feedback.ptemp_co2_ecRad_polyfit;

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
showfeedbacks(1:ix,7) = sum(showfeedbacks(:,[1 2 3 4]),2);

bar(showfeedbacks(1:5,[1 2 3 4 7])')
ylabel('Feedback W/m2/K');
hl = legend('ERA5','MERRA2','THIS WORK','AIRS L3','CLIMCAPS L3','location','south');
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)

%{
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
aslprint([dir0 'feedbackparams_20yrs_UMBC_ERA5_MERRA2_L3.pdf'])
%}

