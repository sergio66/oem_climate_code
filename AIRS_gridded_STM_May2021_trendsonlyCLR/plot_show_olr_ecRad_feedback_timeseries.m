clear showfeedbacks

ix = 1; iNumYears = 05; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback.ptemp_co2_ecRad_polyfit;

ix = 2; iNumYears = 10; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback.ptemp_co2_ecRad_polyfit;

ix = 3; iNumYears = 15; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback.ptemp_co2_ecRad_polyfit;

ix = 4; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback.ptemp_co2_ecRad_polyfit;

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
showfeedbacks(1:ix,7) = sum(showfeedbacks(:,[1 2 3 4]),2);

figure(1); clf
bar(showfeedbacks(1:4,[1 2 3 4 7])')
ylabel('Feedback W/m2/K');
hl = legend('05','10','15','20','location','south');
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)

%{
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
figure(1); aslprint([dir0 'feedbackparams_05_10_15_20yrs.pdf'])
%}

