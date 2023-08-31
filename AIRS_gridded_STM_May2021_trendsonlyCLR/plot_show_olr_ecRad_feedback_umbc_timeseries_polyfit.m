clear showfeedbacks* strfeedbacks

ix = 1; iNumYears = 05; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '05 years ';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_polyfit;

ix = 2; iNumYears = 10; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '10 years ';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_polyfit;

ix = 3; iNumYears = 15; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '15 years ';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_polyfit;

ix = 4; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '20 years ';
showfeedbacks(ix,1) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_polyfit;
showfeedbacks(ix,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_polyfit;
showfeedbacks(ix,3) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_polyfit;
showfeedbacks(ix,4) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_polyfit;
showfeedbacks(ix,5) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_polyfit;
showfeedbacks(ix,6) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_polyfit;

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
ixx = ix;
showfeedbacks(1:ixx,7) = sum(showfeedbacks(1:ixx,[1 2 3 4]),2);

for ix = 1 : 4
  fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},showfeedbacks(ix,[1 2 3 4 7]));
end
trends_paper_show = showfeedbacks(1:ixx,[1 2 3 4 7]);

figure(1); clf
bar(trends_paper_show');
ylabel('Feedback W/m2/K');
hl = legend('05','10','15','20','location','south');
xstr = {'Planck','Lapse','Ozone','Water Vapor','SUM'};
set(gca,'xticklabels',xstr)
xtickangle(45)

%{
dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
figure(1); aslprint([dir0 'feedbackparams_05_10_15_20yrs.pdf'])
%}

