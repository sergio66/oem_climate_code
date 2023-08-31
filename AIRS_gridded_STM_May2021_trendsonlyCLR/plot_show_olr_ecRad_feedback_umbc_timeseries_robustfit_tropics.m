clear showfeedbacks* strfeedbacks

ix = 1; iNumYears = 05; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '05 years ';
%% the value
showfeedbacks_robustfit_tropics(ix,1,1) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_robustfit_tropics(1);  %% (1 == value, 2 == unc)
showfeedbacks_robustfit_tropics(ix,2,1) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,3,1) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,4,1) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,5,1) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,6,1) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_robustfit_tropics(1);
%% the unc
showfeedbacks_robustfit_tropics(ix,1,2) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_robustfit_tropics(2); 
showfeedbacks_robustfit_tropics(ix,2,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,3,2) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,4,2) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,5,2) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,6,2) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_robustfit_tropics(2);

ix = 2; iNumYears = 10; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '10 years ';
%% the value
showfeedbacks_robustfit_tropics(ix,1,1) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,2,1) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,3,1) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,4,1) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,5,1) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,6,1) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_robustfit_tropics(1);
%% the unc
showfeedbacks_robustfit_tropics(ix,1,2) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,2,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,3,2) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,4,2) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,5,2) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,6,2) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_robustfit_tropics(2);

ix = 3; iNumYears = 15; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '15 years ';
%% the value
showfeedbacks_robustfit_tropics(ix,1,1) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,2,1) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,3,1) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,4,1) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,5,1) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,6,1) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_robustfit_tropics(1);
%% the unc
showfeedbacks_robustfit_tropics(ix,1,2) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,2,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,3,2) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,4,2) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,5,2) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,6,2) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_robustfit_tropics(2);

ix = 4; iNumYears = 20; junk = load(['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat']);
strfeedbacks{ix} = '20 years ';
%% the value
showfeedbacks_robustfit_tropics(ix,1,1) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,2,1) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,3,1) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,4,1) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,5,1) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_robustfit_tropics(1);
showfeedbacks_robustfit_tropics(ix,6,1) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_robustfit_tropics(1);
%% the unc
showfeedbacks_robustfit_tropics(ix,1,2) = junk.umbc_spectral_olr.feedback_ecRad.planck_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,2,2) = junk.umbc_spectral_olr.feedback_ecRad.lapse_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,3,2) = junk.umbc_spectral_olr.feedback_ecRad.o3_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,4,2) = junk.umbc_spectral_olr.feedback_ecRad.wv_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,5,2) = junk.umbc_spectral_olr.feedback_ecRad.skt_ecRad_robustfit_tropics(2);
showfeedbacks_robustfit_tropics(ix,6,2) = junk.umbc_spectral_olr.feedback_ecRad.ptemp_co2_ecRad_robustfit_tropics(2);

%% the 6 feedbacks are feedbacks : planck lapse o3 wv skt tz/co2
%% but longwave feedback is um of first 4
ixx = ix;
showfeedbacks_robustfit_tropics(1:ixx,7,1) = sum(squeeze(showfeedbacks_robustfit_tropics(1:ixx,[1 2 3 4],1)),2);
junk = showfeedbacks_robustfit_tropics(1:ixx,[1 2 3 4],2);
junk = sqrt(sum(junk.*junk,2));
showfeedbacks_robustfit_tropics(1:ixx,7,2) = junk;

for ix = 1 : 4
  %fprintf(1,'%s %5.2f %5.2f %5.2f %5.2f    %5.2f \n',strfeedbacks{ix},showfeedbacks_robustfit_tropics(ix,[1 2 3 4 7],1));
  junk = [showfeedbacks_robustfit_tropics(ix,[1 2 3 4 7],1) showfeedbacks_robustfit_tropics(ix,[1 2 3 4 7],2)];
  junk = junk([1 6 2 7 3 8 4 9 5 10]);
  fprintf(1,'%s %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f  %6.3f +/- %6.3f    %6.3f +/- %6.3f \n',strfeedbacks{ix},junk);
end
trends_paper_show = showfeedbacks_robustfit_tropics(1:ixx,[1 2 3 4 7]);

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

