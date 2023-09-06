clear showfeedbacks* strfeedbacks
clear all

figure(1); clf
figure(2); clf
figure(3); clf
figure(4); clf
figure(5); clf

iaComputeWhichFeedback = [0];     %% compute + plot feedbacks only

iNumYears = 20;
read_fileMean17years
h = hMean17years;
p = pMean17years;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('05 years')
a.topts.dataset = 10; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset10_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 05;  %% use CarbonTracker CO2 trends
feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];

fprintf(1,'iNumYears = %2i .. reading in %s \n',iNumYears,strUMBC);
loader = ['load ' strUMBC];
eval(loader);

loader = ['load ' feedbacknameUMBC];
eval(loader)
%umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
 umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',[]    ,[]    ,[]    ,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
saver = ['save ' feedbacknameUMBC ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays'];  %% if you only want to save UMBC
eval(saver)

figure(4); figure(5); figure(6); figure(80); rett(2); %% disp('RET to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('10 years')
a.topts.dataset = 11; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset11_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 10;  %% use CarbonTracker CO2 trends
feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];

fprintf(1,'iNumYears = %2i .. reading in %s \n',iNumYears,strUMBC);
loader = ['load ' strUMBC];
eval(loader);

loader = ['load ' feedbacknameUMBC];
eval(loader)
%umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
 umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',[]    ,[]    ,[]    ,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
saver = ['save ' feedbacknameUMBC ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays'];  %% if you only want to save UMBC
eval(saver)

figure(4); figure(5); figure(6); figure(80); rett(2); %% disp('RET to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('15 years')
a.topts.dataset = 12; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset12_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 15;  %% use CarbonTracker CO2 trends
feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];

fprintf(1,'iNumYears = %2i .. reading in %s \n',iNumYears,strUMBC);
loader = ['load ' strUMBC];
eval(loader);

loader = ['load ' feedbacknameUMBC];
eval(loader)
%umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
 umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',[]    ,[]    ,[]    ,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
saver = ['save ' feedbacknameUMBC ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays'];  %% if you only want to save UMBC
eval(saver)

figure(4); figure(5); figure(6); figure(80); rett(2); %% disp('RET to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('20 years')
a.topts.dataset = 09; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 20;  %% use CarbonTracker CO2 trends
feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];

fprintf(1,'iNumYears = %2i .. reading in %s \n',iNumYears,strUMBC);
loader = ['load ' strUMBC];
eval(loader);

loader = ['load ' feedbacknameUMBC];
eval(loader)
%umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
 umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',[]    ,[]    ,[]    ,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
saver = ['save ' feedbacknameUMBC ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays'];  %% if you only want to save UMBC
iSave = input('Save 20 year stuff??? (default -1) : ');
if length(iSave) == 0
  iSave = -1;
end
if iSave > 0
  eval(saver)
end

figure(4); figure(5); figure(6); figure(80); rett(2); %% disp('RET to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('DONE, now showing results')
%plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfit
%plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfit_smooth
%plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfitsmooth
plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfitsmooth2

