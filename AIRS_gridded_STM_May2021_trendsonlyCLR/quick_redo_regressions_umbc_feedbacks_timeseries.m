clear showfeedbacks* strfeedbacks

iaComputeWhichFeedback = [0];     %% compute + plot feedbacks only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a.topts.dataset = 12; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset12_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 15;  %% use CarbonTracker CO2 trends
feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];

fprintf(1,'iNumYears = %2i .. reading in %s \n',iNumYears,strUMBC);
loader = ['load ' strUMBC];
eval(loader);

read_fileMean17years
h = hMean17years;
p = pMean17years;

loader = ['load ' feedbacknameUMBC];
eval(loader)
umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
saver = ['save ' feedbacknameUMBC ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays'];  %% if you only want to save UMBC
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a.topts.dataset = 11; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset11_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 10;  %% use CarbonTracker CO2 trends
feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];

fprintf(1,'iNumYears = %2i .. reading in %s \n',iNumYears,strUMBC);
loader = ['load ' strUMBC];
eval(loader);

read_fileMean17years
h = hMean17years;
p = pMean17years;

loader = ['load ' feedbacknameUMBC];
eval(loader)
umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
saver = ['save ' feedbacknameUMBC ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays'];  %% if you only want to save UMBC
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error(' ;;; ')

a.topts.dataset = 10; strUMBC = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset10_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat'; iNumYears = 05;  %% use CarbonTracker CO2 trends
feedbacknameUMBC = ['/asl/s1/sergio/JUNK/olr_feedbacks_UMBC_numyears_' num2str(iNumYears,'%02d') '.mat'];

fprintf(1,'iNumYears = %2i .. reading in %s \n',iNumYears,strUMBC);
loader = ['load ' strUMBC];
eval(loader);

read_fileMean17years
h = hMean17years;
p = pMean17years;

loader = ['load ' feedbacknameUMBC];
eval(loader)
umbc_spectral_olr = compute_feedbacks_generic_ecRad(h,p,results,results(:,6)',deltaT,fracWV,fracO3,umbc_spectral_olr,iaComputeWhichFeedback,rlat65,rlon73,+1,'UMBC');
saver = ['save ' feedbacknameUMBC ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays'];  %% if you only want to save UMBC
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_show_olr_ecRad_feedback_umbc_timeseries_globalsstfit

