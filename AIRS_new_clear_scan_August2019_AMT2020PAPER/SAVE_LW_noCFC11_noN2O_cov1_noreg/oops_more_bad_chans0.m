%% this is for kg_ab0
clear all
%load results_nucal_LWMW_ab0.mat
load results_nucal_LWMWSW_ab0_allchans.mat
plot(fairs(g),nanmean(fit_to_rates(g,:)-input_rates(g,:),2),'+-')

r = nanmean(fit_to_rates(g,:)-input_rates(g,:),2);
kgood = find( abs(r) < 0.008 & (fairs(g) < 1080 | fairs(g) > 1145));
kbad0 = setdiff(1:length(g),kgood);

plot(fairs(g),nanmean(fit_to_rates(g,:)-input_rates(g,:),2),'+-',fairs(g(kbad0)),nanmean(fit_to_rates(g(kbad0),:)-input_rates(g(kbad0),:),2),'ro')

load kbad_chans.mat
%kbad = union(kbad,kbad0);

comment = 'see oops_more_bad_chans0.m';
save kbad_chans0.mat kbad kbad0 comment

