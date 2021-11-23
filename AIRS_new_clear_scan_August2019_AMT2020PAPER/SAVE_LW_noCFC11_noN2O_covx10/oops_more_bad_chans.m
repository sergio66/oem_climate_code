%% this is for kg_fixed_ab
clear all
load results_nucal_LWMW.mat
plot(fairs(g),nanmean(fit_to_rates(g,:)-input_rates(g,:),2),'+-')

r = nanmean(fit_to_rates(g,:)-input_rates(g,:),2);
kgood = find( abs(r) < 0.008 & (fairs(g) < 1080 | fairs(g) > 1145));
kbad = setdiff(1:length(g),kgood);

plot(fairs(g),nanmean(fit_to_rates(g,:)-input_rates(g,:),2),'+-',fairs(g(kbad)),nanmean(fit_to_rates(g(kbad),:)-input_rates(g(kbad),:),2),'ro')

comment = 'see oops_more_bad_chans.m';
save kbad_chans.mat kbad comment
