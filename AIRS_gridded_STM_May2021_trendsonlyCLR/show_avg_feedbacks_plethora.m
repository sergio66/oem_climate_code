cosrlat = cos(rlat'*pi/180);

disp(' ')
do_avg_feedback1     %% crude attempt at zonal avg
disp(' ')
do_avg_feedback1cos  %% crude attempt at zonal avg with cosine(lat) wgt
disp(' ')
do_avg_feedback2     %% better attempt at zonal avg
disp(' ')
do_avg_feedback2cos  %% better attempt at zonal avg with cosine(rlat) wgt  BEST

disp(' ')
do_avg_feedback3     %% used in trends paper, try 1 : (polyfit)
disp(' ')
do_avg_feedback4     %% used in trends paper, try 2 : (globalSST division)
disp(' ')
do_avg_feedback5     %% used in trends paper, try 3 : (robustfit division)

%%% Ryan suggested normalizing using dERASST for all, instead of the individual dXSST X=ERA or CMIP6 or UMBC or AIRSL3 
%%% iERAnorm = input('Do you wish to redo the feedback by using only dERA SKT instead of individual d SKT? (-1/default) no (+1) yes : ');
%%% if length(iERAnorm) == 0
%%%   iERAnorm = -1;
%%% end
%%% if iERAnorm > 0
%%%   redo_feedbacks_dERA5ST_dt
%%%   do_avg_feedback2cos_dERA5ST_dt  %% better attempt at zonal avg with cosine(rlat) wgt  BEST
%%% end
