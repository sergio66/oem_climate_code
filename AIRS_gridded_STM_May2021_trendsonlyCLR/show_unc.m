%% see build_cov_matrices.m

%% if cX_wide is tiny, transition is very smooth and gradual ... if cX_wide is large, almost step function transition from trop to strat
%% reg ~ Scov + S_tikonov ~ (sig_q)^-1 + alpha L1'L1       where sig_q is uncertainty in param (eg K^2 for temperature)
%%   if (sig_q)^-1 = 1/(uncertainty)^2 >> alpha ==> covariance regularization
%%      (sig_q)^-1 = 1/(uncertainty)^2 << alpha ==> tikonov    regularization
%%   if sig_q -> 0   then you say you are VERY sure about a-priori ==> do not change ==> delta(param) --> 0
%%      sig_q -> INF then you say you are DO NOT TRUST    a-priori ==>        change ==> delta(param) --> bigly wigly
%%   if alpha -> 0   then you say you are DO NOT TRUST    a-priori ==>        change ==> delta(param) --> bigly wigly
%%      alpha -> INF then you say you are VERY sure about a-priori ==> do not change ==> delta(param) --> 0

disp(' ')
sig_q = min(oem.cov_set(2),oem.cov_set(3)); sig_q = (1./sig_q)^2; alphaa = oem.alpha_temp;
fprintf(1,' unc T      : trop/strat = %8.6f %8.6f K  Tikonov = %8.6f K ==> sigq,alpha = %8.6f %8.6f \n',oem.cov_set(2),oem.cov_set(3),oem.alpha_temp,sig_q,alphaa)
if sig_q/alphaa > 10
  disp('(sig_q)^-1 = 1/(uncertainty)^2 >> alpha ==> covariance regularization for T')
elseif sig_q/alphaa < 0.1
  disp('(sig_q)^-1 = 1/(uncertainty)^2 << alpha ==> tikonov regularization for T')
else
  disp('(sig_q)^-1 = 1/(uncertainty)^2 ~ alpha ==> covariance + tikonov regularization for T')
end

disp(' ')
sig_q = min(oem.cov_set(5),oem.cov_set(6)); sig_q = (1./sig_q)^2; alphaa = oem.alpha_water;
fprintf(1,' unc fracWV : trop/strat = %8.6f %8.6f    Tikonov = %8.6f ==> sigq,alpha = %8.6f %8.6f   \n',oem.cov_set(5),oem.cov_set(6),oem.alpha_water,sig_q,alphaa)
if sig_q/alphaa > 10
  disp('(sig_q)^-1 = 1/(uncertainty)^2 >> alpha ==> covariance regularization for WV')
elseif sig_q/alphaa < 0.1
  disp('(sig_q)^-1 = 1/(uncertainty)^2 << alpha ==> tikonov regularization for WV')
else
  disp('(sig_q)^-1 = 1/(uncertainty)^2 ~ alpha ==> covariance + tikonov regularization for WV')
end

disp(' ')
fprintf(1,' unc fracO3 : trop/strat = %8.6f %8.6f    Tikonov = %8.6f ==> sigq,alpha = %8.6f %8.6f   \n',oem.cov_set(8),oem.cov_set(9),oem.alpha_ozone,sig_q,alphaa)
if sig_q/alphaa > 10
  disp('(sig_q)^-1 = 1/(uncertainty)^2 >> alpha ==> covariance regularization for O3')
elseif sig_q/alphaa < 0.1
  disp('(sig_q)^-1 = 1/(uncertainty)^2 << alpha ==> tikonov regularization for O3')
else
  disp('(sig_q)^-1 = 1/(uncertainty)^2 ~ alpha ==> covariance + tikonov regularization for O3')
end
disp(' ')
