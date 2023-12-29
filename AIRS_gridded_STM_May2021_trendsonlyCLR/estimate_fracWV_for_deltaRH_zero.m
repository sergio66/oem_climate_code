function dfrac_WV = estimate_fracWV_for_deltaRH_zero(dT,T,RH)

%% ocean, so use Isaac Held blog ideas  https://www.gfdl.noaa.gov/blog_held/47-relative-humidity-over-the-oceans/
%% also see Nadir Jeevanjee : Simple Climate Models, arXiv 2018

%{
dfrac_WV = estimate_fracWV_for_deltaRH_zero(0.01*ones(97,4608),p.ptemp(1:97,:));
pcolor(rlat,p.plevs(1:97,3000),squeeze(nanmean(reshape(dfrac_WV,97,72,64),2))); shading interp; set(gca,'ydir','reverse'); colorbar; caxis([-1 +1]*0.015/10)
set(gca,'yscale','linear'); ylim([100 1000])
%}

%% this estimtte delta(WV)/WV needed to keep RH constant if T changes by dT
%% Lo/Rv = K
Lo = 2.5e6;  %%% J/kg
Rv = 461.52; %%% J/kg/K

%moo = exp(Lo/Rv * dT/ppp.stemp(JOBJOBJOB)/ppp.stemp(JOBJOBJOB))-1;
%moo = Lo/Rv * dT/ppp.stemp(JOBJOBJOB)/ppp.stemp(JOBJOBJOB);
%dfrac_WV = moo; %% remember saturation vapor pressure changes at 0.07/K and we want dRH = 0 BUT THIS 0.07 is for 285 K !!!!!!!!!!!!

dfrac_WV = Lo/Rv * dT./T./T;  %% K * K/K/K = []

if nargin == 3
  %% this is adding on Isaac Held/Nadir Jeevanjee arXiv term which is the RH change (0.01/K)
  %% see Bk 48 of my notes
  %% 0.01 is the RH change per K that NK/IH expect both in PBL and in Free Ttrop
  dfrac_WV = dfrac_WV + 0.01./RH .*dT;
end
