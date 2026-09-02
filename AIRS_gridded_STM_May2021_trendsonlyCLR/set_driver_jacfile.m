function [driver,iVersJac,iOldORNew,iXJac,topts] = set_driver_jacfile(driver0,settings,topts0);

driver = driver0;
topts  = topts0;

iOldORNew = -1;

iXJac = settings.iXJac;
%if driver.i16daytimestep > 0
%  iXJac = 0;  %% const geo kcarta jacs, default for trends
%  iXJac = -1; %% const geo saarta jacs,on the fly
%  iXJac = 1;  %% varying geo sarta jacs
%  iXJac = 2;  %% varying geo kcarta jacs, default for anomaly
%  iXJac = 3;  %% should really be Q(X --> 1) but hard to get ERA5 conditions for this!
%end

if driver.i16daytimestep < 0
  set_trends_jacfile
elseif driver.i16daytimestep > 0
  set_anomalies_jacfile
end

disp('in strow_override_defaults_latbins_AIRS_fewlays.m   : finished   [driver,iVersJac,iOldORNew,iXJac,topts] = set_driver_jacfile(driver,settings,topts);')
disp('now onto set_the_jacobians.m')
disp(' ')
