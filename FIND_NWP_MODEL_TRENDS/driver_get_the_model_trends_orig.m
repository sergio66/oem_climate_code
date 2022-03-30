iAorC = input('CMIP6 (-1,default) or AMIP6 (+1) : ');
if length(iAorC) == 0
  iAorC = -1;
end
iJorC = input('Chris Barnet NUCAPS (-1) or Joel Susskind AIRS L2 (+1,default) : ');
if length(iJorC) == 0
  iJorC = +1;
end
if iJorC == +1
  if iAorC < 0
    figure(8); plot_ERA_ERA5_AIRSL3_CMIP6_trends
  else
    figure(8); plot_ERA_ERA5_AIRSL3_AMIP6_trends
    cmip6 = amip6; clear amip6; disp('copying amip6 = cmip6 .. so remember all your plots with title CMIP6 should really be AMIP6')
    airsL3 = climcapsL3; clear climcapsL3; disp('copying climcaps --> L3 .. so remember all your plots with title AIRSL3 should really be CLIMCAPSL3')
  end
elseif iJorC == -1
  if iAorC < 0
    figure(8); plot_ERA_ERA5_AIRSCLIMCAPS_CMIP6_trends
  else
    figure(8); plot_ERA_ERA5_AIRSCLIMCAPS_AMIP6_trends
    cmip6 = amip6;       clear amip6;      disp('copying amip6 --> cmip6 .. so remember all your plots with title CMIP6  should really be AMIP6')
    airsL3 = climcapsL3; clear climcapsL3; disp('copying climcaps --> L3 .. so remember all your plots with title AIRSL3 should really be CLIMCAPSL3')
  end
end
