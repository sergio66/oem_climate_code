%% superceded driver_get_the_model_trends_orig
iJorC = input('(-1) Chris Barnet NUCAPS or (+1, default) Joel Susskind AIRS L2 : ');
if length(iJorC) == 0
  iJorC = +1;
end
if iJorC == -1
  iJorCstr = 'CLIMCAPS';
  disp('remember all your plots with title AIRSL3 should really be CLIMCAPSL3')
else
  iJorCstr = 'AIRSL3';
end

iEorM = input('(1) ERA-I (2) MERRA2 or (5,default) ERA5 : ');
if length(iEorM) == 0
  iEorM = +5;
end
if iEorM == 1
  iEorMstr = 'ERA-I';
  disp('remember all your plots with title ERA5 should really be ERA-I')
elseif iEorM == 2
  iEorMstr = 'MERRA2';
  disp('remember all your plots with title ERA5 should really be MERRA2')
elseif iEorM == 5
  iEorMstr = 'ERA5';
end

iAorC = input('(-1,default) CMIP6  or (+1) AMIP6 : ');
if length(iAorC) == 0
  iAorC = -1;
end
if iAorC == +1
  iAorCstr = 'AMIP6';
  disp('remember all your plots with title CMIP6 should really be AMIP6');
elseif iAorC == -1
  iAorCstr = 'CMIP6';
end

strMODELS = [iJorCstr '_' iEorMstr '_' iAorCstr];
[airsL3,era5,cmip6] = driverchoose_AIRSvsNWPvsXMIP6(iJorC,iEorM,iAorC);

