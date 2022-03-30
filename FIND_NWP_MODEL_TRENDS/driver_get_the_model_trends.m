%% superceded driver_get_the_model_trends_orig
iJorC = input('(-1) Chris Barnet NUCAPS or (+1, default) Joel Susskind AIRS L2 : ');
if length(iJorC) == 0
  iJorC = +1;
end
if iJorC == -1
  disp('remember all your plots with title AIRSL3 should really be CLIMCAPSL3')
end

iEorM = input('(1) ERA-I (2) MERRA2 or (5,default) ERA5 : ');
if length(iEorM) == 0
  iEorM = +5;
end
if iEorM == 1
  disp('remember all your plots with title ERA5 should really be ERA-I')
elseif iEorM == 2
  disp('remember all your plots with title ERA5 should really be MERRA2')
end

iAorC = input('(-1,default) CMIP6  or (+1) AMIP6 : ');
if length(iAorC) == 0
  iAorC = -1;
end
if iAorC == +1
  disp('remember all your plots with title CMIP6 should really be AMIP6');
end

[airsL3,era5,cmip6] = driverchoose_AIRSvsNWPvsXMIP6(iJorC,iEorM,iAorC);

