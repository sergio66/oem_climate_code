% disp('Enter (-7) polar    L/O')
% disp('      (-6) midlat   L/O')
% disp('      (-5) tropical L/O')
% disp('      (-4) polar land    (+4) polar ocean')
% disp('      (-3) midlat land   (+3) midlat ocean')
% disp('      (-2) tropical land (+2) tropical ocean')
% disp('      (-1) land          (+1) ocean');
% disp('      [0,default] ALL trends : ');

if ixAorOorL == +1
  %% ocean only
  ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac < 0.001);
elseif ixAorOorL == -1
  %% land only
  ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1) & landfrac >= 0.001);
elseif ixAorOorL == 0
  %% everything
  ix = find(Airs_Lat >= rlat(jj) & Airs_Lat < rlat(jj+1) & Airs_Lon >= rlon(ii) & Airs_Lon < rlon(ii+1));

elseif ixAorOorL == -7
  %% polar everything
  ix = find(abs(Airs_Lat) >= 60);
elseif ixAorOorL == -6
  %% midlat everything
  ix = find(abs(Airs_Lat) < 60 & abs(Airs_Lat) >= 30);
elseif ixAorOorL == -5
  %% tropical everything
  ix = find(abs(Airs_Lat) < 30);

elseif ixAorOorL == -4
  %% polar land
  ix = find(abs(Airs_Lat) >= 60 & landfrac >= 0.001);
elseif ixAorOorL == -3
  %% midlat land
  ix = find(abs(Airs_Lat) < 60 & abs(Airs_Lat) >= 30 & landfrac >= 0.001);
elseif ixAorOorL == -2
  %% tropical land
  ix = find(abs(Airs_Lat) < 30 & landfrac >= 0.001);

elseif ixAorOorL == +4
  %% polar ocean
  ix = find(abs(Airs_Lat) >= 60 & landfrac < 0.001);
elseif ixAorOorL == +3
  %% midlat ocean
  ix = find(abs(Airs_Lat) < 60 & abs(Airs_Lat) >= 30 & landfrac < 0.001);
elseif ixAorOorL == +2
  %% tropical ocean
  ix = find(abs(Airs_Lat) < 30 & landfrac < 0.001);


end
