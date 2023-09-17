% newz11 = ERA5, newz12 = MERRA2     newz21 = AIRSL3     newz22 = CLIMCAPS  newz11x/newz31 = umbc   newz32 = giss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('SKT trend K/yr  | ERA5    MERRA2    THISWORK    AIRS   CLIMCAPS   GISS');
disp('----------------|-------------------------------------------------------');
junk = [nanmean(newz11)  nanmean(newz12) nanmean(newz31)  nanmean(newz21) nanmean(newz22)  nanmean(newz32)];
fprintf(1,'ALL             |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac == 0);  
junk = [nanmean(newz11(boo))  nanmean(newz12(boo)) nanmean(newz31(boo))  nanmean(newz21(boo)) nanmean(newz22(boo))  nanmean(newz32(boo))];
fprintf(1,'OCEAN           |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac > 0.99);  
junk = [nanmean(newz11(boo))  nanmean(newz12(boo)) nanmean(newz31(boo))  nanmean(newz21(boo)) nanmean(newz22(boo))  nanmean(newz32(boo))];
fprintf(1,'LAND            |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
disp('----------------|-------------------------------------------------------');

disp(' ')
disp('SKT trend K/yr  | ERA5    MERRA2    THISWORK    AIRS   CLIMCAPS   GISS');
disp('area weighted with cos(lat)')
disp('----------------|--------------------------------------------------------');

zoop = 1 : length(YY);
junk = [nansum(newz11(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz31(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(newz21(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'ALL             |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
zoop = tropics;
junk = [nansum(newz11(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz31(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(newz21(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'TROPICS         |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
zoop = midlats;
junk = [nansum(newz11(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz31(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(newz21(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'MIDLATS         |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
zoop = poles;
junk = [nansum(newz11(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz31(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(newz21(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(newz32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'POLAR           |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);

boo = find(landfrac == 0);  
junk = [nansum(newz11(boo).*mu(boo))/nansum(mu(boo))  nansum(newz12(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(newz31(boo).*mu(boo))/nansum(mu(boo))  nansum(newz21(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(newz22(boo).*mu(boo))/nansum(mu(boo))  nansum(newz32(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'OCEAN           |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac > 0.99);  
junk = [nansum(newz11(boo).*mu(boo))/nansum(mu(boo))  nansum(newz12(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(newz31(boo).*mu(boo))/nansum(mu(boo))  nansum(newz21(boo).*mu(boo))/nansum(mu(boo))...
        nansum(newz22(boo).*mu(boo))/nansum(mu(boo))  nansum(newz32(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'LAND            |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
disp('----------------|--------------------------------------------------------');

xnewz11 = newz11 + newz11unc;
xnewz12 = newz12 + newz12unc;
xnewz21 = newz21 + newz21unc;
xnewz22 = newz22 + newz22unc;
xnewz31 = newz11x + newz11xunc;
xnewz32 = newz32 + newz32unc;

disp(' ')
disp('SKT trend K/yr  | ERA5    MERRA2    THISWORK    AIRS   CLIMCAPS   GISS');
disp('SKT with unc    |     ');
disp('area weighted with cos(lat)')
disp('----------------|--------------------------------------------------------');

zoop = 1 : length(YY);
junk = [nansum(xnewz11(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz31(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(xnewz21(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'ALL             |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
zoop = tropics;
junk = [nansum(xnewz11(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz31(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(xnewz21(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'TROPICS         |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
zoop = midlats;
junk = [nansum(xnewz11(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz31(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(xnewz21(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'MIDLATS         |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
zoop = poles;
junk = [nansum(xnewz11(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz31(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(xnewz21(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(xnewz32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'POLAR           |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);

boo = find(landfrac == 0);  
junk = [nansum(xnewz11(boo).*mu(boo))/nansum(mu(boo))  nansum(xnewz12(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(xnewz31(boo).*mu(boo))/nansum(mu(boo))  nansum(xnewz21(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(xnewz22(boo).*mu(boo))/nansum(mu(boo))  nansum(xnewz32(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'OCEAN           |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac > 0.99);  
junk = [nansum(xnewz11(boo).*mu(boo))/nansum(mu(boo))  nansum(xnewz12(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(xnewz31(boo).*mu(boo))/nansum(mu(boo))  nansum(xnewz21(boo).*mu(boo))/nansum(mu(boo))...
        nansum(xnewz22(boo).*mu(boo))/nansum(mu(boo))  nansum(xnewz32(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'LAND            |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
disp('----------------|--------------------------------------------------------');
