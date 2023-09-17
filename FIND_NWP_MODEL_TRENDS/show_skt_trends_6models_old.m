% z11 = AIRS,       z12 = climcaps      z21 = MERRA2        z22 = ERA5            z11x/z31 = umbc      z32 = giss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('SKT trend K/yr  | AIRS   CLIMCAPS   MERRA2      ERA5    UMBC     GISS');
disp('----------------|-------------------------------------------------------');
junk = [nanmean(z11)  nanmean(z12) nanmean(z21)  nanmean(z22) nanmean(z31)  nanmean(z32)];
fprintf(1,'ALL             |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac == 0);  
junk = [nanmean(z11(boo))  nanmean(z12(boo)) nanmean(z21(boo))  nanmean(z22(boo)) nanmean(z31(boo))  nanmean(z32(boo))];
fprintf(1,'OCEAN           |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac > 0.99);  
junk = [nanmean(z11(boo))  nanmean(z12(boo)) nanmean(z21(boo))  nanmean(z22(boo)) nanmean(z31(boo))  nanmean(z32(boo))];
fprintf(1,'LAND            |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
disp('----------------|-------------------------------------------------------');

disp(' ')
disp('SKT trend K/yr  | AIRS   CLIMCAPS   MERRA2      ERA5    UMBC     GISS');
disp('area weighted with cos(lat)')
disp('----------------|--------------------------------------------------------');

zoop = 1 : length(YY);
junk = [nansum(z11(zoop).*mu(zoop))/nansum(mu(zoop))  nansum(z12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(z21(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(z22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(z31(zoop).*mu(zoop))/nansum(mu(zoop))  nansum(z32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'ALL             |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
zoop = tropics;
junk = [nansum(z11(zoop).*mu(zoop))/nansum(mu(zoop))  nansum(z12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(z21(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(z22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(z31(zoop).*mu(zoop))/nansum(mu(zoop))  nansum(z32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'TROPICS         |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
zoop = midlats;
junk = [nansum(z11(zoop).*mu(zoop))/nansum(mu(zoop))  nansum(z12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(z21(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(z22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(z31(zoop).*mu(zoop))/nansum(mu(zoop))  nansum(z32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'MIDLATS         |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
zoop = poles;
junk = [nansum(z11(zoop).*mu(zoop))/nansum(mu(zoop))  nansum(z12(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(z21(zoop).*mu(zoop))/nansum(mu(zoop))  ...
        nansum(z22(zoop).*mu(zoop))/nansum(mu(zoop)) nansum(z31(zoop).*mu(zoop))/nansum(mu(zoop))  nansum(z32(zoop).*mu(zoop))/nansum(mu(zoop))];
fprintf(1,'POLAR           |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);

boo = find(landfrac == 0);  
junk = [nansum(z11(boo).*mu(boo))/nansum(mu(boo))  nansum(z12(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(z21(boo).*mu(boo))/nansum(mu(boo))  nansum(z22(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(z31(boo).*mu(boo))/nansum(mu(boo))  nansum(z32(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'OCEAN           |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac > 0.99);  
junk = [nansum(z11(boo).*mu(boo))/nansum(mu(boo))  nansum(z12(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(z21(boo).*mu(boo))/nansum(mu(boo))  nansum(z22(boo).*mu(boo))/nansum(mu(boo))...
        nansum(z31(boo).*mu(boo))/nansum(mu(boo))  nansum(z32(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'LAND            |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
disp('----------------|--------------------------------------------------------');

xz11 = z11 + z11unc;
xz12 = z12 + z12unc;
xz21 = z21 + z21unc;
xz22 = z22 + z22unc;
xz31 = z11x + z11xunc;
xz32 = z32 + z32unc;

disp(' ')
disp('SKT trend K/yr  | AIRS   CLIMCAPS   MERRA2      ERA5    UMBC     GISS');
disp('SKT with unc    |     ');
disp('area weighted with cos(lat)')
disp('----------------|--------------------------------------------------------');

xzoop = 1 : length(YY);
junk = [nansum(xz11(xzoop).*mu(xzoop))/nansum(mu(xzoop))  nansum(xz12(xzoop).*mu(xzoop))/nansum(mu(xzoop)) nansum(xz21(xzoop).*mu(xzoop))/nansum(mu(xzoop))  ...
        nansum(xz22(xzoop).*mu(xzoop))/nansum(mu(xzoop)) nansum(xz31(xzoop).*mu(xzoop))/nansum(mu(xzoop))  nansum(xz32(xzoop).*mu(xzoop))/nansum(mu(xzoop))];
fprintf(1,'ALL             |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
xzoop = tropics;
junk = [nansum(xz11(xzoop).*mu(xzoop))/nansum(mu(xzoop))  nansum(xz12(xzoop).*mu(xzoop))/nansum(mu(xzoop)) nansum(xz21(xzoop).*mu(xzoop))/nansum(mu(xzoop))  ...
        nansum(xz22(xzoop).*mu(xzoop))/nansum(mu(xzoop)) nansum(xz31(xzoop).*mu(xzoop))/nansum(mu(xzoop))  nansum(xz32(xzoop).*mu(xzoop))/nansum(mu(xzoop))];
fprintf(1,'TROPICS         |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
xzoop = midlats;
junk = [nansum(xz11(xzoop).*mu(xzoop))/nansum(mu(xzoop))  nansum(xz12(xzoop).*mu(xzoop))/nansum(mu(xzoop)) nansum(xz21(xzoop).*mu(xzoop))/nansum(mu(xzoop))  ...
        nansum(xz22(xzoop).*mu(xzoop))/nansum(mu(xzoop)) nansum(xz31(xzoop).*mu(xzoop))/nansum(mu(xzoop))  nansum(xz32(xzoop).*mu(xzoop))/nansum(mu(xzoop))];
fprintf(1,'MIDLATS         |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
xzoop = poles;
junk = [nansum(xz11(xzoop).*mu(xzoop))/nansum(mu(xzoop))  nansum(xz12(xzoop).*mu(xzoop))/nansum(mu(xzoop)) nansum(xz21(xzoop).*mu(xzoop))/nansum(mu(xzoop))  ...
        nansum(xz22(xzoop).*mu(xzoop))/nansum(mu(xzoop)) nansum(xz31(xzoop).*mu(xzoop))/nansum(mu(xzoop))  nansum(xz32(xzoop).*mu(xzoop))/nansum(mu(xzoop))];
fprintf(1,'POLAR           |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);

boo = find(landfrac == 0);  
junk = [nansum(xz11(boo).*mu(boo))/nansum(mu(boo))  nansum(xz12(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(xz21(boo).*mu(boo))/nansum(mu(boo))  nansum(xz22(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(xz31(boo).*mu(boo))/nansum(mu(boo))  nansum(xz32(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'OCEAN           |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
boo = find(landfrac > 0.99);  
junk = [nansum(xz11(boo).*mu(boo))/nansum(mu(boo))  nansum(xz12(boo).*mu(boo))/nansum(mu(boo)) ...
        nansum(xz21(boo).*mu(boo))/nansum(mu(boo))  nansum(xz22(boo).*mu(boo))/nansum(mu(boo))...
        nansum(xz31(boo).*mu(boo))/nansum(mu(boo))  nansum(xz32(boo).*mu(boo))/nansum(mu(boo))];
fprintf(1,'LAND            |  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',junk);
disp('----------------|--------------------------------------------------------');
