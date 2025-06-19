function [] = make_table3(titlestr,newz11,newz12,newz13,newz21,newz22,newz23);

%% this is Table 3 in trend JGR paper
%% see driver_compare_trends_Day_vs_Night.m --> compare_SKT_trends_Day_vs_Night_stats.m --> show_skt_trends_6models_new
addpath /asl/matlib/science/

load llsmap5

do_XX_YY_from_X_Y
[salti,landfrac] =  usgs_deg10_dem(YY,XX);

mu = cos(YY*pi/180);
mu6  =  ones(6,1) * mu;
mu6 = mu6';

newz11 = newz11(:);  %% umbc
newz12 = newz12(:);  %% airsv7
newz13 = newz13(:);  %% climcaps
newz21 = newz21(:);  %% giss
newz22 = newz22(:);  %% era5
newz23 = newz23(:);  %% merra2

figure(30); clf; scatter_coast(XX,YY,100,newz11); title(titlestr); colormap(llsmap5); caxis([-1 +1]*0.151);
iFig = 31;
  figure(iFig); sizefig; ; clf;
  clear plotoptions
  plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
  plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = ['dST/dt : ' titlestr]; plotoptions.cmap = llsmap5;
  plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
  plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
  plotoptions.barstr = ['dSKT/dt [K/yr] ' titlestr];
  plotoptions.smooth = +1; 
  plotoptions.smooth = -1; 
  figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);

tropics = find(abs(YY) <= 30);
midlats = find(abs(YY) <= 60 & abs(YY) > 30);
polar   = find(abs(YY) > 60);
ocean   = find(landfrac == 0);
land    = find(landfrac > 0);

iNatural = -1;
if iNatural > 0
  disp('SKT K yr-1  | AIRS_RT  AIRSv7 CLIMCAPS    GISS      ERA5    MERRA2');  %% this is natural order 11,12,13,21,22,23 of the aslmap_2x3tiledlayout
else
  disp('SKT K yr-1  | AIRS_RT  AIRSv7 CLIMCAPS    ERA5    MERRA2    GISS');    %% this is the table 3 order
end
disp('------------|---------------------------------------------------')
boo = 1:4608;  junk = [newz11(boo) newz12(boo) newz13(boo) newz21(boo) newz22(boo) newz23(boo)]; 
  junk = junk.* mu6(boo,:);
  denom = mu6(boo,:);
  junk = nansum(junk,1)./nansum(denom,1); 
  if iNatural < 0; junk = junk([1 2 3 5 6 4]); end
  fprintf(1,' all        | %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   %6.3f  \n',junk);
disp('------------|---------------------------------------------------')
boo = tropics; junk = [newz11(boo) newz12(boo) newz13(boo) newz21(boo) newz22(boo) newz23(boo)]; 
  junk = junk.* mu6(boo,:);
  denom = mu6(boo,:);
  junk = nansum(junk,1)./nansum(denom,1); 
  if iNatural < 0; junk = junk([1 2 3 5 6 4]); end
  fprintf(1,' tropics    | %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   %6.3f  \n',junk);
boo = midlats; junk = [newz11(boo) newz12(boo) newz13(boo) newz21(boo) newz22(boo) newz23(boo)]; 
  junk = junk.* mu6(boo,:);
  denom = mu6(boo,:);
  junk = nansum(junk,1)./nansum(denom,1); 
  if iNatural < 0; junk = junk([1 2 3 5 6 4]); end
  fprintf(1,' midlats    | %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   %6.3f  \n',junk);
boo = polar;   junk = [newz11(boo) newz12(boo) newz13(boo) newz21(boo) newz22(boo) newz23(boo)]; 
  junk = junk.* mu6(boo,:);
  denom = mu6(boo,:);
  junk = nansum(junk,1)./nansum(denom,1); 
  if iNatural < 0; junk = junk([1 2 3 5 6 4]); end
  fprintf(1,' polar      | %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   %6.3f  \n',junk);
disp('------------|---------------------------------------------------')
boo = ocean;   junk = [newz11(boo) newz12(boo) newz13(boo) newz21(boo) newz22(boo) newz23(boo)]; 
  junk = junk.* mu6(boo,:);
  denom = mu6(boo,:);
  junk = nansum(junk,1)./nansum(denom,1); 
  if iNatural < 0; junk = junk([1 2 3 5 6 4]); end
  fprintf(1,' ocean      | %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   %6.3f  \n',junk);
boo = land;    junk = [newz11(boo) newz12(boo) newz13(boo) newz21(boo) newz22(boo) newz23(boo)]; 
  junk = junk.* mu6(boo,:);
  denom = mu6(boo,:);
  junk = nansum(junk,1)./nansum(denom,1); 
  if iNatural < 0; junk = junk([1 2 3 5 6 4]); end
  fprintf(1,' land       | %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   %6.3f  \n',junk);
disp('------------|---------------------------------------------------')


disp(' ')

  %% correlations against UMBC, except sixth (last) is against ERA5
  [r,chisqr,P] = nanlinearcorrelation(newz11,newz22);    thecorr.ST(1) = r;
  [r,chisqr,P] = nanlinearcorrelation(newz11,newz23);    thecorr.ST(2) = r;
  [r,chisqr,P] = nanlinearcorrelation(newz11,newz12);    thecorr.ST(3) = r;
  [r,chisqr,P] = nanlinearcorrelation(newz11,newz13);    thecorr.ST(4) = r;
  [r,chisqr,P] = nanlinearcorrelation(newz11,newz21);    thecorr.ST(5) = r;
  [r,chisqr,P] = nanlinearcorrelation(newz22,newz21);    thecorr.ST_ERA5_GISS = r;
  printarray(thecorr.ST,'correlations of UMBC with ERA5/MERRA2/AIRS/CLIMCAPS/GISS')

disp(' ')
