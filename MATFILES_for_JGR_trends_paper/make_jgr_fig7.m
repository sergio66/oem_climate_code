load fig7

load llsmap5

  iFig = 1;
    figure(iFig); clf;
    clear plotoptions
    plotoptions.Xstr = ' '; plotoptions.Ystr = ' ';
    plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.cmap = llsmap5;
    plotoptions.str11 = 'AIRS\_RT';     plotoptions.str12 = 'AIRS L3';     plotoptions.str13 = 'CLIMCAPS L3';
    plotoptions.str21 = 'GISS';         plotoptions.str22 = 'ERA5';        plotoptions.str23 = 'MERRA2';
    plotoptions.barstr = 'dSKT/dt [K/yr]';
    plotoptions.smooth = -1;
    %figure(iFig); sizefig; ; clf; aslmap_2x3tiledlayout(newz11,newz12,newz13,newz21,newz22,newz23,iFig,plotoptions);
    figure(iFig); clf; aslmap_2x3tiledlayout(data7.umbc,data7.airsv7,data7.climcaps,data7.giss,data7.era5,data7.merra2,iFig,plotoptions);
