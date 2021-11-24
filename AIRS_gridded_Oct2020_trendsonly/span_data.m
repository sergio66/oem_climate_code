[goodtimes,~] = size(a.all_times);
goodtimes = 1 : goodtimes;
meantime = mean(a.all_times(goodtimes,:),2);

rtime = dtime2tai(meantime);
[yy,mm,dd,hh] = tai2utcSergio(rtime);

iStart = find(yy == yyS & mm == mmS & dd == ddS);
iEnd   = find(yy == yyE & mm == mmE & dd == ddE);

iaList = iStart:iEnd;

a.all_bt_anom = a.all_bt_anom(iaList,:);
a.all_bt_resid = a.all_bt_resid(iaList,:);
a.all_times = a.all_times(iaList,:);

clear goodtime meantime rtime yy mm dd hh iStart iEnd iaList

