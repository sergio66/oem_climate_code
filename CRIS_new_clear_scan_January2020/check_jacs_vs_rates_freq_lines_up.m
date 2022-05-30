addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/CRIS_Hi2Lo/
addpath /umbc/xfs2/strow/asl/packages/ccast/motmsc/time/
addpath /umbc/xfs2/strow/asl/packages/ccast/motmsc/time/

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CLO_Anomaly137_16_12p8/RESULTS/kcarta_085_M_TS_jac_all_5_97_97_97_2235.mat
[fout,jacout] = translate_hi2lo(f,squeeze(M_TS_jac_all(20,:,:)));
i1305 = load('/asl/matlib/cris/ch_std_from1317.mat');
f   = fout.vchan(i1305.ch_std_i);
jac = jacout(i1305.ch_std_i,:);

iLatbin = 20;
iSergioORStrow = +1;
iSergioORStrow = -1;

iSergioORStrow = input('read Sergio(+1/default) or Strow (-1) anomalies : ');
if length(iSergioORStrow) == 0
  iSergioORStrow = +1;
end

if iSergioORStrow == 1
  iiMax = 156; sf = 1;
else
  iiMax = 157; sf = 1000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iSergioORStrow == -1
  fname = ['/home/strow/Work/Cris/Stability/Data/Desc_fits/fit_robs_lat' num2str(iLatbin) '.mat'];
else
  fname = ['ANOM_16dayavg_2012_2019/Daily_Anomalies_2012/sergio_robs_latbin' num2str(iLatbin) '.mat'];
end

strow = load(fname);
plot(f,jac(:,1),f,strow.dbt*sf); xlim([640 840]);                        title(['latbin ' num2str(iLatbin)]); hl = legend('co2 jac','rate dBT/dt','location','best','fontsize',10); rett
plot(f,sum(jac(:,7:100)'),f,strow.dbt*sf); xlim([840 1640]);             title(['latbin ' num2str(iLatbin)]); hl = legend('wv jac','rate dBT/dt','location','best','fontsize',10); rett
plot(f,strow.all_bt_anom(1,:)*sf,f,strow.dbt*sf); xlim([640 840]);       title(['latbin ' num2str(iLatbin)]); hl = legend('anom(t=1,2600 days)','rate dBT/dt','location','best','fontsize',10); rett
plot(f,nanmean(strow.all_bt_anom,1)*sf,f,strow.dbt*sf); xlim([640 840]); title(['latbin ' num2str(iLatbin)]); hl = legend('mean(anom,2600 days)','rate dBT/dt','location','best','fontsize',10); rett

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%oo0 = load('ANOM_16dayavg/latbin_0dayavg_20.mat');
%if iSergioORStrow == -1
%  oo1 = load('ANOM_16dayavgDEBUG/latbin_0dayavg_20.mat');
%else
%  oo1 = load('ANOM_16dayavgDEBUG/sergio_latbin_0dayavg_20.mat');
%end
%plot(std(oo0.avg16_btanom,1)-std(oo1.avg16_btanom,1));   title('156 16 day timesteps : std anom orig-latest'); rett
%plot(mean(oo0.avg16_btanom,1)-mean(oo1.avg16_btanom,1)); title('156 16 day timesteps : mean anom orig-latest'); rett

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

fnamex = ['ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '.mat'];
if iSergioORStrow == -1
  fnamex = ['ANOM_16dayavgDEBUG/latbin_0dayavg_' num2str(iLatbin) '.mat'];
else
  fnamex = ['ANOM_16dayavgDEBUG/sergio_latbin_0dayavg_' num2str(iLatbin) '.mat'];
end
sergio = load(fnamex);
plot(f,nanmean(strow.all_bt_anom,1)*sf,'b.-',f,strow.dbt*sf,'g',f,nanmean(sergio.avg16_btanom,1),'r'); xlim([640 840]); ylim([-0.3 +0.3]); hl = legend('mean 2600 days anom','rate dBT/dt','mean 156 16days anom','location','best','fontsize',10); rett
plot(f,nanmean(strow.all_bt_anom,1)*sf - nanmean(sergio.avg16_btanom,1),'r'); xlim([640 1240]); ylim([-1 +1]*1/10); plotaxis2; title('anom : nanmean(2600 days)-nanmean(16 days)'); rett

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : iiMax;
  boo = (1:16) + (ii-1)*16;
  strow.all_bt_anom_16(ii,:) = nanmean(strow.all_bt_anom(boo,:),1)*sf;
end
plot(f,nanmean(strow.all_bt_anom,1)*sf - nanmean(sergio.avg16_btanom,1),'r'); xlim([640 1640]); ylim([-0.1 +0.1]/10); plotaxis2; title('anom : nanmean(2600 days)-nanmean(16 days)'); rett
plot(f,nanmean(strow.all_bt_anom_16,1) - nanmean(sergio.avg16_btanom,1),'r'); xlim([640 1640]); ylim([-0.1 +0.1]/10); plotaxis2; title('anom : nanmean(2600 days --> 157 16days)-nanmean(16 days)'); rett

plot(f,nanmean(strow.all_bt_anom,1)*sf - nanmean(sergio.avg16_btanom,1),'r'); xlim([640 1640]); ylim([-0.1 +0.1]/10)
hold on
plot(f,nanmean(strow.all_bt_anom_16,1) - nanmean(sergio.avg16_btanom,1),'g'); xlim([640 1640]); ylim([-0.1 +0.1]/10); plotaxis2;
hold off
hl = legend('nanmean(2600 days)-nanmean(16 days)','nanmean(2600 days --> 157 16days)-nanmean(16 days)','location','best','fontsize',10); rett

plot(f,strow.all_bt_anom_16(1:iiMax,:)-sergio.avg16_btanom(1:iiMax,:),'c'); xlim([640 1640]); ylim([-1 +1]*4); plotaxis2;
title('(2600 days --> 157 16days)-157 days'); rett

%%% plot(f,strow.all_bt_anom_16(1,:),'b.-',f,sergio.avg16_btanom(1,:),'r'); xlim([640 1640]);  plotaxis2; ylim([-1 +1]*0.5); title('timestep 1'); hl = legend('2600 days --> 157 16days','157 16 days','location','best','fontsize',10); rett

[yyx,mmx,ddx,hhx] = tai2utcSergio(dtime2tai(strow.all_times(:,1)));
[yyx,mmx,ddx,hhx] = tai2utcSergio(dtime2tai(nanmean(strow.all_times,2)));
days_since_2012x = change2days(yyx,mmx,ddx,2012);

for ii = 1 : 10
  figure(1); clf; plot(f,strow.all_bt_anom_16(ii,:),'b.-',f,sergio.avg16_btanom(ii,:),'r'); xlim([640 1640]);  plotaxis2; ylim([-1 +1]*0.5);   title(['loopy loop ' num2str(ii)])
  figure(2); clf; plot(f,strow.all_bt_anom_16(ii,:)-sergio.avg16_btanom(ii,:),'r'); xlim([640 1640]);  plotaxis2;   title(['loopy loop ' num2str(ii)])
  rett
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is just to check times, no need to do
return
  %% see clust_convert_cris_strowrates2oemrates_anomalyDEBUG.m
  iPlot = +1;

  a = strow;
  [mmlen,nnlen] = size(a.all_times);
  one2mmlen = 1:mmlen;

  [rowi,coli] = find(isinf(a.all_bt_anom));
  [rown,coln] = find(isnan(a.all_bt_anom));
  if (length(coli) > 0) | (length(coln) > 0) | (length(rowi) > 0) | (length(rown) > 0)
    disp('found bad all_bt_anom')
    whos rowi coli rown coln
  end

  [rowi,coli] = find(isinf(a.all_times));
  [rown,coln] = find(isnan(a.all_times));
  if (length(coli) > 0) | (length(coln) > 0) | (length(rowi) > 0) | (length(rown) > 0)
    disp('found bad all_times')
    whos rowi coli rown coln
  end

  rad400 = a.all_bt_anom(:,400);
  nan400 = find(isnan(rad400));
  inf400 = find(isinf(rad400));
  if length(nan400) > 0 | length(inf400) > 0
    disp('oh oh found bad rads for ch400 : ')
    inf400
    nan400
    %error('bad rad')
  end

  [rowi,coli] = find(isinf(a.all_times(:,1:1300)));
  [rown,coln] = find(isnan(a.all_times(:,1:1300)));

  [rowi,coli] = find(isinf(a.all_times(:,100)));
  [rown,coln] = find(isnan(a.all_times(:,100)));

  badrow = union(rowi,rown);
  goodtimes = setdiff(1:mmlen,badrow);
  meantime = nanmean(a.all_times(goodtimes,:),2);
  stdtime  = nanstd(a.all_times(goodtimes,:)');

  if iPlot > 0
    figure(1)
    errorbar(goodtimes,meantime,stdtime);
  end

  %% figure out start/stop UTC dates
  rtime = dtime2tai(meantime);
  [yy,mm,dd,hh] = tai2utcSergio(rtime);

  %% figure out full time data
  days_since_2012 = change2days(yy,mm,dd,2012);
  [C,IA,IC] = unique(days_since_2012);
  if length(C) < length(days_since_2012)
    disp('whoops : some days are not unique, fixing!')
    days_since_2012 = days_since_2012(IA);
    rtime     = rtime(IA);
    goodtimes = goodtimes(IA);
    yy = yy(IA);
    mm = mm(IA);
    dd = dd(IA);
    hh = hh(IA);
  end

  iStart = find(yy == 2012 & mm == 5 & dd >= 01); iStart = iStart(1);
  if iSergioORStrow == -1
    %iEnd   = find(yy == 2018 & mm == 4 & dd <= 30); iEnd = iEnd(end);
    %iEnd   = find(yy == 2019 & mm == 4 & dd <= 30); iEnd = iEnd(end);
    iEnd   = find(yy == 2019 & mm == 7 & dd <= 30); iEnd = iEnd(end);
  else
    iEnd   = find(yy == 2019 & mm == 3 & dd <= 30); iEnd = iEnd(end);
  end

  %% this is for iLatbin 1,2
  if length(iStart) == 0
    iStart = 1;
  end
  if length(iEnd) == 0
    iEnd = length(yy);
  end

  dStart = days_since_2012(iStart);
  dEnd   = days_since_2012(iEnd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
