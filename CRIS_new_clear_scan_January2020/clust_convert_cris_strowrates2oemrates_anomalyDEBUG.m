addpath  /asl/matlib/time
addpath /home/sergio/MATLABCODE/TIME

iPlot = -1;
iDebug = -1;

i_obsORcal = -1; %% cal
i_obsORcal = +1; %% obs

iAvgNumDays = 180;
iAvgNumDays = 000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i1231 = 731;

iLS = 1; iLE = 40;

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));  %% latbins 1-40
%JOB = 16
%JOB = 20
%JOB = 30

iLS = JOB; iLE = JOB;

%%% iLS = 1; iLE = 4;

for iLatbin = iLS : iLE
  clear avg16_rtime avg_doy_since2012 avg16_btanom
  if i_obsORcal > 0
    fname = ['/home/strow/Work/Cris/Stability/Data/Desc_fits/fit_robs_lat' num2str(iLatbin) '.mat'];     sf = 1000;   %% I think there are problems here in strow stuff, since he forgot x1000
    fname = ['ANOM_16dayavg_2012_2019/Daily_Anomalies_2012/sergio_robs_latbin' num2str(iLatbin) '.mat']; sf = 1;
  else
    fname = ['/home/strow/Work/Cris/Stability/Data/Desc_fits/fit_rclr_lat' num2str(iLatbin) '.mat'];     sf = 1000;   %% I think there are problems here in strow stuff, since he forgot x1000
    fname = ['ANOM_16dayavg_2012_2019/Daily_Anomalies_2012/sergio_rclr_latbin' num2str(iLatbin) '.mat']; sf = 1;      
  end
  a = load(fname);

  %% figure out "bad" data
  fprintf(1,'latbin %2i fname %s \n',iLatbin,fname);

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

  %iEnd   = find(yy == 2018 & mm == 4 & dd <= 30); iEnd = iEnd(end);
  %iEnd   = find(yy == 2019 & mm == 4 & dd <= 30); iEnd = iEnd(end);
  %iEnd   = find(yy == 2019 & mm == 7 & dd <= 30); iEnd = iEnd(end);

  iEnd   = find(yy == 2019 & mm == 3 & dd <= 30); iEnd = iEnd(end);

  %% this is for iLatbin 1,2
  if length(iStart) == 0
    iStart = 1;
  end
  if length(iEnd) == 0
    iEnd = length(yy);
  end

  dStart = days_since_2012(iStart);
  dEnd   = days_since_2012(iEnd);

  ok = iStart : iEnd;
    rtime     = rtime(ok);
    goodtimes = goodtimes(ok);
    yy = yy(ok);
    mm = mm(ok);
    dd = dd(ok);
    hh = hh(ok);
    days_since_2012 = days_since_2012(ok);

  %% figure out "existing days" (iA) and "missing days" (iB)
  iaAllDays = dStart : dEnd;
  rtimeAllDays  = interp1(days_since_2012,rtime,iaAllDays);
  [Y,iA,founddays]   = intersect(iaAllDays,days_since_2012);

  missingdays = setdiff(iaAllDays,days_since_2012);
  [Y,iB,notfounddays]   = intersect(iaAllDays,missingdays);

  %% loop over channels, first making "whole" time series by interpolation, then do 16 day avgs
  fprintf(1,'--> latbin %2i, smoothing channels \n',iLatbin)
  for iChan = 1 : 1305
    clear data* ysmooth
    if mod(iChan,100) == 0
      fprintf(1,'.');
    end
    data1  = ones(size(iaAllDays)) * nan;   %% orig, since 2020
    data2 = ones(size(iaAllDays)) * nan;    %% new, May 2022

    data1(iA)  = a.all_bt_anom(goodtimes,iChan);
    data2(iA) = data1(iA);

    data1(iB)  = interp1(days_since_2012,data1(iA),missingdays);
    data2(iB) = nan;

    ooer = find(isnan(data2(iA)));
    oo1 = find(isnan(data1));
    oo2 = find(isnan(data2));

    data = data1;   %% orig
    data = data2;   %% new May 2022

    %% smooth the data
    if iAvgNumDays > 0
      ysmooth = smooth(data,iAvgNumDays);
    else
      ysmooth = data;
    end

    if (iDebug > 0 | iChan == i1231) & iPlot > 0
      figure(1); clf
      plot(iaAllDays,data)
      plot(iaAllDays,data,days_since_2012,a.all_bt_anom(goodtimes,iChan),'r.');
      plot(iaAllDays,data,iaAllDays,ysmooth,'k'); 
      title(num2str(iLatbin))
      pause(0.1);
    end

    %% then find units of 16 days
    ia16 = 1 : 16 : length(ysmooth);
    for iDay16 = 1 : length(ia16)-1
      ind = ia16(iDay16):ia16(iDay16+1);
      ind = ia16(iDay16):ia16(iDay16+1)-1;
      avg16_rtime(iDay16) = nanmean(rtimeAllDays(ind));
      avg_doy_since2012(iDay16) = nanmean(iaAllDays(ind));
      avg16_btanom(iDay16,iChan) = nanmean(ysmooth(ind)) * sf;
    end     %% loop over 16 days
     
    if (iDebug > 0 | iChan == i1231) & iPlot > 0
      figure(1); clf
      plot(iaAllDays,data,iaAllDays,ysmooth,'k',avg_doy_since2012,avg16_btanom(:,iChan),'r')
      title(num2str(iLatbin))
      pause(0.1)
    end

    [yy16,mm16,dd16,hh16] = tai2utcSergio(avg16_rtime);
  
  end       %% loop over channels iChan

  fprintf(1,' saving\n');

  avg16_btanom = real(avg16_btanom);
  %if i_obsORcal > 0
  %  avg16_btanom = 1000 * avg16_btanom;
  %end

  comment = 'see clust_convert_cris_strowrates2oemrates_anomalyDEBUG.m';

  if length(avg16_rtime) > 157
    avg16_rtime = avg16_rtime(1:157);
    avg_doy_since2012  =  avg_doy_since2012(1:157);
    avg16_btanom = avg16_btanom(1:157,:);
  end

  if i_obsORcal > 0
    fsave = ['ANOM_16dayavgDEBUG/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];        %% redoing strow stuff, carefully
    fsave = ['ANOM_16dayavgDEBUG/sergio_latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat']; %% doing sergio
    saver = ['save ' fsave ' avg16_rtime avg_doy_since2012 avg16_btanom comment'];
  elseif i_obsORcal < 0
    fsave = ['ANOM_16dayavgDEBUG/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];         %% redoing strow stuff, carefully
    fsave = ['ANOM_16dayavgDEBUG/sergio_latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];  %% doing sergio
    saver = ['save ' fsave ' avg16_rtime avg_doy_since2012 avg16_btanom comment'];
  end

  if ~exist(fsave)
    eval(saver)
  else
    fprintf(1,'%s already exists \n',fsave)
  end

  if iPlot > 0
    figure(2); clf
    pcolor(2012+(avg_doy_since2012)/365,1:1305,avg16_btanom'); colormap jet; colorbar; shading interp;
    caxis([-2 +2]);
    xlabel('days since 1/1/2012'); ylabel('Channel ID')
    title(num2str(iLatbin))

    figure(3); clf
      nd = length(days_since_2012);
      plot(days_since_2012/365+2012,1000*(a.all_bt_anom(1:nd,227)-a.all_bt_anom(1:nd,230)),'g')
      hold
      plot(2012+avg_doy_since2012/365,avg16_btanom(:,227)-avg16_btanom(:,230),'k','linewidth',2)
      hold
      title('BT 790-792')
      grid

    pause(0.1);
  end
end         %% loop over latbins  iLatbin

[mmbad,nnbad] = find(isnan(avg16_btanom) | isinf(avg16_btanom)); if length(mmbad) > 0; plot(mmbad,nnbad,'.'); end
whos mmbad nnbad avg16_btanom
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now plot all results

clear all

i_obsORcal = +1; %% obs
i_obsORcal = -1; %% cal

iLatbin = 20;
  if i_obsORcal > 0  
    fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
  else
    fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
  end
  load(fin)
save_days = avg_doy_since2012;
save_dat_1231 = zeros(length(save_days),40);

for iLatbin = 1 : 40
  if i_obsORcal > 0
    fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
  else
    fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
  end

  load(fin)

  figure(2); clf
  pcolor(2012+(avg_doy_since2012)/365,1:1305,avg16_btanom'); colormap jet; colorbar; shading interp;
  caxis([-2 +2]);
  xlabel('days since 1/1/2012'); ylabel('Channel ID')
  title(num2str(iLatbin))

  save_dat_day(:,iLatbin) = [avg_doy_since2012(1) avg_doy_since2012(end) length(avg_doy_since2012)];

  if length(avg_doy_since2012) == length(save_days)
    save_dat_1231(:,iLatbin) = avg16_btanom(:,i1231);
  else
    for dd = 1 : length(avg_doy_since2012)
      moo = abs(avg_doy_since2012(dd)-save_days);
      moo = find(moo == min(moo));
      save_dat_1231(moo,iLatbin) = avg16_btanom(dd,i1231);
    end
  end
  pause(0.1);
end

pcolor(2012+save_days/365,1:40,save_dat_1231'); caxis([-1 +1]); colorbar; title('ANOM BT1231');
