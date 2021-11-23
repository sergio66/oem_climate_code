addpath  /asl/matlib/time
addpath /home/sergio/MATLABCODE/TIME

iPlot = -1;
iDebug = -1;

i_obsORcal = -1; %% cal
i_obsORcal = +1; %% obs

iAvgNumDays = 180;
iAvgNumDays = 000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iLS = 1; iLE = 40;
JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
iLS = JOB; iLE = JOB;

%%% iLS = 1; iLE = 4;

for iLatbin = iLS : iLE
  clear avg16_rtime avg_doy_since2002 avg16_btanom
  if i_obsORcal > 0
    fname = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/fit_nucal_robs_lat' num2str(iLatbin) '.mat'];
  else
    fname = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/fit_rclr_lat' num2str(iLatbin) '.mat'];
  end
  a = load(fname);

  %% figure out "bad" data
  fprintf(1,'latbin %2i fname %s \n',iLatbin,fname);

  [mmlen,nnlen] = size(a.all_times);
  one2mmlen = 1:mmlen;
  [rowi,coli] = find(isinf(a.all_times));
  [rown,coln] = find(isnan(a.all_times));
  badrow = union(rowi,rown);
  goodtimes = setdiff(1:mmlen,badrow);
  meantime = mean(a.all_times(goodtimes,:),2);
  stdtime  = std(a.all_times(goodtimes,:)');

  if iPlot > 0
    figure(1)
    errorbar(goodtimes,meantime,stdtime);
  end

  %% figure out start/stop UTC dates
  rtime = dtime2tai(meantime);
  [yy,mm,dd,hh] = tai2utcSergio(rtime);

  %% figure out full time data
  days_since_2002 = change2days(yy,mm,dd,2002);
  [C,IA,IC] = unique(days_since_2002);
  if length(C) < length(days_since_2002)
    disp('whoops : some days are not unique, fixing!')
    days_since_2002 = days_since_2002(IA);
    rtime     = rtime(IA);
    goodtimes = goodtimes(IA);
    yy = yy(IA);
    mm = mm(IA);
    dd = dd(IA);
    hh = hh(IA);
  end

  iStart = find(yy == 2002 & mm == 9 & dd == 01);
  iEnd   = find(yy == 2018 & mm == 8 & dd == 30);
  %% this is for iLatbin 1,2
  if length(iStart) == 0
    iStart = 1;
  end
  if length(iEnd) == 0
    iEnd = length(yy);
  end

  dStart = days_since_2002(iStart);
  dEnd   = days_since_2002(iEnd);
  
  %% figure out "existing days" (iA) and "missing days" (iB)
  iaAllDays = dStart : dEnd;
  rtimeAllDays  = interp1(days_since_2002,rtime,iaAllDays);
  [Y,iA,founddays]   = intersect(iaAllDays,days_since_2002);

  missingdays = setdiff(iaAllDays,days_since_2002);
  [Y,iB,notfounddays]   = intersect(iaAllDays,missingdays);

  %% loop over channels, first making "whole" time series by interpolation, then do 16 day avgs
  fprintf(1,'--> latbin %2i, smoothing channels \n',iLatbin)
  for iChan = 1 : 2645
    clear data ysmooth
    if mod(iChan,100) == 0
      fprintf(1,'.');
    end
    data = ones(size(iaAllDays)) * -99;
    data(iA) = a.all_bt_anom(goodtimes,iChan);
    data(iB) = interp1(days_since_2002,data(iA),missingdays);

    %% smooth the data
    if iAvgNumDays > 0
      ysmooth = smooth(data,iAvgNumDays);
    else
      ysmooth = data;
    end

    if (iDebug > 0 | iChan == 1520) & iPlot > 0
      figure(1); clf
      plot(iaAllDays,data)
      plot(iaAllDays,data,days_since_2002,a.all_bt_anom(goodtimes,iChan),'r.');
      plot(iaAllDays,data,iaAllDays,ysmooth,'k'); 
      title(num2str(iLatbin))
      pause(0.1);
    end

    %% then find units of 16 days
    ia16 = 1 : 16 : length(ysmooth);
    for iDay16 = 1 : length(ia16)-1
      ind = ia16(iDay16):ia16(iDay16+1);
      avg16_rtime(iDay16) = mean(rtimeAllDays(ind));
      avg_doy_since2002(iDay16) = mean(iaAllDays(ind));
      avg16_btanom(iDay16,iChan) = mean(ysmooth(ind));
    end     %% loop over 16 days
     
    if (iDebug > 0 | iChan == 1520) & iPlot > 0
      figure(1); clf
      plot(iaAllDays,data,iaAllDays,ysmooth,'k',avg_doy_since2002,avg16_btanom(:,iChan),'r')
      title(num2str(iLatbin))
      pause(0.1)
    end

    [yy16,mm16,dd16,hh16] = tai2utcSergio(avg16_rtime);
  
  end       %% loop over channels iChan
  fprintf(1,' saving\n');

  avg16_btanom = real(avg16_btanom);
  if i_obsORcal > 0
    saver = ['save ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat avg16_rtime avg_doy_since2002 avg16_btanom'];
  elseif i_obsORcal < 0
    saver = ['save ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat avg16_rtime avg_doy_since2002 avg16_btanom'];
  end

  eval(saver)

  if iPlot > 0
    figure(2); clf
    pcolor(2002+(avg_doy_since2002)/365,1:2645,avg16_btanom'); colormap jet; colorbar; shading interp;
    caxis([-2 +2]);
    xlabel('days since 1/1/2002'); ylabel('Channel ID')
    title(num2str(iLatbin))
    pause(0.1);
  end

end         %% loop over latbins  iLatbin

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
save_days = avg_doy_since2002;
save_dat_1231 = zeros(length(save_days),40);

for iLatbin = 1 : 40
  if i_obsORcal > 0
    fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
  else
    fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
  end

  load(fin)

  figure(2); clf
  pcolor(2002+(avg_doy_since2002)/365,1:2645,avg16_btanom'); colormap jet; colorbar; shading interp;
  caxis([-2 +2]);
  xlabel('days since 1/1/2002'); ylabel('Channel ID')
  title(num2str(iLatbin))

  save_dat_day(:,iLatbin) = [avg_doy_since2002(1) avg_doy_since2002(end) length(avg_doy_since2002)];

  if length(avg_doy_since2002) == length(save_days)
    save_dat_1231(:,iLatbin) = avg16_btanom(:,1520);
  else
    for dd = 1 : length(avg_doy_since2002)
      moo = abs(avg_doy_since2002(dd)-save_days);
      moo = find(moo == min(moo));
      save_dat_1231(moo,iLatbin) = avg16_btanom(dd,1520);
    end
  end
  pause(0.1);
end

pcolor(2002+save_days/365,1:40,save_dat_1231'); caxis([-1 +1]); colorbar; title('ANOM BT1231');
