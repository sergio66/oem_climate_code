%% this is AIRS

clear iCnt obs cal fairs fcris
iCnt = 0;
for iLatbin = 11:30
  iCnt = iCnt + 1;
  fnamex = ['../AIRS_new_clear_scan_August2019_AMT2020PAPER/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '.mat'];
  fnamey = ['../AIRS_new_clear_scan_August2019_AMT2020PAPER/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '_cal.mat'];
  fnamez = ['../AIRS_new_clear_scan_August2019_AMT2020PAPER/ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_2012/sergio_latbin' num2str(iLatbin) '.mat'];

  sergiox = load(fnamex);
  sergioy = load(fnamey);
  obs(iCnt,:,:) = sergiox.avg16_btanom;
  cal(iCnt,:,:) = sergioy.avg16_btanom;
end

fairs = instr_chans2645;
  iC1231 = find(fairs >= 1231,1);
  iC791 = find(fairs >= 791,1);
  iC792 = find(fairs >= 792,1);
  iC791 = find(fairs >= 791-0.5,1);
  iC792 = find(fairs >= 792-0.5,1);
  iC730 = find(fairs >= 730,1);

  iC1300 = find(fairs >= 1300,1);
  iC1310 = find(fairs >= 1310,1);

  iC1419 = find(fairs >= 1419,1);
  iC1619 = find(fairs >= 1619,1);

rtime = sergiox.avg16_rtime;
[yy,mm,dd,hh] = tai2utcSergio(rtime);
days_since_2012 = change2days(yy,mm,dd,2012);
tdays = yy + (mm-1)/12 + (dd-1)/30/12;

%obsmean = rad2bt(fairs,squeeze(nanmean(obs,1))');
%calmean = rad2bt(fairs,squeeze(nanmean(cal,1))');
obsmean = squeeze(nanmean(obs,1))';
calmean = squeeze(nanmean(cal,1))';

figure(1)
plot(tdays,obsmean(iC1231,:)-calmean(iC1231,:)); title('AIRS L1C rad 1231 obs-cal')
plot(tdays,obsmean(iC730,:)-calmean(iC730,:)); title('AIRS L1C rad 730 obs-cal')
plot(tdays,(obsmean(iC791,:)-obsmean(iC792,:)) - (calmean(iC791,:)-calmean(iC792,:))); title('AIRS L1C rad (791-792) obs-cal')

disp('ret to continue'); pause
for ii = 1 : iC792 + 30
  figure(1); plot(tdays,obsmean(ii,:),tdays,calmean(ii,:)); title(['AIRS L1C rad ' num2str(fairs(ii)) '  (b) obs (r) cal'])
  figure(2); plot(tdays,obsmean(ii,:)-calmean(ii,:)); title(['AIRS L1C rad ' num2str(fairs(ii)) '  obs-cal'])
  disp('ret to continue'); pause
  pause(0.1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this is CRIS
iSergioORStrow = input('read Sergio(+1/default) or Strow (-1) anomalies : ');
if length(iSergioORStrow) == 0
  iSergioORStrow = +1;
end

if iSergioORStrow == 1
  iiMax = 156; sf = 1;
else
  iiMax = 157; sf = 1000;
end

clear iCnt obs cal fairs fcris
iCnt = 0;
for iLatbin = 11:30
  iCnt = iCnt + 1;
  if iSergioORStrow == -1
    %% this is me mucking around
    fnamex = ['ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '.mat'];
    fnamey = ['ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '_cal.mat'];
  else
    fnamex = ['ANOM_16dayavgDEBUG/sergio_latbin_0dayavg_' num2str(iLatbin) '.mat'];
    fnamey = ['ANOM_16dayavgDEBUG/sergio_latbin_0dayavg_' num2str(iLatbin) '_cal.mat'];
  end
  sergiox = load(fnamex);
  sergioy = load(fnamey);
  obs(iCnt,:,:) = sergiox.avg16_btanom;
  cal(iCnt,:,:) = sergioy.avg16_btanom;
end

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/cris_ichan_vchan_nsr_fsr.mat
usethese = find(nsr.ichan <= 1305);
fcris = nsr.vchan(usethese);
fairs = fcris;
  iC1231 = find(fairs >= 1231,1);
  iC791 = find(fairs >= 791,1);
  iC792 = find(fairs >= 792,1);
  iC791 = find(fairs >= 791-0.5,1);
  iC792 = find(fairs >= 792-0.5,1);
  iC730 = find(fairs >= 730,1);

  iC1300 = find(fairs >= 1300,1);
  iC1310 = find(fairs >= 1310,1);

  iC1419 = find(fairs >= 1419,1);
  iC1619 = find(fairs >= 1619,1);

rtime = sergiox.avg16_rtime;
[yy,mm,dd,hh] = tai2utcSergio(rtime);
days_since_2012 = change2days(yy,mm,dd,2012);
tdays = yy + (mm-1)/12 + (dd-1)/30/12;

%obsmean = rad2bt(fcris,squeeze(nanmean(obs,1))');
%calmean = rad2bt(fcris,squeeze(nanmean(cal,1))');
obsmean = squeeze(nanmean(obs,1))';
calmean = squeeze(nanmean(cal,1))';

figure(1)
plot(tdays,obsmean(iC1231,:)-calmean(iC1231,:)); title('CRIS NSR rad 1231 obs-cal')
plot(tdays,obsmean(iC730,:)-calmean(iC730,:)); title('CRIS NSR rad 730 obs-cal')
plot(tdays,(obsmean(iC791,:)-obsmean(iC792,:)) - (calmean(iC791,:)-calmean(iC792,:))); title('CRIS NSR rad (791-792) obs-cal')

disp('ret to continue'); pause
for ii = 1 : iC792 + 30
  figure(1); plot(tdays,obsmean(ii,:),tdays,calmean(ii,:)); title(['CRIS NSR rad ' num2str(fcris(ii)) '  (b) obs (r) cal'])
  figure(2); plot(tdays,obsmean(ii,:)-calmean(ii,:)); title(['CRIS NSR rad ' num2str(fcris(ii)) '  obs-cal'])
  disp('ret to continue'); pause
  pause(0.1)
end

for ii = iC1300:iC1310
  figure(1); plot(tdays,obsmean(ii,:),tdays,calmean(ii,:)); title(['CRIS NSR rad ' num2str(fcris(ii)) '  (b) obs (r) cal'])
  figure(2); plot(tdays,obsmean(ii,:)-calmean(ii,:)); title(['CRIS NSR rad ' num2str(fcris(ii)) '  obs-cal'])
  disp('ret to continue'); pause
  pause(0.1)
end

for ii = iC1419:iC1619
  figure(1); plot(tdays,obsmean(ii,:),tdays,calmean(ii,:)); title(['CRIS NSR rad ' num2str(fcris(ii)) '  (b) obs (r) cal'])
  figure(2); plot(tdays,obsmean(ii,:)-calmean(ii,:)); title(['CRIS NSR rad ' num2str(fcris(ii)) '  obs-cal'])
  disp('ret to continue'); pause
  pause(0.1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
