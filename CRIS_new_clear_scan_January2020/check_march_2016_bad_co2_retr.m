iLatbin = 20;

iSergioORStrow = input('read Sergio(+1/default) or Strow (-1) anomalies : ');
if length(iSergioORStrow) == 0
  iSergioORStrow = +1;
end

if iSergioORStrow == 1
  iiMax = 156; sf = 1;
else
  iiMax = 157; sf = 1000;
end

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

junk = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CLO_Anomaly137_16_12p8/RESULTS/kcarta_085_M_TS_jac_all_5_97_97_97_2235.mat','f');
[fout,jacout] = translate_hi2lo(junk.f,zeros(size(junk.f)));
i1305 = load('/asl/matlib/cris/ch_std_from1317.mat');
fcris   = fout.vchan(i1305.ch_std_i);

[yyC,mmC,ddC] = tai2utcSergio(sergiox.avg16_rtime(88)); %% this gives real bad CO2
figure(1); plot(fcris,sergiox.avg16_btanom(88,:),'c.-',fcris,sergioy.avg16_btanom(88,:),'m.-'); hl = legend('CRIS obs','CRIS cal','location','best','fontsize',10);
title(['CRIS anom : ' num2str(yyC) '/' num2str(mmC,'%02d') '/' num2str(ddC,'%02d') ])
xlim([min(fcris) max(fcris)])
xlim([min(fcris) 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fairs = instr_chans2645;

fnamex = ['../AIRS_new_clear_scan_August2019_AMT2020PAPER/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '.mat'];
fnamey = ['../AIRS_new_clear_scan_August2019_AMT2020PAPER/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '_cal.mat'];
airssergiox = load(fnamex);
airssergioy = load(fnamey);

fnamez = ['../AIRS_new_clear_scan_August2019_AMT2020PAPER/ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_2012/sergio_latbin' num2str(iLatbin) '.mat'];
airssergioz = load(fnamez);
for ii = 1 : 1813/16
  ind = (1:16) + (ii-1)*16;
  airssergioz.all_bt_anom_16day(ii,:) = nanmean(airssergioz.all_bt_anom(ind,:),1);
  airssergioz.all_times_16day(ii)     = nanmean(nanmean(airssergioz.all_times(ind,:),1));
end
[yyz,mmz,ddz,hhz] = tai2utcSergio(dtime2tai(airssergioz.all_times_16day));

[yyA,mmA,ddA] = tai2utcSergio(airssergiox.avg16_rtime(309)); %% this gives fine CO2
figure(2); plot(fairs,airssergiox.avg16_btanom(309,:),'b.-',fairs,airssergioy.avg16_btanom(309,:),'r.-',fairs,airssergioz.all_bt_anom_16day(88,:),'g'); hl = legend('AIRS obs','AIRS cal','AIRS obs since 2012','location','best','fontsize',10);
title(['AIRSanom : ' num2str(yyA) '/' num2str(mmA,'%02d') '/' num2str(ddA,'%02d') ])
xlim([min(fairs) max(fairs)])
xlim([min(fairs) 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
plot(fcris,sergiox.avg16_btanom(88,:),'c.-',fcris,sergioy.avg16_btanom(88,:),'m.-',fairs,airssergiox.avg16_btanom(309,:),'b.-',fairs,airssergioy.avg16_btanom(309,:),'r.-',fairs,airssergioz.all_bt_anom_16day(88,:),'g');
xlim([min(fairs) max(fairs)])
xlim([min(fairs) 840])
hl = legend('CRIS obs','CRIS cal','AIRS obs','AIRS cal','AIRS obs since 2012','location','best','fontsize',10);

plot(fcris,sergiox.avg16_btanom(88,:),'c.-',fairs,airssergiox.avg16_btanom(309,:),'b.-',fairs,airssergioz.all_bt_anom_16day(88,:),'g')
xlim([min(fairs) max(fairs)])
xlim([min(fairs) 840])
hl = legend('CRIS obs','AIRS obs','AIRS obs since 2012','location','best','fontsize',10);
