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

junk = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CLO_Anomaly137_16_12p8/RESULTS/kcarta_085_M_TS_jac_all_5_97_97_97_2235.mat');
qrenorm = junk.qrenorm;
joutCO2 = squeeze(junk.M_TS_jac_all(20,:,1));
joutST  = squeeze(junk.M_TS_jac_all(20,:,6));
joutWV  = squeeze(junk.M_TS_jac_all(20,:,007:103));
joutT   = squeeze(junk.M_TS_jac_all(20,:,104:200));
junk = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CLO_Anomaly137_16_12p8/RESULTS/kcarta_085_M_TS_jac_all_5_97_97_97_2235.mat','f');
[fout,jacout] = translate_hi2lo(junk.f,zeros(size(junk.f)));
[fout,jacoutCO2] = translate_hi2lo(junk.f,joutCO2'); plot(fout.vchan,jacoutCO2);
[fout,jacoutST] = translate_hi2lo(junk.f,joutST');   plot(fout.vchan,jacoutST);
[fout,jacoutWV] = translate_hi2lo(junk.f,joutWV);    plot(fout.vchan,sum(jacoutWV'));
[fout,jacoutT] = translate_hi2lo(junk.f,joutT);      plot(fout.vchan,sum(jacoutT'));
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

iA = 309; %% this is March 2016 for AIRS 2002
iC = 88;  %% this is March 2016 for CrIS 2021, AIRS 2012

iOffSet = input('enter offset to the date plotted (March 2016 is default date, offset = 0 default, remeber +/-23 = +/- 1 year) : ');
if length(iOffSet) == 0
  iOffSet = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
junkAIRS = load(['../AIRS_new_clear_scan_August2019_AMT2020PAPER/OutputAnomaly_OBS/20/anomtest_timestep' num2str(iA+iOffSet) '.mat']);
plot(instr_chans2645,junkAIRS.rateset.rates)
fprintf(' <<<< AIRS retrieval >>>> for timestep %3i = %8.4f ppm while xb(1) was %8.6f ppm \n',iA+iOffSet,junkAIRS.oem.finalrates(1),junkAIRS.oem.xb(1));

junkCrIS = load(['OutputAnomaly_OBS/20/anomtest_timestep' num2str(iC) '.mat']);
plot(fcris,junkCrIS.rateset.rates)
fprintf(' <<<< CrIS retrieval >>>> for timestep %3i = %8.4f ppm while xb(1) was %8.6f ppm \n',iC,junkCrIS.oem.finalrates(1),junkCrIS.oem.xb(1));

junk = [junkAIRS.oem.xb(1) junkCrIS.oem.xb(1) junkAIRS.oem.xb(1)-junkCrIS.oem.xb(1) junkAIRS.oem.finalrates(1)  junkCrIS.oem.finalrates(1)   junkAIRS.oem.finalrates(1)-junkCrIS.oem.finalrates(1)];
fprintf(1,' so << AIRS - CRIS >> start CO2 = %8.6f - %8.6f = %8.6f ppm     while retrieved =  %8.6f - %8.6f = %8.6f ppm \n',junk);
plot(fairs,junkAIRS.rateset.rates,fcris,junkCrIS.rateset.rates); hl = legend('AIRS','CRIS','location','best');

fprintf(1,'reading in latbin 20 AIRS 365 steps (23/year) : ');
for ii = 1 : 365
  if mod(ii,23) == 0
    fprintf(1,'.');
  end
  junkAIRS = load(['../AIRS_new_clear_scan_August2019_AMT2020PAPER/OutputAnomaly_OBS/20/anomtest_timestep' num2str(ii) '.mat']);
  %junkAIRS = load(['../AIRS_new_clear_scan_August2019_AMT2020PAPER/SAVE_LW_noCFC11_noN2O//OutputAnomaly_OBS/20/anomtest_timestep' num2str(ii) '.mat']);
  %plot(instr_chans2645,junkAIRS.rateset.rates); title(num2str(ii)); pause(0.1);
  airs0_CO2(ii) = junkAIRS.oem.xb(1);
  airsF_CO2(ii) = junkAIRS.oem.finalrates(1);
end
fprintf(1,'\n');
miaow = load('../AIRS_new_clear_scan_August2019_AMT2020PAPER/SAVE_LW_noCFC11_noN2O/anomaly_0dayavg_results.mat');
[yyA,mmA,ddA] = tai2utcSergio(airssergiox.avg16_rtime);
ttA = yyA + (mmA-1)/12 + (ddA-1)/12/30;

plot(ttA,airs0_CO2,ttA,airsF_CO2);
line([2016.25 2016.25],[-20 +50],'color','k','linewidth',2);
xlim([min(ttA) max(ttA)]);
hl = legend('Start CO2','Retr CO2','location','best');

fprintf(1,'reading in latbin 20 CRIS 157 steps (23/year) : ');
for ii = 1 : 157
  if mod(ii,23) == 0
    fprintf(1,'.');
  end
  junkCrIS = load(['OutputAnomaly_OBS/20/anomtest_timestep' num2str(ii) '.mat']);
  %plot(fcris,junkCrIS.rateset.rates); title(num2str(ii)); pause(0.1);
  cris0_CO2(ii) = junkCrIS.oem.xb(1);
  crisF_CO2(ii) = junkCrIS.oem.finalrates(1);
end
fprintf(1,'\n');
[yyC,mmC,ddC] = tai2utcSergio(sergiox.avg16_rtime);
ttC = yyC + (mmC-1)/12 + (ddC-1)/12/30;

plot(ttC,cris0_CO2,ttC,crisF_CO2);
line([2016.25 2016.25],[-10 +20],'color','k','linewidth',2);
xlim([min(ttC) max(ttC)]); 
hl = legend('Start CO2','Retr CO2','location','best'); 

plot(ttA,airs0_CO2,'b',ttA,airsF_CO2,'r',ttC,cris0_CO2,'c',ttC,crisF_CO2,'m');
line([2016.25 2016.25],[-20 +50],'color','k','linewidth',2);
xlim([min(ttA) max(ttA)]);
hl = legend('AIRS Start CO2','AIRS Retr CO2','CRIS Start CO2','CRIS Retr CO2','location','best','fontsize',10); ylabel('CO2 ppmv'); xlabel('time')

%plot(ttA,airs0_CO2,'b',ttA,airsF_CO2,'r',ttA,miaow.co2(20,:))
error(';skjdg')
%%%%%%%%%%%%%%%%%%%%%%%%%
  
fnamez = ['../AIRS_new_clear_scan_August2019_AMT2020PAPER/ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_2012/sergio_latbin' num2str(iLatbin) '.mat'];
airssergioz = load(fnamez);
for ii = 1 : 1813/16
  ind = (1:16) + (ii-1)*16;
  airssergioz.all_bt_anom_16day(ii,:) = nanmean(airssergioz.all_bt_anom(ind,:),1);
  airssergioz.all_times_16day(ii)     = nanmean(nanmean(airssergioz.all_times(ind,:),1));
end
[yyz,mmz,ddz,hhz] = tai2utcSergio(dtime2tai(airssergioz.all_times_16day));

[yyA,mmA,ddA] = tai2utcSergio(airssergiox.avg16_rtime(iA+iOffSet)); %% this gives fine CO2
figure(2); plot(fairs,airssergiox.avg16_btanom(iA+iOffSet,:),'b.-',fairs,airssergioy.avg16_btanom(iA+iOffSet,:),'r.-',fairs,airssergioz.all_bt_anom_16day(iC+iOffSet,:),'g'); 
  hl = legend('AIRS obs','AIRS cal','AIRS obs since 2012','location','best','fontsize',10);
title(['AIRSanom : ' num2str(yyA) '/' num2str(mmA,'%02d') '/' num2str(ddA,'%02d') ])
xlim([min(fairs) max(fairs)])
xlim([min(fairs) 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
plot(fcris,sergiox.avg16_btanom(iC+iOffSet,:),'c.-',fcris,sergioy.avg16_btanom(iC+iOffSet,:),'m.-',fairs,airssergiox.avg16_btanom(iA+iOffSet,:),'b.-',fairs,airssergioy.avg16_btanom(iA+iOffSet,:),'r.-',fairs,airssergioz.all_bt_anom_16day(iC+iOffSet,:),'g');
xlim([min(fairs) max(fairs)])
xlim([min(fairs) 840])
hl = legend('CRIS obs','CRIS cal','AIRS obs','AIRS cal','AIRS obs since 2012','location','best','fontsize',10);

plot(fcris,sergiox.avg16_btanom(iC+iOffSet,:),'b.-',fairs,airssergiox.avg16_btanom(iA+iOffSet,:),'r.-',fairs,airssergioz.all_bt_anom_16day(iC+iOffSet,:),'k')
xlim([min(fairs) max(fairs)])
xlim([min(fairs) 840])
hl = legend('CRIS obs','AIRS obs','AIRS obs since 2012','location','best','fontsize',10);

figure(4); plot(fout.vchan,jacoutCO2*10,'r',fcris,sergiox.avg16_btanom(iC+iOffSet,:),'b.-'); hl = legend('CO2 jac','CRIS anom obs','location','best','fontsize',10); 
title(['CRIS anom : ' num2str(yyC) '/' num2str(mmC,'%02d') '/' num2str(ddC,'%02d') ])
xlim([640 840])

figure(5); plot(fout.vchan,-1 + -15*jacoutCO2 ,'r',fcris,sergioy.avg16_btanom(iC+iOffSet,:),'b.-'); plotaxis2; hl = legend('CO2 jac','CRIS anom cal','location','best','fontsize',10); 

simple = -200*jacoutST + -15*jacoutCO2 + (sum(jacoutT'))';
simple = -20*jacoutST + -15*jacoutCO2 + -100*(sum(jacoutT'))';
simple = -20*jacoutST + -15*jacoutCO2 + -100*(sum(jacoutT'))' + 10*(sum(jacoutWV'))';
simple = +15*jacoutST + -15*jacoutCO2 + -100*(sum(jacoutT'))' + 10*(sum(jacoutWV'))';
         [+15*qrenorm(6) -15*qrenorm(1) -100*qrenorm(150) +10*qrenorm(10)]
simple = +1.14/qrenorm(6)*jacoutST + -1/qrenorm(1)*jacoutCO2 + +1/qrenorm(150)*(sum(jacoutT'))' + 0.25/qrenorm(10)*(sum(jacoutWV'))';

figure(5); plot(fout.vchan,simple,'r',fcris,sergioy.avg16_btanom(iC+iOffSet,:),'b.-'); plotaxis2; hl = legend('CO2 jac','CRIS anom cal','location','best','fontsize',10); 
title(['jac scale factor = -15 CRIS anom : ' num2str(yyC) '/' num2str(mmC,'%02d') '/' num2str(ddC,'%02d') ])
xlim([640 840])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6); clf; 
  boo1 = load('individual_prof_convolved_kcarta_crisHI_crisMED_1_hamming.mat');
  plot(boo1.lo_fcris,rad2bt(boo1.lo_fcris,boo1.lo_rcris_all),'b','linewidth',2); hold on
  boo2 = load('individual_prof_convolved_kcarta_crisHI_crisMED_1_sinc.mat');
  hold on
  plot(boo2.lo_fcris,rad2bt(boo2.lo_fcris,boo2.lo_rcris_all),'r'); xlim([640 840])
  hold off
hl = legend('hamming','sinc','location','best'); title('tropical profile BT');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fac = 2.2/400;
figure(7); clf; 
  boo1 = load('individual_prof_convolved_kcarta_crisHI_crisMED_1_coljac_hamming.mat');
  plot(boo1.lo_fcris,boo1.lo_rcris_all(:,1)*fac,'b','linewidth',2); hold on
  boo2 = load('individual_prof_convolved_kcarta_crisHI_crisMED_1_coljac_sinc.mat');
  hold on
  plot(boo2.lo_fcris,boo2.lo_rcris_all(:,1)*fac,'r'); xlim([640 840])
  hold off
hl = legend('hamming','sinc','location','best'); title('CO2 jac');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8); clf; 
  boo1 = load('individual_prof_convolved_kcarta_crisHI_crisMED_1_jac_hamming_100layers.mat');
  plot(boo1.lo_fcris,sum(boo1.lo_rcris_all(:,1:97),2)*fac,'b','linewidth',2); hold on
  boo2 = load('individual_prof_convolved_kcarta_crisHI_crisMED_1_jac_sinc_100layers.mat');
  hold on
  plot(boo2.lo_fcris,sum(boo2.lo_rcris_all(:,1:97),2)*fac,'r'); xlim([640 840])
  hold off
hl = legend('hamming','sinc','location','best'); title('CO2 jac');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(9)
finX = ['/home/strow/Work/Cris/Stability/Data/Desc/statlat' num2str(20) '.mat'];
aX = load(finX);
load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/cris_ichan_vchan_nsr_fsr.mat
plot(nsr.vchan,rad2bt(nsr.vchan,squeeze(aX.robs(1,5,:))),nsr.vchan,rad2bt(nsr.vchan,squeeze(aX.rclr(1,5,:)))); xlim([640 840])
plot(nsr.vchan,rad2bt(nsr.vchan,squeeze(nanmean(aX.robs(1,:,:),2))),nsr.vchan,rad2bt(nsr.vchan,squeeze(nanmean(aX.rclr(1,:,:),2)))); xlim([640 840])
