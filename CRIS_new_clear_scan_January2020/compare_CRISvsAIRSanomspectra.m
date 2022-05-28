function [cris,airs] = compare_CRISvsAIRSanomspectra(iLatbin,iObsORCal,iTimeStepCRIS,iTimeStepAIRS)

%% function [cris,airs] = compare_CRISvsAIRSanomspectra(iLatbin,iObsORCal,iTimeStepCRIS,iTimeStepAIRS)
%% iTimeStepAIRS is optional, code tries to figure out which iTimeStepCRIS = iTimeStepAIRS

addpath /home/sergio/MATLABCODE/CRIS_Hi2Lo/

if iObsORCal > 0
  resultsC = load('../CRIS_new_clear_scan_January2020/anomaly_0dayavg_results.mat');
  resultsA = load('../AIRS_new_clear_scan_August2019/SAVE_LW_noCFC11_noN2O/anomaly_0dayavg_results.mat');

  booC = load(['../CRIS_new_clear_scan_January2020/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '.mat']);
  booA = load(['../AIRS_new_clear_scan_August2019/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '.mat']);
else
  resultsC = load('../CRIS_new_clear_scan_January2020/anomaly_0dayavg_results_cal.mat');
  resultsA = load('../AIRS_new_clear_scan_August2019/SAVE_LW_noCFC11_noN2O/anomaly_0dayavg_results_cal.mat');

  booC = load(['../CRIS_new_clear_scan_January2020/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '_cal.mat']);
  booA = load(['../AIRS_new_clear_scan_August2019/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '_cal.mat']);
end

latbins = resultsC.latbins;
trop = find(abs(latbins) <= 30)
dt = (2012.3-2002.75);
figure(1); plot(resultsA.okdates,nanmean(resultsA.co2(trop,:)),'b',resultsC.okdates,dt*2.2+nanmean(resultsC.co2(trop,:)),'r','linewidth',2); title('co2'); grid
figure(2); plot(resultsA.okdates,nanmean(resultsA.n2o(trop,:)),'b',resultsC.okdates,dt*1.0+nanmean(resultsC.n2o(trop,:)),'r','linewidth',2); title('n2o'); grid
figure(3); plot(resultsA.okdates,nanmean(resultsA.ch4(trop,:)),'b',resultsC.okdates,dt*4.0+nanmean(resultsC.ch4(trop,:)),'r','linewidth',2); title('ch4'); grid
figure(4); plot(resultsA.okdates,nanmean(resultsA.stemp(trop,:)),'b',resultsC.okdates,dt*0.0+nanmean(resultsC.stemp(trop,:)),'r','linewidth',2); title('stemp'); grid
disp('ret'); pause

load f1305.mat;
cris.f1305 = f1305;

load ../AIRS_new_clear_scan_August2019/f2645.mat
airs.f2645 = f2645;

airstime = booA.avg_doy_since2002/365 + 2002;
cristime = booC.avg_doy_since2012/365 + 2012;
oo = find(airstime >= cristime(1)); 
figure(1); 
plot(cris.f1305,std(booC.avg16_btanom),'b',airs.f2645,std(booA.avg16_btanom(oo,:)),'r')
title('std dev of the anomalies over all time since 2012')

addpath /home/sergio/MATLABCODE/TIME
[yy,mm,dd,hh] = tai2utcSergio(booC.avg16_rtime(iTimeStepCRIS));
fprintf(1,'iTimeStepCRIS = %3i corresponds to %4i/%2i/%2i \n',iTimeStepCRIS,yy,mm,dd)

if nargin == 1 | nargin == 2 
  error('need at least3 args : function [cris,airs] = compare_CRISvsAIRSanomspectra(iLatbin,iObsORCal,iTimeStepCRIS,iTimeStepAIRS)');
end

if nargin == 3
  boo = booC.avg16_rtime(iTimeStepCRIS);
  iTimeStepAIRS = find(booA.avg16_rtime >= boo,1);
  [yyA,mmA,ddA,hhA] = tai2utcSergio(booA.avg16_rtime(iTimeStepAIRS));
  fprintf(1,'iTimeStepCRIS = %3i corresponds to iTimeStepAIRS = %3i around %4i/%2i/%2i \n',iTimeStepCRIS,iTimeStepAIRS,yyA,mmA,ddA);
else
  [yyA,mmA,ddA,hhA] = tai2utcSergio(booA.avg16_rtime(iTimeStepAIRS));
  fprintf(1,'iTimeStepCRIS = %3i while iTimeStepAIRS = %3i around %4i/%2i/%2i \n',iTimeStepCRIS,iTimeStepAIRS,yyA,mmA,ddA);
end

cris.anomspectra = squeeze(booC.avg16_btanom(iTimeStepCRIS,:));

airs.anomspectra = squeeze(booA.avg16_btanom(iTimeStepAIRS,:));

figure(2)
plot(airs.f2645,airs.anomspectra,'b',cris.f1305,cris.anomspectra,'r'); 
str = ['Anomaly : iTimeStep=' num2str(iTimeStepCRIS) ' around ' num2str(yy) '/' num2str(mm,'%02d') '/' num2str(dd,'%02d')];
hl = legend('AIRS','CRIS','location','best');
title(str);
ax = axis; ax(1) = 640; ax(2) = 1800; axis(ax); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iObsORCal > 0
  dirA = '../AIRS_new_clear_scan_August2019/SAVE_LW_noCFC11_noN2O_covx10/OutputAnomaly_OBS/';
  dirA = '../AIRS_new_clear_scan_August2019/SAVE_LW_noCFC11_noN2O/OutputAnomaly_OBS/';
  dirC = '../CRIS_new_clear_scan_January2020/OutputAnomaly_OBS/';
else
  dirA = '../AIRS_new_clear_scan_August2019/SAVE_LW_noCFC11_noN2O_covx10/OutputAnomaly_CAL/';
  dirA = '../AIRS_new_clear_scan_August2019/SAVE_LW_noCFC11_noN2O/OutputAnomaly_CAL/';
  dirC = '../CRIS_new_clear_scan_January2020/OutputAnomaly_CAL/';
end

fretrievA = [dirA '/' num2str(iLatbin,'%02d') '/anomtest_timestep' num2str(iTimeStepAIRS) '.mat'];
fretrievC = [dirC '/' num2str(iLatbin,'%02d') '/anomtest_timestep' num2str(iTimeStepCRIS) '.mat'];

retA = load(fretrievA);
retC = load(fretrievC);

airsfit = retA.oem.fit; airsret = retA.oem.finalrates; airsretsig = retA.oem.finalsigs;
crisfit = retC.oem.fit; crisret = retC.oem.finalrates; crisretsig = retC.oem.finalsigs;

figure(3)
plot(airs.f2645,airs.anomspectra,'b',cris.f1305,cris.anomspectra,'r',airs.f2645,airsfit,'c',cris.f1305,crisfit,'m'); 
if iObsORCal > 0
  str = ['OBS Anomaly : iTimeStep=' num2str(iTimeStepCRIS) ' ' num2str(yy) '/' num2str(mm,'%02d') '/' num2str(dd,'%02d')];
else
  str = ['CAL Anomaly : iTimeStep=' num2str(iTimeStepCRIS) ' ' num2str(yy) '/' num2str(mm,'%02d') '/' num2str(dd,'%02d')];
end
hl = legend('AIRS','CRIS','AIRSfit','CRISfit','location','best');
title(str);
ax = axis; ax(1) = 640; ax(2) = 1800; axis(ax); grid

iW = (1:20)+ 0*20 + 6;
iT = (1:20)+ 1*20 + 6;
iO = (1:20)+ 2*20 + 6;

if iObsORCal > 0
  %% rescale back to 2102 values
  airsret(1) = airsret(1)-2.2*(2012.3-2002.75);
  airsret(2) = airsret(2)-1.0*(2012.3-2002.75);
  airsret(3) = airsret(3)-4.0*(2012.3-2002.75);
  airsret(4) = airsret(4)+2.0*(2012.3-2002.75);
  airsret(5) = airsret(5)+2.0*(2012.3-2002.75);
end

disp('   var              AIRSadj               CRIS');
disp('------------------------------------------------------------');
ii=1; fprintf(1,' CO2         %8.4f +/- %8.4f   %8.4f +/- %8.4f \n',airsret(ii),airsretsig(ii),crisret(ii),crisretsig(ii))
ii=2; fprintf(1,' N2O         %8.4f +/- %8.4f   %8.4f +/- %8.4f \n',airsret(ii),airsretsig(ii),crisret(ii),crisretsig(ii))
ii=3; fprintf(1,' CH4         %8.4f +/- %8.4f   %8.4f +/- %8.4f \n',airsret(ii),airsretsig(ii),crisret(ii),crisretsig(ii))
ii=4; fprintf(1,' CFC11       %8.4f +/- %8.4f   %8.4f +/- %8.4f \n',airsret(ii),airsretsig(ii),crisret(ii),crisretsig(ii))
ii=5; fprintf(1,' CFC12       %8.4f +/- %8.4f   %8.4f +/- %8.4f \n',airsret(ii),airsretsig(ii),crisret(ii),crisretsig(ii))
ii=6; fprintf(1,' stemp/unadj %8.4f +/- %8.4f   %8.4f +/- %8.4f \n',airsret(ii),airsretsig(ii),crisret(ii),crisretsig(ii))

eraTz = load('era_ptempanom.mat'); eraTz = squeeze(eraTz.era_ptempanom(iLatbin,:,iTimeStepCRIS));
eraWV = load('era_gas_1anom.mat'); eraWV = squeeze(eraWV.era_gas_1anom(iLatbin,:,iTimeStepCRIS));
eraO3 = load('era_gas_3anom.mat'); eraO3 = squeeze(eraO3.era_gas_3anom(iLatbin,:,iTimeStepCRIS));

AeraTz = load('../AIRS_new_clear_scan_August2019/era_ptempanom.mat'); AeraTz = squeeze(AeraTz.era_ptempanom(iLatbin,:,iTimeStepAIRS));
%AeraWV = load('../AIRS_new_clear_scan_August2019/era_gas_1anom.mat'); AeraWV = squeeze(AeraWV.era_gas_1anom(iLatbin,:,iTimeStepAIRS));
AeraO3 = load('../AIRS_new_clear_scan_August2019/era_gas_3anom.mat'); AeraO3 = squeeze(AeraO3.era_gas_3anom(iLatbin,:,iTimeStepAIRS));

figure(4);
subplot(131); plot(airsret(iW),1:20,'b',crisret(iW),1:20,'r',eraWV,1:20,'k'); set(gca,'ydir','reverse'); title('WV'); grid
subplot(132); plot(airsret(iT),1:20,'b',crisret(iT),1:20,'r',eraTz,1:20,'k',AeraTz,1:20,'kx-'); set(gca,'ydir','reverse'); title('Tz'); grid
subplot(133); plot(airsret(iO),1:20,'b',crisret(iO),1:20,'r',eraO3,1:20,'k'); set(gca,'ydir','reverse'); title('O3'); grid

error(';lksg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk = num2str(iTimeStepCRIS,'%03d');
newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CLO_Anomaly137_16_12p8/RESULTS_FiniteDiff_Try3/'];       %% no seasonal, redone
newjacname = [newjacname '/kcarta_' junk '_tracegas_finitediff_5_2235_V4.mat']; %% Aug 21,         better CO2/CH4/N2O/CFC11/CFC12 prof, fixed 2002/09 yay
crisjac = load(newjacname);
f = crisjac.f;
co2cris = squeeze(crisjac.tracegas(iLatbin,:,1))';
[foutC,co2jacoutC] = translate_hi2lo(f,co2cris);

junk = num2str(iTimeStepCRIS,'%03d');
newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CLO_Anomaly137_16_12p8/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2235.mat']; 
crisjac = load(newjacname);
f = crisjac.f;
crisjac = squeeze(crisjac.M_TS_jac_all(iLatbin,:,:));
[foutC,crisjac] = translate_hi2lo(f,crisjac);
foutC = foutC.vchan';

junk = num2str(iTimeStepAIRS,'%03d');
newjacname = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16_12p8/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2645.mat'];
airsjac = load(newjacname);
foutA = airsjac.f;
airsjac = squeeze(airsjac.M_TS_jac_all(iLatbin,:,:));

%whos fout* *jac
iW = (1:97)+ 0*97 + 6;
iT = (1:97)+ 1*97 + 6;
iO = (1:97)+ 2*97 + 6;

figure(2);
  plot(foutA,sum(airsjac(:,iW),2),'b',foutC,sum(crisjac(:,iW),2),'r');
  str = ['WV jac : iTimeStep=' num2str(iTimeStepCRIS) ' around ' num2str(yy) '/' num2str(mm,'%02d') '/' num2str(dd,'%02d')];
  hl = legend('AIRS','CRIS','location','best');
  title(str);
  ax = axis; ax(1) = 640; ax(2) = 1800; axis(ax); grid

figure(3);
  plot(foutA,sum(airsjac(:,iT),2),'b',foutC,sum(crisjac(:,iT),2),'r');
  str = ['Tz jac : iTimeStep=' num2str(iTimeStepCRIS) ' around ' num2str(yy) '/' num2str(mm,'%02d') '/' num2str(dd,'%02d')];
  hl = legend('AIRS','CRIS','location','best');
  title(str);
  ax = axis; ax(1) = 640; ax(2) = 1800; axis(ax); grid

figure(4);
  plot(foutA,sum(airsjac(:,iO),2),'b',foutC,sum(crisjac(:,iO),2),'r');
  str = ['O3 jac : iTimeStep=' num2str(iTimeStepCRIS) ' around ' num2str(yy) '/' num2str(mm,'%02d') '/' num2str(dd,'%02d')];
  hl = legend('AIRS','CRIS','location','best');
  title(str);
  ax = axis; ax(1) = 640; ax(2) = 1800; axis(ax); grid

figure(5)
  iCO2 = 1;
  plot(foutA,airsjac(:,iCO2),'b',foutC,crisjac(:,iCO2),'r',foutC,co2jacoutC,'k');
  str = ['CO2 jac : iTimeStep=' num2str(iTimeStepCRIS) ' around ' num2str(yy) '/' num2str(mm,'%02d') '/' num2str(dd,'%02d')];
  hl = legend('AIRS','CRIS','CRIS FINITE','location','best');
  title(str);
  ax = axis; ax(1) = 640; ax(2) = 1800; axis(ax); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
