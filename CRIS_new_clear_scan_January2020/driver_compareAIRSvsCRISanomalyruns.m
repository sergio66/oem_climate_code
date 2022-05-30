% see clust_convert_cris_strowrates2oemrates_anomaly.m hich makes eg junk = load(['ANOM_16dayavg/latbin_0dayavg_' num2str(iii) '.mat']);
load f1305.mat
junk = load('/home/strow/Work/Cris/Stability/Data/Desc_fits/fit_robs_lat20.mat');
iC791 = find(f1305 >= 791-0.5,1);
iC792 = find(f1305 >= 792-0.5,1);
plot(f1305,junk.dbt,f1305([iC791 iC792]),junk.dbt([iC791 iC792]),'ro'); xlim([790 795])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('seeing CRIS anom runs for diff chans')
load f1305.mat
iLatbin = 20;
  booC_obs = load(['../CRIS_new_clear_scan_January2020/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '.mat']);
  booC_cal = load(['../CRIS_new_clear_scan_January2020/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '_cal.mat']);
i709 = find(f1305 >= 709,1);
i791 = find(f1305 >= 791.5,1);
i900 = find(f1305 >= 900,1);
i1420 = find(f1305 >= 1420,1);
iaChan = [i709 i791 i900 i1420];
figure(1)
plot(booC_obs.avg_doy_since2012,booC_obs.avg16_btanom(:,i709),'b',booC_obs.avg_doy_since2012,booC_cal.avg16_btanom(:,i709),'c')
plot(2012 + booC_obs.avg_doy_since2012/365.25,booC_obs.avg16_btanom(:,i709),'b',2012 + booC_obs.avg_doy_since2012/365.25,booC_cal.avg16_btanom(:,i709),'c')
plot(2012 + booC_obs.avg_doy_since2012/365.25,booC_obs.avg16_btanom(:,i709)-booC_cal.avg16_btanom(:,i709),'c')
plot(2012 + booC_obs.avg_doy_since2012/365.25,booC_obs.avg16_btanom(:,iaChan)-booC_cal.avg16_btanom(:,iaChan),'linewidth',2)
  hl = legend('709','791','900','1420','location','best'); grid; title('CRIS Latbin 20 : ObsAnom-CalAnom')

disp('seeing AIRS anom runs for diff chans')
load ../AIRS_new_clear_scan_August2019_AMT2020PAPER/f2645.mat
iLatbin = 20;
  booA_obs = load(['../AIRS_new_clear_scan_August2019_AMT2020PAPER/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '.mat']);
  booA_cal = load(['../AIRS_new_clear_scan_August2019_AMT2020PAPER/ANOM_16dayavg/latbin_0dayavg_' num2str(iLatbin) '_cal.mat']);
i709 = find(f2645 >= 709,1);
i791 = find(f2645 >= 791.5,1);
i900 = find(f2645 >= 900,1);
i1420 = find(f2645 >= 1420,1);
iaChan = [i709 i791 i900 i1420];
figure(2)
plot(booA_obs.avg_doy_since2002,booA_obs.avg16_btanom(:,i709),'b',booA_obs.avg_doy_since2002,booA_cal.avg16_btanom(:,i709),'c')
plot(2002 + booA_obs.avg_doy_since2002/365.25,booA_obs.avg16_btanom(:,i709),'b',2002 + booA_obs.avg_doy_since2002/365.25,booA_cal.avg16_btanom(:,i709),'c')
plot(2002 + booA_obs.avg_doy_since2002/365.25,booA_obs.avg16_btanom(:,i709)-booA_cal.avg16_btanom(:,i709),'c')
plot(2002 + booA_obs.avg_doy_since2002/365.25,booA_obs.avg16_btanom(:,iaChan)-booA_cal.avg16_btanom(:,iaChan),'linewidth',2)
  hl = legend('709','791','900','1420','location','best'); grid; title('AIRS Latbin 20 : ObsAnom-CalAnom')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret'); pause
file1C = 'SAVE_LW_noCFC11_Feb14_2020/anomaly_0dayavg_results.mat';
file2C = 'SAVE_LW_noCFC11_Feb14_2020/anomaly_0dayavg_cal_results.mat';
fileC  = 'SAVE_LW_noCFC11_Feb14_2020/anomaly_0dayavg_results_spectra.mat';

file1A = '../AIRS_new_clear_scan_August2019_AMT2020PAPER/SAVE_LW_noCFC11/anomaly_0dayavg_results.mat';
file2A = '../AIRS_new_clear_scan_August2019_AMT2020PAPER/SAVE_LW_noCFC11/anomaly_0dayavg_cal_results.mat';
fileA  = '../AIRS_new_clear_scan_August2019_AMT2020PAPER/SAVE_LW_noCFC11/anomaly_0dayavg_results_spectra.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 6
  figure(ii); clf
end

[A,C] = compare_anomaly_runs2datasets(file1A,file2A,file1C,file2C);
btA = load(fileA);
btC = load(fileC);

plot_compareAIRSvsCRISanomalyruns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/TIME
[h,ha,pairs,pa] = rtpread('../../oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Aug20_2019_Clr/Desc/16DayAvgNoS/latbin20_16day_avg.rp.rtp');
[h,ha,pcris,pa] = rtpread('../../oem_pkg_run_sergio_AuxJacs/MakeProfs/CLO_LAT40_avg_made_Dec2019_Clr/Desc/16DayAvgNoS/latbin20_16day_avg.rp.rtp');

[yyA,mmA,ddA,hhA] = tai2utcSergio(pairs.rtime);
[yyC,mmC,ddC,hhC] = tai2utcSergio(pcris.rtime);

daysA = change2days(yyA,mmA,ddA,2002);
daysC = change2days(yyC,mmC,ddC,2002);

figure(9); plot(daysC,pcris.stemp,'b',daysA,pairs.stemp,'r'); title('Stemp latbin20 (b)CRIS (r)AIRS')
  xlabel('days since 2002')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iCheckJac = -1;
%% this is VERY slow as it reads in 157 jacs for CRIS and 365 jacs for AIRS, about 10 minutes
if iCheckJac > 0
  plot_compareAIRSvsCRISanomalyjacs
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
latbinsx = equal_area_spherical_bands(20);
latbinsx = meanvaluebin(latbinsx);
trop = find(abs(latbinsx) <= 30);
junkkk = zeros(157,1305);
for iii = 1 : length(trop)
  junk = load(['ANOM_16dayavg/latbin_0dayavg_' num2str(trop(ii)) '.mat']);
  junkkk = junkkk + junk.avg16_btanom;
end
junkkk = junkkk/length(trop);

iii = 20;  junk = load(['ANOM_16dayavg/latbin_0dayavg_' num2str(iii) '.mat']);
figure(10); clf; pcolor(C.okdates,f1305,junk.avg16_btanom'); shading flat; colorbar; colormap jet; title('rtp CrIS : <btAnom>'); ylabel('f1305'); xlabel('time')
  disp('blanks are where hi2lo puts nan for guard channels etc')
figure(11); clf; plot(C.okdates,junk.avg16_btanom(:,iC1231));                           title('rtp CrIS : <btAnom>'); ylabel('btnanom 1231 cm-1');  plotaxis2;
figure(12); clf; plot(C.okdates,junk.avg16_btanom(:,iC791)-junk.avg16_btanom(:,iC792)); title('rtp CrIS : <btAnom>'); ylabel('btnanom 791-792 cm-1'); plotaxis2;
[mmbad,nnbad] = find(isnan(junk.avg16_btanom)); unique(nnbad)

figure(10); clf; pcolor(C.okdates,f1305,junkkk'); shading flat; colorbar; colormap jet; title('rtp CrIS : <btAnom>'); ylabel('f1305'); xlabel('time')
  disp('blanks are where hi2lo puts nan for guard channels etc')
figure(11); clf; plot(C.okdates,junkkk(:,iC1231));                           title('rtp CrIS : <btAnom>'); ylabel('btnanom 1231 cm-1');  plotaxis2;
figure(12); clf; plot(C.okdates,junkkk(:,iC791)-junkkk(:,iC792)); title('rtp CrIS : <btAnom>'); ylabel('btnanom 791-792 cm-1'); plotaxis2;
[mmbad,nnbad] = find(isnan(junkkk)); unique(nnbad)
