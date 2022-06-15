function A = compare_anomaly_runs(file1,file2,file3);

%% file1 = 'anomaly_0dayavg_results.mat'; file2 = 'anomaly_0dayavg_cal_results.mat'; compare_anomaly_runs(file1,file2);

for ii=1:8; figure(ii); clf; end

%{

crisL1c = load('/home/strow/Work/Cris/Stability/Data/Desc_fits/fit_robs_lat20.mat');
airsL1c = load('/home/strow/Work/Airs/Stability/Data/Desc_fits/fit_robs_lat20.mat');

i1231airs = 1520;
i1231cris = 731;

[yyA,mmA,hhA,ddA] = tai2utcSergio(dtime2tai(airsL1c.all_times(:,i1231airs)));
[yyC,mmC,hhC,ddC] = tai2utcSergio(dtime2tai(crisL1c.all_times(:,i1231cris)));

airs_dayssince2002 = change2days(yyA,mmA,ddA,2002);
cris_dayssince2002 = change2days(yyC,mmC,ddC,2002);

figure(1)
plot(airs_dayssince2002,airsL1c.all_bt_anom(:,i1231airs),'b',cris_dayssince2002,crisL1c.all_bt_anom(:,i1231cris)*1000,'r')
plot(airs_dayssince2002,smooth(airsL1c.all_bt_anom(:,i1231airs),180),'b',cris_dayssince2002,smooth(crisL1c.all_bt_anom(:,i1231cris)*1000,180),'r')
plot(2002+airs_dayssince2002/365,smooth(airsL1c.all_bt_anom(:,i1231airs),180),'b',2002+cris_dayssince2002/365,smooth(crisL1c.all_bt_anom(:,i1231cris)*1000,180),'r','linewidth',2)
   axis([2012 2020 0 +2]); grid; hl = legend('AIRS','CRIS','location','best'); title('BT1231 anomaly smoothed 180 days')

%%%%%%%%%%%%%%%%%%%%%%%%%
f2378 = instr_chans; f2378([445 449])

hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
f = vchan2834;
load sarta_chans_for_l1c.mat
fairs = f(ichan);

xyz = load('f1305.mat');
fcris = xyz.f1305;

iA1 = find(fairs >= f2378(445),1); iA2 = find(fairs >= f2378(449),1);
iC1 = find(fcris >= f2378(445),1); iC2 = find(fcris >= f2378(449),1);

figure(2)
plot(2002+airs_dayssince2002/365,smooth(airsL1c.all_bt_anom(:,iA1)-airsL1c.all_bt_anom(:,iA2),180),'b',2002+cris_dayssince2002/365,smooth(crisL1c.all_bt_anom(:,iC1)*1000-crisL1c.all_bt_anom(:,iC2)*1000,180)+0.5,'r','linewidth',2)
   axis([2012 2020 -1 +1]); grid; hl = legend('AIRS','CRIS','location','best'); title('BT 790-791 anomaly smoothed 180 days')

%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = 'anomaly_0dayavg_results.mat'; f2 = 'anomaly_0dayavg_cal_results.mat'; cris = compare_anomaly_runs(f1,f2);

f1 = '../AIRS_new_clear_scan_August2019/SAVE_BESTRUNv1/anomaly_0dayavg_results_save_bestrunv1.mat'; 
f2 = '../AIRS_new_clear_scan_August2019/SAVE_BESTRUNv1/anomaly_0dayavg_cal_results_save_bestrunv1.mat'; airs = compare_anomaly_runs(f1,f2);

figure(1); plot(cris.okdates,cris.co2_1,airs.okdates,airs.co2_1-(2012.333-2002.75)*2.2,'linewidth',2); hl = legend('CRIS','AIRS','location','best'); title('OBS'); grid
figure(2); plot(cris.okdates,cris.co2,airs.okdates,airs.co2-(2012.333-2002.75)*2.2,'linewidth',2); hl = legend('CRIS','AIRS','location','best'); title('OBS-CAL'); grid

figure(1); plot(cris.okdates,cris.ch4_1,airs.okdates,airs.ch4_1-(2012.333-2002.75)*5,'linewidth',2); hl = legend('CRIS','AIRS','location','best'); title('OBS'); grid
figure(2); plot(cris.okdates,cris.ch4,airs.okdates,airs.ch4-(2012.333-2002.75)*5,'linewidth',2); hl = legend('CRIS','AIRS','location','best'); title('OBS-CAL'); grid

figure(1); plot(cris.okdates,cris.stemp,airs.okdates,airs.stemp,'linewidth',2); hl = legend('CRIS','AIRS','location','best'); title('OBS-CAL'); grid
figure(2); plot(cris.okdates,cris.stemp_1,airs.okdates,airs.stemp_1,'linewidth',2); hl = legend('CRIS','AIRS','location','best'); title('OBS'); grid


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = struct;

iMaxTimeSteps   = 157;

iMaxTimeSteps   = input('Enter iMaxTimeSteps (157 for CriS <default> , 365 for AIRS or -223 (9/8 yerars x 365/16) = -223 for AIRS to start roughtlu at 2012/05) : ');
if length(iMaxTimeSteps) == 0
   iMaxTimeSteps = 157;
end

iMaxTimeStepsM1 = abs(iMaxTimeSteps)-1;

if nargin < 2
  error('need at least 2 input files, at most 3')
end
%if nargin == 0
%  file1 = 'anomaly_0dayavg_results.mat'; file2 = 'anomaly_0dayavg_cal_results.mat'; 
%end

iSmooth = input('Enter how many timesteps to smooth over (each timestep = 16 days = 0.0438 yr; default = 11 timesteps = 0.5 yr) : ');
if length(iSmooth) == 0
  iSmooth = 1;
end

a1 = load(file1);
a2 = load(file2); 
if nargin == 3
  a3 = load(file3); 
end

addpath /home/sergio/MATLABCODE/PLOTTER
iaLat = equal_area_spherical_bands(20);
iaLat = 0.5 * (iaLat(1:end-1) + iaLat(2:end));

xtropics = input('Enter range (-3) SP (-2) SML (-1) ST (+1) NT (+2) NML (+3) NP      (0) tropics (-90) SH (+90) NH (-60) SML and ST (+60) NML and NT : ');
if length(xtropics) == 0
  xtropics = 0;
end
if xtropics == -3
  disp('subbing for Southern Polar')
  tropics = find(iaLat <= -60);  
elseif xtropics == +3
  disp('subbing for Northern Polar')
  tropics = find(iaLat >= +60);  
elseif xtropics == -2
  disp('subbing for Southern MidLats')
  tropics = find(iaLat > -60 & iaLat <= -30);  
elseif xtropics == +2
  disp('subbing for Norththern MidLats')
  tropics = find(iaLat > +30 & iaLat <= +60);  
elseif xtropics == -1
  disp('subbing for Southern Tropics')
  tropics = find(iaLat > -30 & iaLat <= -0);  
elseif xtropics == +1
  disp('subbing for Northern Tropics')
  tropics = find(iaLat > +0 & iaLat <= +30);  
elseif xtropics == 0
  disp('subbing for ALL Tropics')
  tropics = find(iaLat > -30 & iaLat <= +30);  
elseif xtropics == -90
  disp('subbing for Southern Hemisphere')
  tropics = find(iaLat <= +0);  
elseif xtropics == +90
  disp('subbing for Northern Hemisphere')
  tropics = find(iaLat >= +0);  
elseif xtropics == -60
  disp('subbing for Southern (ML and T)')
  tropics = find(iaLat >= -60 & iaLat < 0);  
elseif xtropics == +60
  disp('subbing for Northern (ML and T)')
  tropics = find(iaLat >= 0 & iaLat < 60);  
else
  error('error in your choice')
end

fprintf(1,'your zonal avg choice contains %2i latbins \n',length(tropics));

if nargin == 2
  %{
  useful for eg 
  f2 = 'anomaly_0dayavg_resultsXnoloop_unc0p0025.mat';
  f3 = 'anomaly_0dayavg_cal_resultsXnoloop_unc0p0025.mat';
  compare_anomaly_runs(f2,f3)
  %}
  iYes = input('pretend a3 = a1-a2 (-1/+1) ? ');
  if length(iYes) == 0
    iYes = 1;
  end
  if iYes > 0
    a3.stemp = a1.stemp - a2.stemp;
    a3.co2 = a1.co2 - a2.co2;
    a3.n2o = a1.n2o - a2.n2o;
    a3.ch4 = a1.ch4 - a2.ch4;
    a3.cfc11 = a1.cfc11 - a2.cfc11;
    a3.cfc12 = a1.cfc12 - a2.cfc12;

    if ~isfield(a1,'topts')
      a1.topts.nloop = -9999;
    end
    if ~isfield(a2,'topts')
      a2.topts.nloop = -9999;
    end

    if ~isfield(a1,'topts')
      a1.topts.nloop = -9999;
    end
    if ~isfield(a2,'topts')
      a2.topts.nloop = -9999;
    end

    if isfield(a1,'topts') & isfield(a2,'topts')
      if isfield(a1.topts,'nloop') & isfield(a2.topts,'nloop')
        a3.topts.nloop = min(a1.topts.nloop,a2.topts.nloop);
      else
        a3.topts.nloop = -9999;
      end
    else
      a3.topts.nloop = -9999;
    end
    a3.okdates = a2.okdates;
  end
end

if exist('a3')
  [mmbad,nnbad] = find(abs(a1.stemp) > 5 | abs(a2.stemp) > 5 |abs(a3.stemp) > 5 | abs(a1.co2 > 50) | abs(a2.co2 > 50) | abs(a3.co2 > 50)); 
end

[mmbad,nnbad] = find(abs(a1.stemp) > 5);  %% what was used in Sep 11, 2019
if length(mmbad) > 0
  for zz = 1 : length(mmbad)
    fprintf(1,'%3i %3i \n',mmbad(zz),nnbad(zz))
  end
end
%% to debug, if length(mmbad) > 0, but I did not find anything wrong with input data
iDebug = -1;

if length(mmbad) > 100 & iDebug < 0
  disp('oooh LOTS of bad points, really should try to debug')
  mmu = unique(mmbad);
  mmn = unique(nnbad);
  whos *bad* mmu mmn
  figure(1); plot(1:365,a1.stemp); pause(0.1)
  figure(2); pcolor(a1.stemp); caxis([-5 +5]); shading flat; colorbar; hold on; plot(nnbad,mmbad,'ko','linewidth',2); hold off;
  disp('wll be setting "bad" stemp and co2 and ch4 to NAN : ret'); pause

  a1.stemp(mmbad,nnbad) = NaN;
  a1.co2(mmbad,nnbad) = NaN;
  a1.n2o(mmbad,nnbad) = NaN;
  a1.ch4(mmbad,nnbad) = NaN;
  a1.cfc11(mmbad,nnbad) = NaN;
  a1.cfc12(mmbad,nnbad) = NaN;

  a2.stemp(mmbad,nnbad) = NaN;
  a2.co2(mmbad,nnbad) = NaN;
  a2.n2o(mmbad,nnbad) = NaN;
  a2.ch4(mmbad,nnbad) = NaN;
  a2.cfc11(mmbad,nnbad) = NaN;
  a2.cfc12(mmbad,nnbad) = NaN;

  a3.stemp(mmbad,nnbad) = NaN;
  a3.co2(mmbad,nnbad) = NaN;
  a3.n2o(mmbad,nnbad) = NaN;
  a3.ch4(mmbad,nnbad) = NaN;
  a3.cfc11(mmbad,nnbad) = NaN;
  a3.cfc12(mmbad,nnbad) = NaN;

elseif length(mmbad) > 100 & iDebug > 0
  for zz = 1 : length(mmbad)
    fprintf(1,'%3i %3i \n',mmbad(zz),nnbad(zz))
  end

  whos *bad*
  unique(mmbad)
  unique(nnbad)
  figure(1); plot(1:365,a1.stemp); pause(0.1)
  figure(2); pcolor(a1.stemp); caxis([-5 +5]); shading flat; colorbar; hold on; plot(nnbad,mmbad,'ko','linewidth',2); hold off;

  anom_noise = load('btn_avg.mat');
  btn_avg = anom_noise.btn_avg;
  uunn = unique(nnbad);
  for tt = 1 : length(uunn)
    figure(3); semilogy(1:2645,squeeze(btn_avg(tropics,uunn(tt),:)),'b',1:2645,mean(squeeze(btn_avg(tropics,uunn(tt),:))),'k'); grid 
    title(['All da tropical noise at timestep ' num2str(uunn(tt))]); pause;
    boo = squeeze(btn_avg(tropics,uunn(tt),:));
    ohoh = find(isnan(boo) | isinf(boo));
    if length(ohoh) > 0
      ohoh
    end
  end

  uunn = unique(nnbad);
  uumm = unique(mmbad);
  for zz = 1 : length(uumm)
    driver.rateset.datafile = ['ANOM_16dayavg/latbin_0dayavg_' num2str(uumm(zz)) '.mat'];
    boo = load(driver.rateset.datafile);
    for tt = 1 : length(uunn)
      %figure(3); semilogy(1:2645,boo.avg16_btanom(uunn(tt),:),'b',1:2645,mean(squeeze(btn_avg(tropics,uunn(tt),:))),'k'); grid 
      %title(['All da tropical signal (b) and noise (r) at timestep ' num2str(uunn(tt))]); pause;
      figure(3); plot(1:2645,boo.avg16_btanom(uunn(tt),:),'b'); grid 
      ohoh = find(isnan(boo.avg16_btanom(uunn(tt),:)) | isinf(boo.avg16_btanom(uunn(tt),:)));
      if length(ohoh) > 0
        ohoh
      end
      title(['signal (b) at latbin ' num2str(uumm(zz)) ' timestep ' num2str(uunn(tt))]); pause;
    end
  end
end

if iMaxTimeSteps > 0
  usethese = (1:iMaxTimeSteps);
else
  usethese = (abs(iMaxTimeSteps):length(a1.okdates));
end

if nargin == 2 & iYes < 0
  junk = [nansum(nansum(a1.co2(tropics,usethese)-a2.co2(tropics,usethese))) nansum(nansum(a1.stemp(tropics,usethese)-a2.stemp(tropics,usethese)))];
  fprintf(1,'nansum(diff co2 and stemp) = %8.6f %8.6f \n',junk)
  fprintf(1,'co2   at timestep iMaxTimeStepsM1, latbin 21 = %8.6f %8.6f \n',[a1.co2(21,iMaxTimeStepsM1) a2.co2(21,iMaxTimeStepsM1)])
  fprintf(1,'stemp at timestep iMaxTimeStepsM1, latbin 21 = %8.6f %8.6f \n',[a1.stemp(21,iMaxTimeStepsM1) a2.stemp(21,iMaxTimeStepsM1)])
  fprintf(1,'nloop at timestep iMaxTimeStepsM1, latbin 21 = %8.6f %8.6f \n',[a1.topts.nloop a2.topts.nloop])

  figure(1); plot(a1.okdates(usethese),smooth(nanmean(a1.stemp(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.stemp(tropics,usethese)),iSmooth),'r.-'); plotaxis2; title('stemp');
  figure(2); plot(a1.okdates(usethese),smooth(nanmean(a1.co2(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.co2(tropics,usethese)),iSmooth),'r.-'); plotaxis2; title('co2');
  figure(3); plot(a1.okdates(usethese),smooth(nanmean(a1.n2o(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.n2o(tropics,usethese)),iSmooth),'r.-'); plotaxis2; title('n2o');
  figure(4); plot(a1.okdates(usethese),smooth(nanmean(a1.ch4(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.ch4(tropics,usethese)),iSmooth),'r.-'); plotaxis2; title('ch4');
  figure(5); plot(a1.okdates(usethese),smooth(nanmean(a1.cfc11(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.cfc11(tropics,usethese)),iSmooth),'r.-'); plotaxis2; title('cfc11');
  figure(6); plot(a1.okdates(usethese),smooth(nanmean(a1.cfc12(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.cfc12(tropics,usethese)),iSmooth),'r.-'); plotaxis2; title('cfc12');

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.stemp(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.stemp(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed stemp rates = %8.6f %8.6f \n',P1(1),P2(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.co2(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.co2(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed co2 rates = %8.6f %8.6f \n',P1(1),P2(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.n2o(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.n2o(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed n2o rates = %8.6f %8.6f \n',P1(1),P2(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.ch4(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.ch4(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed ch4 rates = %8.6f %8.6f \n',P1(1),P2(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.cfc11(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.cfc11(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed cfc11 rates = %8.6f %8.6f \n',P1(1),P2(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.cfc12(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.cfc12(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed cfc12 rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.stemp(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.stemp(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed stemp rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.co2(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.co2(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed co2 rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.n2o(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.n2o(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed n2o rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.ch4(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.ch4(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed ch4 rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.cfc11(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.cfc11(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed cfc11 rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.cfc12(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.cfc12(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed cfc12 rates = %8.6f %8.6f \n',P1(1),P2(1))

  for ii = 1 : 4
    figure(ii); ax = axis; ax(1) = min(a1.okdates); ax(2) = max(a1.okdates); axis(ax); grid on;
  end

  figure(5); plot(a1.okdates(usethese),a1.stemp(tropics,usethese)-a2.stemp(tropics,usethese)); plotaxis2; title('obs-cal retr : difference in tropical stemp'); grid
  figure(6); plot(a1.okdates(usethese),a1.co2(tropics,usethese)-a2.co2(tropics,usethese)); plotaxis2; title('obs-cal retr : difference in tropical co2'); grid

  figure(5); plot(a1.okdates(usethese),smooth(mean(a1.stemp(tropics,usethese)-a2.stemp(tropics,usethese)),iSmooth),'b','linewidth',2);
    plotaxis2; title('obs-cal retr : difference in tropical stemp'); grid

    figure(6); plot(a1.okdates(usethese),smooth(mean(a1.co2(tropics,usethese)-a2.co2(tropics,usethese)),iSmooth),'b','linewidth',2)
    plotaxis2; title('obs-cal retr : difference in tropical co2'); grid

else
  junk = [nansum(nansum(a1.co2(tropics,usethese)-a2.co2(tropics,usethese))) nansum(nansum(a1.stemp(tropics,usethese)-a2.stemp(tropics,usethese)))];
  fprintf(1,'nansum(diff 1,2 co2 and stemp) = %8.6f %8.6f \n',junk(1:2))
  junk = [nansum(nansum(a1.co2(tropics,usethese)-a3.co2(tropics,usethese))) nansum(nansum(a1.stemp(tropics,usethese)-a3.stemp(tropics,usethese)))];
  fprintf(1,'nansum(diff 1,3 co2 and stemp) = %8.6f %8.6f \n',junk(1:2))
  junk = [nansum(nansum(a2.co2(tropics,usethese)-a3.co2(tropics,usethese))) nansum(nansum(a2.stemp(tropics,usethese)-a3.stemp(tropics,usethese)))];
  fprintf(1,'nansum(diff 2,3 co2 and stemp) = %8.6f %8.6f \n',junk(1:2))

  fprintf(1,'co2   at timestep iMaxTimeStepsM1, latbin 21 = %8.6f %8.6f %8.6f \n',[a1.co2(21,iMaxTimeStepsM1) a2.co2(21,iMaxTimeStepsM1) a3.co2(21,iMaxTimeStepsM1)])
  fprintf(1,'stemp at timestep iMaxTimeStepsM1, latbin 21 = %8.6f %8.6f %8.6f \n',[a1.stemp(21,iMaxTimeStepsM1) a2.stemp(21,iMaxTimeStepsM1) a3.stemp(21,iMaxTimeStepsM1)])
  %fprintf(1,'nloop at timestep iMaxTimeStepsM1, latbin 21 = %8.6f %8.6f %8.6f \n',[a1.topts.nloop a2.topts.nloop a3.topts.nloop])

  figure(1); plot(a1.okdates(usethese),smooth(nanmean(a1.stemp(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.stemp(tropics,usethese)),iSmooth),'g.-',...
                  a1.okdates(usethese),smooth(nanmean(a3.stemp(tropics,usethese)),iSmooth),'kx-'); plotaxis2; title('stemp : ignore black curve');
  figure(2); plot(a1.okdates(usethese),smooth(nanmean(a1.co2(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.co2(tropics,usethese)),iSmooth),'g.-',...
                a3.okdates(usethese),smooth(nanmean(a3.co2(tropics,usethese)),iSmooth),'r.-'); plotaxis2; title('co2');
  figure(3); plot(a1.okdates(usethese),smooth(nanmean(a1.ch4(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.ch4(tropics,usethese)),iSmooth),'g.-',...
                a3.okdates(usethese),smooth(nanmean(a3.ch4(tropics,usethese)),iSmooth),'r.-'); plotaxis2; title('ch4');
  figure(4); plot(a1.okdates(usethese),smooth(nanmean(a1.n2o(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.n2o(tropics,usethese)),iSmooth),'g.-',...
                a3.okdates(usethese),smooth(nanmean(a3.n2o(tropics,usethese)),iSmooth),'r.-'); plotaxis2; title('n2o');
  figure(5); plot(a1.okdates(usethese),smooth(nanmean(a1.cfc11(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.cfc11(tropics,usethese)),iSmooth),'g.-',...
                a3.okdates(usethese),smooth(nanmean(a3.cfc11(tropics,usethese)),iSmooth),'r.-'); plotaxis2; title('cfc11');
  figure(6); plot(a1.okdates(usethese),smooth(nanmean(a1.cfc12(tropics,usethese)),iSmooth),'b.-',a1.okdates(usethese),smooth(nanmean(a2.cfc12(tropics,usethese)),iSmooth),'g.-',...
                a3.okdates(usethese),smooth(nanmean(a3.cfc12(tropics,usethese)),iSmooth),'r.-'); plotaxis2; title('cfc12');

  P1 = polyfit(a2.okdates(usethese)',smooth(nanmean(a1.stemp(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates(usethese)',smooth(nanmean(a2.stemp(tropics,usethese)),iSmooth),1);
  P3 = polyfit(a3.okdates(usethese)',smooth(nanmean(a3.stemp(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed stemp rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates(usethese)',smooth(nanmean(a1.co2(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates(usethese)',smooth(nanmean(a2.co2(tropics,usethese)),iSmooth),1);
  P3 = polyfit(a3.okdates(usethese)',smooth(nanmean(a3.co2(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed co2 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates(usethese)',smooth(nanmean(a1.n2o(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates(usethese)',smooth(nanmean(a2.n2o(tropics,usethese)),iSmooth),1);
  P3 = polyfit(a3.okdates(usethese)',smooth(nanmean(a3.n2o(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed n2o rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates(usethese)',smooth(nanmean(a1.ch4(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates(usethese)',smooth(nanmean(a2.ch4(tropics,usethese)),iSmooth),1);
  P3 = polyfit(a3.okdates(usethese)',smooth(nanmean(a3.ch4(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed ch4 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates(usethese)',smooth(nanmean(a1.cfc11(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates(usethese)',smooth(nanmean(a2.cfc11(tropics,usethese)),iSmooth),1);
  P3 = polyfit(a3.okdates(usethese)',smooth(nanmean(a3.cfc11(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed cfc11 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates(usethese)',smooth(nanmean(a1.cfc12(tropics,usethese)),iSmooth),1);
  P2 = polyfit(a2.okdates(usethese)',smooth(nanmean(a2.cfc12(tropics,usethese)),iSmooth),1);
  P3 = polyfit(a3.okdates(usethese)',smooth(nanmean(a3.cfc12(tropics,usethese)),iSmooth),1);
  fprintf(1,'  method1 smoothed cfc12 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  %%%%%%%%%%%%%%%%%%%%%%%%%
  A.okdates = a1.okdates;
  A.stemp = smooth(nanmean(a3.stemp(tropics,usethese)),iSmooth);
  A.co2   = smooth(nanmean(a3.co2(tropics,usethese)),iSmooth);
  A.n2o   = smooth(nanmean(a3.n2o(tropics,usethese)),iSmooth);
  A.ch4   = smooth(nanmean(a3.ch4(tropics,usethese)),iSmooth);
  A.cfc11   = smooth(nanmean(a3.cfc11(tropics,usethese)),iSmooth);
  A.cfc12   = smooth(nanmean(a3.cfc12(tropics,usethese)),iSmooth);

  A.stemp_1 = smooth(nanmean(a1.stemp(tropics,usethese)),iSmooth);
  A.co2_1   = smooth(nanmean(a1.co2(tropics,usethese)),iSmooth);
  A.ch4_1   = smooth(nanmean(a1.ch4(tropics,usethese)),iSmooth);
  A.stemp_2 = smooth(nanmean(a2.stemp(tropics,usethese)),iSmooth);
  A.co2_2   = smooth(nanmean(a2.co2(tropics,usethese)),iSmooth);
  A.ch4_2   = smooth(nanmean(a2.ch4(tropics,usethese)),iSmooth);
  %%%%%%%%%%%%%%%%%%%%%%%%%

  ind = 3 : length(a1.okdates)-2;
  ind = 3 : length(a1.okdates(usethese))-2;

  gah = smooth(nanmean(a1.stemp(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a2.stemp(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a3.stemp(tropics,usethese)),iSmooth); P3 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  fprintf(1,'  Method 2 smoothed stemp rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  gah = smooth(nanmean(a1.co2(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a2.co2(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a3.co2(tropics,usethese)),iSmooth); P3 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  fprintf(1,'  Method 2 smoothed co2 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  gah = smooth(nanmean(a1.n2o(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a2.n2o(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a3.n2o(tropics,usethese)),iSmooth); P3 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  fprintf(1,'  Method 2 smoothed n2o rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  gah = smooth(nanmean(a1.ch4(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a2.ch4(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a3.ch4(tropics,usethese)),iSmooth); P3 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  fprintf(1,'  Method 2 smoothed ch4 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  gah = smooth(nanmean(a1.cfc11(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a2.cfc11(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a3.cfc11(tropics,usethese)),iSmooth); P3 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  fprintf(1,'  Method 2 smoothed cfc11 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  gah = smooth(nanmean(a1.cfc12(tropics,usethese)),iSmooth); P1 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a2.cfc12(tropics,usethese)),iSmooth); P2 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  gah = smooth(nanmean(a3.cfc12(tropics,usethese)),iSmooth); P3 = polyfit(a1.okdates(usethese(ind))',gah(ind),1);
  fprintf(1,'  Method 2 smoothed cfc12 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  for ii = 1 : 6
    figure(ii); ax = axis; ax(1) = min(a1.okdates(usethese)); ax(2) = max(a1.okdates(usethese)); axis(ax); grid on;
  end

  figure(7); plot(a1.okdates(usethese),a1.stemp(tropics(1),usethese)-a2.stemp(tropics(1),usethese),'b',...
                  a1.okdates(usethese),a1.stemp(tropics(1),usethese)-a3.stemp(tropics(1),usethese),'r',...
                  a1.okdates(usethese),a1.stemp(tropics,usethese)-a2.stemp(tropics,usethese),'b',a1.okdates(usethese),a1.stemp(tropics,usethese)-a3.stemp(tropics,usethese),'r'); 
    plotaxis2; title('obs-cal retr : difference in tropical stemp'); hl = legend('1-2','1-3','location','best'); grid

  figure(8); plot(a1.okdates(usethese),a1.co2(tropics(1),usethese)-a2.co2(tropics(1),usethese),'b',...
                  a1.okdates(usethese),a1.co2(tropics(1),usethese)-a3.co2(tropics(1),usethese),'r',...
                  a1.okdates(usethese),a1.co2(tropics,usethese)-a2.co2(tropics,usethese),'b',a1.okdates(usethese),a1.co2(tropics,usethese)-a3.co2(tropics,usethese),'r'); 
    plotaxis2; title('obs-cal retr : difference in tropical co2'); hl = legend('1-2','1-3','location','best'); grid

    figure(7); plot(a1.okdates(usethese),smooth(mean(a1.stemp(tropics,usethese)-a2.stemp(tropics,usethese)),iSmooth),'b',...
  		    a1.okdates(usethese),smooth(mean(a1.stemp(tropics,usethese)-a3.stemp(tropics,usethese)),iSmooth),'r','linewidth',2)
    plotaxis2; title('obs-cal retr : difference in tropical stemp'); hl = legend('1-2','1-3','location','best'); grid

    figure(8); plot(a1.okdates(usethese),smooth(mean(a1.co2(tropics,usethese)-a2.co2(tropics,usethese)),iSmooth),'b',...
   	            a1.okdates(usethese),smooth(mean(a1.co2(tropics,usethese)-a3.co2(tropics,usethese)),iSmooth),'r','linewidth',2)
    plotaxis2; title('obs-cal retr : difference in tropical co2'); hl = legend('1-2','1-3','location','best'); grid

    figure(8); plot(a1.okdates(usethese),smooth(mean(a1.co2(tropics,usethese)-a2.co2(tropics,usethese)),iSmooth),'r',...
   	            a1.okdates(usethese),smooth(mean(a1.co2(tropics,usethese)),iSmooth),'b',a1.okdates(usethese),smooth(mean(a2.co2(tropics,usethese)),iSmooth),'c','linewidth',2)
    plotaxis2; title('obs-cal retr : difference in tropical co2'); hl = legend('1-2','1','2','location','best'); grid

end

disp(' ')
