function compare_anomaly_runs(file1,file2,file3);

%{
(1) old and new sarta with finite diff strow jacs
f1a = 'ANOM_RESULTS_TEST1/anomaly_0dayavg_results_oldchans.mat';
f1b = 'ANOM_RESULTS_TEST1/anomaly_0dayavg_cal_results_oldchans.mat';
compare_anomaly_runs(f1a,f1b)
  Method 2 smoothed stemp rates = 0.015246 0.008036 0.007211
  Method 2 smoothed co2 rates = 1.692519 -0.225477 1.917996
  Method 2 smoothed n2o rates = 0.605222 -0.184678 0.789901
  Method 2 smoothed ch4 rates = 3.903583 -1.076112 4.979695

f2a = 'ANOM_RESULTS_TEST1/anomaly_0dayavg_results_newchans.mat';
f2b = 'ANOM_RESULTS_TEST1/anomaly_0dayavg_cal_results_newchans.mat';
compare_anomaly_runs(f2a,f2b)
  Method 2 smoothed stemp rates = 0.015419 0.009667 0.005751
  Method 2 smoothed co2 rates = 1.703375 -0.249410 1.952784
  Method 2 smoothed n2o rates = 0.591044 -0.208141 0.799185
  Method 2 smoothed ch4 rates = 3.999539 -1.203400 5.202939
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
  error('need at least 2 input files, at most 3')
end

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

[mmbad,nnbad] = find(abs(a1.stemp) > 5 | abs(a2.stemp) > 5 |abs(a3.stemp) > 5 | abs(a1.co2 > 50) | abs(a2.co2 > 50) | abs(a3.co2 > 50)); 
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

if nargin == 2 & iYes < 0
  junk = [sum(sum(a1.co2(tropics,:)-a2.co2(tropics,:))) sum(sum(a1.stemp(tropics,:)-a2.stemp(tropics,:)))];
  fprintf(1,'sum(diff co2 and stemp) = %8.6f %8.6f \n',junk)
  fprintf(1,'co2   at timestep 364, latbin 21 = %8.6f %8.6f \n',[a1.co2(21,364) a2.co2(21,364)])
  fprintf(1,'stemp at timestep 364, latbin 21 = %8.6f %8.6f \n',[a1.stemp(21,364) a2.stemp(21,364)])
  fprintf(1,'nloop at timestep 364, latbin 21 = %8.6f %8.6f \n',[a1.topts.nloop a2.topts.nloop])

  figure(1); plot(a1.okdates,smooth(nanmean(a1.stemp(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.stemp(tropics,:)),iSmooth),'r.-'); title('stemp');
  figure(2); plot(a1.okdates,smooth(nanmean(a1.co2(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.co2(tropics,:)),iSmooth),'r.-'); title('co2');
  figure(3); plot(a1.okdates,smooth(nanmean(a1.n2o(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.n2o(tropics,:)),iSmooth),'r.-'); title('n2o');
  figure(4); plot(a1.okdates,smooth(nanmean(a1.ch4(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.ch4(tropics,:)),iSmooth),'r.-'); title('ch4');
  figure(5); plot(a1.okdates,smooth(nanmean(a1.cfc11(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.cfc11(tropics,:)),iSmooth),'r.-'); title('cfc11');
  figure(6); plot(a1.okdates,smooth(nanmean(a1.cfc12(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.cfc12(tropics,:)),iSmooth),'r.-'); title('cfc12');

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.stemp(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.stemp(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed stemp rates = %8.6f %8.6f \n',P1(1),P2(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.co2(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.co2(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed co2 rates = %8.6f %8.6f \n',P1(1),P2(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.n2o(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.n2o(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed n2o rates = %8.6f %8.6f \n',P1(1),P2(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.ch4(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.ch4(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed ch4 rates = %8.6f %8.6f \n',P1(1),P2(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.cfc11(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.cfc11(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed cfc11 rates = %8.6f %8.6f \n',P1(1),P2(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.cfc12(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.cfc12(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed cfc12 rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.stemp(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.stemp(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed stemp rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.co2(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.co2(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed co2 rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.n2o(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.n2o(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed n2o rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.ch4(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.ch4(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed ch4 rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.cfc11(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.cfc11(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed cfc11 rates = %8.6f %8.6f \n',P1(1),P2(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.cfc12(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.cfc12(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed cfc12 rates = %8.6f %8.6f \n',P1(1),P2(1))

  for ii = 1 : 4
    figure(ii); ax = axis; ax(1) = min(a1.okdates); ax(2) = max(a1.okdates); axis(ax); grid on;
  end

  figure(5); plot(a1.okdates,a1.stemp(tropics,:)-a2.stemp(tropics,:)); title('difference in tropical stemp'); grid
  figure(6); plot(a1.okdates,a1.co2(tropics,:)-a2.co2(tropics,:)); title('difference in tropical co2'); grid

  figure(5); plot(a1.okdates,smooth(mean(a1.stemp(tropics,:)-a2.stemp(tropics,:)),iSmooth),'b','linewidth',2);
    title('difference in tropical stemp'); grid

    figure(6); plot(a1.okdates,smooth(mean(a1.co2(tropics,:)-a2.co2(tropics,:)),iSmooth),'b','linewidth',2)
    title('difference in tropical co2'); grid

else
  junk = [sum(sum(a1.co2(tropics,:)-a2.co2(tropics,:))) sum(sum(a1.stemp(tropics,:)-a2.stemp(tropics,:)))];
  fprintf(1,'sum(diff 1,2 co2 and stemp) = %8.6f %8.6f \n',junk(1:2))
  junk = [sum(sum(a1.co2(tropics,:)-a3.co2(tropics,:))) sum(sum(a1.stemp(tropics,:)-a3.stemp(tropics,:)))];
  fprintf(1,'sum(diff 1,3 co2 and stemp) = %8.6f %8.6f \n',junk(1:2))
  junk = [sum(sum(a2.co2(tropics,:)-a3.co2(tropics,:))) sum(sum(a2.stemp(tropics,:)-a3.stemp(tropics,:)))];
  fprintf(1,'sum(diff 2,3 co2 and stemp) = %8.6f %8.6f \n',junk(1:2))

  fprintf(1,'co2   at timestep 364, latbin 21 = %8.6f %8.6f %8.6f \n',[a1.co2(21,364) a2.co2(21,364) a3.co2(21,364)])
  fprintf(1,'stemp at timestep 364, latbin 21 = %8.6f %8.6f %8.6f \n',[a1.stemp(21,364) a2.stemp(21,364) a3.stemp(21,364)])
  fprintf(1,'nloop at timestep 364, latbin 21 = %8.6f %8.6f %8.6f \n',[a1.topts.nloop a2.topts.nloop a3.topts.nloop])

  figure(1); plot(a1.okdates,smooth(nanmean(a1.stemp(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.stemp(tropics,:)),iSmooth),'g.-',...
                  a1.okdates,smooth(nanmean(a3.stemp(tropics,:)),iSmooth),'kx-'); title('stemp : ignore black curve');
  figure(2); plot(a1.okdates,smooth(nanmean(a1.co2(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.co2(tropics,:)),iSmooth),'g.-',...
                a3.okdates,smooth(nanmean(a3.co2(tropics,:)),iSmooth),'r.-'); title('co2');
  figure(3); plot(a1.okdates,smooth(nanmean(a1.ch4(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.ch4(tropics,:)),iSmooth),'g.-',...
                a3.okdates,smooth(nanmean(a3.ch4(tropics,:)),iSmooth),'r.-'); title('ch4');
  figure(4); plot(a1.okdates,smooth(nanmean(a1.n2o(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.n2o(tropics,:)),iSmooth),'g.-',...
                a3.okdates,smooth(nanmean(a3.n2o(tropics,:)),iSmooth),'r.-'); title('n2o');
  figure(5); plot(a1.okdates,smooth(nanmean(a1.cfc11(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.cfc11(tropics,:)),iSmooth),'g.-',...
                a3.okdates,smooth(nanmean(a3.cfc11(tropics,:)),iSmooth),'r.-'); title('cfc11');
  figure(6); plot(a1.okdates,smooth(nanmean(a1.cfc12(tropics,:)),iSmooth),'b.-',a1.okdates,smooth(nanmean(a2.cfc12(tropics,:)),iSmooth),'g.-',...
                a3.okdates,smooth(nanmean(a3.cfc12(tropics,:)),iSmooth),'r.-'); title('cfc12');

  P1 = polyfit(a2.okdates',smooth(nanmean(a1.stemp(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.stemp(tropics,:)),iSmooth),1);
  P3 = polyfit(a3.okdates',smooth(nanmean(a3.stemp(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed stemp rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.co2(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.co2(tropics,:)),iSmooth),1);
  P3 = polyfit(a3.okdates',smooth(nanmean(a3.co2(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed co2 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.n2o(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.n2o(tropics,:)),iSmooth),1);
  P3 = polyfit(a3.okdates',smooth(nanmean(a3.n2o(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed n2o rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.ch4(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.ch4(tropics,:)),iSmooth),1);
  P3 = polyfit(a3.okdates',smooth(nanmean(a3.ch4(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed ch4 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.cfc11(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.cfc11(tropics,:)),iSmooth),1);
  P3 = polyfit(a3.okdates',smooth(nanmean(a3.cfc11(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed cfc11 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  P1 = polyfit(a1.okdates',smooth(nanmean(a1.cfc12(tropics,:)),iSmooth),1);
  P2 = polyfit(a2.okdates',smooth(nanmean(a2.cfc12(tropics,:)),iSmooth),1);
  P3 = polyfit(a3.okdates',smooth(nanmean(a3.cfc12(tropics,:)),iSmooth),1);
  fprintf(1,'  method1 smoothed cfc12 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.stemp(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.stemp(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a3.stemp(tropics,:)),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed stemp rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.co2(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.co2(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a3.co2(tropics,:)),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed co2 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.n2o(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.n2o(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a3.n2o(tropics,:)),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed n2o rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.ch4(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.ch4(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a3.ch4(tropics,:)),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed ch4 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.cfc11(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.cfc11(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a3.cfc11(tropics,:)),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed cfc11 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  ind = 3 : length(a1.okdates)-2;
  gah = smooth(nanmean(a1.cfc12(tropics,:)),iSmooth); P1 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a2.cfc12(tropics,:)),iSmooth); P2 = polyfit(a1.okdates(ind)',gah(ind),1);
  gah = smooth(nanmean(a3.cfc12(tropics,:)),iSmooth); P3 = polyfit(a1.okdates(ind)',gah(ind),1);
  fprintf(1,'  Method 2 smoothed cfc12 rates = %8.6f %8.6f %8.6f \n',P1(1),P2(1),P3(1))

  for ii = 1 : 6
    figure(ii); ax = axis; ax(1) = min(a1.okdates); ax(2) = max(a1.okdates); axis(ax); grid on;
  end

  figure(7); plot(a1.okdates,a1.stemp(tropics(1),:)-a2.stemp(tropics(1),:),'b',...
                  a1.okdates,a1.stemp(tropics(1),:)-a3.stemp(tropics(1),:),'r',...
                  a1.okdates,a1.stemp(tropics,:)-a2.stemp(tropics,:),'b',a1.okdates,a1.stemp(tropics,:)-a3.stemp(tropics,:),'r'); 
    title('difference in tropical stemp'); hl = legend('1-2','1-3','location','best'); grid

  figure(8); plot(a1.okdates,a1.co2(tropics(1),:)-a2.co2(tropics(1),:),'b',...
                  a1.okdates,a1.co2(tropics(1),:)-a3.co2(tropics(1),:),'r',...
                  a1.okdates,a1.co2(tropics,:)-a2.co2(tropics,:),'b',a1.okdates,a1.co2(tropics,:)-a3.co2(tropics,:),'r'); 
    title('difference in tropical co2'); hl = legend('1-2','1-3','location','best'); grid

    figure(7); plot(a1.okdates,smooth(mean(a1.stemp(tropics,:)-a2.stemp(tropics,:)),iSmooth),'b',...
  		    a1.okdates,smooth(mean(a1.stemp(tropics,:)-a3.stemp(tropics,:)),iSmooth),'r','linewidth',2)
    title('difference in tropical stemp'); hl = legend('1-2','1-3','location','best'); grid

    figure(8); plot(a1.okdates,smooth(mean(a1.co2(tropics,:)-a2.co2(tropics,:)),iSmooth),'b',...
   	            a1.okdates,smooth(mean(a1.co2(tropics,:)-a3.co2(tropics,:)),iSmooth),'r','linewidth',2)
    title('difference in tropical co2'); hl = legend('1-2','1-3','location','best'); grid

end

disp(' ')
