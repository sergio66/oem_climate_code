addpath0

% addpath /home/sergio/MATLABCODE
% addpath /home/sergio/MATLABCODE/TROPOPAUSE
% addpath /home/sergio/MATLABCODE/COLORMAP
% addpath /home/sergio/MATLABCODE/COLORMAP/LLS
% addpath /home/sergio/MATLABCODE/PLOTTER
% addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
% addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD
% addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
% addpath /home/sergio/MATLABCODE/matlib/science/
% addpath /home/sergio/MATLABCODE/NANROUTINES/
% addpath /home/sergio/MATLABCODE/SHOWSTATS/
% addpath /asl/matlib/aslutil
% addpath /asl/matlib/h4tools
% addpath /asl/matlib/maps
% addpath /home/sergio/MATLABCODE/TIME

addpath /home/sergio/git/matlabcode/COLORMAP/COLORBREWER/cbrewer/cbrewer/
addpath /home/sergio/git/matlabcode/NANROUTINES/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%^

cmap = cbrewer('div', 'RdBu', 128);
cmap = cbrewer('div', 'RdBu', 16);
cmap(cmap < 0) = 0; cmap(cmap > 1) = 1;
cmap = flipud(cmap);
jett = jet(128); jett(1,:) = 1;

f2645 = instr_chans2645;
i0900 = find(f2645 >= 900,1);
i0723 = find(f2645 >= 723,1);
i1305 = find(f2645 >= 1305,1);
i1419 = find(f2645 >= 1419,1);

set_anomaly_info    %% also in clust_run_retrieval_setlatbin_AIRS_loop_lonbin.m

if strfind(anomalydatafile,'globalavg_and_tropics') | strfind(anomalydatafile,'globalavg_and_TWPlat35') 
  iStartOffset = 3;  % first two are global and tropic
else
  iStartOffset = 2;
end
iStartOffset = input('Enter iStartOffset (1 for global, sometimes 2 for TRP, 3-24 for latbins, sometimes 2-65 for latbins) : ');

disp('  ')
disp('make sure you do this before starting Matlab, if you want to run ecRad!!! module load netCDF-Fortran/4.4.4-intel-2018b');
disp('make sure you do this before starting Matlab, if you want to run ecRad!!! module load netCDF-Fortran/4.4.4-intel-2018b');
disp('make sure you do this before starting Matlab, if you want to run ecRad!!! module load netCDF-Fortran/4.4.4-intel-2018b');
disp('  ')

load llsmap5
%if length(llsmap5) == 64
%  llsmap5 = llsmap5(2:end,:);
%end

llsmap5NAN = [[0.0 0.0 0.0]; llsmap5];
llsmap5_0 = llsmap5;
llsmap5   = llsmap5NAN;

disp(' ')
%disp('saved a version as   save -v7.3 /asl/s1/sergio/JUNK/gather_tileCLRday.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_startwithERA5trends.mat')
%disp('saved a version as   save -v7.3 /asl/s1/sergio/JUNK/gather_tileCLRday.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_startwithERA5trends.mat')
%disp('saved a version as   save -v7.3 /asl/s1/sergio/JUNK/gather_tileCLRday.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight.mat or /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_startwithERA5trends.mat')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%iAK = input('do AvgKernels ???? (-1/+1 = yes = default) ');
iAK = -1;
if length(iAK) == 0
  iAK = +1;
end
if iAK > 0
  %% see Plotutils/plot_retrieval_latbins_fewlays
  iNumYears = input('  Enter number of years from 2002-X (so we can load in appropriate ERA5 trend file !!!! [default = 20] ');  
  if length(iNumYears) == 0
    iNumYears = 20;
  end
  junk = ['../FIND_NWP_MODEL_TRENDS/MEAN_PROFILES/ERA5_atm_data_2002_09_to_' num2str(2002 + iNumYears) '_08_trends_desc.mat'];
  junk = ['../FIND_NWP_MODEL_TRENDS/MEAN_PROFILES/ERA5_atm_N_cld_data_2002_09_to_' num2str(2002 + iNumYears) '_08_trends_desc_surf.mat'];  
  junk = ['../FIND_NWP_MODEL_TRENDS/MEAN_PROFILES/ERA5_atm_N_cld_data_2002_09_to_' num2str(2002 + iNumYears) '_08_trends_desc.mat'];
  era5 = load(junk);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%iJunk = input('clear 50 figs???? (-1 = no = default) ');
%if length(iJunk) == 0
%  iJunk = -1;
%end
iJunk = +1;
if iJunk > 0
  for ii = 1 : 30
    figure(ii); clf; colormap jet
  end
end

iNumLay = 20;
iNumLay = input('Enter number of expected fat layers (20 [5 thick AIRS layers], having good luck with (default 49 [2 thick AIRS layers]) : ');
if length(iNumLay) == 0
  iNumLay = 20;
  iNumLay = 49;
end
iNavgLay = floor(100/iNumLay);
iWVind = 07:26; iWVind = (1:iNumLay)+6+0*iNumLay;
iTind  = 27:46; iTind  = (1:iNumLay)+6+1*iNumLay;
iO3ind = 47:66; iO3ind = (1:iNumLay)+6+2*iNumLay;

iNorD = input('night (+1) or day (-1)  [+1 = night = default] ');
if length(iNorD) == 0
  iNorD = +1;
end
iDorA = iNorD;

iOCBset = input('obs cal or bias (0,+1,-1) [default 0] : ');
if length(iOCBset) == 0
  iOCBset = 0;
end

data_anom = load(anomalydatafile);
if strfind(anomalydatafile,'_tile')
  rlat = [data_anom.LatBin];
else
  rlat = [0; meanvaluebin(data_anom.newLatGrid)];
  rlat = [meanvaluebin(data_anom.newLatGrid)];
end

daysSince2002 = change2days(data_anom.yy,data_anom.mm,data_anom.dd,2002);
yymm = data_anom.yy + (data_anom.mm-1)/12 + (data_anom.dd-1)/30/12;

iNumAnomData = iNumAnomTimeSteps * iNumAnomTiles;

fnamelastloaded = 'none';
if ~exist('iaFound')
  clear results*
  iaFound = zeros(1,iNumAnomData);
  existfname = zeros(1,iNumAnomData);

  save_cov_set.cov_set = nan(13,iNumAnomData);
  save_cov_set.fmat = nan(6,iNumAnomData);

  cdofs = nan(iNumAnomData,66);
  lencdofs = nan(1,iNumAnomData);
  thedofs = nan(1,iNumAnomData);

  results = nan(iNumAnomData,6);
  resultsWV = nan(iNumAnomData,iNumLay);
  resultsT  = nan(iNumAnomData,iNumLay);
  resultsO3 = nan(iNumAnomData,iNumLay);

  resultsunc = nan(iNumAnomData,6);
  resultsWVunc = nan(iNumAnomData,iNumLay);
  resultsTunc  = nan(iNumAnomData,iNumLay);
  resultsO3unc = nan(iNumAnomData,iNumLay);

  spectral_deltan00 = nan(2645,iNumAnomData);
  rates             = nan(2645,iNumAnomData);
  fits              = nan(2645,iNumAnomData);
  componentfits     = nan(5,2645,iNumAnomData);
  nedt              = nan(2645,iNumAnomData);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% can do something like      watch "ls -lt OutputAnomaly/Quantile03/*.mat | wc -l"   
%% can do something like      watch "ls -lt OutputAnomaly/Quantile03/*.mat | wc -l"   
%% can do something like      watch "ls -lt OutputAnomaly/Quantile03/*.mat | wc -l"   

JOB = -1;
%get_anomaly_processors  %% this sets mapAnomData_to_processor
[jobjunk,mapAnomData_to_processor] = get_anomaly_processors(ia_OorC_DataSet_Quantile,iNumAnomTimeSteps,iNumAnomTiles,iNumAnomJobsPerProc,anomalydatafile,-1);

driver.i16daytimestep = 2;
topts.ocb_set = 2;
driver.NorD = iNorD;
iQuantile = 1;

iaaFound = zeros(iNumAnomTimeSteps,iNumAnomTiles);
iCnt = 0;
for jj = 1 : (iNumAnomTiles)
  for ii = 1 : iNumAnomTimeSteps
    iCnt = iCnt + 1;
    iInd = iCnt;
    set_anom_outfilename    
    fname = [zanom_outdir '/Quantile' num2str(iQuantile,'%02d') '/test' num2str(iCnt) '.mat']; %% stored here before before July 2021, and fornew test comparisons
    if exist(fname) > 0
      iaaFound(ii,jj) = +1;
    end
  end
end
fprintf(1,'found %6i of the expected %6i files \n',sum(iaaFound(:)),iNumAnomTimeSteps*(iNumAnomTiles))
printarray([1:iNumAnomTiles; sum(iaaFound,1)]',['checking how many of the ' num2str(iNumAnomTimeSteps) ' have been made for the ' num2str(iNumAnomTiles) ' tiles']);
figure(1); clf; imagesc(iaaFound'); colorbar; title('Jobs found'); xlabel('iNumAnomTimeSteps'); ylabel('iNumAnomTiles'); 
disp('ret to continue to reading in the files'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% iWarning = 0;
% clear fname
% 
% iCnt = 0;
% for ii = 1 : iNumAnomTimeSteps : iNumAnomTimeSteps * iNumAnomTiles
%   iInd = ii;
%   set_anom_outfilename    
%   fname = [zanom_outdir '/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021, and fornew test comparisons
% 
%   if exist(fname) > 0
%     iCnt = iCnt + 1;
%     junkdir = dir(fname);
%     fprintf(1,'%s %s \n',[junkdir.folder     fname],junkdir.date);
%   else
%     fprintf(1,'%s DNE \n',fname);  
%   end
% end
% fprintf(1,'found %4i of %4i (subset) files \n',iCnt,iNumAnomTiles);
% if iCnt == 0
%   fprintf(1,'with iOCBset = %2i dataset = %2i iNorD = %2i seem to have found nothing nada zilch in %s \n',iOCBset,dataset,iNorD,fname )
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iJunk = input('correct dates/names etc etc???  Proceed or quit (+1 default/-1) : ');
if length(iJunk) == 0
  iJunk = +1;
end
if iJunk < 0
  return
end

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN = plevs(1:100)-plevs(2:101);
playsD = log(plevs(1:100)./plevs(2:101));
plays = playsN./playsD;
plays = flipud(plays);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear pavg

iLoopXY = +1; %% loop over lonbins, plot zonal avg once after all 64 lats,72 lons tried to be read
iLoopXY = -1; %% loop over latbins, plot zonal avg each time after all 64 lats tried to be read, so 72 plots DEFAULT

iDoAgain = +1;
iWarning = 0;
if iLoopXY > 0
  anomaly_read_loop_latbins
else
  anomaly_read_loop_lonbins
end

a = load(fnamelastloaded);
fprintf(1,'last loaded file %s has xb(1:6) = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f \n',fnamelastloaded,a.oem.xb(1:6))

bah = nanmean(reshape(rtime,iNumAnomTimeSteps,iNumAnomTiles),2);
[yyy,mmm,ddd] = tai2utcSergio(bah);
yymmx = yyy + (mmm-1)/12 + (ddd-1)/12/30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk_globalavg_rawdata = data_anom.btavgAnomFinal(:,1:iNumAnomTimeSteps);
junk_globalavg_spectral_deltan00 = spectral_deltan00(:,1:iNumAnomTimeSteps);
for ii = 1 : 2645
  P = nanpolyfit(1:iNumAnomTimeSteps,junk_globalavg_rawdata(ii,:),1);
  trend_globalavg_raw(ii) = P(1)*365/16;
  P = nanpolyfit(1:iNumAnomTimeSteps,junk_globalavg_spectral_deltan00(ii,:),1);
  trend_globalavg(ii) = P(1)*365/16;
end

figure(2); clf; plot(f2645,trend_globalavg_raw,'b',f2645,trend_globalavg,'r'); title('Global Avg Quick dBT/dt'); plotaxis2; xlim([645 1620])
  hl = legend('Raw from data','after CO2/CH4/N2O removed','location','best','fontsize',10);

figure(3); clf; plot(f2645,trend_globalavg_raw*20,'b',f2645,trend_globalavg*20,'r',...
                      f2645,data_anom.btavgAnomFinal(:,1),f2645,data_anom.btavgAnomFinal(:,iNumAnomTimeSteps)); title('Global Avg Quick dBT/dt'); plotaxis2; xlim([645 1620])
  hl = legend('20 x (Raw data)','20 x (after CO2/CH4/N2O removed)','TimeStep 1','TimeStep Nyears','location','best','fontsize',10);

figure(4); clf; 
  plot(yymm,smooth(junk_globalavg_spectral_deltan00(i0723,:),23),'r.-',...
       yymm,smooth(junk_globalavg_spectral_deltan00(i0900,:),23),'gx-',...
       yymm,smooth(junk_globalavg_spectral_deltan00(i1305,:),23),'k.-',...
       yymm,smooth(junk_globalavg_spectral_deltan00(i1419,:),23),'bo-',...
       yymm,smooth(junk_globalavg_rawdata(i0723,:),23),'r--',...
       yymm,smooth(junk_globalavg_rawdata(i0900,:),23),'y--',...
       yymm,smooth(junk_globalavg_rawdata(i1305,:),23),'k--',...
       yymm,smooth(junk_globalavg_rawdata(i1419,:),23),'c--',...
       'linewidth',2)
plotaxis2; hl = legend('723 CO2 (700 mb)','900 Window','1305 CH4','1419 WV (300 mb)','location','best','fontsize',8); xlim([2002 2025]);
title('Fitted (obs-cal) data in thick \newline raw data in dashes')

figure(5); clf
  plot(yymm,save_cov_set.xb_trace(1:6,1:iNumAnomTimeSteps),'linewidth',2); hl = legend('CO2','N2O','CH4','CFC11','CFC12','ST','location','best');
  set(gca,'fontsize',10); title('xb(trace gas + ST)');
figure(5); clf
  plot(yymm,save_cov_set.xf_trace(1:6,1:iNumAnomTimeSteps),'linewidth',2); hl = legend('CO2','N2O','CH4','CFC11','CFC12','ST','location','best');
  set(gca,'fontsize',10); title('xfinal(trace gas + ST)');
figure(5); clf
  wah = save_cov_set.xf_trace(1:6,1:iNumAnomTimeSteps);
  plot(yymm,wah(1,:),'b',yymm,wah(2,:),'m',yymm,wah(3,:),'g',yymm,wah(4,:),'y',yymm,wah(5,:),'k',yymm,wah(6,:),'r','linewidth',2); hold on; 
  wah = save_cov_set.xb_trace(1:6,1:iNumAnomTimeSteps);
  plot(yymm,wah(1,:),'b--',yymm,wah(2,:),'m--',yymm,wah(3,:),'g--',yymm,wah(4,:),'y--',yymm,wah(5,:),'k--',yymm,wah(6,:),'r--','linewidth',2); hold on; 
  hl = legend('CO2','N2O','CH4','CFC11','CFC12','ST','location','best');
  set(gca,'fontsize',10); title('xb(trace gas + ST)');

figure(6); clf
  plot(yymm,save_cov_set.xb_wvz(:,1:iNumAnomTimeSteps),'linewidth',2); set(gca,'fontsize',10); title('xb(WV)');
figure(6); clf
  plot(yymm,save_cov_set.xb_tzz(:,1:iNumAnomTimeSteps),'linewidth',2); set(gca,'fontsize',10); title('xb(TZ)');
figure(6); clf
  plot(yymm,save_cov_set.xb_ozz(:,1:iNumAnomTimeSteps),'linewidth',2); set(gca,'fontsize',10); title('xb(OZ)');

figure(7); clf; 
  rates_global = rates(:,1:iNumAnomTimeSteps);
  fits_global  = fits(:,1:iNumAnomTimeSteps);
  plot(f2645,nanmean(rates_global'-fits_global'),f2645,nanmean(rates_global')); plotaxis2; xlim([645 1620])
    hl = legend('rates-fits','rates','location','best','fontsize',10); title('Global set')

rates_global_sum = zeros(size(rates_global));
fits_global_sum  = zeros(size(rates_global));
for ii = iStartOffset : iNumAnomTiles
  ind = (1:iNumAnomTimeSteps) + (ii-1)*iNumAnomTimeSteps;
  coslat = cos(rlat(ii-iStartOffset+1)*pi/180);
  rates_global_sum = coslat * rates(:,ind) + rates_global_sum;
  fits_global_sum  = coslat * fits(:,ind) + fits_global_sum;
end
rates_global_sum = rates_global_sum/sum(cos(rlat*pi/180));
fits_global_sum = fits_global_sum/sum(cos(rlat*pi/180));
figure(8); clf
  plot(f2645,nanmean(rates_global_sum'-fits_global_sum'),f2645,nanmean(rates_global_sum')); plotaxis2; xlim([645 1620])
    hl = legend('rates-fits','rates','location','best','fontsize',10); title('Sum over latbins')
  
figure(9); clf
  plot(f2645,nanmean(rates_global'-rates_global_sum'),f2645,nanmean(rates_global')); plotaxis2; xlim([645 1620])
    hl = legend('rates- rates sum','rates','location','best','fontsize',10); title('Sum over latbins')

disp('ret to continue to plot_anomalies_ALL'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

if strfind(anomalydatafile,'_tile')
  plot_anomalies_1
else
  simple = -1;
  plot_anomalies_All
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ../FIND_NWP_MODEL_TRENDS/
look_at_anomalies_computeERA5_monthly_trends_desc_or_asc_64fast(40,1);

figure(46); ax = axis; figure(36); axis(ax)  %% WV 10-25 km, 2022 - 2025
figure(42); ax = axis; figure(35); axis(ax)  %% TZ 10-25 km, 2022 - 2025

