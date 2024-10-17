addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TROPOPAUSE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/matlib/science/
addpath /home/sergio/MATLABCODE/NANROUTINES/
addpath /home/sergio/MATLABCODE/SHOWSTATS/
addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /asl/matlib/maps
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/COLORMAP/COLORBREWER/cbrewer/cbrewer/

cmap = cbrewer('div', 'RdBu', 128);
cmap = cbrewer('div', 'RdBu', 16);
cmap(cmap < 0) = 0; cmap(cmap > 1) = 1;
cmap = flipud(cmap);
jett = jet(128); jett(1,:) = 1;

f2645 = instr_chans2645;
i0900 = find(f2645 >= 900,1);
i0723 = find(f2645 >= 723,1);
i0667 = find(f2645 >= 667.2,1);
i1305 = find(f2645 >= 1305,1);
i1419 = find(f2645 >= 1419,1);
i1400 = find(f2645 >= 1400,1);
i1598 = find(f2645 >= 1598,1);

set_anomaly_info

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
  junk = ['../FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_' num2str(2002 + iNumYears) '_08_trends_desc.mat'];
  junk = ['../FIND_NWP_MODEL_TRENDS/ERA5_atm_N_cld_data_2002_09_to_' num2str(2002 + iNumYears) '_08_trends_desc.mat'];
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

disp('quants for dataset 1-8        = [0 0.01 0.02 0.03 0.04 0.05 0.10 0.25 0.50 0.75 0.9 0.95 0.96 0.97 0.98 0.99 1.00]');
disp('quants for dataset 9,10,11,12 = [0.50 0.80 0.90 0.95 0.97 1.00]');

dataset = 9;
if length(dataset) == 0
  dataset = 9;
end

if dataset == 1 | dataset == -1
  iNumYears = 18;
elseif dataset == 2
  iNumYears = 19;
elseif abs(dataset) == 3
  iNumYears = 19;
elseif dataset == 4
  iNumYears = 19;
elseif dataset == 5
  iNumYears = 12;
elseif dataset == 6
  iNumYears = 07;
elseif dataset == 7
  iNumYears = 20;
elseif dataset == 8
  iNumYears = 07;
elseif dataset == 9
  iNumYears = 20;
elseif dataset == 10
  iNumYears = 05;
elseif dataset == 11
  iNumYears = 10;
elseif dataset == 12
  iNumYears = 15;
end

if dataset == -3
  iQuantile = 00;
elseif dataset >= 9
  iQuantile = 05; 
  if iOCBset == 0
    iQuantile = input('Dataset = 9,10,11,12, iOCBset = 0 (obs)  ==> Which quantile 1..5   [3 = Default] : ');
    if length(iQuantile) == 0
      iQuantile = 3;
    end
  elseif iOCBset == 1
    iQuantile = input('Dataset = 9,10,11,12, iOCBset = 1 (cal)  ==> Which quantile 1..16   [16 = Default] : ');
    if length(iQuantile) == 0
      iQuantile = 16;
    end
  else
    iOCBset
    error('huh??? iOCBset = 0,1 only!!!')
  end
elseif dataset ~= 3 & dataset < 9
  iQuantile = 16;  %% AIRS STM 2021, hottest
  iQuantile = 08;  %% 
  iQuantile = input('Which quantile -1 for extremes [(1--16) (99 for orig Q16, done for AIRS STM)]   [16 = Default] : ');
  if length(iQuantile) == 0
    iQuantile = 16;
  end
else
  iQuantile = [];
end

data_anom = load(anomalydatafile);
if strfind(anomalydatafile,'_tile')
  rlat = [data_anom.LatBin];
else
  rlat = [0 meanvaluebin(data_anom.newLatGrid)];
  rlat = [meanvaluebin(data_anom.newLatGrid)];
end

daysSince2002 = change2days(data_anom.yy,data_anom.mm,data_anom.dd,2002);
yymm = data_anom.yy + (data_anom.mm-1)/12 + (data_anom.dd-1)/30/12;

iNumAnomTiles = 1;  %% global only
iNumAnomTiles = 2;  %% global + tropics only
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

iaaFound = zeros(iNumAnomTimeSteps,iNumAnomTiles);
iCnt = 0;
for jj = 1 : iNumAnomTiles
  for ii = 1 : iNumAnomTimeSteps
    iCnt = iCnt + 1;
    if iNorD > 0
      fname = ['OutputAnomaly/Quantile' num2str(iQuantile,'%02d') '/test' num2str(iCnt) '.mat']; %% stored here before before July 2021, and fornew test comparisons
    elseif iNorD < 0
      fname = ['OutputAnomaly_Day/Quantile' num2str(iQuantile,'%02d') '/test' num2str(iCnt) '.mat']; %% stored here before before July 2021
    end
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

iWarning = 0;
clear fname

iCnt = 0;
for ii = 1 : iNumAnomTimeSteps : iNumAnomTimeSteps * iNumAnomTiles
  if iNorD > 0
    fname = ['OutputAnomaly/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021, and fornew test comparisons
  elseif iNorD < 0
    fname = ['OutputAnomaly_Day/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021
  end

  if exist(fname) > 0
    iCnt = iCnt + 1;
    junkdir = dir(fname);
    fprintf(1,'%s %s \n',[junkdir.folder     fname],junkdir.date);
  else
    fprintf(1,'%s DNE \n',fname);  
  end
end
fprintf(1,'found %4i of %4i (subset) files \n',iCnt,iNumAnomTiles);
if iCnt == 0
  fprintf(1,'with iOCBset = %2i dataset = %2i iNorD = %2i seem to have found nothing nada zilch in %s \n',iOCBset,dataset,iNorD,fname )
end

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

iDoAgain = +1;
fprintf(1,'will loop through %3i timeteps for %2i latbins .. the timesteps will mark off as "o" for 100 and "." for 10 \n',iNumAnomTimeSteps,iNumAnomTiles)

while iDoAgain > 0
  for ii = 1 : iNumAnomTimeSteps * iNumAnomTiles
    if iNorD > 0
      fname = ['OutputAnomaly/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021, and fornew test comparisons
    elseif iNorD < 0
      fname = ['OutputAnomaly_Day/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021
    end
  
    existfname(ii) = exist(fname);
    i10sec = -1;
    if exist(fname)
      %% datenum : A serial date number represents the whole and fractional number of days from a fixed, preset date (January 0, 0000) in the proleptic ISO calendar.
      moo = dir(fname);
      rightnow = datenum(datetime('now'));
      if (rightnow-moo.datenum)*24*60*60 > 10
        i10sec = 1;
      end
    end
    %fprintf(1,'%5i    %s \n',existfname(ii),fname)
    if exist(fname) > 0 & iaFound(ii) == 0 & i10sec == 1
      fnamelastloaded = fname;
      iExist = +1;
      loader = ['load ' fname];
      eval(loader);
      iaFound(ii) = +1;
      
      results(ii,1:6)    = oem.finalrates(1:6);
      resultsunc(ii,1:6) = oem.finalsigs(1:6);
      [mmn,nn] = size(oem.ak_water);

      if length(jacobian.wvjaclays_used) == iNumLay
        for iii = 1 : length(jacobian.wvjaclays_used)
          %iavg = jacobian.wvjaclays_used{iNumLay}-6;
          iavg = jacobian.wvjaclays_used{iii}-jacobian.wvjaclays_offset;
          pavg(iii) = mean(plays(iavg));
        end
      end

      save_cov_set.cov_set(:,ii)   = oem.cov_set;
      save_cov_set.fmat(:,ii)      = sqrt(diag(oem.fmat));
      save_cov_set.reg_type        = oem.reg_type;
      save_cov_set.xb_traceI(:,ii) = jacobian.scalar_i([1 length(jacobian.scalar_i)]); %% Co2/N2O/CH4/CFC11/CFC11/stemp etc  indices (1,6)
      save_cov_set.xb_wvzI(:,ii)   = jacobian.water_i([1 length(jacobian.water_i)]);   %% top/bottom indices
      save_cov_set.xb_tzzI(:,ii)   = jacobian.temp_i([1 length(jacobian.temp_i)]);     %% top/bottom indices
      save_cov_set.xb_ozzI(:,ii)   = jacobian.ozone_i([1 length(jacobian.ozone_i)]);   %% top/bottom indices
      save_cov_set.xb_trace(:,ii)  = oem.xb(1:6);                                             %% Co2/N2O/CH4/CFC11/CFC11/stemp etc  init values (1,6)
      save_cov_set.xb_wvz(:,ii)    = oem.xb(jacobian.water_i([1 length(jacobian.water_i)]));  %% top/bottom init values
      save_cov_set.xb_tzz(:,ii)    = oem.xb(jacobian.temp_i([1 length(jacobian.temp_i)]));    %% top/bottom init values
      save_cov_set.xb_ozz(:,ii)    = oem.xb(jacobian.ozone_i([1 length(jacobian.ozone_i)]));  %% top/bottom init values
      save_cov_set.xf_trace(:,ii)  = oem.finalrates(1:6);                                             %% Co2/N2O/CH4/CFC11/CFC11/stemp etc  final values (1,6)
      save_cov_set.xf_wvz(:,ii)    = oem.finalrates(jacobian.water_i([1 length(jacobian.water_i)]));  %% top/bottom final values
      save_cov_set.xf_tzz(:,ii)    = oem.finalrates(jacobian.temp_i([1 length(jacobian.temp_i)]));    %% top/bottom final values
      save_cov_set.xf_ozz(:,ii)    = oem.finalrates(jacobian.ozone_i([1 length(jacobian.ozone_i)]));  %% top/bottom final values

      nlays_straight_from_results(ii) = nn;
      nn0 = min(nn,iNumLay);

      rtime(ii)     = anomalyinfo.rtime;
      xb(ii,1:6)    = oem.xb(1:6);
      xbWV(ii,1:nn) = oem.xb((1:nn)+6+nn*0);
      xbT(ii,1:nn)  = oem.xb((1:nn)+6+nn*1);
      xbO3(ii,1:nn) = oem.xb((1:nn)+6+nn*2);

      if nn0 == iNumLay
        thedofs(ii) = oem.dofs;
        lencdofs(ii) = length(oem.cdofs);
        cdofs(ii,1:length(oem.cdofs)) = oem.cdofs;
        
        resultsWV(ii,1:nn) = oem.finalrates((1:nn)+6+nn*0);
        resultsT(ii,1:nn)  = oem.finalrates((1:nn)+6+nn*1);
        resultsO3(ii,1:nn) = oem.finalrates((1:nn)+6+nn*2);
        resultsWVunc(ii,1:nn) = oem.finalsigs((1:nn)+6+nn*0);
        resultsTunc(ii,1:nn)  = oem.finalsigs((1:nn)+6+nn*1);
        resultsO3unc(ii,1:nn) = oem.finalsigs((1:nn)+6+nn*2);

      else
        iWarning = iWarning + 1;
        iaWarning(iWarning) = ii;
  
        thedofs(ii) = oem.dofs;
        lencdofs(ii) = length(oem.cdofs);
        cdofs(ii,1:length(oem.cdofs)) = oem.cdofs;
  
        resultsWV(ii,:) = NaN;
        resultsT(ii,:)  = NaN;
        resultsO3(ii,:) = NaN;
        resultsWVunc(ii,:) = NaN;
        resultsTunc(ii,:)  = NaN;
        resultsO3unc(ii,:) = NaN;
  
        wah = oem.finalrates((1:nn)+6+nn*0);   resultsWV(ii,1:nn0) = wah(1:nn0);
        wah = oem.finalrates((1:nn)+6+nn*1);   resultsT(ii,1:nn0)  = wah(1:nn0);
        wah = oem.finalrates((1:nn)+6+nn*2);   resultsO3(ii,1:nn0) = wah(1:nn0);
  
        wah = oem.finalsigs((1:nn)+6+nn*0);   resultsWVunc(ii,1:nn0) = wah(1:nn0);
        wah = oem.finalsigs((1:nn)+6+nn*1);   resultsTunc(ii,1:nn0)  = wah(1:nn0);
        wah = oem.finalsigs((1:nn)+6+nn*2);   resultsO3unc(ii,1:nn0) = wah(1:nn0);
      end

      %%%%%%%%%%%%%%%%%%%%%%%%% DO AK %%%%%%%%%%%%%%%%%%%%%%%%% DO AK %%%%%%%%%%%%%%%%%%%%%%%%%
      if iAK > 0
        if ~exist('waterrate_akF_era5')
           waterrate_akF_era5 = nan(size(resultsWV));
           o3rate_akF_era5 = nan(size(resultsWV));
           temprrate_akF_era5 = nan(size(resultsWV));

           waterrate_ak0_era5 = nan(size(resultsWV));
           o3rate_ak0_era5 = nan(size(resultsWV));
           temprrate_ak0_era5 = nan(size(resultsWV));

           mean_ak_wv = nan(size(resultsWV));;
           mean_ak_T  = nan(size(resultsWV));;
           mean_ak_o3 = nan(size(resultsWV));;
        end
  
        %figure(2); plot(oem.ak_water',pjunk20,'c',max(oem.ak_water'),pjunk20,'rx-',mean(oem.ak_water'),pjunk20,'bx-'); set(gca,'ydir','reverse'); ylim([0.1 1000])
        %figure(3); plot(oem.ak_ozone',pjunk20,'c',max(oem.ak_ozone'),pjunk20,'rx-',mean(oem.ak_ozone'),pjunk20,'bx-'); set(gca,'ydir','reverse'); ylim([0.1 1000])
        %figure(4); plot(oem.ak_temp',pjunk20,'c',max(oem.ak_temp'),pjunk20,'rx-',mean(oem.ak_temp'),pjunk20,'bx-'); set(gca,'ydir','reverse'); ylim([0.1 1000])

        clear waterrate_ak0 waterrate_ak1 o3rate_ak0 o3rate_ak1 temprate_ak0 temprate_ak1
        ix = ii;
        waterrate_ak0 = ones(ix,1)*era5.trend_gas_1(1:100,ix)';
          for iii = 1 : length(jacobian.wvjaclays_used)
            junk = jacobian.wvjaclays_used{iii}-6;
            waterrate_ak1(:,iii) = mean(waterrate_ak0(:,junk)');
          end
          ak = oem.ak_water;
          mean_ak_wv(ix,1:length(ak)) = max(ak');
          waterrate_ak0_era5(ix,1:length(ak)) = (waterrate_ak1(ix,:)')';
          waterrate_akF_era5(ix,1:length(ak)) = (ak * waterrate_ak1(ix,:)')';
        o3rate_ak0 = ones(ix,1)*era5.trend_gas_3(1:100,ix)';
          for iii = 1 : length(jacobian.wvjaclays_used)
            junk = jacobian.wvjaclays_used{iii}-6;
            o3rate_ak1(:,iii) = mean(o3rate_ak0(:,junk)');
          end
          ak = oem.ak_ozone;
          mean_ak_o3(ix,1:length(ak)) = max(ak');
          o3rate_ak0_era5(ix,1:length(ak)) = (o3rate_ak1(ix,:)')';
          o3rate_akF_era5(ix,1:length(ak)) = (ak * o3rate_ak1(ix,:)')';
        temprate_ak0 = ones(ix,1)*era5.trend_ptemp(1:100,ix)';
          for iii = 1 : length(jacobian.wvjaclays_used)
            junk = jacobian.wvjaclays_used{iii}-6;
            temprate_ak1(:,iii) = mean(temprate_ak0(:,junk)');
          end
          ak = oem.ak_temp;
          mean_ak_T(ix,1:length(ak)) = max(ak');
          temprate_ak0_era5(ix,1:length(ak)) = (temprate_ak1(ix,:)')';
          temprate_akF_era5(ix,1:length(ak)) = (ak * temprate_ak1(ix,:)')';
      end   %% if iAK > 0
      %%%%%%%%%%%%%%%%%%%%%%%%% DO AK %%%%%%%%%%%%%%%%%%%%%%%%% DO AK %%%%%%%%%%%%%%%%%%%%%%%%%

      junknoise  = nan(2645,1);
      junknoise2 = nan(2645,1);
      junknoise(jacobian.chanset)  = sqrt(diag(oem.se));
      junknoise2(jacobian.chanset) = rateset.unc_rates(jacobian.chanset);
      junknoise2                   = rateset.unc_rates;

      %% depending on the lag1 correction, maybe junknoise and junknoise2 are different
      if sum([junknoise(jacobian.chanset)-junknoise2(jacobian.chanset)]) > eps
        fprintf(1,'WARNING fov %4i we have junknoise and junknoise2 differing \n',ii)
      end
  
      if isfield(oem,'spectral_deltan00')
        spectral_deltan00(:,ii) = oem.spectral_deltan00;
      end

      rates(:,ii)           = rateset.rates;
      fits(:,ii)            = oem.fit';
      if isfield(oem,'fitXcomponents')
        componentfits(1,:,ii) = oem.fitXcomponents(1,:);   %% trace gases CO2/N2O/CH4
        componentfits(2,:,ii) = oem.fitXcomponents(2,:);   %% ST
        componentfits(3,:,ii) = oem.fitXcomponents(3,:);   %% WV(z)
        componentfits(4,:,ii) = oem.fitXcomponents(4,:);   %% T(z)
        componentfits(5,:,ii) = oem.fitXcomponents(5,:);   %% O3(z)
      end

      nedt(:,ii)  = junknoise;  %% will have lots of NaNs since based on oem.se which only used selected channels
      nedt(:,ii)  = junknoise2; %% should be nicely filled in

    elseif exist(fname) > 0 & iaFound(ii) == 1
      %% do nothing, things are cool
    else
      iExist = -1;
      iaFound(ii) = 0;
      results(ii,1:6) = NaN;
      resultsWV(ii,:) = NaN;
      resultsT(ii,:)  = NaN;
      resultsO3(ii,:) = NaN;
  
      resultsunc(ii,1:6) = NaN;
      resultsWVunc(ii,:) = NaN;
      resultsTunc(ii,:)  = NaN;
      resultsO3unc(ii,:) = NaN;
  
      rates(:,ii) = NaN;
      fits(:,ii) = NaN;
      nedt(:,ii) = NaN;
    end

    if mod(ii,iNumAnomTimeSteps) == 0
      fprintf(1,'+ fat latbin   %2i of %2i \n',ii/iNumAnomTimeSteps,iNumAnomTiles);
    elseif iExist == 1
      if mod(ii,100) == 0
        fprintf(1,'o');
      elseif mod(ii,10) == 0
        fprintf(1,'.');
      end
    elseif iExist == -1
      fprintf(1,' ');
    end
  end
  
  fprintf(1,'\n');
  resultsunc = real(resultsunc);
  resultsWVunc = real(resultsWVunc);
  resultsTunc  = real(resultsTunc);
  resultsO3unc = real(resultsO3unc);
  
  display('Trying to look at atributes of the last file that should have been read in .... ')
  lser = ['!ls -lt ' fnamelastloaded]; eval(lser)
  
  fprintf(1,'found %4i of %4i \n',sum(iaFound),iNumAnomData)
  iDoAgain = -1;


  if sum(iaFound) < iNumAnomData
    if strfind(anomalydatafile,'_tile')
      plot_anomalies_1
    else
      simple = +1;
      plot_anomalies_global_trp
    end

    junk_globalavg_rawdata = data_anom.btavgAnomFinal(:,1:iNumAnomTimeSteps);
    junk_globalavg_spectral_deltan00 = spectral_deltan00(:,1:iNumAnomTimeSteps);
    figure(2); clf; 
    plot(yymm,smooth(junk_globalavg_spectral_deltan00(i0723,:),23),'r.-',...
        yymm,smooth(junk_globalavg_rawdata(i0723,:),23),'r--',...
       'linewidth',2); plotaxis2;
    title('Fitted BT0723 data in thick \newline raw data in dashes, CO2 only'); pause(0.1)

    figure(3); clf
    hold off; wah = save_cov_set.xf_trace(1:6,1:iNumAnomTimeSteps);
      plot(yymm,wah(1,:),'b',yymm,wah(2,:),'m',yymm,wah(3,:),'g',yymm,wah(4,:),'y',yymm,wah(5,:),'k',yymm,wah(6,:),'r','linewidth',2); hold on; 
    hold on; wah = save_cov_set.xb_trace(1:6,1:iNumAnomTimeSteps);
      plot(yymm,wah(1,:),'b--',yymm,wah(2,:),'m--',yymm,wah(3,:),'g--',yymm,wah(4,:),'y--',yymm,wah(5,:),'k--',yymm,wah(6,:),'r--','linewidth',2); hold on; 
      hl = legend('CO2','N2O','CH4','CFC11','CFC12','ST','location','best');
      set(gca,'fontsize',10); title('xb(trace gas + ST)');
    hold off

    iDoAgain = input('read in remaining files (-1/+1 Default) : '); 
    if length(iDoAgain) == 0
      iDoAgain = +1;
    end
  end
end

a = load(fnamelastloaded);
fprintf(1,'last loaded file %s has xb(1:6) = %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f \n',fnamelastloaded,a.oem.xb(1:6))

bah = nanmean(reshape(rtime,iNumAnomTimeSteps,iNumAnomTiles),2);
[yyy,mmm,ddd] = tai2utcSergio(bah);
yymmx = yyy + (mmm-1)/12 + (ddd-1)/12/30;

%%%%%%%%%%%%%%%%%%%%%%%%%

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
title('Fitted data in thick \newline raw data in dashes')

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

figure(5); clf
  subplot(211); 
  wah = save_cov_set.xf_trace(1:6,1:iNumAnomTimeSteps);
  plot(yymm,wah(1,:),'b',yymm,wah(2,:),'m',yymm,wah(3,:),'g',yymm,wah(6,:),'r','linewidth',2); hold on; 
  wah = save_cov_set.xb_trace(1:6,1:iNumAnomTimeSteps);
  plot(yymm,wah(1,:),'b--',yymm,wah(2,:),'m--',yymm,wah(3,:),'g--',yymm,wah(6,:),'r--','linewidth',2); hold on; 
  hl = legend('CO2','N2O','CH4','ST','location','best');
  set(gca,'fontsize',10); title('xb(trace gas + ST)');
figure(5);
  subplot(212); 
  wah = save_cov_set.xf_trace(1:6,1:iNumAnomTimeSteps);
  plot(yymm,wah(4,:),'y',yymm,wah(5,:),'k','linewidth',2); hold on; 
  wah = save_cov_set.xb_trace(1:6,1:iNumAnomTimeSteps);
  plot(yymm,wah(4,:),'y--',yymm,wah(5,:),'k--','linewidth',2); hold on; 
  hl = legend('CFC11','CFC12','location','best');
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

figure(8); clf; 
  rates_tropics = rates(:,(1:iNumAnomTimeSteps)+iNumAnomTimeSteps);
  fits_tropics  = fits(:,(1:iNumAnomTimeSteps)+iNumAnomTimeSteps);
  plot(f2645,nanmean(rates_tropics'-fits_tropics'),f2645,nanmean(rates_tropics')); plotaxis2; xlim([645 1620])
    hl = legend('rates-fits','rates','location','best','fontsize',10); title('Tropical set')

figure(9); clf;
plot(yymm,smooth(rates_global(i0667,:),23),yymm,smooth(rates_global(i0723,:),23),yymm,smooth(rates_global(i1598,:),23),yymm,smooth(rates_global(i1419,:),23),'linewidth',2); plotaxis2;
  legend('0667 S','0723 T','1598 T','1419 S','location','best','fontsize',10); title('Global average')

figure(10); clf;
plot(yymm,smooth(rates_tropics(i0667,:),23),yymm,smooth(rates_tropics(i0723,:),23),yymm,smooth(rates_tropics(i1598,:),23),yymm,smooth(rates_tropics(i1419,:),23),'linewidth',2); plotaxis2;
  legend('0667','0723','1598 T','1419 S','location','best','fontsize',10); title('Tropics average')

%[c,lags]=xcorr(rates_tropics(i0667,:),rates_tropics(i0723,:));  % compute cross correlation; keep lags vector
%[~,iLag]=max(c(find(lags==0):end));                             % find the max in one-sided
%s3=circshift(rates_tropics(i0723,:),[0 iLag]);                  % correct for the shift

disp('ret to continue to plot_anomalies_global_trp'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%
plot_anomalies_global_trp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath ../FIND_NWP_MODEL_TRENDS/
look_at_anomalies_computeERA5_monthly_trends_desc_or_asc_64fast(20,1);

figure(26); ax = axis; cx = caxis; figure(17); axis(ax); caxis(cx); colormap(cmap); %% WV 10-25 km, 2022 - 2025
figure(22); ax = axis; cx = caxis; figure(16); axis(ax); caxis(cx); colormap(cmap); %% Tz 10-25 km, 2022 - 2025

disp('hit RET to continue to axis = [2015 2025 0 45]'); pause
figure(26); axis([2015 2025 0 25])
figure(17); axis([2015 2025 0 25])
figure(22); axis([2015 2025 0 40])
figure(16); axis([2015 2025 0 40])

