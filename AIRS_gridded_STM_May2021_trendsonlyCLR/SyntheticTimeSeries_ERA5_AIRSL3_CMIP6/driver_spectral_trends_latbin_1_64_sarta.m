addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil/
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/

system_slurm_stats
JOB = str2num(getenv('SLURM_ARRAY_TASK_ID')); %JOB = 1 -- 64
if length(JOB) == 0
  JOB = 32;
end

iStartYY = 2018; iNumYears = 04;

% ls -lt /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/ERA5/simulate64binsERA5_15*.rp.rtp
% -rw-rw-r-- 1 sergio pi_strow 316109325 Oct 29 15:35 /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/ERA5/simulate64binsERA5_15_2002_09_2022_08.rp.rtp
% -rw-rw-r-- 1 sergio pi_strow 110653581 Sep 19 07:01 /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/ERA5/simulate64binsERA5_15_2014_09_2021_08.rp.rtp
% -rw-rw-r-- 1 sergio pi_strow 300305037 Jul 19 15:42 /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/ERA5/simulate64binsERA5_15.rp.rtp
if ~exist('iStartYY')
  iiStartYY = 2002;
end;

if ~exist('iNumYears')
  iNumYears = 19;
  iNumYears = 20;
end

if ~exist('iModel')
  iModel = 2;   %% MERRA2
  iModel = 1;   %% UMBC
  iModel = -3;  %% CLIMCAPS
  iModel = +3;  %% AIRS v7
  iModel = 5;   %% ERA5
  iModel = 51;  %% ERA5 Const TraceGas
end

iNorD = -1; %% day
iNorD = +1; %% DEFAULT, night

iConstORVary = -1;   %% vary the CO2 N2O CH4 as function of time  %% DEFAULT
iConstORVary = +1;   %% constant CO2 N2O CH4 as function of time  %% to replace driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2_constracegas.m

if iNorD == +1
  dndir = 'NIGHTorAVG';
else
  dndir = 'DAY';
end

if length(JOB) == 0
  JOB = 32;
end

workdir = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/';
JOB_YYS_YYE_str = [num2str(JOB) '_' num2str(iStartYY) '_09_' num2str(iStartYY+iNumYears) '_08.rp.rtp'];
STS = 'SimulateTimeSeries/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iModel == 1
  dirout = ['UMBC_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(JOB,'%02d')  '/'];
elseif iModel == 5
  dirout = ['ERA5_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(JOB,'%02d')  '/'];
elseif iModel == 51
  dirout = ['ERA5_CONST_G_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(JOB,'%02d')  '/'];
elseif iModel == 2
  dirout = ['MERRA2_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(JOB,'%02d')  '/'];
elseif iModel == +3
  dirout = ['AIRSL3_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(JOB,'%02d')  '/'];
elseif iModel == -3
  dirout = ['CLIMCAPSL3_SARTA_SPECTRAL_RATES/KCARTA_latbin' num2str(JOB,'%02d')  '/'];
else
  error('cannot assign dirout')
end

if ~exist(dirout)
  mker = ['!mkdir ' dirout];
  eval(mker);
end

if iNorD  == 1
  fnameOUT = [dirout '/sarta_spectral_trends_latbin' num2str(JOB,'%02d')     '_' num2str(iStartYY) '_09_' num2str(iStartYY+iNumYears) '_08.mat'];
else
  fnameOUT = [dirout '/sarta_spectral_trends_asc_latbin' num2str(JOB,'%02d') '_' num2str(iStartYY) '_09_' num2str(iStartYY+iNumYears) '_08.mat'];
end

if ~exist(fnameOUT)
  fprintf(1,'will save to %s \n',fnameOUT);
else
  fprintf(1,'%s already exists \n',fnameOUT)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iModel == 1
  disp('UMBC')
  % if iNumYears == 07
  %   frp = ['STS/' dndir '/ERA5/simulate64binsERA5_' num2str(JOB) '_2014_09_2021_08.rp.rtp'];  %% OCO2
  % elseif iNumYears == 19
  %   frp = ['STS/' dndir '/ERA5/simulate64binsERA5_' num2str(JOB) '.rp.rtp'];                  %% GENERIC, 19 year
  %   %%% OOPSY I THINIK I DELETED THESE SO RERUN IF NEEDED, and then 2002_09_2021_08 will now be added to name
  % elseif iNumYears == 20
  %   frp = ['STS/' dndir '/UMBC/simulate64binsUMBC_' num2str(JOB) '_2002_09_2022_08.rp.rtp'];  %% 20 year
  % end
  frp = ['STS/' dndir '/UMBC/simulate64binsUMBC_' JOB_YYS_YYE_str];  %% N year

elseif iModel == 5
  disp('ERA5')
  % if iNumYears == 07
  %   frp = ['STS/' dndir '/ERA5/simulate64binsERA5_' num2str(JOB) '_2014_09_2021_08.rp.rtp'];  %% OCO2
  % elseif iNumYears == 19
  %   frp = ['STS/' dndir '/ERA5/simulate64binsERA5_' num2str(JOB) '.rp.rtp'];                  %% GENERIC, 19 year
  %   %%% OOPSY I THINIK I DELETED THESE SO RERUN IF NEEDED, and then 2002_09_2021_08 will now be added to name
  % elseif iNumYears == 20
  %   frp = ['STS/' dndir '/ERA5/simulate64binsERA5_' num2str(JOB) '_2002_09_2022_08.rp.rtp'];  %% 20 year
  % end
  frp = ['STS/' dndir '/ERA5/simulate64binsERA5_' JOB_YYS_YYE_str];  %% 20 year

elseif iModel == 51
  disp('ERA5_ConstTG Const Tragce Gas')
  % if iNumYears == 07
  %   frp = ['STS/' dndir '/ERA5_ConstG/simulate64binsERA5_' num2str(JOB) '_2014_09_2021_08.rp.rtp'];  %% OCO2
  % elseif iNumYears == 19
  %   frp = ['STS/' dndir '/ERA5_ConstG/simulate64binsERA5_' num2str(JOB) '.rp.rtp'];                  %% GENERIC, 19 year
  %   %%% OOPSY I THINIK I DELETED THESE SO RERUN IF NEEDED, and then 2002_09_2021_08 will now be added to name
  % elseif iNumYears == 20
  %   frp = ['STS/' dndir '/ERA5_ConstG/simulate64binsERA5_' num2str(JOB) '_2002_09_2022_08.rp.rtp'];  %% 20 year
  % end
  frp = ['STS/' dndir '/ERA5_ConstG/simulate64binsERA5_' JOB_YYS_YYE_str];  %% 20 year

elseif iModel == 2
  disp('MERRA2')
  % if iNumYears == 19
  %   frp = ['STS/' dndir '/MERRA2/simulate64binsMERRA2_' num2str(JOB) '.rp.rtp'];                  %% GENERIC, 19 year
  %   %%% OOPSY I THINIK I DELETED THESE SO RERUN IF NEEDED, and then 2002_09_2021_08 will now be added to name
  % elseif iNumYears == 20
  %   frp = ['STS/' dndir '/MERRA2/simulate64binsMERRA2_' num2str(JOB) '_2002_09_2022_08.rp.rtp'];  %% 20 year
  % end
  frp = ['STS/' dndir '/MERRA2/simulate64binsMERRA2_' JOB_YYS_YYE_str];

elseif iModel == 3
  disp('AIRS L3 ')
  % if iNumYears == 19
  %   frp = ['STS/' dndir '/AIRSL3/simulate64binsAIRSL3_' num2str(JOB) '.rp.rtp'];                  %% GENERIC, 19 year
  %   %%% OOPSY I THINIK I DELETED THESE SO RERUN IF NEEDED, and then 2002_09_2021_08 will now be added to name
  % elseif iNumYears == 20
  %   frp = ['STS/' dndir '/AIRSL3/simulate64binsAIRSL3_' num2str(JOB) '_2002_09_2022_08.rp.rtp'];  %% 20 year
  % end
  frp = ['STS/' dndir '/AIRSL3/simulate64binsAIRSL3_' JOB_YYS_YYE_str];  %% 20 year

elseif iModel == -3
  disp('AIRS CLIMCAPSL3 ')
  % if iNumYears == 19
  %   frp = ['STS/' dndir '/AIRSCLIMCAPSL3/simulate64binsAIRSCLIMCAPSL3_' num2str(JOB) '.rp.rtp'];    %% GENERIC, 19 year
  %   %%% OOPSY I THINIK I DELETED THESE SO RERUN IF NEEDED, and then 2002_09_2021_08 will now be added to name
  % elseif iNumYears == 20
  %   frp = ['STS/' dndir '/CLIMCAPSL3/simulate64binsCLIM_' num2str(JOB) '_2002_09_2022_08.rp.rtp'];  %% 20 year
  % end
  frp = ['STS/' dndir '/CLIMCAPSL3/simulate64binsCLIM_' JOB_YYS_YYE_str];  %% 20 year

else
  error('iModel = 2 (MERRA2) or 5 (ERA5) or +/-3 (AIRS v7/CLIMCAPS) ')
end

fprintf(1,'JOB = %2i iNumYears = %2i frp = %s \n',JOB,iNumYears,frp);

[h,ha,pppp,pa] = rtpread(frp);

numdatapts = iNumYears * 12 * 72; %% eg 19 years --> 16416 points

[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2ppmv(h,pppp,1:length(pppp.stemp),2);

bad_lonbins = [];
if isfield(pppp,'verybad')
  bad_lonbins = find(pppp.verybad > 0);
  if length(bad_lonbins) > 0
    disp('the are the verybad lonbins, so you may not want to believe the spectral or geophysical trends from this set of routines')
    bad_lonbins = unique(pppp.lonbin(bad_lonbins));
  end
end

for xx = 1 : 72
  fprintf(1,'JOB ( = latbin) = %2i       : Lonbin %2i of 72 : doing 2645 chans : ',JOB,xx)
  ind = (1:72:numdatapts);
  ind = ind + (xx-1);
  boo = find(ind <= iNumYears*12*72);
  ind = ind(boo); 
  plot(pppp.rlon(ind),'o-'); title(num2str(xx)); pause(0.1);

  [yy,mm,dd,hh] = tai2utcSergio(pppp.rtime(ind));
  dayOFtime = change2days(yy,mm,dd,2002);

  raaData = rad2bt(h.vchan,pppp.rcalc(:,ind));
  for iii = 1 : 2645
    if mod(iii,1000) == 0
      fprintf(1,'+');
    elseif mod(iii,100) == 0
      fprintf(1,'.');
    end
    data = double(raaData(iii,:));
    zoo = find(isfinite(data));
    if length(zoo) > 20
      [junk err] = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
      thesave.xtrend(iii,xx)    = junk(2);
      thesave.xtrendErr(iii,xx) = err.se(2);
      thesave.const(iii,xx)     = junk(1);
    else
      thesave.xtrend(iii,xx)    = NaN;
      thesave.xtrendErr(iii,xx) = NaN;
      thesave.const(iii,xx)     = NaN;
    end
  end
  fprintf(1,' 1231 cm-1 trend (AIRS L1C Ch 1520) = %8.6f \n',thesave.xtrend(1520,xx));
end

comment = 'see MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/driver_spectral_trends_latbin_1_64_sarta.m';
% if iNumYears == 07
%   saver = ['save ' dirout '/sarta_spectral_trends_latbin' num2str(JOB,'%02d') '_2014_09_2021_08.mat thesave comment'];
% elseif iNumYears == 19
%   saver = ['save ' dirout '/sarta_spectral_trends_latbin' num2str(JOB,'%02d') '.mat thesave comment'];
%   saver = ['save ' dirout '/sarta_spectral_trends_latbin' num2str(JOB,'%02d') '_2002_09_2021_08.mat thesave comment'];
% elseif iNumYears == 20
%   saver = ['save ' dirout '/sarta_spectral_trends_latbin' num2str(JOB,'%02d') '_2002_09_2022_08.mat thesave comment'];
% end


%% STS/NIGHTorAVG/ERA5_ConstG/reconstruct_era5_spectra_geo_rlat31_2018_09_2022_08.mat
if ~exist(fnameOUT)
  fprintf(1,'saving to %s \n',fnameOUT);
  saver = ['save '  fnameOUT ' thesave comment bad_lonbins'];
  eval(saver)
else
  fprintf(1,'%s already exists \n',fnameOUT)
end
