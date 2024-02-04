addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE
%addpath ../../../FIND_TRENDS/
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME

disp('Note cluster job needs to be run from AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6')
disp('Note cluster job needs to be run from AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6')
disp('Note cluster job needs to be run from AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6')

if length(strfind(pwd,'SyntheticTimeSeries_ERA5_AIRSL3_CMIP6')) == 0
  error('you need to cd to MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/SyntheticTimeSeries_ERA5_AIRSL3_CMIP6')
end

system_slurm_stats

%% kleenslurm; sbatch  --exclude=cnode[204,225,260,267] --array=1-64 sergio_matlab_jobB.sbatch 10

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 64 for the 64 latbins
if length(JOB) == 0
  JOB = 3;
end

iNorD = +1; %% default, night
iNorD = -1; %% day
iNumYears = 20;

YMStart = [2015 01];  YMEnd = [2021 12];  %% OCO2
YMStart = [2014 09];  YMEnd = [2021 08];  %% OCO2
YMStart = [2002 09];  YMEnd = [2021 08];  %% 19 years
YMStart = [2002 09];  YMEnd = [2007 08];  %% 05 years
YMStart = [2002 09];  YMEnd = [2012 08];  %% 10 years
YMStart = [2002 09];  YMEnd = [2017 08];  %% 15 years
YMStart = [2002 09];  YMEnd = [2022 08];  %% 20 years

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('llsmap5.mat');

iOops = +1;
if iOops < 0
  % load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat
  %%%% kaboom = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat','h','p');  OOPS DELETERD DAMN
  h = kaboom.h;
  p = kaboom.p;
  clear kaboom;
    
  RH000 = layeramt2RH(h,p);
  
  pAIRSCLIMCAPSL3 = p;
  
  %% function airsChoice  = getdata_AIRSL3vsCLIMCAPSL3(iA,iNorD,iAorOrL,iNumYears);
  %%  iA = +1 (AIRS L3)
  %%       -1 CLIMCAPS
  %%       +3 CESM3
  climcapsL3trend = getdata_AIRSL3vsCLIMCAPSL3(-1,iNorD,0,iNumYears);
  
  pAIRSCLIMCAPSL3.stemp          = pAIRSCLIMCAPSL3.stemp          + reshape(climcapsL3trend.thestats64x72.stemprate,1,4608);
  pAIRSCLIMCAPSL3.ptemp(1:100,:) = pAIRSCLIMCAPSL3.ptemp(1:100,:) + reshape(permute(climcapsL3trend.thestats64x72.ptemprate,[3 1 2]),100,4608);
  pAIRSCLIMCAPSL3.gas_1(35:100,:) = pAIRSCLIMCAPSL3.gas_1(35:100,:).*(1 + reshape(permute(climcapsL3trend.thestats64x72.waterrate,[3 1 2]),66,4608));
  pAIRSCLIMCAPSL3.gas_3(35:100,:) = pAIRSCLIMCAPSL3.gas_3(35:100,:).*(1 + reshape(permute(climcapsL3trend.thestats64x72.ozonerate,[3 1 2]),66,4608));
  %{
  pAIRSCLIMCAPSL3.stemp          = pAIRSCLIMCAPSL3.stemp          + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
  pAIRSCLIMCAPSL3.ptemp(1:100,:) = pAIRSCLIMCAPSL3.ptemp(1:100,:) + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp;
  pAIRSCLIMCAPSL3.gas_1(1:100,:) = pAIRSCLIMCAPSL3.gas_1(1:100,:).*(1 + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1);
  pAIRSCLIMCAPSL3.gas_3(1:100,:) = pAIRSCLIMCAPSL3.gas_3(1:100,:).*(1 + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_3);
  %}
  RHAIRSCLIMCAPSL3 = layeramt2RH(h,pAIRSCLIMCAPSL3);
  
  RHAIRSCLIMCAPSL3rate = RHAIRSCLIMCAPSL3 - RH000;
  zonalRHAIRSCLIMCAPSL3rate = reshape(RHAIRSCLIMCAPSL3rate,100,72,64);
  zonalRHAIRSCLIMCAPSL3rate = squeeze(nanmean(zonalRHAIRSCLIMCAPSL3rate,2));
  
  kaboom = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat','rlat');
  rlat = kaboom.rlat; 
  
  zonalrlat = rlat;
  zonalplays = p.plays(1:100,3000);
  figure(1); clf; pcolor(zonalrlat,zonalplays,zonalRHAIRSCLIMCAPSL3rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
    set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('reconstruct Trate,WVrate \newline -> RH rate')
  
  TAIRSCLIMCAPSL3rate = permute(climcapsL3trend.thestats64x72.ptemprate,[3 1 2]);
  zonalTAIRSCLIMCAPSL3rate = TAIRSCLIMCAPSL3rate;
  zonalTAIRSCLIMCAPSL3rate = squeeze(nanmean(zonalTAIRSCLIMCAPSL3rate,2));
  figure(2); clf; pcolor(zonalrlat,zonalplays,zonalTAIRSCLIMCAPSL3rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.15 +0.15]); 
    set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('reconstruct Trate')
  
  WVAIRSCLIMCAPSL3rate = permute(climcapsL3trend.thestats64x72.waterrate,[3 1 2]);
  zonalWVAIRSCLIMCAPSL3rate = WVAIRSCLIMCAPSL3rate;
  zonalWVAIRSCLIMCAPSL3rate = squeeze(nanmean(zonalWVAIRSCLIMCAPSL3rate,2));
  figure(3); clf; pcolor(zonalrlat,zonalplays(35:100),zonalWVAIRSCLIMCAPSL3rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.15 +0.15]/10); 
    set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('reconstruct WV rate')

else
  kaboom = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_startwithERA5trends.mat','h','p','rlat'); 
  h = kaboom.h;
  p = kaboom.p;

  zonalrlat = kaboom.rlat;   
  zonalplays = p.plays(1:100,3000);

  zonalRHAIRSCLIMCAPSL3rate = [];
  zonalTAIRSCLIMCAPSL3rate = [];

end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see  FIND_NWP_MODEL_TRENDS/driver_computeERA5_monthly_trends.m  and do_the_AIRSL3_trends.m
disp('if you get silly messages like "YM timeperiod  = 2002/ 9 --> 2022/ 8 needs 240 of 228 timesteps" then check this >>>>>>>>')
disp('if you get silly messages like "YM timeperiod  = 2002/ 9 --> 2022/ 8 needs 240 of 228 timesteps" then check this >>>>>>>>')

%% NOTE THIS IS JUST HOW MUCH DATA YOU HAVE, AND IS DIFFERENT THAN YMStart,YMEnd where you set ACTUAL startYY/stopYY for trends
iYS = 2002; iYE = 2021;
iYS = 2002; iYE = 2022;

%% see  FIND_NWP_MODEL_TRENDS/driver_computeAIRSCLIMCAPSL3_monthly_trends.m  and do_the_AIRSCLIMCAPSL3_trends.m
%% see  FIND_NWP_MODEL_TRENDS/driver_compute_AIRS_CLIMCAPS_trends_desc_or_asc.m or FIND_NWP_MODEL_TRENDS/driver_compute_AIRS_CLIMCAPS_trends_desc_or_ascNOQuestioN.m  

if iNorD > 0
  %airsclimcapsl3_64x72 = load('/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_Sept2002_Aug2021_19yr_desc.mat');
  airsclimcapsl3_64x72 = load('/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_Sept2002_Aug2022_20yr_desc.mat');
else
  airsclimcapsl3_64x72 = load('/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_Sept2002_Aug2022_20yr_asc.mat');
end

[numtimesteps0] = length(airsclimcapsl3_64x72.days);
numtimesteps = numtimesteps0;
fprintf(1,'driver_check_WV_T_RH_AIRSCLIMCAPSL3_geo_and_spectral_rates2.m : numtimesteps = %3i \n',numtimesteps)
iNumYears = numtimesteps/12;

rlat = load('latB64.mat'); rlat = 0.5*(rlat.latB2(1:end-1)+rlat.latB2(2:end));
rlon = (1:72); rlon = -177.5 + (rlon-1)*5;

iNumYears = iYE-iYS;

yy = []; mm = []; dd = [];
for ii = iYS : iYE
  clear yyx mmx ddx
  if ii == iYS
    inum = 4;
    yyx(1:inum) = ii;
    mmx = 9:12;
    ddx = ones(size(mmx)) * 15;
%  elseif ii == 2021
%    %inum = 7;
%    %yyx(1:inum) = ii;
%    %mmx = 1 : 7;
%    inum = 8;
%    yyx(1:inum) = ii;
%    mmx = 1 : 8;
%    ddx = ones(size(mmx)) * 15;
  elseif ii == iYE
    inum = 8;
    yyx(1:inum) = ii;
    mmx = 1 : 8;
    ddx = ones(size(mmx)) * 15;
  else
    inum = 12;
    yyx(1:inum) = ii;
    mmx = 1:12;
    ddx = ones(size(mmx)) * 15;
  end
  fprintf(1,'%4i %2i \n',[ii inum])
  yy = [yy yyx];
  mm = [mm mmx];
  dd = [dd ddx];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

daysSince2002 = change2days(yy,mm,dd,2002);
whos daysSince2002

daysSince2002Start = change2days(YMStart(1),YMStart(2),15,2002);
daysSince2002End   = change2days(YMEnd(1),  YMEnd(2),  15,2002);

iNumYears = YMEnd(1) - YMStart(1);

usethese = find(daysSince2002  >= daysSince2002Start & daysSince2002 <= daysSince2002End);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,' YM timeperiod  = %4i/%2i --> %4i/%2i needs %3i of %3i timesteps \n',YMStart,YMEnd,length(usethese),numtimesteps)
yy = yy(usethese);
mm = mm(usethese);
dd = dd(usethese);
 
numtimesteps = length(yy);
airsclimcapsl3_64x72.days                = airsclimcapsl3_64x72.days(usethese);
airsclimcapsl3_64x72.save64x72_olr       = airsclimcapsl3_64x72.save64x72_olr(:,:,usethese);
airsclimcapsl3_64x72.save64x72_stemp     = airsclimcapsl3_64x72.save64x72_stemp(:,:,usethese);
airsclimcapsl3_64x72.save64x72_Q         = airsclimcapsl3_64x72.save64x72_Q(:,:,:,usethese);
airsclimcapsl3_64x72.save64x72_RH        = airsclimcapsl3_64x72.save64x72_RH(:,:,:,usethese);
airsclimcapsl3_64x72.save64x72_T         = airsclimcapsl3_64x72.save64x72_T(:,:,:,usethese);

airsclimcapsl3_64x72.save64x72_cld_frac  = airsclimcapsl3_64x72.save64x72_cld_frac(:,:,:,usethese);
airsclimcapsl3_64x72.save64x72_cld_pres  = airsclimcapsl3_64x72.save64x72_cld_pres(:,:,:,usethese);
airsclimcapsl3_64x72.save64x72_RHSurf    = airsclimcapsl3_64x72.save64x72_RHSurf(:,:,usethese);
airsclimcapsl3_64x72.save64x72_TWetSurf  = airsclimcapsl3_64x72.save64x72_TWetSurf(:,:,usethese);
airsclimcapsl3_64x72.save64x72_clrolr    = airsclimcapsl3_64x72.save64x72_clrolr(:,:,usethese);
airsclimcapsl3_64x72.save64x72_iceT      = airsclimcapsl3_64x72.save64x72_iceT(:,:,usethese);
airsclimcapsl3_64x72.save64x72_ice_od    = airsclimcapsl3_64x72.save64x72_ice_od(:,:,usethese);
airsclimcapsl3_64x72.save64x72_icesze    = airsclimcapsl3_64x72.save64x72_icesze(:,:,usethese);
airsclimcapsl3_64x72.save64x72_liq_water = airsclimcapsl3_64x72.save64x72_liq_water(:,:,usethese);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rtime = utc2taiSergio(yy,mm,dd,ones(size(yy))*12.0);
co2ppm = 370 + 2.2*((yy+mm/12)-2002);

%% see ~/MATLABCODE/CRODGERS_FAST_CLOUD/driver_stage2_ESRL_set_CO2_CH4_N2O.m
time_so_far = (yy-2000) + ((mm-1)+1)/12;
co2ppm = 368 + 2.1*time_so_far;
n2oppm = 315  + (332-315)/(2020-2000)*time_so_far; n2oppm = n2oppm/1000;
ch4ppm = 1.75 + (1.875-1.750)/(2020-2000)*time_so_far;

iConstORVary = -1;
if iConstORVary < 0
  %% see ~/MATLABCODE/CRODGERS_FAST_CLOUD/driver_stage2_ESRL_set_CO2_CH4_N2O.m
  co2ppm = 368 + 2.1*time_so_far;
  n2oppm = 315  + (332-315)/(2020-2000)*time_so_far; n2oppm = n2oppm/1000;
  ch4ppm = 1.75 + (1.875-1.750)/(2020-2000)*time_so_far;
else
  co2ppm = 385 * ones(size(co2ppm));
  n2oppm = 323 * ones(size(co2ppm))/1000;
  ch4ppm = 1813 * ones(size(co2ppm))/1000;
end

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';;

if iNorD > 0
  dirout = '../../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries';
  dirout = 'SimulateTimeSeries/CLIMCAPSL3/';
  dirout = 'SimulateTimeSeries/NIGHTorAVG/CLIMCAPSL3/';
  dirout = 'STS/NIGHTorAVG/CLIMCAPSL3/';
else
  dirout = 'SimulateTimeSeries/DAY/CLIMCAPSL3/';
  dirout = 'STS/DAY/CLIMCAPSL3/';
end

co2ppm_t = [];
n2oppm_t = [];
ch4ppm_t = [];

nemis_t = [];
emis_t  = [];
efreq_t = [];
rho_t   = [];

for ii = JOB

  h72 = struct;
  p72 = struct;
  
  h72 = h;
  h72.ptype = 0;
  h72.pfields = 1;

%% https://docserver.gesdisc.eosdis.nasa.gov/public/project/Sounder/CLIMCAPS_V2_L3_README.pdf
%% spec_hum   orbit_pass,    mass fraction of water   kg / kg
%%            air_pres_h2o,  vapor in moist air
%%            lat, lon
%% rel_hum    orbit_pass,    relative humidity over   unitless
%%            air_pres_h2o,  equilibrium phase
%%            lat, lon

  h72.ngas = 2;
  h72.gunit = [20 12]';  %%  g/kg and VMR
  h72.gunit = [21 12]';  %% kg/kg and VMR
  h72.glist = [ 1 3 ]';

  h72.ngas = 4;
  h72.gunit = [20 12 12 12]';  %%  g/kg and VMR  
  h72.gunit = [21 12 12 12]';  %% kg/kg and VMR
  h72.glist = [ 1  3  5  6]';

  h72.ngas = 1;
  h72.gunit = [20]';  %%  g/kg and VMR
  h72.gunit = [21]';  %% kg/kg and VMR
  h72.glist = [ 1]';
  
  p72.rtime = [];
  p72.co2ppm = [];
  for iii = 1 : numtimesteps
    co2ppm_t   = [co2ppm_t ones(1,72)*co2ppm(iii)];
    n2oppm_t   = [n2oppm_t ones(1,72)*n2oppm(iii)];
    ch4ppm_t   = [ch4ppm_t ones(1,72)*ch4ppm(iii)];

    p72.rtime  = [p72.rtime ones(1,72)*rtime(iii)];
    p72.co2ppm = [p72.co2ppm ones(1,72)*co2ppm(iii)];
  end
  
  iNlev = 24;
  iNlev  = length(airsclimcapsl3_64x72.Tlevs);
  iNlevQ = length(airsclimcapsl3_64x72.Qlevs);
  diffTQ = iNlev - iNlevQ;
  plevsnwp  = airsclimcapsl3_64x72.Tlevs/100;
  plevsnwpQ = airsclimcapsl3_64x72.Qlevs/100;
  p72.nlevs = ones(size(p72.rtime)) * iNlev;
  p72.plevs = airsclimcapsl3_64x72.Tlevs/100 * ones(1,72*numtimesteps);
  
  p72.rlat = rlat(ii) * ones(1,72*numtimesteps);  p72.rlat = p72.rlat(:)';
  p72.rlon = rlon' * ones(1,numtimesteps);        p72.rlon = p72.rlon(:)';
  p72.plat = p72.rlat;
  p72.plon = p72.rlon;

  junk = airsclimcapsl3_64x72.save64x72_stemp;     junk = permute(junk,[3 2 1]);    junk = squeeze(junk(:,:,ii));     junk = junk'; p72.stemp = reshape(junk,1,72*numtimesteps);
  junk = airsclimcapsl3_64x72.save64x72_T;         junk = permute(junk,[4 3 2 1]);  junk = squeeze(junk(:,:,:,ii));   junk = permute(junk,[2 3 1]); p72.ptemp = reshape(junk,iNlev,72*numtimesteps);
  junk = airsclimcapsl3_64x72.save64x72_Q;         junk = permute(junk,[4 3 2 1]);  junk = squeeze(junk(:,:,:,ii));   junk = permute(junk,[2 3 1]); p72.gas_1 = reshape(junk,iNlevQ,72*numtimesteps);
  %junk = airsclimcapsl3_64x72.save64x72_O3;       junk = permute(junk,[4 3 2 1]);  junk = squeeze(junk(:,:,:,ii));   junk = permute(junk,[2 3 1]); p72.gas_3 = reshape(junk,iNlev,72*numtimesteps);
  %junk = airsclimcapsl3_64x72.save64x72_CO;       junk = permute(junk,[4 3 2 1]);  junk = squeeze(junk(:,:,:,ii));   junk = permute(junk,[2 3 1]); p72.gas_5 = reshape(junk,iNlev,72*numtimesteps);
  %junk = airsclimcapsl3_64x72.save64x72_CH4;      junk = permute(junk,[4 3 2 1]);  junk = squeeze(junk(:,:,:,ii));   junk = permute(junk,[2 3 1]); p72.gas_6 = reshape(junk,iNlev,72*numtimesteps);
  junk = airsclimcapsl3_64x72.save64x72_RH*100;    junk = permute(junk,[4 3 2 1]);  junk = squeeze(junk(:,:,:,ii));   junk = permute(junk,[2 3 1]); p72.rh    = reshape(junk,iNlevQ,72*numtimesteps);

  junk = p72.gas_1; p72.gas_1 = zeros(size(p72.ptemp)); p72.gas_1(diffTQ+1:iNlev,:) = junk; for jj = 1 : diffTQ; frac = (jj-1)/diffTQ; p72.gas_1(jj,:) = junk(1,:) * frac; end;
  junk = p72.rh;    p72.rh = zeros(size(p72.ptemp));    p72.rh(diffTQ+1:iNlev,:) = junk;    for jj = 1 : diffTQ; frac = (jj-1)/diffTQ; p72.rh(jj,:) = junk(1,:) * frac; end;
  
%  p72.gas_1(isnan(p72.gas_1)) = 0.0;

  junk = airsclimcapsl3_64x72.save64x72_T;      miaowT = squeeze(junk(JOB,36,:,1));  figure(1); semilogy(miaowT,plevsnwp,p72.ptemp(:,1),p72.plevs(:,1),'r.-'); set(gca,'ydir','reverse'); ylim([0.01 1000])
  junk = airsclimcapsl3_64x72.save64x72_Q;      miaowQ = squeeze(junk(JOB,36,:,1));  figure(2); loglog(miaowQ,plevsnwpQ,p72.gas_1(:,1),p72.plevs(:,1),'r.-');  set(gca,'ydir','reverse'); ylim([0.01 1000])
  junk = airsclimcapsl3_64x72.save64x72_Q;      miaowQ = squeeze(junk(JOB,36,:,1));  figure(2); loglog(miaowQ,plevsnwpQ,p72.gas_1(:,1),p72.plevs(:,1),'r.-');  set(gca,'ydir','reverse'); ylim([0.01 1000])
  junk = airsclimcapsl3_64x72.save64x72_RH*100; miaowRH = squeeze(junk(JOB,36,:,1)); figure(3); semilogy(miaowRH,plevsnwpQ,p72.rh(:,1),p72.plevs(:,1),'r.-');   set(gca,'ydir','reverse'); ylim([0.01 1000])

%  p72.scanang = zeros(size(p72.stemp));
%  p72.satzen = zeros(size(p72.stemp));
  p72.zobs = 705000 * ones(size(p72.stemp));
  p72.scanang = ones(size(p72.stemp)) * 22;
  p72.satzen = vaconv(p72.scanang, p72.zobs, zeros(size(p72.zobs)));
  p72.solzen = 150 * ones(size(p72.stemp));
  %p72.spres = 1000 * ones(size(p72.stemp));
  %p72.salti = 0 * ones(size(p72.stemp));
  p72.spres = reshape(p.spres,72,64); p72.spres = p72.spres(:,ii) * ones(1,1*numtimesteps); p72.spres = p72.spres(:)';
  p72.salti = reshape(p.salti,72,64); p72.salti = p72.salti(:,ii) * ones(1,1*numtimesteps); p72.salti = p72.salti(:)';

  iCheck = 50; %% CLIMCAPS has 100 levels for T, 66 for WV
  p72.verybad = zeros(size(p72.salti));
  p72.lonbin = zeros(size(p72.salti));
  for jjj = 1 : 72
    lonbinx = (jjj-1) + (1:72:length(p72.stemp));
    p72.lonbin(lonbinx) = jjj;
  end
  verybad = find( isnan(p72.ptemp(iCheck,:)) | isnan(p72.gas_1(iCheck,:)) | isnan(p72.rh(iCheck,:)) );
  p72.verybad(verybad) = 1;
  if length(verybad) > 0  
    for jjj = 1 : length(verybad)
      bah = verybad(jjj);
      badah = (1:72:12*iNumYears*72) + (mod(bah,72)-1);
      [Y,I1,I2] = intersect(bah,badah);
      if I2 > 1
        p72.stemp(bah)   = p72.stemp(badah(I2-1));
        p72.ptemp(:,bah) = p72.ptemp(:,badah(I2-1));
        p72.gas_1(:,bah) = p72.gas_1(:,badah(I2-1));
%        p72.gas_3(:,bah) = p72.gas_3(:,badah(I2-1));
        p72.rh(:,bah)    = p72.rh(:,badah(I2-1));
      else
        p72.stemp(bah)   = p72.stemp(badah(I2+1));
        p72.ptemp(:,bah) = p72.ptemp(:,badah(I2+1));
        p72.gas_1(:,bah) = p72.gas_1(:,badah(I2+1));
%        p72.gas_3(:,bah) = p72.gas_3(:,badah(I2+1));
        p72.rh(:,bah)    = p72.rh(:,badah(I2+1));
      end
    end    
  end

  for jjj = 1 : length(p72.stemp)
    junklevs = p72.plevs(:,jjj);
    junkT    = p72.ptemp(:,jjj);
    junkQ    = p72.gas_1(:,jjj);
    junkrh   = p72.rh(:,jjj);
    moo = find(junklevs >= p72.spres(jjj),1);
    if length(moo) == 0
      moo = length(junklevs);
    end
    p72.nlevs(jjj) = moo;

    bad = find(junkT(1:moo) < 10 | isnan(junkT(1:moo)));
    good = find(junkT(1:moo) > 150);
    if length(bad) > 0
      junkT(bad) = interp1(log(junklevs(good)),junkT(good),log(junklevs(bad)),[],'extrap');
      p72.ptemp(bad,jjj) = junkT(bad);
    end

    bad = find(junkQ(1:moo) < 0 | isnan(junkQ(1:moo)));
    good = find(junkQ(1:moo) > 0);
    if length(bad) > 0
      junkQ(bad) = interp1(log(junklevs(good)),junkQ(good),log(junklevs(bad)),[],'extrap');
      haha = find(junkQ(bad) < 0);
      junkQ(bad(haha)) = junkQ(good(end));
      p72.gas_1(bad,jjj) = junkQ(bad);
    end

    bad = find(junkrh(1:moo) < 0 | isnan(junkrh(1:moo)));
    good = find(junkrh(1:moo) > 0);
    if length(bad) > 0
      junkrh(bad) = interp1(log(junklevs(good)),junkrh(good),log(junklevs(bad)),[],'extrap');
      haha = find(junkrh(bad) < 0);
      junkrh(bad(haha)) = junkrh(good(end));
      p72.rh(bad,jjj) = junkrh(bad);
    end
  end
    
  pcolor(rlon,1:numtimesteps,reshape(p72.stemp,72,numtimesteps)'); colormap jet; colorbar
  pcolor(rlon,1:numtimesteps,reshape(p72.spres,72,numtimesteps)'); colormap jet; colorbar
  pcolor(rlon,1:numtimesteps,reshape(p72.salti,72,numtimesteps)'); colormap jet; colorbar
  
  p72.cngwat = zeros(size(p72.stemp));
  p72.cngwat2 = zeros(size(p72.stemp));
  p72.cfrac = zeros(size(p72.stemp));
  p72.cfrac2 = zeros(size(p72.stemp));
  p72.cfrac12 = zeros(size(p72.stemp));
  
  p72.nemis = ones(size(p72.stemp)) * 2;
  p72.efreq(1,:) = ones(size(p72.stemp)) * 200;
  p72.efreq(2,:) = ones(size(p72.stemp)) * 3200;
  p72.emis(1,:) = ones(size(p72.stemp)) * 0.98;
  p72.emis(2,:) = ones(size(p72.stemp)) * 0.98;
  p72.rho = (1-p72.emis)/pi;
  
  fstr = ['_' num2str(YMStart(1),'%04d') '_' num2str(YMStart(2),'%02d') '_' num2str(YMEnd(1),'%04d') '_' num2str(YMEnd(2),'%02d')];
  %fip = [dirout '/simulate64binsAIRSCLIMCAPSL3_' num2str(ii) fstr '.ip.rtp'];
  %fop = [dirout '/simulate64binsAIRSCLIMCAPSL3_' num2str(ii) fstr '.op.rtp'];
  %frp = [dirout '/simulate64binsAIRSCLIMCAPSL3_' num2str(ii) fstr '.rp.rtp'];
  fip = [dirout '/simulate64binsCLIM_' num2str(ii) fstr '.ip.rtp'];
  fop = [dirout '/simulate64binsCLIM_' num2str(ii) fstr '.op.rtp'];
  frp = [dirout '/simulate64binsCLIM_' num2str(ii) fstr '.rp.rtp'];

  klayers_sarta_check_WV_T_RH_geo_and_spectral_rates2

  iType = 4;
  plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2

  pwd
  disp(' ')
  disp('output in eg           SimulateTimeSeries/ NIGHTorAVG   or    DAY /CLIMCAPSL3/reconstruct_climcapsL3_spectra_geo_rlat[01-64]_2002_09_2022_08.mat')
  disp('         also makes eg SimulateTimeSeries/ NIGHTorAVG   or    DAY /CLIMCAPSL3/simulate64binsAIRSCLIMCAPSL3_[01-64]_2002_09_2022_08.[i/o/r]p.rtp');
  disp('then run driver_spectral_trends_latbin_1_64_sarta.m using   sbatch  --array=1-64 sergio_matlab_jobB.sbatch 12; make sure you set iModel correctly')
  disp(' ')

end
