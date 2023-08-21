addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE
%addpath ../../../FIND_TRENDS/
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME

system_slurm_stats

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 64 for the 64 latbins
%JOB = 40

kaboom = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat','h','p','nwp_spectral_trends_cmip6_era5_airsL3_umbc');
h = kaboom.h;
p = kaboom.p;
nwp_spectral_trends_cmip6_era5_airsL3_umbc = kaboom.nwp_spectral_trends_cmip6_era5_airsL3_umbc;
clear kaboom;

load('llsmap5.mat');

RH000 = layeramt2RH(h,p);

pAIRSL3 = p;
%{
airsL3trend = getdata_AIRSL3vsCLIMCAPSL3(+1,1,0,20);  %% unfortunately at 12, 24 levels
pAIRSL3.stemp          = pAIRSL3.stemp          + reshape(airsL3trend.thestats64x72.stemprate,1,4608);
pAIRSL3.ptemp(1:100,:) = pAIRSL3.ptemp(1:100,:) + reshape(permute(airsL3trend.thestats64x72.ptemprate,[3 1 2]),100,4608);
pAIRSL3.gas_1(35:100,:) = pAIRSL3.gas_1(35:100,:).*(1 + reshape(permute(airsL3trend.thestats64x72.waterrate,[3 1 2]),66,4608));
pAIRSL3.gas_3(35:100,:) = pAIRSL3.gas_3(35:100,:).*(1 + reshape(permute(airsL3trend.thestats64x72.ozonerate,[3 1 2]),66,4608));
%}
pAIRSL3.stemp          = pAIRSL3.stemp          + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
pAIRSL3.ptemp(1:100,:) = pAIRSL3.ptemp(1:100,:) + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp;
pAIRSL3.gas_1(1:100,:) = pAIRSL3.gas_1(1:100,:).*(1 + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1);
pAIRSL3.gas_3(1:100,:) = pAIRSL3.gas_3(1:100,:).*(1 + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_3);
RHAIRSL3 = layeramt2RH(h,pAIRSL3);

RHAIRSL3rate = RHAIRSL3 - RH000;
zonalRHAIRSL3rate = reshape(RHAIRSL3rate,100,72,64);
zonalRHAIRSL3rate = squeeze(nanmean(zonalRHAIRSL3rate,2));

kaboom = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat','rlat');
rlat = kaboom.rlat; 

zonalrlat = rlat;
zonalplays = p.plays(1:100,3000);
figure(1); clf; pcolor(zonalrlat,zonalplays,zonalRHAIRSL3rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('reconstruct Trate,WVrate \newline -> RH rate')

TAIRSL3rate = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp;
zonalTAIRSL3rate = reshape(TAIRSL3rate,100,72,64);
zonalTAIRSL3rate = squeeze(nanmean(zonalTAIRSL3rate,2));
figure(2); clf; pcolor(zonalrlat,zonalplays,zonalTAIRSL3rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.15 +0.15]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('reconstruct Trate')

WVAIRSL3rate = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1;
zonalWVAIRSL3rate = reshape(WVAIRSL3rate,100,72,64);
zonalWVAIRSL3rate = squeeze(nanmean(zonalWVAIRSL3rate,2));
figure(3); clf; pcolor(zonalrlat,zonalplays,zonalWVAIRSL3rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.15 +0.15]/10); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('reconstruct WVrate')

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

%% see  FIND_NWP_MODEL_TRENDS/driver_computeAIRSL3_monthly_trends.m  and do_the_AIRSL3_trends.m
%airsl3_64x72 = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_Sept2002_Jul2021_19yr_desc.mat');
%airsl3_64x72 = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_Sept2002_Aug2021_19yr_desc.mat');
airsl3_64x72 = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_Sept2002_Aug2022_20yr_desc.mat');

[numtimesteps0] = length(airsl3_64x72.days);
numtimesteps = numtimesteps0;
fprintf(1,'driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m : numtimesteps = %3i \n',numtimesteps)

rlat = load('latB64.mat'); rlat = 0.5*(rlat.latB2(1:end-1)+rlat.latB2(2:end));
rlon = (1:72); rlon = -177.5 + (rlon-1)*5;

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

YMStart = [2015 01];  YMEnd = [2021 12];  %% OCO2
YMStart = [2014 09];  YMEnd = [2021 08];  %% OCO2
YMStart = [2002 09];  YMEnd = [2021 08];  %% 19 years

YMStart = [2002 09];  YMEnd = [2022 08];  %% 20 years
YMStart = [2002 09];  YMEnd = [2017 08];  %% 15 years
YMStart = [2002 09];  YMEnd = [2012 08];  %% 10 years
YMStart = [2002 09];  YMEnd = [2007 08];  %% 05 years

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
airsl3_64x72.days                = airsl3_64x72.days(usethese);
airsl3_64x72.save64x72_olr       = airsl3_64x72.save64x72_olr(:,:,usethese);
airsl3_64x72.save64x72_stemp     = airsl3_64x72.save64x72_stemp(:,:,usethese);
airsl3_64x72.save64x72_CH4       = airsl3_64x72.save64x72_CH4(:,:,:,usethese);
airsl3_64x72.save64x72_CO        = airsl3_64x72.save64x72_CO(:,:,:,usethese);
airsl3_64x72.save64x72_O3        = airsl3_64x72.save64x72_O3(:,:,:,usethese);
airsl3_64x72.save64x72_Q         = airsl3_64x72.save64x72_Q(:,:,:,usethese);
airsl3_64x72.save64x72_RH        = airsl3_64x72.save64x72_RH(:,:,:,usethese);
airsl3_64x72.save64x72_T         = airsl3_64x72.save64x72_T(:,:,:,usethese);
airsl3_64x72.save64x72_cld_frac  = airsl3_64x72.save64x72_cld_frac(:,:,:,usethese);
airsl3_64x72.save64x72_cld_pres  = airsl3_64x72.save64x72_cld_pres(:,:,:,usethese);
airsl3_64x72.save64x72_RHSurf    = airsl3_64x72.save64x72_RHSurf(:,:,usethese);
airsl3_64x72.save64x72_TWetSurf  = airsl3_64x72.save64x72_TWetSurf(:,:,usethese);
airsl3_64x72.save64x72_clrolr    = airsl3_64x72.save64x72_clrolr(:,:,usethese);
airsl3_64x72.save64x72_iceT      = airsl3_64x72.save64x72_iceT(:,:,usethese);
airsl3_64x72.save64x72_ice_od    = airsl3_64x72.save64x72_ice_od(:,:,usethese);
airsl3_64x72.save64x72_icesze    = airsl3_64x72.save64x72_icesze(:,:,usethese);
airsl3_64x72.save64x72_liq_water = airsl3_64x72.save64x72_liq_water(:,:,usethese);

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

dirout = '../../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries';
dirout = 'SimulateTimeSeries/AIRSL3/';

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

  h72.ngas = 2;
  h72.gunit = [20 12]';  %% g/kg and VMR
  h72.glist = [ 1 3 ]';

%% https://docserver.gesdisc.eosdis.nasa.gov/public/project/AIRS/V7_L3_Product_User_Guide.pdf
%% pg 16, 
%% H2O_MMR_Lyr 32-bit floating point H2OPressureLay(12) Water vapor mass mixing ratio
%%                                                      averaged over each of standard
%%                                                      pressure layers (g/kg dry air)
%% H2O_MMR 32-bit floating point  H2OPressureLev(12)    Water vapor mass mixing ratio at 
%%                                                      standard pressure levels (g/kg dryair)
%% 
%  h72.ngas = 4;
%  h72.gunit = [20 12 12 12]';  %% g/kg and VMR
%  h72.glist = [ 1  3  5  6]';
  
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
  iNlev  = length(airsl3_64x72.Tlevs);
  iNlevQ = length(airsl3_64x72.Qlevs);
  diffTQ = iNlev - iNlevQ;
  plevsnwp  = airsl3_64x72.Tlevs;
  plevsnwpQ = airsl3_64x72.Qlevs;
  p72.nlevs = ones(size(p72.rtime)) * iNlev;
  p72.plevs = airsl3_64x72.Tlevs' * ones(1,72*numtimesteps);
  
  p72.rlat = rlat(ii) * ones(1,72*numtimesteps);  p72.rlat = p72.rlat(:)';
  p72.rlon = rlon' * ones(1,numtimesteps);        p72.rlon = p72.rlon(:)';
  p72.plat = p72.rlat;
  p72.plon = p72.rlon;
  
  junk = airsl3_64x72.save64x72_stemp;     junk = permute(junk,[3 2 1]);   junk = squeeze(junk(:,:,ii));    junk = junk'; p72.stemp = reshape(junk,1,72*numtimesteps);
  junk = airsl3_64x72.save64x72_T;         junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.ptemp = reshape(junk,iNlev,72*numtimesteps);
  junk = airsl3_64x72.save64x72_Q;         junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.gas_1 = reshape(junk,iNlev/2,72*numtimesteps);
  junk = airsl3_64x72.save64x72_O3;        junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.gas_3 = reshape(junk,iNlev,72*numtimesteps);
%  junk = airsl3_64x72.save64x72_CO;        junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.gas_5 = reshape(junk,iNlev,72*numtimesteps);
%  junk = airsl3_64x72.save64x72_CH4;       junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.gas_6 = reshape(junk,iNlev,72*numtimesteps);
  junk = airsl3_64x72.save64x72_RH;        junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.rh    = reshape(junk,iNlev/2,72*numtimesteps);

  junk = p72.gas_1; p72.gas_1 = zeros(size(p72.ptemp)); p72.gas_1(1:iNlev/2,:) = junk; for jj = iNlev/2+1 : iNlev; frac = (iNlev-jj+1)/(iNlev/2+1); p72.gas_1(jj,:) = junk(iNlev/2,:) * frac; end;
  junk = p72.rh;    p72.rh = zeros(size(p72.ptemp));    p72.rh(1:iNlev/2,:) = junk;    for jj = iNlev/2+1 : iNlev; frac = (iNlev-jj+1)/(iNlev/2+1); p72.rh(jj,:) = junk(iNlev/2,:) * frac; end;

%  junk = p72.gas_1; p72.gas_1 = zeros(size(p72.ptemp)); p72.gas_1(diffTQ+1:iNlev,:) = junk; for jj = 1 : diffTQ; frac = (jj-1)/diffTQ; p72.gas_1(jj,:) = junk(1,:) * frac; end;
%  junk = p72.rh;    p72.rh = zeros(size(p72.ptemp));    p72.rh(diffTQ+1:iNlev,:) = junk;    for jj = 1 : diffTQ; frac = (jj-1)/diffTQ; p72.rh(jj,:) = junk(1,:) * frac; end;

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

  p72.plevs = flipud(p72.plevs);
  p72.ptemp = flipud(p72.ptemp);
  p72.gas_1 = flipud(p72.gas_1);
  p72.gas_3 = flipud(p72.gas_3);
%  p72.gas_5 = flipud(p72.gas_5);
%  p72.gas_6 = flipud(p72.gas_6);
  p72.rh    = flipud(p72.rh);

  p72.verybad = zeros(size(p72.salti));
  p72.lonbin = zeros(size(p72.salti));
  for jjj = 1 : 72
    lonbinx = (jjj-1) + (1:72:length(p72.stemp));
    p72.lonbin(lonbinx) = jjj;
  end
  iCheck = 12; %% AIRS L3 has 24 levels for T, 12 for WV
  verybad1 = find(isnan(p72.ptemp(iCheck,:)) | isnan(p72.gas_1(iCheck,:)) | isnan(p72.rh(iCheck,:)) | isnan(p72.gas_3(iCheck,:)))
  iCheck = 18; %% AIRS L3 has 24 levels for T, 12 for WV, but we have ""augmented" the latter
  verybad2 = find(isnan(p72.ptemp(iCheck,:)) | isnan(p72.gas_1(iCheck,:)) | isnan(p72.rh(iCheck,:)) | isnan(p72.gas_3(iCheck,:)))
  verybad = union(verybad1,verybad2);
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
        p72.gas_3(:,bah) = p72.gas_3(:,badah(I2-1));
        p72.rh(:,bah)    = p72.rh(:,badah(I2-1));
      else
        p72.stemp(bah)   = p72.stemp(badah(I2+1));
        p72.ptemp(:,bah) = p72.ptemp(:,badah(I2+1));
        p72.gas_1(:,bah) = p72.gas_1(:,badah(I2+1));
        p72.gas_3(:,bah) = p72.gas_3(:,badah(I2+1));
        p72.rh(:,bah)    = p72.rh(:,badah(I2+1));
      end
    end    
  end

  for jjj = 1 : length(p72.stemp)
    junklevs = p72.plevs(:,jjj);
    junkT    = p72.ptemp(:,jjj);
    junkQ    = p72.gas_1(:,jjj);
    junkZ    = p72.gas_3(:,jjj);
%    junk5    = p72.gas_5(:,jjj);
%    junk6    = p72.gas_6(:,jjj);
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

    bad = find(junkZ(1:moo) < 0 | isnan(junkZ(1:moo)));
    good = find(junkZ(1:moo) > 0);
    if length(bad) > 0
      junkZ(bad) = interp1(log(junklevs(good)),junkZ(good),log(junklevs(bad)),[],'extrap');
      haha = find(junkZ(bad) < 0);
      junkZ(bad(haha)) = junkZ(good(end));
      p72.gas_3(bad,jjj) = junkZ(bad);      
    end

%     bad = find(junk5(1:moo) < 0 | isnan(junk5(1:moo)));
%     good = find(junk5(1:moo) > 0);
%     if length(bad) > 0
%       junk5(bad) = interp1(log(junklevs(good)),junk5(good),log(junklevs(bad)),[],'extrap');
%       haha = find(junk5(bad) < 0);
%       junk5(bad(haha)) = junk5(good(end));
%       p72.gas_5(bad,jjj) = junk5(bad);      
%     end
% 
%     bad = find(junk6(1:moo) < 0 | isnan(junk6(1:moo)));
%     good = find(junk6(1:moo) > 0);
%     if length(bad) > 0
%       junk6(bad) = interp1(log(junklevs(good)),junk6(good),log(junklevs(bad)),[],'extrap');
%       haha = find(junk6(bad) < 0);
%       junk6(bad(haha)) = junk6(good(end));
%       p72.gas_6(bad,jjj) = junk6(bad);      
%     end

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
  fip = [dirout 'simulate64binsAIRSL3_' num2str(ii) fstr '.ip.rtp'];
  fop = [dirout 'simulate64binsAIRSL3_' num2str(ii) fstr '.op.rtp'];
  frp = [dirout 'simulate64binsAIRSL3_' num2str(ii) fstr '.rp.rtp'];

  klayers_sarta_check_WV_T_RH_geo_and_spectral_rates2

  iType = 3;
  plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2

  pwd
  disp(' ')
  disp('output in eg           SimulateTimeSeries/AIRSL3/reconstruct_AIRSL3_spectra_geo_rlat[01-64]_2002_09_2022_08.mat')
  disp('         also makes eg SimulateTimeSeries/AIRSL3/simulate64binsAIRSL3_[01-64]_2002_09_2022_08.[i/o/r]p.rtp');
  disp('then run driver_spectral_trends_latbin_1_64_sarta.m using   sbatch  --array=1-64 sergio_matlab_jobB.sbatch 5')
  disp(' ')

end
