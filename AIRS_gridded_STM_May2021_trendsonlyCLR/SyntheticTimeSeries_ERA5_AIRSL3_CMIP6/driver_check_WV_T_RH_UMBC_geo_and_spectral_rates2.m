%% this tries to loop over the 64 zonal bins using the cluster
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

addpath  /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/Strow_humidity/convert_humidity/

system_slurm_stats

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 64 for the 64 latbins
%JOB = 31

%load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat
wah = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat','h');
h = wah.h;
wah = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat','p');
p = wah.p;

load('llsmap5.mat');

RH000 = layeramt2RH(h,p);

pUMBC = p; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hkcarta_emis,~,kcarta_emis,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp');
ind = (1:72) + (JOB-1)*72;
[hkcarta_emism,kcarta_emis] = subset_rtp(hkcarta_emis,kcarta_emis,[],[],ind);

%% see  FIND_NWP_MODEL_TRENDS/driver_computeERA5_monthly_trends.m  and do_the_AIRSL3_trends.m
disp('if you get silly messages like "YM timeperiod  = 2002/ 9 --> 2022/ 8 needs 240 of 228 timesteps" then check this >>>>>>>>')
disp('if you get silly messages like "YM timeperiod  = 2002/ 9 --> 2022/ 8 needs 240 of 228 timesteps" then check this >>>>>>>>')

iYS = 2002; iYE = 2021;
iYS = 2002; iYE = 2022;

%% see  FIND_NWP_MODEL_TRENDS/driver_computeERA5_monthly_trends.m  and do_the_AIRSL3_trends.m
%% this will have to be modified as per my UMBC retrievals!!!!!!!
%era5_64x72 = load('../../FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_07_desc.mat');
%era5_64x72 = load('../../FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_08_desc.mat');
era5_64x72 = load('../../FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2022_08_desc.mat');

umbc = load('/asl/s1/sergio/JUNK/test9_guessstartWV_Vers1_march22_2023.mat','results','resultsO3','resultsT','resultsWV','pjunk20');
plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
plays = flipud(plevs2plays(plevs));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[numtimesteps0,~] = size(era5_64x72.all.mmw);
numtimesteps = numtimesteps0;
fprintf(1,'driver_check_WV_T_RH_UMBC_geo_and_spectral_rates2.m : numtimesteps = %3i \n',numtimesteps)

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

daysSince2002Start = change2days(YMStart(1),YMStart(2),15,2002);
daysSince2002End   = change2days(YMEnd(1),  YMEnd(2),  15,2002);

usethese = find(daysSince2002  >= daysSince2002Start & daysSince2002 <= daysSince2002End);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,' YM timeperiod  = %4i/%2i --> %4i/%2i needs %3i of %3i timesteps \n',YMStart,YMEnd,length(usethese),numtimesteps)
yy = yy(usethese);
mm = mm(usethese);
dd = dd(usethese);

numtimesteps = length(yy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rtime = utc2taiSergio(yy,mm,dd,ones(size(yy))*12.0);
time_so_far = (yy-2000) + ((mm-1)+1)/12;

co2ppm = 370 + 2.2*((yy+mm/12)-2002);

%% see ~/MATLABCODE/CRODGERS_FAST_CLOUD/driver_stage2_ESRL_set_CO2_CH4_N2O.m
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
dirout = 'SimulateTimeSeries/UMBC/';
if iConstORVary == +1
  dirout = 'SimulateTimeSeries/UMBC_ConstTracegas/';
end

time_slope = time_so_far-mean(time_so_far);

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
  h72.gunit = [21 21]';  %% g/g
  h72.glist = [ 1 3]';
  
  p72.rtime = [];
  p72.co2ppm = [];
  for iii = 1 : numtimesteps
    nemis_t  = [nemis_t  kcarta_emis.nemis];
    efreq_t  = [efreq_t  kcarta_emis.efreq];
    emis_t   = [emis_t   kcarta_emis.emis];
    rho_t    = [rho_t    kcarta_emis.rho];

    co2ppm_t = [co2ppm_t ones(1,72)*co2ppm(iii)];
    n2oppm_t = [n2oppm_t ones(1,72)*n2oppm(iii)];
    ch4ppm_t = [ch4ppm_t ones(1,72)*ch4ppm(iii)];

    p72.rtime  = [p72.rtime ones(1,72)*rtime(iii)];
    p72.co2ppm = [p72.co2ppm ones(1,72)*co2ppm(iii)];
  end
  
  iNlev = 37;
  plevsnwp  = squeeze(era5_64x72.all.nwp_plevs(1,:,3000))';
  p72.nlevs = ones(size(p72.rtime)) * iNlev;
  p72.plevs = squeeze(era5_64x72.all.nwp_plevs(1,:,3000))' * ones(1,72*numtimesteps);
  
  p72.rlat  = rlat(ii) * ones(1,72*numtimesteps); p72.rlat = p72.rlat(:)';
  p72.rlon = rlon' * ones(1,numtimesteps);        p72.rlon = p72.rlon(:)';
  p72.plat = p72.rlat;
  p72.plon = p72.rlon;

  junk = era5_64x72.all.stemp;                  junk = reshape( junk,numtimesteps,72,64);    junk = squeeze( junk(:,:,ii));  junk =  junk';  rara = mean(junk,2);
  xjunk = time_slope' * (umbc.results(:,6)');  xjunk = reshape(xjunk,numtimesteps,72,64);   xjunk = squeeze(xjunk(:,:,ii)); xjunk = xjunk';
  bad = find(isnan(xjunk)); xjunk(bad) = 0;
  for lonlon = 1 : 72
    yjunk(lonlon,:) = rara(lonlon) + xjunk(lonlon,:);
  end
  p72.stemp = reshape(yjunk,1,72*numtimesteps);

  junk = era5_64x72.all.nwp_plevs;               junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); eraplevs = squeeze(nanmean(junk,1));

  junk = era5_64x72.all.nwp_ptemp;               junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); rara = squeeze(nanmean(junk,1));
  trend = umbc.resultsT';                        trend = reshape(trend,49,72,64);               trend = squeeze(trend(:,:,ii)); 
  for kkkk = 1 : 72
    xtrend(:,kkkk) = interp1(log(umbc.pjunk20),trend(:,kkkk),log(eraplevs(:,kkkk)),[],'extrap');
  end
  bad = find(isnan(xtrend)); xtrend(bad) = 0;
  xjunk = zeros([size(rara) length(time_slope)]);
  for kkkk = 1 : 37
    for llll = 1 : 72
      xjunk(kkkk,llll,:) = rara(kkkk,llll) + xtrend(kkkk,llll)*time_slope;
    end
  end
  wah = squeeze(xjunk(:,1,:));
  figure(1); pcolor(1:72,p72.plevs(:,1),xtrend); colorbar; shading interp; colormap(usa2); colorbar; set(gca,'ydir','reverse')
  figure(2); pcolor(1:72,p72.plevs(:,1),rara); colorbar; shading interp; colormap(jet); colorbar; set(gca,'ydir','reverse')
  figure(3); pcolor(1:240,p72.plevs(:,1),wah); colorbar; shading interp; colormap(jet); colorbar; set(gca,'ydir','reverse')
  figure(4); pcolor(1:240,p72.plevs(:,1),nanmean(wah')'-wah); colorbar; shading interp; colormap(usa2); colorbar; set(gca,'ydir','reverse')
  p72.ptemp = reshape(xjunk,iNlev,72*numtimesteps);

  junk = era5_64x72.all.nwp_gas_1;               junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); rara = squeeze(nanmean(junk,1));
  trend = umbc.resultsWV';                       trend = reshape(trend,49,72,64);               trend = squeeze(trend(:,:,ii)); 
  for kkkk = 1 : 72
    xtrend(:,kkkk) = interp1(log(umbc.pjunk20),trend(:,kkkk),log(eraplevs(:,kkkk)),[],'extrap');
  end
  bad = find(isnan(xtrend)); xtrend(bad) = 0;
  xjunk = zeros([size(rara) length(time_slope)]);
  for kkkk = 1 : 37
    for llll = 1 : 72
      xjunk(kkkk,llll,:) = rara(kkkk,llll) * (1 + xtrend(kkkk,llll)*time_slope);
    end
  end
  wah = squeeze(xjunk(:,1,:));
  figure(1); pcolor(1:72,p72.plevs(:,1),xtrend); colorbar; shading interp; colormap(usa2); colorbar; set(gca,'ydir','reverse')
  figure(2); pcolor(1:72,p72.plevs(:,1),rara); colorbar; shading interp; colormap(jet); colorbar; set(gca,'ydir','reverse')
  figure(3); pcolor(1:240,p72.plevs(:,1),wah); colorbar; shading interp; colormap(jet); colorbar; set(gca,'ydir','reverse')
  figure(4); pcolor(1:240,p72.plevs(:,1),nanmean(wah')'-wah); colorbar; shading interp; colormap(usa2); colorbar; set(gca,'ydir','reverse')
  p72.gas_1 = reshape(xjunk,iNlev,72*numtimesteps);

  junk = era5_64x72.all.nwp_gas_3;               junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); rara = squeeze(nanmean(junk,1));
  trend = umbc.resultsO3';                       trend = reshape(trend,49,72,64);               trend = squeeze(trend(:,:,ii)); 
  for kkkk = 1 : 72
    xtrend(:,kkkk) = interp1(log(umbc.pjunk20),trend(:,kkkk),log(eraplevs(:,kkkk)),[],'extrap');
  end
  bad = find(isnan(xtrend)); xtrend(bad) = 0;
  xjunk = zeros([size(rara) length(time_slope)]);
  for kkkk = 1 : 37
    for llll = 1 : 72
      xjunk(kkkk,llll,:) = rara(kkkk,llll) * (1 + xtrend(kkkk,llll)*time_slope);
    end
  end
  wah = squeeze(xjunk(:,1,:));
  figure(1); pcolor(1:72,p72.plevs(:,1),xtrend); colorbar; shading interp; colormap(usa2); colorbar; set(gca,'ydir','reverse')
  figure(2); pcolor(1:72,p72.plevs(:,1),rara); colorbar; shading interp; colormap(jet); colorbar; set(gca,'ydir','reverse')
  figure(3); pcolor(1:240,p72.plevs(:,1),wah); colorbar; shading interp; colormap(jet); colorbar; set(gca,'ydir','reverse')
  figure(4); pcolor(1:240,p72.plevs(:,1),nanmean(wah')'-wah); colorbar; shading interp; colormap(usa2); colorbar; set(gca,'ydir','reverse')
  p72.gas_3 = reshape(xjunk,iNlev,72*numtimesteps);

% function out = convert_humidity (P, T, in, type_in, type_out, method_es, tol_T)
% pressure in Pa (not hPa nor mbar)
% temperature in K (not in degree Celsius)
% humidity:
% - partial pressure of water vapor in Pa (not hPa nor mbar)
% - specific humidity in kg/kg (not g/kg)
% - mixing ratio in kg/kg (not g/kg)
% - relative humidity in percent
% - dew point temperature in K (not degree Celsius)
% - virtual temperature in K (not degree Celsius)
  p72.rh = convert_humidity(p72.plevs*100,p72.ptemp,p72.gas_1,'specific humidity','relative humidity');
  
%  p72.scanang = zeros(size(p72.stemp));
%  p72.satzen = zeros(size(p72.stemp));
  p72.zobs = 705000 * ones(size(p72.stemp));
  p72.scanang = ones(size(p72.stemp)) * 22;
  p72.satzen = vaconv(p72.scanang, p72.zobs, zeros(size(p72.zobs)));
  p72.satzen = ones(size(p72.stemp)) * 24; %%% to match what is in home/sergio/KCARTA/WORK/RUN/RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp which was used to make jacs
  p72.solzen = 150 * ones(size(p72.stemp));
  %p72.spres = 1000 * ones(size(p72.stemp));
  %p72.salti = 0 * ones(size(p72.stemp));
  p72.spres = reshape(p.spres,72,64); p72.spres = p72.spres(:,ii) * ones(1,1*numtimesteps); p72.spres = p72.spres(:)';
  p72.salti = reshape(p.salti,72,64); p72.salti = p72.salti(:,ii) * ones(1,1*numtimesteps); p72.salti = p72.salti(:)';

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

  p72.nemis = nemis_t;
  p72.efreq = efreq_t;
  p72.emis  = emis_t;
  p72.rho   = rho_t;
    
  fstr = ['_' num2str(YMStart(1),'%04d') '_' num2str(YMStart(2),'%02d') '_' num2str(YMEnd(1),'%04d') '_' num2str(YMEnd(2),'%02d')];
  fip = [dirout 'simulate64binsUMBC_' num2str(ii) fstr '.ip.rtp'];
  fop = [dirout 'simulate64binsUMBC_' num2str(ii) fstr '.op.rtp'];
  frp = [dirout 'simulate64binsUMBC_' num2str(ii) fstr '.rp.rtp'];

  klayers_sarta_check_WV_T_RH_geo_and_spectral_rates2

  iType = 1;
  plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2

  pwd
  disp(' ')
  disp('output in eg           SimulateTimeSeries/UMBC/reconstruct_umbc_spectra_geo_rlat[01-64]_2002_09_2022_08.mat')
  disp('         also makes eg SimulateTimeSeries/UMBC/simulate64binsUMBC_[01-64]_2002_09_2022_08.[i/o/r]p.rtp');
  disp('then run driver_spectral_trends_latbin_1_64_sarta.m using   sbatch  --array=1-64 sergio_matlab_jobB.sbatch 12')
  disp(' ')

end