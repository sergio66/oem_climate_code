%%% look_at_olr_trends.m
clear f*

%% GULP
%umbc_day_file   = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2E_day.mat';
%umbc_night_file = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2E.mat';

iV = input('Enter version of UMBC file (0/DEFAULT) for what is in test.pdf (dRHdt = 0 ... trends paper) (1) for (Held/Jeevanjee) with d lnP/dSKT = 0.02 (constant) (2) for (Held/Jeevanjee) with zonally varying d lnP/dSKT from IPCC 2007 : ');
if length(iV) == 0
  iV = 0;
end

if iV == 0
  %% SEQN
  %% dRH = 0 (Sergio)
  umbc_day_file    = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2D_day.mat';  
  umbc_night_file  = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2D.mat';
  
elseif iV == 1
  %% dRH = f(dP/dt) (Held/Jeevanjee) with d lnP/dSKT = 0.02 (constant)
  umbc_day_file    = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_day.mat';
  umbc_night_file  = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F.mat';
  
elseif iV == 2
  %% dRH = f(dP/dt) (Held/Jeevanjee) with zonally varying d lnP/dSKT from IPCC 2007
  umbc_day_file    = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_vary_dlnPdST_day.mat';
  umbc_night_file  = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_vary_dlnPdST.mat';

end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

era5_day_file    = 'ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc.mat';
era5_night_file  = 'ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc.mat';

airsL3_day_file   = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_asc.mat';
airsL3_night_file = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat';

climcapsL3_day_file   = '/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_asc.mat';
climcapsL3_night_file = '/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat';

merra2_file = 'MERRA2_atm_data_2002_09_to_2022_08_trends_desc.mat';
giss_file   = '/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ChrisHTrends/giss_trends_2002_2022.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fUMBC_day   = load(umbc_day_file,'results','resultsunc'); 
%fUMBC_night = load(umbc_night_file,'results','resultsunc'); 

fERA5_night = load(era5_night_file,'trend_stemp');
fERA5_day   = load(era5_day_file,'trend_stemp');

fAIRSL3_day   = load(airsL3_day_file,'thestats64x72');
fAIRSL3_night = load(airsL3_night_file,'thestats64x72');

fCLIMCAPSL3_day = load(climcapsL3_day_file,'thestats64x72');
fCLIMCAPSL3_night = load(climcapsL3_night_file,'thestats64x72');

fMERRA2 = load(merra2_file,'trend_stemp');
fGISS = load(giss_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk = load(era5_night_file,'trend_RH');     fERA5_night.RHtrend = junk.trend_RH;
junk = load(era5_night_file,'trend_RHSurf'); fERA5_night.RHsurftrend = junk.trend_RHSurf;
junk = load(era5_night_file,'trend_mmw');    fERA5_night.mmwtrend = junk.trend_mmw;
junk = load(era5_night_file,'trend_ptemp');  fERA5_night.ptemptrend = junk.trend_ptemp;
junk = load(era5_night_file,'trend_stemp');  fERA5_night.stemptrend = junk.trend_stemp;
junk = load(era5_night_file,'trend_gas_1');  fERA5_night.gas_1trend = junk.trend_gas_1;
pjunk = p; 
  pjunk.stemp = pjunk.stemp + fERA5_night.stemptrend;
  pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + fERA5_night.ptemptrend;
  pjunk.gas_1(1:100,:) = pjunk.gas_1(1:100,:) .* (1 + fERA5_night.gas_1trend(1:100,:));
mmwX = mmwater_rtp(h,pjunk);  
RHX  = layeramt2RH(h,pjunk);
QX = 1000 * layers2gg(h,pjunk,1:4608,1);
[mm,nn] = size(QX);
QX = Lv * QX + cp * pjunk.ptemp(1:mm,:) + g * p.palts(1:mm,:);
fERA5_night.MSEtrend                 = QX - Q0;
fERA5_night.mmwtrend_deux            = mmwX - mmw0;

junk = load(era5_day_file,'trend_RH');       fERA5_day.RHtrend = junk.trend_RH;
junk = load(era5_day_file,'trend_RHSurf');   fERA5_day.RHsurftrend = junk.trend_RHSurf;
junk = load(era5_day_file,'trend_mmw');      fERA5_day.mmwtrend = junk.trend_mmw;
junk = load(era5_day_file,'trend_ptemp');    fERA5_day.ptemptrend = junk.trend_ptemp;
junk = load(era5_day_file,'trend_stemp');    fERA5_day.stemptrend = junk.trend_stemp;
junk = load(era5_day_file,'trend_gas_1');    fERA5_day.gas_1trend = junk.trend_gas_1;
pjunk = p; 
  pjunk.stemp = pjunk.stemp + fERA5_day.stemptrend;
  pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + fERA5_day.ptemptrend;
  pjunk.gas_1(1:100,:) = pjunk.gas_1(1:100,:) .* (1 + fERA5_day.gas_1trend(1:100,:));
mmwX = mmwater_rtp(h,pjunk);  
RHX  = layeramt2RH(h,pjunk);
QX = 1000 * layers2gg(h,pjunk,1:4608,1);
[mm,nn] = size(QX);
QX = Lv * QX + cp * pjunk.ptemp(1:mm,:) + g * p.palts(1:mm,:);
fERA5_day.MSEtrend                 = QX - Q0;
fERA5_day.mmwtrend_deux            = mmwX - mmw0;

junk = load(merra2_file,'trend_RH');         fMERRA2.RHtrend = junk.trend_RH;
junk = load(merra2_file,'trend_RHSurf');     fMERRA2.RHsurftrend = junk.trend_RHSurf;
junk = load(merra2_file,'trend_mmw');        fMERRA2.mmwtrend = junk.trend_mmw;
junk = load(merra2_file,'trend_ptemp');      fMERRA2.ptemptrend = junk.trend_ptemp;
junk = load(merra2_file,'trend_stemp');      fMERRA2.stemptrend = junk.trend_stemp;
junk = load(merra2_file,'trend_gas_1');      fMERRA2.gas_1trend = junk.trend_gas_1;
pjunk = p; 
  pjunk.stemp = pjunk.stemp + fMERRA2.stemptrend;
  pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + fMERRA2.ptemptrend;
  pjunk.gas_1(1:100,:) = pjunk.gas_1(1:100,:) .* (1 + fMERRA2.gas_1trend(1:100,:));
mmwX = mmwater_rtp(h,pjunk);  
RHX  = layeramt2RH(h,pjunk);
QX = 1000 * layers2gg(h,pjunk,1:4608,1);
[mm,nn] = size(QX);
QX = Lv * QX + cp * pjunk.ptemp(1:mm,:) + g * p.palts(1:mm,:);
fMERRA2.MSEtrend                 = QX - Q0;
fMERRA2.mmwtrend_deux            = mmwX - mmw0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk = load(climcapsL3_night_file);
junk.thestats64x72.RHSurfrate = junk.thestats64x72.RHSurfrate * 100;
junk.thestats64x72.RHSurfrate(isnan(junk.thestats64x72.RHSurfrate)) = 0;
junk.thestats64x72.RHrate(isnan(junk.thestats64x72.RHrate)) = 0;
junk.thestats64x72.stemprate(isnan(junk.thestats64x72.stemprate)) = 0;
junk.thestats64x72.stempratestd(isnan(junk.thestats64x72.stempratestd)) = 0;
junk.thestats64x72.ptemprate(isnan(junk.thestats64x72.ptemprate)) = 0;
junk.thestats64x72.ptempratestd(isnan(junk.thestats64x72.ptempratestd)) = 0;
junk.thestats64x72.waterrate(isnan(junk.thestats64x72.waterrate)) = 0;
junk.thestats64x72.waterratestd(isnan(junk.thestats64x72.waterratestd)) = 0;
pjunk = p; 
  pjunk.stemp = pjunk.stemp + reshape(junk.thestats64x72.stemprate,1,4608);
  pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + (reshape(junk.thestats64x72.ptemprate,4608,100)');
  pjunk.gas_1((1:066)+34,:) = pjunk.gas_1((1:066)+34,:) .* (1 + (reshape(junk.thestats64x72.waterrate,4608,066)'));
mmwX = mmwater_rtp(h,pjunk);  
QX = 1000 * layers2gg(h,pjunk,1:4608,1);
[mm,nn] = size(QX);
QX = Lv * QX + cp * pjunk.ptemp(1:mm,:) + g * p.palts(1:mm,:);

fCLIMCAPSL3_night.MSEtrend                 = QX - Q0;
fCLIMCAPSL3_night.mmwtrend                 = mmwX - mmw0;
fCLIMCAPSL3_night.stemptrend               = reshape(junk.thestats64x72.stemprate,1,4608);
fCLIMCAPSL3_night.ptemptrend               = (reshape(junk.thestats64x72.ptemprate,4608,100)');
fCLIMCAPSL3_night.gas_1trend((1:066)+34,:) = (reshape(junk.thestats64x72.waterrate,4608,066)');
fCLIMCAPSL3_night.RHtrend((1:066)+34,:)    = (reshape(junk.thestats64x72.RHrate,4608,066)')*100;
fCLIMCAPSL3_night.RHsurftrend              = junk.thestats64x72.RHSurfrate;
for junkll = 1 : 4608
  nlays = p.nlevs(junkll)-1;
  fCLIMCAPSL3_night.RHsurftrend(junkll) = 0.5*(fCLIMCAPSL3_night.RHtrend(nlays,junkll)+fCLIMCAPSL3_night.RHtrend(nlays-1,junkll));
end

junk = load(climcapsL3_day_file);
junk.thestats64x72.RHSurfrate = junk.thestats64x72.RHSurfrate * 100;
junk.thestats64x72.RHSurfrate(isnan(junk.thestats64x72.RHSurfrate)) = 0;
junk.thestats64x72.RHrate(isnan(junk.thestats64x72.RHrate)) = 0;
junk.thestats64x72.stemprate(isnan(junk.thestats64x72.stemprate)) = 0;
junk.thestats64x72.stempratestd(isnan(junk.thestats64x72.stempratestd)) = 0;
junk.thestats64x72.ptemprate(isnan(junk.thestats64x72.ptemprate)) = 0;
junk.thestats64x72.ptempratestd(isnan(junk.thestats64x72.ptempratestd)) = 0;
junk.thestats64x72.waterrate(isnan(junk.thestats64x72.waterrate)) = 0;
junk.thestats64x72.waterratestd(isnan(junk.thestats64x72.waterratestd)) = 0;
pjunk = p; 
  pjunk.stemp = pjunk.stemp + reshape(junk.thestats64x72.stemprate,1,4608);
  pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + (reshape(junk.thestats64x72.ptemprate,4608,100)');
  pjunk.gas_1((1:066)+34,:) = pjunk.gas_1((1:066)+34,:) .* (1 + (reshape(junk.thestats64x72.waterrate,4608,066)'));
mmwX = mmwater_rtp(h,pjunk);  
QX = 1000 * layers2gg(h,pjunk,1:4608,1);
[mm,nn] = size(QX);
QX = Lv * QX + cp * pjunk.ptemp(1:mm,:) + g * p.palts(1:mm,:);

fCLIMCAPSL3_day.MSEtrend                 = QX - Q0;
fCLIMCAPSL3_day.mmwtrend                 = mmwX - mmw0;
fCLIMCAPSL3_day.stemptrend               = reshape(junk.thestats64x72.stemprate,1,4608);
fCLIMCAPSL3_day.ptemptrend               = (reshape(junk.thestats64x72.ptemprate,4608,100)');
fCLIMCAPSL3_day.gas_1trend((1:066)+34,:) = (reshape(junk.thestats64x72.waterrate,4608,066)');
fCLIMCAPSL3_day.RHtrend((1:066)+34,:)    = (reshape(junk.thestats64x72.RHrate,4608,066)')*100;
fCLIMCAPSL3_day.RHsurftrend              = junk.thestats64x72.RHSurfrate;
for junkll = 1 : 4608
  nlays = p.nlevs(junkll)-1;
  fCLIMCAPSL3_day.RHsurftrend(junkll) = 0.5*(fCLIMCAPSL3_day.RHtrend(nlays,junkll)+fCLIMCAPSL3_day.RHtrend(nlays-1,junkll));
end

clear junk pjunk junk1 junk2 mmwX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk = load(airsL3_night_file);
junk.thestats64x72.RHSurfrate(isnan(junk.thestats64x72.RHSurfrate)) = 0;
junk.thestats64x72.RHrate(isnan(junk.thestats64x72.RHrate)) = 0;
junk.thestats64x72.stemprate(isnan(junk.thestats64x72.stemprate)) = 0;
junk.thestats64x72.stempratestd(isnan(junk.thestats64x72.stempratestd)) = 0;
junk.thestats64x72.ptemprate(isnan(junk.thestats64x72.ptemprate)) = 0;
junk.thestats64x72.ptempratestd(isnan(junk.thestats64x72.ptempratestd)) = 0;
junk.thestats64x72.waterrate(isnan(junk.thestats64x72.waterrate)) = 0;
junk.thestats64x72.waterratestd(isnan(junk.thestats64x72.waterratestd)) = 0;
for ii = 1 : 72
  for jj = 1 : 64
    junk1 = squeeze(junk.thestats64x72.waterrate(ii,jj,:));
    junk2 = interp1(log(junk.Qlevs),junk1,log(plays100),[],'extrap');
    junk.thestats64x72.waterrate100(ii,jj,:) = junk2;

    junk1 = squeeze(junk.thestats64x72.RHrate(ii,jj,:));
    junk2 = interp1(log(junk.Qlevs),junk1,log(plays100),[],'extrap');
    junk.thestats64x72.RHrate100(ii,jj,:) = junk2;

    junk1 = squeeze(junk.thestats64x72.ptemprate(ii,jj,:));
    junk2 = interp1(log(junk.Tlevs),junk1,log(plays100),[],'extrap');
    junk.thestats64x72.ptemprate100(ii,jj,:) = junk2;
  end
end
pjunk = p; 
  pjunk.stemp = pjunk.stemp + reshape(junk.thestats64x72.stemprate,1,4608);
  pjunk.ptemp(1:100,:)  = pjunk.ptemp(1:100,:) + (reshape(junk.thestats64x72.ptemprate100,4608,100)');
  pjunk.gas_1(1:100,:) = pjunk.gas_1(1:100,:) .* (1 + (reshape(junk.thestats64x72.waterrate100,4608,100)'));
mmwX = mmwater_rtp(h,pjunk);
QX = 1000 * layers2gg(h,pjunk,1:4608,1);
[mm,nn] = size(QX);
QX = Lv * QX + cp * pjunk.ptemp(1:mm,:) + g * p.palts(1:mm,:);

fAIRSL3_night.MSEtrend    = QX - Q0;  
fAIRSL3_night.mmwtrend    = mmwX - mmw0;
fAIRSL3_night.ptemptrend  = (reshape(junk.thestats64x72.ptemprate100,4608,100)');
fAIRSL3_night.gas_1trend  = (reshape(junk.thestats64x72.waterrate100,4608,100)');
fAIRSL3_night.stemptrend  = reshape(junk.thestats64x72.stemprate,1,4608);
fAIRSL3_night.RHtrend     = (reshape(junk.thestats64x72.RHrate100,4608,100)');
fAIRSL3_night.RHsurftrend = junk.thestats64x72.RHSurfrate;

junk = load(airsL3_day_file);
junk.thestats64x72.RHSurfrate(isnan(junk.thestats64x72.RHSurfrate)) = 0;
junk.thestats64x72.RHrate(isnan(junk.thestats64x72.RHrate)) = 0;
junk.thestats64x72.stemprate(isnan(junk.thestats64x72.stemprate)) = 0;
junk.thestats64x72.stempratestd(isnan(junk.thestats64x72.stempratestd)) = 0;
junk.thestats64x72.ptemprate(isnan(junk.thestats64x72.ptemprate)) = 0;
junk.thestats64x72.ptempratestd(isnan(junk.thestats64x72.ptempratestd)) = 0;
junk.thestats64x72.waterrate(isnan(junk.thestats64x72.waterrate)) = 0;
junk.thestats64x72.waterratestd(isnan(junk.thestats64x72.waterratestd)) = 0;
for ii = 1 : 72
  for jj = 1 : 64
    junk1 = squeeze(junk.thestats64x72.waterrate(ii,jj,:));
    junk2 = interp1(log(junk.Qlevs),junk1,log(plays100),[],'extrap');
    junk.thestats64x72.waterrate100(ii,jj,:) = junk2;

    junk1 = squeeze(junk.thestats64x72.RHrate(ii,jj,:));
    junk2 = interp1(log(junk.Qlevs),junk1,log(plays100),[],'extrap');
    junk.thestats64x72.RHrate100(ii,jj,:) = junk2;

    junk1 = squeeze(junk.thestats64x72.ptemprate(ii,jj,:));
    junk2 = interp1(log(junk.Tlevs),junk1,log(plays100),[],'extrap');
    junk.thestats64x72.ptemprate100(ii,jj,:) = junk2;
  end
end
pjunk = p; 
  pjunk.stemp = pjunk.stemp + reshape(junk.thestats64x72.stemprate,1,4608);
  pjunk.ptemp(1:100,:)  = pjunk.ptemp(1:100,:) + (reshape(junk.thestats64x72.ptemprate100,4608,100)');
  pjunk.gas_1(1:100,:) = pjunk.gas_1(1:100,:) .* (1 + (reshape(junk.thestats64x72.waterrate100,4608,100)'));
mmwX = mmwater_rtp(h,pjunk);  
QX = 1000 * layers2gg(h,pjunk,1:4608,1);
[mm,nn] = size(QX);
QX = Lv * QX + cp * pjunk.ptemp(1:mm,:) + g * p.palts(1:mm,:);

fAIRSL3_day.MSEtrend    = QX - Q0;  
fAIRSL3_day.mmwtrend    = mmwX - mmw0;
fAIRSL3_day.stemptrend  = reshape(junk.thestats64x72.stemprate,1,4608);
fAIRSL3_day.ptemptrend  = (reshape(junk.thestats64x72.ptemprate100,4608,100)');
fAIRSL3_day.gas_1trend  = (reshape(junk.thestats64x72.waterrate100,4608,100)');
fAIRSL3_day.RHtrend     = (reshape(junk.thestats64x72.RHrate100,4608,100)');
fAIRSL3_day.RHsurftrend = junk.thestats64x72.RHSurfrate;

clear junk pjunk junk1 junk2 mmwX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fUMBC_day   = load(umbc_day_file,   'results','resultsT','resultsWV','resultsunc','resultsTunc','resultsWVunc','pavg'); %% SEQN
fUMBC_night = load(umbc_night_file','results','resultsT','resultsWV','resultsunc','resultsTunc','resultsWVunc','pavg'); %% SEQN

for ii = 1 : 4608
  junk1 = fUMBC_day.resultsT(ii,:);
  junk2 = interp1(log(fUMBC_day.pavg),junk1,log(plays100),[],'extrap');
  junk.resultsT100(:,ii) = junk2;

  junk1 = fUMBC_day.resultsWV(ii,:);
  junk2 = interp1(log(fUMBC_day.pavg),junk1,log(plays100),[],'extrap');
  junk.resultsWV100(:,ii) = junk2;

  junk1 = fUMBC_day.resultsTunc(ii,:);
  junk2 = interp1(log(fUMBC_day.pavg),junk1,log(plays100),[],'extrap');
  junk.resultsTunc100(:,ii) = junk2;

  junk1 = fUMBC_day.resultsWVunc(ii,:);
  junk2 = interp1(log(fUMBC_day.pavg),junk1,log(plays100),[],'extrap');
  junk.resultsWVunc100(:,ii) = junk2;
end
pjunk = p; 
  pjunk.stemp = pjunk.stemp + fUMBC_day.results(:,6);
  pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + junk.resultsT100(1:100,:);
  pjunk.gas_1(1:100,:) = pjunk.gas_1(1:100,:) .* (1 + junk.resultsWV100(1:100,:));
mmwX = mmwater_rtp(h,pjunk);  
RHX  = layeramt2RH(h,pjunk);
QX = 1000 * layers2gg(h,pjunk,1:4608,1);
[mm,nn] = size(QX);
QX = Lv * QX + cp * pjunk.ptemp(1:mm,:) + g * p.palts(1:mm,:);

fUMBC_day.MSEtrend = QX - Q0;
fUMBC_day.mmwtrend = mmwX - mmw0;
fUMBC_day.ptemptrend = junk.resultsT100(1:100,:);
fUMBC_day.gas_1trend = junk.resultsWV100(1:100,:);
fUMBC_day.RHtrend    = RHX - RH0;
pjunk = p; 
  pjunk.stemp = pjunk.stemp + fUMBC_day.results(:,6) + fUMBC_day.resultsunc(:,6);
  pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + junk.resultsT100(1:100,:) + junk.resultsTunc100(1:100,:);
  pjunk.gas_1(1:100,:) = pjunk.gas_1(1:100,:) .* (1 + junk.resultsWV100(1:100,:) + junk.resultsWVunc100(1:100,:));
mmwX2 = mmwater_rtp(h,pjunk);  
RHX2  = layeramt2RH(h,pjunk);
fUMBC_day.mmwtrendunc = abs(mmwX - mmwX2);
fUMBC_day.RHtrendunc  = abs(RHX - RHX2);
fUMBC_day.ptemptrendunc = junk.resultsT100(1:100,:);
fUMBC_day.gas_1trendunc = junk.resultsWV100(1:100,:);
for junkll = 1 : 4608
  nlays = p.nlevs(junkll)-1;
  fUMBC_day.RHsurftrend(junkll)    = fUMBC_day.RHtrend(nlays,junkll);
  fUMBC_day.RHsurftrendunc(junkll) = fUMBC_day.RHtrendunc(nlays,junkll);

  fUMBC_day.RHsurftrend(junkll)    = 0.5*(fUMBC_day.RHtrend(nlays,junkll)+fUMBC_day.RHtrend(nlays-1,junkll));
  fUMBC_day.RHsurftrendunc(junkll) = 0.5*(fUMBC_day.RHtrendunc(nlays,junkll)+fUMBC_day.RHtrendunc(nlays-1,junkll));
end

%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 4608
  junk1 = fUMBC_night.resultsT(ii,:);
  junk2 = interp1(log(fUMBC_night.pavg),junk1,log(plays100),[],'extrap');
  junk.resultsT100(:,ii) = junk2;

  junk1 = fUMBC_night.resultsWV(ii,:);
  junk2 = interp1(log(fUMBC_night.pavg),junk1,log(plays100),[],'extrap');
  junk.resultsWV100(:,ii) = junk2;

  junk1 = fUMBC_night.resultsTunc(ii,:);
  junk2 = interp1(log(fUMBC_night.pavg),junk1,log(plays100),[],'extrap');
  junk.resultsTunc100(:,ii) = junk2;

  junk1 = fUMBC_night.resultsWVunc(ii,:);
  junk2 = interp1(log(fUMBC_night.pavg),junk1,log(plays100),[],'extrap');
  junk.resultsWVunc100(:,ii) = junk2;
end
pjunk = p; 
  pjunk.stemp = pjunk.stemp + fUMBC_night.results(:,6);
  pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + junk.resultsT100(1:100,:);
  pjunk.gas_1(1:100,:) = pjunk.gas_1(1:100,:) .* (1 + junk.resultsWV100(1:100,:));
mmwX = mmwater_rtp(h,pjunk);  
RHX  = layeramt2RH(h,pjunk);
QX = 1000 * layers2gg(h,pjunk,1:4608,1);
[mm,nn] = size(QX);
QX = Lv * QX + cp * pjunk.ptemp(1:mm,:) + g * p.palts(1:mm,:);

fUMBC_night.MSEtrend = QX - Q0;
fUMBC_night.mmwtrend = mmwX - mmw0;
fUMBC_night.ptemptrend = junk.resultsT100(1:100,:);
fUMBC_night.gas_1trend = junk.resultsWV100(1:100,:);
fUMBC_night.RHtrend    = RHX - RH0;
pjunk = p; 
  pjunk.stemp = pjunk.stemp + fUMBC_night.results(:,6) + fUMBC_night.resultsunc(:,6);
  pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + junk.resultsT100(1:100,:) + junk.resultsTunc100(1:100,:);
  pjunk.gas_1(1:100,:) = pjunk.gas_1(1:100,:) .* (1 + junk.resultsWV100(1:100,:) + junk.resultsWVunc100(1:100,:));
mmwX2 = mmwater_rtp(h,pjunk);  
RHX2  = layeramt2RH(h,pjunk);
fUMBC_night.mmwtrendunc = abs(mmwX - mmwX2);
fUMBC_night.RHtrendunc  = abs(RHX - RHX2);
fUMBC_night.ptemptrendunc = junk.resultsT100(1:100,:);
fUMBC_night.gas_1trendunc = junk.resultsWV100(1:100,:);
for junkll = 1 : 4608
  nlays = p.nlevs(junkll)-1;
  fUMBC_night.RHsurftrend(junkll)    = fUMBC_night.RHtrend(nlays,junkll);
  fUMBC_night.RHsurftrendunc(junkll) = fUMBC_night.RHtrendunc(nlays,junkll);

  fUMBC_night.RHsurftrend(junkll)    = 0.5*(fUMBC_night.RHtrend(nlays,junkll)+fUMBC_night.RHtrend(nlays-1,junkll));
  fUMBC_night.RHsurftrendunc(junkll) = 0.5*(fUMBC_night.RHtrendunc(nlays,junkll)+fUMBC_night.RHtrendunc(nlays-1,junkll));

end

clear junk* mmwX* pjunk
