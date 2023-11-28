addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTMISC/

iMin = 0;
if exist('summary_mmw_stemp_corr_umbc_era5.mat');
  load summary_mmw_stemp_corr_umbc_era5.mat
  [iMin,~] = size(thesummary.umbc_mean(:,1))
  iMin
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii = 0; 
ii = ii + 1; 
  strUMBC{ii} = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q16_newERA5_2021jacs_startwith0_50fatlayers_ERA5calcs_spectraltrends.mat'; comment{ii} = 'GULP use SPECTRAL ERA5 to test code';
ii = ii + 1; 
  strUMBC{ii} = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2.mat';         comment{ii} = 'GULP this is default GULP used at Sounder STM Oct 2023';
ii = ii + 1; 
  strUMBC{ii} = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_MLS.mat';     comment{ii} = 'GULP this is default GULP + MLS';
ii = ii + 1; 
  strUMBC{ii} = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2.mat';   
  comment{ii} = 'SEQN : correct channels, test2A should be same as test2 but slope at S.pole maybe 3 (for test2) and 2(for test2A) : so slope probably -2/31';
ii = ii + 1; 
  strUMBC{ii} = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2A.mat';   
  comment{ii} = 'SEQN : correct channels, test2A should be same as test2 but slope at S.pole maybe 3 (for test2) and 2(for test2A) : so slope probably -2/31';
ii = ii + 1; 
  strUMBC{ii} = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test3_xb_WV_000000.mat'; comment{ii} = 'SEQN, with xb(WV) = 0';
ii = ii + 1; 
  strUMBC{ii} = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test3_xb_WV_000000.mat'; comment{ii} = 'GULP, with xb(WV) = 0';
ii = ii + 1; 
  strUMBC{ii} = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2B.mat';             comment{ii} = 'SEQN: same as test2A but has iAdjLowerAtmWVfrac = 0.0625, slope = -1/31';
ii = ii + 1; 
  strUMBC{ii} = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2C.mat';             comment{ii} = 'SEQN: same as test2B but has iAdjLowerAtmWVfrac = 1.0000, slope = -1/31';
ii = ii + 1; 
  strUMBC{ii} = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2D.mat';             comment{ii} = 'SEQN: same as test2C but has iAdjLowerAtmWVfrac = 1.0000, slope = 0';
ii = ii + 1; 
  strUMBC{ii} = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2E.mat';             comment{ii} = 'GULP: same as test2D with iAdjLowerAtmWVfrac = 1.0000, slope = 0';

xstrstr = {'GULP ER5 test','GULP SounderSTM10/23','GULP with MLS','SEQN correct channels test2','SEQN correct channels test2A','SEQN xb(WV)=0','GULP xb(WV)=0',...
           'SEQN test3,iAdj=0.0625,slope=-1/31','SEQN test3,iAdj=1,slope=-1/31','SEQN test3,iAdj=0.0625,slope=0','GULP test3,iAdj=0.0625,slope=0'};
if length(xstrstr) ~= length(strUMBC)  
  [ii length(strUMBC) length(comment) length(xstrstr)]
  error('need to have length(xstrstr) == length(strUMBC)')
end
                          
for ii = 1 : length(strUMBC) 
  boo{ii} = load(strUMBC{ii},'topts');
end
for ii = 1 : length(strUMBC)
  fprintf(1,'conparing settings for topts=2 DEFAULT GULF    vs  topts=%2i "%s"\n',ii,comment{ii});
  compare_two_structures(boo{2}.topts,boo{ii}.topts);
  disp('ret to continue'); pause
end

for ii = iMin+1 : length(strUMBC) 
  fprintf(1,'%2i of %2i : comment %s \n using strUMBC = %s \n',ii,length(strUMBC),comment{ii},strUMBC{ii})
  [era5,merra2,airsL3,climcapsL3,umbc,thecorr,amp,saverates_rlat_pres,summary_olr_mmw_stats] = driver_show_AIRSV7_L3_vs_CLIMCAPS_vs_MERRA2_vs_ERA5_trends(strUMBC{ii},20,+1,comment{ii});

  thesummary.umbc_era5_corr(ii,:) = summary_olr_mmw_stats.umbc.allXchi;
  thesummary.umbc_era5_mean(ii,:) = summary_olr_mmw_stats.umbc.allXmean;
  thesummary.umbc_era5_std(ii,:)  = summary_olr_mmw_stats.umbc.allXstd;
  thesummary.umbc_era5_frac(ii,:) = summary_olr_mmw_stats.umbc.allX_frac_neg0pos;
  thesummary.umbc_era5_corr_mmw(ii,:) = summary_olr_mmw_stats.corr_mmw;

  ifoo = 0;
  thesummary.umbc_mean(ii,1) = summary_olr_mmw_stats.umbc.stemp(1);
  thesummary.umbc_mean(ii,2) = summary_olr_mmw_stats.umbc.mmw(1);
  thesummary.umbc_mean(ii,3) = summary_olr_mmw_stats.umbc.mmwpc(1);
  thesummary.era5_mean(ii,1) = summary_olr_mmw_stats.era5.stemp(1);
  thesummary.era5_mean(ii,2) = summary_olr_mmw_stats.era5.mmw(1);
  thesummary.era5_mean(ii,3) = summary_olr_mmw_stats.era5.mmwpc(1);
  thesummary.merra2_mean(ii,1) = summary_olr_mmw_stats.merra2.stemp(1);
  thesummary.merra2_mean(ii,2) = summary_olr_mmw_stats.merra2.mmw(1);
  thesummary.merra2_mean(ii,3) = summary_olr_mmw_stats.merra2.mmwpc(1); 

  thesummary.umbc_std(ii,1) = summary_olr_mmw_stats.umbc.stemp(2);
  thesummary.umbc_std(ii,2) = summary_olr_mmw_stats.umbc.mmw(2);
  thesummary.umbc_std(ii,3) = summary_olr_mmw_stats.umbc.mmwpc(2);
  thesummary.era5_std(ii,1) = summary_olr_mmw_stats.era5.stemp(2);
  thesummary.era5_std(ii,2) = summary_olr_mmw_stats.era5.mmw(2);
  thesummary.era5_std(ii,3) = summary_olr_mmw_stats.era5.mmwpc(2);
  thesummary.merra2_std(ii,1) = summary_olr_mmw_stats.merra2.stemp(2);
  thesummary.merra2_std(ii,2) = summary_olr_mmw_stats.merra2.mmw(2);
  thesummary.merra2_std(ii,3) = summary_olr_mmw_stats.merra2.mmwpc(2); 

  ifoo = 3;
  thesummary.umbc_mean(ii,ifoo+1) = summary_olr_mmw_stats.umbc.tropical_stemp(1);
  thesummary.umbc_mean(ii,ifoo+2) = summary_olr_mmw_stats.umbc.tropical_mmw(1);
  thesummary.umbc_mean(ii,ifoo+3) = summary_olr_mmw_stats.umbc.tropical_mmwpc(1);
  thesummary.era5_mean(ii,ifoo+1) = summary_olr_mmw_stats.era5.tropical_stemp(1);
  thesummary.era5_mean(ii,ifoo+2) = summary_olr_mmw_stats.era5.tropical_mmw(1);
  thesummary.era5_mean(ii,ifoo+3) = summary_olr_mmw_stats.era5.tropical_mmwpc(1);
  thesummary.merra2_mean(ii,ifoo+1) = summary_olr_mmw_stats.merra2.tropical_stemp(1);
  thesummary.merra2_mean(ii,ifoo+2) = summary_olr_mmw_stats.merra2.tropical_mmw(1);
  thesummary.merra2_mean(ii,ifoo+3) = summary_olr_mmw_stats.merra2.tropical_mmwpc(1);

  thesummary.umbc_std(ii,ifoo+1) = summary_olr_mmw_stats.umbc.tropical_stemp(2);
  thesummary.umbc_std(ii,ifoo+2) = summary_olr_mmw_stats.umbc.tropical_mmw(2);
  thesummary.umbc_std(ii,ifoo+3) = summary_olr_mmw_stats.umbc.tropical_mmwpc(2);
  thesummary.era5_std(ii,ifoo+1) = summary_olr_mmw_stats.era5.tropical_stemp(2);
  thesummary.era5_std(ii,ifoo+2) = summary_olr_mmw_stats.era5.tropical_mmw(2);
  thesummary.era5_std(ii,ifoo+3) = summary_olr_mmw_stats.era5.tropical_mmwpc(2);
  thesummary.merra2_std(ii,ifoo+1) = summary_olr_mmw_stats.merra2.tropical_stemp(2);
  thesummary.merra2_std(ii,ifoo+2) = summary_olr_mmw_stats.merra2.tropical_mmw(2);
  thesummary.merra2_std(ii,ifoo+3) = summary_olr_mmw_stats.merra2.tropical_mmwpc(2);

  ifoo = 6;
  thesummary.umbc_mean(ii,ifoo+1) = summary_olr_mmw_stats.umbc.midlat_stemp(1);
  thesummary.umbc_mean(ii,ifoo+2) = summary_olr_mmw_stats.umbc.midlat_mmw(1);
  thesummary.umbc_mean(ii,ifoo+3) = summary_olr_mmw_stats.umbc.midlat_mmwpc(1);
  thesummary.era5_mean(ii,ifoo+1) = summary_olr_mmw_stats.era5.midlat_stemp(1);
  thesummary.era5_mean(ii,ifoo+2) = summary_olr_mmw_stats.era5.midlat_mmw(1);
  thesummary.era5_mean(ii,ifoo+3) = summary_olr_mmw_stats.era5.midlat_mmwpc(1);
  thesummary.merra2_mean(ii,ifoo+1) = summary_olr_mmw_stats.merra2.midlat_stemp(1);
  thesummary.merra2_mean(ii,ifoo+2) = summary_olr_mmw_stats.merra2.midlat_mmw(1);
  thesummary.merra2_mean(ii,ifoo+3) = summary_olr_mmw_stats.merra2.midlat_mmwpc(1);

  thesummary.umbc_std(ii,ifoo+1) = summary_olr_mmw_stats.umbc.midlat_stemp(2);
  thesummary.umbc_std(ii,ifoo+2) = summary_olr_mmw_stats.umbc.midlat_mmw(2);
  thesummary.umbc_std(ii,ifoo+3) = summary_olr_mmw_stats.umbc.midlat_mmwpc(2);
  thesummary.era5_std(ii,ifoo+1) = summary_olr_mmw_stats.era5.midlat_stemp(2);
  thesummary.era5_std(ii,ifoo+2) = summary_olr_mmw_stats.era5.midlat_mmw(2);
  thesummary.era5_std(ii,ifoo+3) = summary_olr_mmw_stats.era5.midlat_mmwpc(2);
  thesummary.merra2_std(ii,ifoo+1) = summary_olr_mmw_stats.merra2.midlat_stemp(2);
  thesummary.merra2_std(ii,ifoo+2) = summary_olr_mmw_stats.merra2.midlat_mmw(2);
  thesummary.merra2_std(ii,ifoo+3) = summary_olr_mmw_stats.merra2.midlat_mmwpc(2);

  ifoo = 9;
  thesummary.umbc_mean(ii,ifoo+1) = summary_olr_mmw_stats.umbc.polar_stemp(1);
  thesummary.umbc_mean(ii,ifoo+2) = summary_olr_mmw_stats.umbc.polar_mmw(1);
  thesummary.umbc_mean(ii,ifoo+3) = summary_olr_mmw_stats.umbc.polar_mmwpc(1);
  thesummary.era5_mean(ii,ifoo+1) = summary_olr_mmw_stats.era5.polar_stemp(1);
  thesummary.era5_mean(ii,ifoo+2) = summary_olr_mmw_stats.era5.polar_mmw(1);
  thesummary.era5_mean(ii,ifoo+3) = summary_olr_mmw_stats.era5.polar_mmwpc(1);
  thesummary.merra2_mean(ii,ifoo+1) = summary_olr_mmw_stats.merra2.polar_stemp(1);
  thesummary.merra2_mean(ii,ifoo+2) = summary_olr_mmw_stats.merra2.polar_mmw(1);
  thesummary.merra2_mean(ii,ifoo+3) = summary_olr_mmw_stats.merra2.polar_mmwpc(1);

  thesummary.umbc_std(ii,ifoo+1) = summary_olr_mmw_stats.umbc.polar_stemp(2);
  thesummary.umbc_std(ii,ifoo+2) = summary_olr_mmw_stats.umbc.polar_mmw(2);
  thesummary.umbc_std(ii,ifoo+3) = summary_olr_mmw_stats.umbc.polar_mmwpc(2);
  thesummary.era5_std(ii,ifoo+1) = summary_olr_mmw_stats.era5.polar_stemp(2);
  thesummary.era5_std(ii,ifoo+2) = summary_olr_mmw_stats.era5.polar_mmw(2);
  thesummary.era5_std(ii,ifoo+3) = summary_olr_mmw_stats.era5.polar_mmwpc(2);
  thesummary.merra2_std(ii,ifoo+1) = summary_olr_mmw_stats.merra2.polar_stemp(2);
  thesummary.merra2_std(ii,ifoo+2) = summary_olr_mmw_stats.merra2.polar_mmw(2);
  thesummary.merra2_std(ii,ifoo+3) = summary_olr_mmw_stats.merra2.polar_mmwpc(2);

  ifoo = 12;
  thesummary.umbc_mean(ii,ifoo+1) = summary_olr_mmw_stats.umbc.ocean_stemp(1);
  thesummary.umbc_mean(ii,ifoo+2) = summary_olr_mmw_stats.umbc.ocean_mmw(1);
  thesummary.umbc_mean(ii,ifoo+3) = summary_olr_mmw_stats.umbc.ocean_mmwpc(1);
  thesummary.era5_mean(ii,ifoo+1) = summary_olr_mmw_stats.era5.ocean_stemp(1);
  thesummary.era5_mean(ii,ifoo+2) = summary_olr_mmw_stats.era5.ocean_mmw(1);
  thesummary.era5_mean(ii,ifoo+3) = summary_olr_mmw_stats.era5.ocean_mmwpc(1);
  thesummary.merra2_mean(ii,ifoo+1) = summary_olr_mmw_stats.merra2.ocean_stemp(1);
  thesummary.merra2_mean(ii,ifoo+2) = summary_olr_mmw_stats.merra2.ocean_mmw(1);
  thesummary.merra2_mean(ii,ifoo+3) = summary_olr_mmw_stats.merra2.ocean_mmwpc(1);

  thesummary.umbc_std(ii,ifoo+1) = summary_olr_mmw_stats.umbc.ocean_stemp(2);
  thesummary.umbc_std(ii,ifoo+2) = summary_olr_mmw_stats.umbc.ocean_mmw(2);
  thesummary.umbc_std(ii,ifoo+3) = summary_olr_mmw_stats.umbc.ocean_mmwpc(2);
  thesummary.era5_std(ii,ifoo+1) = summary_olr_mmw_stats.era5.ocean_stemp(2);
  thesummary.era5_std(ii,ifoo+2) = summary_olr_mmw_stats.era5.ocean_mmw(2);
  thesummary.era5_std(ii,ifoo+3) = summary_olr_mmw_stats.era5.ocean_mmwpc(2);
  thesummary.merra2_std(ii,ifoo+1) = summary_olr_mmw_stats.merra2.ocean_stemp(2);
  thesummary.merra2_std(ii,ifoo+2) = summary_olr_mmw_stats.merra2.ocean_mmw(2);
  thesummary.merra2_std(ii,ifoo+3) = summary_olr_mmw_stats.merra2.ocean_mmwpc(2);

  ifoo = 15;
  thesummary.umbc_mean(ii,ifoo+1) = summary_olr_mmw_stats.umbc.land_stemp(1);
  thesummary.umbc_mean(ii,ifoo+2) = summary_olr_mmw_stats.umbc.land_mmw(1);
  thesummary.umbc_mean(ii,ifoo+3) = summary_olr_mmw_stats.umbc.land_mmwpc(1);
  thesummary.era5_mean(ii,ifoo+1) = summary_olr_mmw_stats.era5.land_stemp(1);
  thesummary.era5_mean(ii,ifoo+2) = summary_olr_mmw_stats.era5.land_mmw(1);
  thesummary.era5_mean(ii,ifoo+3) = summary_olr_mmw_stats.era5.land_mmwpc(1);
  thesummary.merra2_mean(ii,ifoo+1) = summary_olr_mmw_stats.merra2.land_stemp(1);
  thesummary.merra2_mean(ii,ifoo+2) = summary_olr_mmw_stats.merra2.land_mmw(1);
  thesummary.merra2_mean(ii,ifoo+3) = summary_olr_mmw_stats.merra2.land_mmwpc(1);

  thesummary.umbc_std(ii,ifoo+1) = summary_olr_mmw_stats.umbc.land_stemp(2);
  thesummary.umbc_std(ii,ifoo+2) = summary_olr_mmw_stats.umbc.land_mmw(2);
  thesummary.umbc_std(ii,ifoo+3) = summary_olr_mmw_stats.umbc.land_mmwpc(2);
  thesummary.era5_std(ii,ifoo+1) = summary_olr_mmw_stats.era5.land_stemp(2);
  thesummary.era5_std(ii,ifoo+2) = summary_olr_mmw_stats.era5.land_mmw(2);
  thesummary.era5_std(ii,ifoo+3) = summary_olr_mmw_stats.era5.land_mmwpc(2);
  thesummary.merra2_std(ii,ifoo+1) = summary_olr_mmw_stats.merra2.land_stemp(2);
  thesummary.merra2_std(ii,ifoo+2) = summary_olr_mmw_stats.merra2.land_mmw(2);
  thesummary.merra2_std(ii,ifoo+3) = summary_olr_mmw_stats.merra2.land_mmwpc(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

  %% allX_BLAH(iY,iX,:) will have size 9 x 5 x 4 
  %%     == [200/500/800 RH/WV/T] X ['ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK'] x [GLOBAL TRP MIDLAT POLAR]
  %%   first index  (iY)   1 .. 9 is  1:3=200 mb RH/WV/T   4:6=500 mb RH/WV/T   7:9=500 mb RH/WV/T
  %%   second index (iX)   1 .. 5 is  'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK' if iBiasWRT_ERA5orUMBC > 0, if iBiasWRT_ERA5orUMBC < 0 then swap UMBC, ERA5
  %%   third index  (iG)   1 .. 4 is  global,tropical,midlat,polar

  %% iY = 1,2,3  200 mb    4,5,6 500 mb    7,8,9 800 mb             = 9 total   [RH,WVfrac,T][RH,WVfrac,T][RH,WVfrac,T]
  %% iX = 1,2,3,4,5 === all, tropics,midlats,midlats+tropics,poles  = 5 total
  %% whos allXchi = 9 x 5 x 4                                       correlate ERA5 with [AIRS L3, CLIMCAPS, MERRA2, UMBC]
  %%   first index  (iY)   1 .. 9 is  1:3=200 mb RH/WV/T   4:6=500 mb RH/WV/T   7:9=500 mb RH/WV/T
  %%   second index (iX)   1 .. 5 is  'ERA5','AIRS L3','CLIMCAPS L3','MERRA2','THIS WORK' if iBiasWRT_ERA5orUMBC > 0, if iBiasWRT_ERA5orUMBC < 0 then swap UMBC, ERA5
  %%   third index  (iG)   1 .. 4 is  global,tropical,midlat,polar
%}

len = length(strUMBC);
if iMin < len
  if ~exist('era5X')
    [~,~,p,~] = rtpread('summary_atm_N_cld_20years_all_lat_all_lon_2002_2022_monthlyERA5.ip.rtp');
    plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
    plays = flipud(plevs2plays(plevs));
    i800 = find(plays >= 800,1);  i500 = find(plays >= 500,1); i200 = find(plays >= 200,1); 
    commentA = 'era5X : pressures levels at 200, 500, 800 mb save from era5 rates';
    commentB = 'allX_BLAH(iY,iX,:) will have size 9 x 5 x 4 ===== [200/500/800 RH/WV/T] X [''ERA5'',''AIRS L3'',''CLIMCAPS L3'',''MERRA2'',''THIS WORK''] x [GLOBAL TRP MIDLAT POLAR]';
    era5X.rlat = p.rlat; 
    era5X.rhrate(1,:) = era5.RHrate(i200,:); era5X.rhrate(2,:) = era5.RHrate(i500,:);era5X.rhrate(3,:) = era5.RHrate(i800,:);
    era5X.waterrate(1,:) = era5.waterrate(i200,:); era5X.waterrate(2,:) = era5.waterrate(i500,:);era5X.waterrate(3,:) = era5.waterrate(i800,:);
    era5X.ptemprate(1,:) = era5.ptemprate(i200,:); era5X.ptemprate(2,:) = era5.ptemprate(i500,:);era5X.ptemprate(3,:) = era5.ptemprate(i800,:);
  end

  save summary_mmw_stemp_corr_umbc_era5.mat thesummary strUMBC era5X comment*

  for ii = 1 : 3 
    [m,s,m0,s0] = weighted_mean_stddev(era5X.rhrate(ii,:),cos(era5X.rlat*pi/180));    cosine_avg_era5rh(ii) = m;
    [m,s,m0,s0] = weighted_mean_stddev(era5X.waterrate(ii,:),cos(era5X.rlat*pi/180)); cosine_avg_era5water(ii) = m;
    [m,s,m0,s0] = weighted_mean_stddev(era5X.ptemprate(ii,:),cos(era5X.rlat*pi/180)); cosine_avg_era5ptemp(ii) = m;
  end

  figure(32); figure(31); figure(33); figure(35); figure(37); figure(39)

else
  fprintf(1,'already have saved off all the important data for the %2i different mat files .. will plot summary\n',len)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iBest = 10; tiled_loop_driver_summary_plots
iBest = 10; tiled_loop_driver_summary_plots_errorbars

