addpath /home/sergio/MATLABCODE/COLORMAP/LLS

load llsmap5

%% compare effects of chaging land emiss
%% get_umbc_day_night_name.m
umbc_file1N = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjacV2.mat';
umbc_file2N = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_removeemisstrends.mat';

umbc_file1D = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_dayV2.mat';
umbc_file2D = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2F_const_dlnPdST_clrjac_day_removeemisstrends.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see prep_colWV_T_WV_trends_Day_vs_Night.m
era5_day_file    = 'ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc.mat';
era5_night_file  = 'ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc.mat';

airsL3_day_file   = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_asc.mat';
airsL3_night_file = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat';

climcapsL3_day_file   = '/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_asc.mat';
climcapsL3_night_file = '/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat';

merra2_file = 'MERRA2_atm_data_2002_09_to_2022_08_trends_desc.mat';
giss_file   = '/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ChrisHTrends/giss_trends_2002_2022.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oldN = load(umbc_file1N,'results');  oldD = load(umbc_file1D,'results');
newN = load(umbc_file2N,'results');  newD = load(umbc_file2D,'results');

era5N = load(era5_night_file);          era5D = load(era5_day_file);
climL3N = load(climcapsL3_night_file);  climL3D = load(climcapsL3_day_file);
airsL3N = load(airsL3_night_file);      airsL3D = load(airsL3_day_file);

giss = load('ChrisHTrends/giss_trends_2002_2022.mat');
merra2 = load(merra2_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,ha,p,pa] = rtpread('summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.op.rtp');
[salti,landfrac] =  usgs_deg10_dem(p.rlat,p.rlon);
p.landfrac = landfrac;
do_XX_YY_from_X_Y
land = find(p.landfrac == 1);
ocean = find(p.landfrac == 0);

iNorD = input('Enter (+1/default) for night (-1) for day (-0) for D/N -9999 to stop : ');
if length(iNorD) == 0
  iNorD = 1;
end
while iNorD > -10
  if iNorD == 0
    fD = 1; fN = 1;
  elseif iNorD == +1
    fD = 0; fN = 2;
  elseif iNorD == -1
    fD = 2; fN = 0;
  end
  
  old.results = (fN*oldN.results + fD*oldD.results)/2;
  new.results = (fN*newN.results + fD*newD.results)/2;
  era5.trend_stemp = (fN*era5N.trend_stemp + fD*era5D.trend_stemp)/2;
  airsL3.thestats64x72.stemprate = (fN*airsL3N.thestats64x72.stemprate + fD*airsL3D.thestats64x72.stemprate)/2;
  climL3.thestats64x72.stemprate = (fN*climL3N.thestats64x72.stemprate + fD*climL3D.thestats64x72.stemprate)/2;  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ii = 0; 
  ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(old.results(:,6),72,64)',1),                 [-90 +90],[-180 +180]); title('dSKT/dt : d/dt emiss = 0');
  ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(new.results(:,6),72,64)',1),                 [-90 +90],[-180 +180]); title('dSKT/dt : d/dt emiss = nonzero')
  ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape(old.results(:,6)-new.results(:,6),72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : 0 - nonzero');
  ii = ii + 1; figure(ii); clf; aslmap(ii,rlat65,rlon73,smoothn(reshape((old.results(:,6)-new.results(:,6))./(old.results(:,6)+eps),72,64)',1),[-90 +90],[-180 +180]); title('dSKT/dt : fraction change 0-nonzero');
  
  figure(1); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(2); caxis([-1 +1]*0.15); colormap(llsmap5)
  figure(3); caxis([-1 +1]*0.15/3); colormap(llsmap5)
  figure(4); caxis([-1 +1]); colormap(llsmap5)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  figure(5); clf
  old.landresults = old.results(:,6);        old.landresults(ocean) = NaN;
  new.landresults = new.results(:,6);        new.landresults(ocean) = NaN;
  era5.landresults = era5.trend_stemp;       era5.landresults(ocean) = NaN;
  merra2.landresults = merra2.trend_stemp;   merra2.landresults(ocean) = NaN;
  airsL3.landresults = airsL3.thestats64x72.stemprate(:); airsL3.landresults(ocean) = NaN;
  climL3.landresults = climL3.thestats64x72.stemprate(:); climL3.landresults(ocean) = NaN;
  giss.landresults = giss.giss_trend4608(:); giss.landresults(ocean) = NaN;
  plot(rlat,nanmean(reshape(old.landresults,72,64),1),'g',rlat,nanmean(reshape(new.landresults,72,64),1),'y',...
       rlat,nanmean(reshape(giss.landresults,72,64),1),'k',....
       rlat,nanmean(reshape(airsL3.landresults,72,64),1),'b',rlat,nanmean(reshape(climL3.landresults,72,64),1),'c',...
       rlat,nanmean(reshape(era5.landresults,72,64),1),'g',rlat,nanmean(reshape(merra2.landresults,72,64),1),'m',...
       'linewidth',2)
  plotaxis2;  hl = legend('constant emiss','changing emiss','GISS','AIRSL3','CLIMCAPSL3','ERA5','MERRA2','location','best','fontsize',8); 
  ylabel('dSKT/dt K/yr'); xlabel('Latitude'); title('dSKT/dt land')
  
  figure(6); clf
  old.oceanresults = old.results(:,6);        old.oceanresults(land) = NaN;
  new.oceanresults = new.results(:,6);        new.oceanresults(land) = NaN;
  era5.oceanresults = era5.trend_stemp;       era5.oceanresults(land) = NaN;
  merra2.oceanresults = merra2.trend_stemp;   merra2.oceanresults(land) = NaN;
  airsL3.oceanresults = airsL3.thestats64x72.stemprate(:); airsL3.oceanresults(land) = NaN;
  climL3.oceanresults = climL3.thestats64x72.stemprate(:); climL3.oceanresults(land) = NaN;
  giss.oceanresults = giss.giss_trend4608(:); giss.oceanresults(land) = NaN;
  plot(rlat,nanmean(reshape(old.oceanresults,72,64),1),'g',rlat,nanmean(reshape(new.oceanresults,72,64),1),'y',...
       rlat,nanmean(reshape(giss.oceanresults,72,64),1),'k',....
       rlat,nanmean(reshape(airsL3.oceanresults,72,64),1),'b',rlat,nanmean(reshape(climL3.oceanresults,72,64),1),'c',...
       rlat,nanmean(reshape(era5.oceanresults,72,64),1),'g',rlat,nanmean(reshape(merra2.oceanresults,72,64),1),'m',...
       'linewidth',2)
  plotaxis2;  hl = legend('constant emiss','changing emiss','GISS','AIRSL3','CLIMCAPSL3','ERA5','MERRA2','location','best','fontsize',8); 
  ylabel('dSKT/dt K/yr'); xlabel('Latitude'); title('dSKT/dt ocean')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(5); clf
  old.oceanresults = old.results(:,6);        old.oceanresults(land) = NaN;
  new.oceanresults = new.results(:,6);        new.oceanresults(land) = NaN;
  era5.oceanresults = era5.trend_stemp;       era5.oceanresults(land) = NaN;
  merra2.oceanresults = merra2.trend_stemp;   merra2.oceanresults(land) = NaN;
  giss.oceanresults = giss.giss_trend4608(:); giss.oceanresults(land) = NaN;
  airsL3.oceanresults = airsL3.thestats64x72.stemprate(:); airsL3.oceanresults(land) = NaN;
  climL3.oceanresults = climL3.thestats64x72.stemprate(:); climL3.oceanresults(land) = NaN;
  plot(rlat,nanmean(reshape(old.oceanresults,72,64),1),'b',rlat,nanmean(reshape(new.oceanresults,72,64),1),'r',...
       rlat,nanmean(reshape(era5.oceanresults,72,64),1),'g',rlat,nanmean(reshape(giss.oceanresults,72,64),1),'k','linewidth',2)
  plotaxis2;  hl = legend('constant emiss','changing emiss','ERA5','GISS','location','best','fontsize',10); 
  ylabel('dSKT/dt K/yr'); xlabel('Latitude'); title('dSKT/dt ocean')
  
  figure(6); clf
  old.landresults = old.results(:,6);        old.landresults(ocean) = NaN;
  new.landresults = new.results(:,6);        new.landresults(ocean) = NaN;
  era5.landresults = era5.trend_stemp;       era5.landresults(ocean) = NaN;
  merra2.landresults = merra2.trend_stemp;   merra2.landresults(ocean) = NaN;
  giss.landresults = giss.giss_trend4608(:); giss.landresults(ocean) = NaN;
  airsL3.landresults = airsL3.thestats64x72.stemprate(:); airsL3.landresults(ocean) = NaN;
  climL3.landresults = climL3.thestats64x72.stemprate(:); climL3.landresults(ocean) = NaN;
  plot(rlat,nanmean(reshape(old.landresults,72,64),1),'b',rlat,nanmean(reshape(new.landresults,72,64),1),'r',...
       rlat,nanmean(reshape(era5.landresults,72,64),1),'g',rlat,nanmean(reshape(giss.landresults,72,64),1),'k','linewidth',2)
  plotaxis2;  hl = legend('constant emiss','changing emiss','ERA5','GISS','location','best','fontsize',10); 
  ylabel('dSKT/dt K/yr'); xlabel('Latitude'); title('dSKT/dt land')
   
  iNorD = input('Enter (+1/default) for night (-1) for day (-0) for D/N -9999 to stop : ');
  if length(iNorD) == 0
    iNorD = 1;
  end
end

