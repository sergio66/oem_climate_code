%% see How closely do changes in surface and column water vapor follow Clausius-Clapeyron scaling in climate-change simulations?
%% P A O’Gorman, C J Muller, https://core.ac.uk/download/pdf/4426849.pdf, or saved in PDF/change_of_RH_with_stemp_GCM_PGorman.pdf
%% P A O'Gorman and C J Muller 2010 Environ. Res. Lett. 5 025207 DOI 10.1088/1748-9326/5/2/025207

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools

load llsmap5

plays100 = load('ERA5_atm_data_2002_09_to_2022_08_trends_desc.mat','trend_plays');
plays100 = plays100.trend_plays;

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2;
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

[h,ha,p,pa] = rtpread('summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.op.rtp');
mmw0 = mmwater_rtp(h,p);

for ii = 1 : 8;
  figure(ii); clf
end

%iX = input('fastgrib (+1/default) or slow calc (-1) : ');
%if length(iX) == 0
%  iX = 1;
%end
iX = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iX > 0
  load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat
else
  load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug2022_20yr_desc.mat
end

%{
figure(1); 
pcolor(thestats64x72.lats,Tlevs,squeeze(nanmean(thestats64x72.ptemprate,1))'); colormap(llsmap5); shading interp; colorbar;
set(gca,'ydir','reverse'); caxis([-1 +1]*0.15); title('AIRS L3 20 years dT/dt'); set(gca,'yscale','log'); ylim([10 1000])

figure(2)
pcolor(thestats64x72.lats,Qlevs,squeeze(nanmean(thestats64x72.waterrate,1))'); colormap(llsmap5); shading interp; colorbar; 
set(gca,'ydir','reverse'); caxis([-1 +1]*0.015); title('AIRS L3 20 years dWVfrac/dt'); ylim([100 1000])

figure(3)
pcolor(thestats64x72.lats,Qlevs,squeeze(nanmean(thestats64x72.RHrate,1))'); colormap(llsmap5); shading interp; colorbar; 
set(gca,'ydir','reverse'); caxis([-1 +1]*0.5); title('AIRS L3 20 years dRH/dt'); ylim([100 1000])
%}

z11 = thestats64x72.stemprate(:);
z11 = z11';

for ii = 1 : 72
  for jj = 1 : 64
    junk = squeeze(thestats64x72.RHrate(ii,jj,:));
    airsL3.RHrate(ii,jj,:) = interp1(log10(Qlevs),junk,log10(plays100),[],'extrap');
    junk = squeeze(thestats64x72.waterrate(ii,jj,:));
    airsL3.waterrate(ii,jj,:) = interp1(log10(Qlevs),junk,log10(plays100),[],'extrap');
    junk = squeeze(thestats64x72.ptemprate(ii,jj,:));
    airsL3.ptemprate(ii,jj,:) = interp1(log10(Tlevs),junk,log10(plays100),[],'extrap');
  end
end

%{
booAIRSL3 = thestats64x72.trend_mmw;
booAIRSL3 = booAIRSL3./mmw0 * 100;
zonal_st = nanmean(reshape(trend_stemp,72,64),1);
booAIRSL3   = nanmean(reshape(booAIRSL3,72,64),1);
booAIRSL3 = booAIRSL3 ./ zonal_st;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load /asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat

%{
figure(4); 
pcolor(thestats64x72.lats,Tlevs/100,squeeze(nanmean(thestats64x72.ptemprate,1))'); colormap(llsmap5); shading interp; colorbar;
set(gca,'ydir','reverse'); caxis([-1 +1]*0.15); title('CLIMCAPS 20 years dT/dt'); set(gca,'yscale','log'); ylim([10 1000])

figure(5)
pcolor(thestats64x72.lats,Qlevs/100,squeeze(nanmean(thestats64x72.waterrate,1))'); colormap(llsmap5); shading interp; colorbar; 
set(gca,'ydir','reverse'); caxis([-1 +1]*0.015); title('CLIMCAPS 20 years dWVfrac/dt'); ylim([100 1000])

figure(6)
pcolor(thestats64x72.lats,Qlevs/100,100*squeeze(nanmean(thestats64x72.RHrate,1))'); colormap(llsmap5); shading interp; colorbar; 
set(gca,'ydir','reverse'); caxis([-1 +1]*0.5); title('CLIMCAPS 20 years dRH/dt'); ylim([100 1000])
%}

z12 = thestats64x72.stemprate(:);
z12 = z12';

for ii = 1 : 72
  for jj = 1 : 64
    junk = squeeze(thestats64x72.RHrate(ii,jj,:));
    climcapsL3.RHrate(ii,jj,:) = interp1(log10(Qlevs/100),100*junk,log10(plays100),[],'extrap');
    junk = squeeze(thestats64x72.waterrate(ii,jj,:));
    climcapsL3.waterrate(ii,jj,:) = interp1(log10(Qlevs/100),junk,log10(plays100),[],'extrap');
    junk = squeeze(thestats64x72.ptemprate(ii,jj,:));
    climcapsL3.ptemprate(ii,jj,:) = interp1(log10(Tlevs/100),junk,log10(plays100),[],'extrap');
  end
end

%{
booCLIMCAPSL3 = thestats64x72.trend_mmw;
booCLIMCAPSL3 = booCLIMCAPSL3./mmw0 * 100;
zonal_st = nanmean(reshape(trend_stemp,72,64),1);
booCLIMCAPSL3   = nanmean(reshape(booCLIMCAPSL3,72,64),1);
booCLIMCAPSL3 = booCLIMCAPSL3 ./ zonal_st;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load MERRA2_atm_data_2002_09_to_2022_08_trends_desc.mat

%{
figure(7); 
pcolor(trend_rlat,trend_plays,trend_ptemp); colormap(llsmap5); shading interp; colorbar;
pcolor(trend_rlat64,trend_plays,squeeze(nanmean(reshape(trend_ptemp,100,72,64),2))); colormap(llsmap5); shading interp; colorbar;
set(gca,'ydir','reverse'); caxis([-1 +1]*0.15); title('MERRA2 20 years dT/dt'); set(gca,'yscale','log'); ylim([10 1000])

figure(8)
pcolor(trend_rlat,trend_plays,trend_gas_1); colormap(llsmap5); shading interp; colorbar; 
pcolor(trend_rlat64,trend_plays,squeeze(nanmean(reshape(trend_gas_1,100,72,64),2))); colormap(llsmap5); shading interp; colorbar;
set(gca,'ydir','reverse'); caxis([-1 +1]*0.015); title('MERRA2 20 years dWVfrac/dt'); ylim([100 1000])

figure(9)
pcolor(trend_rlat,trend_plays,trend_RH); colormap(llsmap5); shading interp; colorbar; 
pcolor(trend_rlat64,trend_plays,squeeze(nanmean(reshape(trend_RH,100,72,64),2))); colormap(llsmap5); shading interp; colorbar;
set(gca,'ydir','reverse'); caxis([-1 +1]*0.5); title('MERRA2 20 years dRH/dt'); ylim([100 1000])
%}

z21 = trend_stemp;
merra2.RHrate = trend_RH;
merra2.waterrate = trend_gas_1;
merra2.ptemprate = trend_ptemp;

iFig = 14;
iFig = 1;
figure(iFig); clf
[h,ha,p,pa] = rtpread('summary_19years_all_lat_all_lon_2002_2021_monthlyMERRA2.op.rtp');
mmw0 = mmwater_rtp(h,p);
booMERRA2 = trend_mmw;
booMERRA2 = booMERRA2./mmw0 * 100;
zonal_st = nanmean(reshape(trend_stemp,72,64),1);
booMERRA2   = nanmean(reshape(booMERRA2,72,64),1);
booMERRA2 = booMERRA2 ./ zonal_st;
plot(trend_rlat64,booMERRA2,trend_rlat64,smooth(booMERRA2,10),'linewidth',2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load ERA5_atm_data_2002_09_to_2022_08_trends_desc.mat

%{
figure(10); 
pcolor(trend_rlat,trend_plays,trend_ptemp); colormap(llsmap5); shading interp; colorbar;
pcolor(trend_rlat64,trend_plays,squeeze(nanmean(reshape(trend_ptemp,100,72,64),2))); colormap(llsmap5); shading interp; colorbar;
set(gca,'ydir','reverse'); caxis([-1 +1]*0.15); title('ERA5 20 years dT/dt'); set(gca,'yscale','log'); ylim([10 1000])

figure(11)
pcolor(trend_rlat,trend_plays,trend_gas_1); colormap(llsmap5); shading interp; colorbar; 
pcolor(trend_rlat64,trend_plays,squeeze(nanmean(reshape(trend_gas_1,100,72,64),2))); colormap(llsmap5); shading interp; colorbar;
set(gca,'ydir','reverse'); caxis([-1 +1]*0.015); title('ERA5 20 years dWVfrac/dt'); ylim([100 1000])

figure(12)
pcolor(trend_rlat,trend_plays,trend_RH); colormap(llsmap5); shading interp; colorbar; 
pcolor(trend_rlat64,trend_plays,squeeze(nanmean(reshape(trend_RH,100,72,64),2))); colormap(llsmap5); shading interp; colorbar;
set(gca,'ydir','reverse'); caxis([-1 +1]*0.5); title('ERA5 20 years dRH/dt'); ylim([100 1000])
%}

z22 = trend_stemp;
era5.RHrate = trend_RH;
era5.waterrate = trend_gas_1;
era5.ptemprate = trend_ptemp;

figure(iFig); clf
[h,ha,p,pa] = rtpread('summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.op.rtp');
mmw0 = mmwater_rtp(h,p);
booERA5 = trend_mmw;
booERA5 = booERA5./mmw0 * 100;
zonal_st = nanmean(reshape(trend_stemp,72,64),1);
booERA5   = nanmean(reshape(booERA5,72,64),1);
booERA5 = booERA5 ./ zonal_st;
plot(trend_rlat64,booERA5,trend_rlat64,smooth(booERA5,10),'linewidth',2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(iFig); clf
plot(trend_rlat64,booMERRA2,'b',trend_rlat64,smooth(booMERRA2,10),'c',trend_rlat64,booERA5,'r',trend_rlat64,smooth(booERA5,10),'m','linewidth',2); 
plot(trend_rlat64,smooth(booMERRA2,10),'b',trend_rlat64,smooth(booERA5,10),'r','linewidth',2); 
plotaxis2; hl = legend('MERRA2','ERA5','location','best'); ylabel('d%mmw/dST'); xlabel('Latitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = iFig + 1;
figure(iFig); clf; clear plotoptions
maskLF = ones(size(z22));
plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'AIRS L3';
plotoptions.str12 = 'CLIMCAPS L3';
plotoptions.str21 = 'MERRA2';
plotoptions.str22 = 'ERA5';
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = +1;
aslmap_4tiledlayout(z11,z12,z21,z22,iFig,plotoptions);

clear plotoptions
plotoptions.plotcolors = llsmap5;
plotoptions.str11 = 'AIRS L3';
plotoptions.str12 = 'CLIMCAPS L3';
plotoptions.str21 = 'ERA5';
plotoptions.str22 = 'MERRA2';
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
plotoptions.yLinearOrLog = +1;
plotoptions.yReverseDir  = +1;
plotoptions.yLimits = [100 1000];

iFig = iFig + 1;
figure(iFig); clf; 
plotoptions.maintitle = 'dRH/dt'; 
plotoptions.cx = [-1 +1]*0.5; 
miaow11 = squeeze(nanmean(permute(airsL3.RHrate,[3 1 2]),2));
miaow12 = squeeze(nanmean(permute(climcapsL3.RHrate,[3 1 2]),2));
miaow21 = squeeze(nanmean(reshape(era5.RHrate,100,72,64),2));
miaow22 = squeeze(nanmean(reshape(merra2.RHrate,100,72,64),2));
profile_plots_4tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

iFig = iFig + 1;
figure(iFig); clf; 
plotoptions.maintitle = 'dWVfrac/dt'; 
plotoptions.cx = [-1 +1]*0.015; 
miaow11 = squeeze(nanmean(permute(airsL3.waterrate,[3 1 2]),2));
miaow12 = squeeze(nanmean(permute(climcapsL3.waterrate,[3 1 2]),2));
miaow21 = squeeze(nanmean(reshape(era5.waterrate,100,72,64),2));
miaow22 = squeeze(nanmean(reshape(merra2.waterrate,100,72,64),2));
profile_plots_4tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

iFig = iFig + 1;
figure(iFig); clf; 
plotoptions.maintitle = 'dT/dt'; 
plotoptions.cx = [-1 +1]*0.15; 
plotoptions.yLinearOrLog = -1;
plotoptions.yReverseDir  = +1;
plotoptions.yLimits = [10 1000];
miaow11 = squeeze(nanmean(permute(airsL3.ptemprate,[3 1 2]),2));
miaow12 = squeeze(nanmean(permute(climcapsL3.ptemprate,[3 1 2]),2));
miaow21 = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
miaow22 = squeeze(nanmean(reshape(merra2.ptemprate,100,72,64),2));
profile_plots_4tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iUMBC = input('Load UMBC (-1 default no /  +1 yes) : ');
if length(iUMBC) == 0
  iUMBC = -1;
end
if iUMBC  > 0
  strUMBC = input('Enter UMBC results fname eg ''/asl/s1/sergio/JUNK/test5_march9_2023.mat'' : ');
  umbc = load(strUMBC,'deltaRH','deltaT','fracWV','results');

  if isfield(umbc,'fracWV')
    umbc.fracWV = umbc.fracWV(1:100,:);
    umbc.deltaT = umbc.deltaT(1:100,:);
  else
    umbc = load(strUMBC,'resultsT','resultsWV','pjunk20','results','deltaRH');
    for ii = 1 : 4608
      junk = umbc.resultsWV(ii,:);
      umbc.fracWV(ii,:) = interp1(log10(umbc.pjunk20),junk,log10(plays100),[],'extrap')';
      junk = umbc.resultsT(ii,:);
      umbc.deltaT(ii,:) = interp1(log10(umbc.pjunk20),junk,log10(plays100),[],'extrap')';
    end
    umbc.fracWV = umbc.fracWV';
    umbc.deltaT = umbc.deltaT';
  end

  umbc.stemp = umbc.results(:,6)';

  iFig = iFig + 1;
  figure(iFig); clf; clear plotoptions
  maskLF = ones(size(z22));
  plotoptions.cx = [-1 +1]*0.15; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'UMBC';
  plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'MERRA2';
  plotoptions.str22 = 'ERA5';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  z11x = umbc.stemp;
  aslmap_4tiledlayout(z11x,z12,z21,z22,iFig,plotoptions);
  
  clear plotoptions
  plotoptions.plotcolors = llsmap5;
  plotoptions.str11 = 'UMBC';
  plotoptions.str12 = 'CLIMCAPS L3';
  plotoptions.str21 = 'ERA5';
  plotoptions.str22 = 'MERRA2';
  plotoptions.xstr = ' ';        plotoptions.ystr = ' ';
  plotoptions.yLinearOrLog = +1;
  plotoptions.yReverseDir  = +1;
  
  iFig = iFig + 1;
  figure(iFig); clf; 
  plotoptions.maintitle = 'dRH/dt'; 
  plotoptions.cx = [-1 +1]*0.5; 
  plotoptions.yLimits = [100 1000];
  miaow11 = squeeze(nanmean(reshape(umbc.deltaRH,100,72,64),2));
  miaow12 = squeeze(nanmean(permute(climcapsL3.RHrate,[3 1 2]),2));
  miaow21 = squeeze(nanmean(reshape(era5.RHrate,100,72,64),2));
  miaow22 = squeeze(nanmean(reshape(merra2.RHrate,100,72,64),2));
  profile_plots_4tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);
  
  iFig = iFig + 1;
  figure(iFig); clf; 
  plotoptions.maintitle = 'dWVfrac/dt'; 
  plotoptions.cx = [-1 +1]*0.015; 
  plotoptions.yLimits = [100 1000];
  miaow11 = squeeze(nanmean(reshape(umbc.fracWV,100,72,64),2));
  miaow12 = squeeze(nanmean(permute(climcapsL3.waterrate,[3 1 2]),2));
  miaow21 = squeeze(nanmean(reshape(era5.waterrate,100,72,64),2));
  miaow22 = squeeze(nanmean(reshape(merra2.waterrate,100,72,64),2));
  profile_plots_4tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);
  
  iFig = iFig + 1;
  figure(iFig); clf; 
  plotoptions.maintitle = 'dT/dt'; 
  plotoptions.cx = [-1 +1]*0.15; 
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits = [10 1000];
  miaow11 = squeeze(nanmean(reshape(umbc.deltaT,100,72,64),2));
  miaow12 = squeeze(nanmean(permute(climcapsL3.ptemprate,[3 1 2]),2));
  miaow21 = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
  miaow22 = squeeze(nanmean(reshape(merra2.ptemprate,100,72,64),2));
  profile_plots_4tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);

  iFig = iFig + 1;
  figure(iFig); clf; 
  plotoptions.maintitle = 'dT/dt'; 
  plotoptions.cx = [-1 +1]*0.15; 
  plotoptions.yLinearOrLog = -1;
  plotoptions.yReverseDir  = +1;
  plotoptions.yLimits = [0.1 100];
  miaow11 = squeeze(nanmean(reshape(umbc.deltaT,100,72,64),2));
  miaow12 = squeeze(nanmean(permute(climcapsL3.ptemprate,[3 1 2]),2));
  miaow21 = squeeze(nanmean(reshape(era5.ptemprate,100,72,64),2));
  miaow22 = squeeze(nanmean(reshape(merra2.ptemprate,100,72,64),2));
  profile_plots_4tiledlayout(trend_rlat64,plays100,miaow11,miaow12,miaow21,miaow22,iFig,plotoptions);
end

