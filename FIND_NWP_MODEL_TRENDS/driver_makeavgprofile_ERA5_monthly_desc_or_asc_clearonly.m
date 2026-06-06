% addpath /asl/matlib/h4tools
% addpath /home/sergio/MATLABCODE/TIME
% addpath /home/sergio/MATLABCODE/PLOTTER
% addpath /asl/matlib/aslutil/

addpath /home/sergio/git/rtpmake/CLUST_RTPMAKE/COMMON_SETTINGS
addpath0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see eg ~/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_AIRS_STM_Oct2020_allstarts_Jan20XY/call_save_split_apart_rtp_howard_bins.m
%% see eg ~/git/oem_climate_jacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_AIRS_STM_Oct2020_allstarts_Jan20XY/call_save_split_apart_rtp_howard_bins.m
fRTP = '/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/summary_17years_all_lat_all_lon_2002_2019.rtp';
if exist(fRTP)
  [hERAI,ha,pERAI,pa] = rtpread(fRTP)
else
  if ~exist('better_angles')
    better_surface = load('/home/sergio/git/oem_climate_jacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/tile_avg_salti_lf_emis.mat');  
    better_angles = load('/home/sergio/git/oem_climate_jacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/asc_desc_solzen_time_525_64x72.mat');
  end
end  

%topts.klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
%topts.sarta = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';
set_path_to_execs
topts.klayers = klayers;
topts.sarta   = sartaCld;

iStartY = 2002; iStartYM = 09; iStopY = 2021; iStopYM = 08; 
iStartY = 2002; iStartYM = 09; iStopY = 2014; iStopYM = 08; %% Joao wants AMIP/CMIP
iStartY = 2012; iStartYM = 05; iStopY = 2019; iStopYM = 04; %% CRIS NSR
iStartY = 2002; iStartYM = 09; iStopY = 2022; iStopYM = 08;
iStartY = 2002; iStartYM = 09; iStopY = 2025; iStopYM = 08; %% AIRS

%% likely made by /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_computeERA5_monthly_trends_desc_or_asc.m
if ~exist('all') & iStopY == 2021
  load ERA5_atm_data_2002_09_to_2021_08_desc.mat
elseif ~exist('all') & iStopY == 2022
  load ERA5_atm_data_2002_09_to_2022_08_desc.mat
elseif ~exist('pall') & iStopY == 2025
  load MEAN_PROFILES/ERA5_atm_N_cld_data_2002_09_to_2025_08_desc.mat
end
%% likely made by /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_computeERA5_monthly_trends_desc_or_asc.m

%% before it ws "all" not it is "pall"
booS = find(pall.yy == iStartY & pall.mm == iStartYM);
booE = find(pall.yy == iStopY & pall.mm == iStopYM);
ind = booS : booE;

if round(length(ind)/12) ~= length(ind)/12
  error('need full years')
end

%% these are actually DESC
if length(ind)/12 == 19 & iStartY == 2002 & iStopY == 2021
  foutNyearaverageIP = ['summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.ip.rtp'];
  foutNyearaverageOP = ['summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.op.rtp'];
  foutNyearaverageRP = ['summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp'];
else
  %% generic 2002/09 - 20XY/08
  foutNyearaverageIP = ['summary_' num2str(length(ind)/12,'%02d') 'years_all_lat_all_lon_' num2str(iStartY) '_' num2str(iStopY) '_monthlyERA5.ip.rtp'];
  foutNyearaverageOP = ['summary_' num2str(length(ind)/12,'%02d') 'years_all_lat_all_lon_' num2str(iStartY) '_' num2str(iStopY) '_monthlyERA5.op.rtp'];
  foutNyearaverageRP = ['summary_' num2str(length(ind)/12,'%02d') 'years_all_lat_all_lon_' num2str(iStartY) '_' num2str(iStopY) '_monthlyERA5.rp.rtp'];
end

if ~exist(foutNyearaverageIP)
  pNyearaverage = [];
  for ii = 1 : 4608
    yavg.stemp(ii) = nanmean(pall.stemp(ind,ii));
    yavg.ptemp(:,ii) = nanmean(squeeze(pall.nwp_ptemp(ind,:,ii)),1);
    yavg.gas_1(:,ii) = nanmean(squeeze(pall.nwp_gas_1(ind,:,ii)),1);
    yavg.gas_3(:,ii) = nanmean(squeeze(pall.nwp_gas_3(ind,:,ii)),1);
    yavg.plevs(:,ii) = nanmean(squeeze(pall.nwp_plevs(ind,:,ii)),1);
    yavg.rlon(ii)  = pall.rlon(ii);
    yavg.rlat(ii)  = pall.rlat(ii);
  end

  if exist('pERAI')
    error('why use 17 year data?????')
    yavg.spres   = pERAI.spres;
    yavg.solzen  = pERAI.solzen;
    yavg.satzen  = pERAI.satzen;
    yavg.rtime   = pERAI.rtime;
    yavg.zobs    = pERAI.zobs;
    yavg.salti   = pERAI.salti;
    yavg.scanang = pERAI.scanang;
    yavg.nemis   = pERAI.nemis;
    yavg.emis    = pERAI.emis;
    yavg.efreq   = pERAI.efreq;
    yavg.rho     = pERAI.rho;
  else
    [mmjunk,nnjunk,oojunk] = size(better_angles.thedata.solzen_desc);
    yavg.spres    = reshape(nanmean(pall.spres,1),1,4608);
    yavg.zobs     = ones(1,4608) * 705000;    
    yavg.solzen   = reshape(nanmean(reshape(better_angles.thedata.solzen_desc,4608,oojunk),2),1,4608);
    yavg.satzen   = reshape(nanmean(reshape(better_angles.thedata.satzen_desc,4608,oojunk),2),1,4608);
    yavg.rtime    = reshape(nanmean(reshape(better_angles.thedata.rtime_desc,4608,oojunk),2),1,4608);
    yavg.scanang  = saconv(yavg.satzen,yavg.zobs);
    yavg.salti    = better_surface.salti;
    yavg.landfrac = better_surface.landfrac;
    yavg.nemis    = better_surface.nemis;
    yavg.emis     = better_surface.emis;
    yavg.efreq    = better_surface.efreq;
    yavg.rho      = better_surface.rho;
  end
  
  tS = utc2taiSergio(2002,09,01,00);
  tE = utc2taiSergio(2021,08,31,23.99);
  
  [yy,mm,dd,hh] = tai2utcSergio(yavg.rtime(ind));
  tS = utc2taiSergio(iStartY,iStartYM,01,00);
  tE = utc2taiSergio(iStopY, iStopYM ,31,23.99);

  [yyX,mmX,ddX,hhX] = tai2utcSergio((tS+tE)/2);
  yavg.rtime = ones(size(yavg.rtime)) * (tS+tE)/2;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  yavg.co2ppm = yyX + (mmX-1)/12; yavg.co2ppm = ones(size(yavg.rtime)) * (yavg.co2ppm-2002)*2.2 + 370;
  
    addpath /home/sergio/MATLABCODE/ESRL_TRACE_GAS
    disp('putting in CO2,N2O.CH4 ESRL')
    co2ppm = read_trace_gas(double(yavg.rlat),double(yavg.rtime),2);
    ch4ppm = read_trace_gas(double(yavg.rlat),double(yavg.rtime),6)/1000;
  
    n2oslope = (332-315)/(2020-2000); %% https://www.esrl.noaa.gov/gmd/hats/combined/N2O.html
    [xyy,xmm,xdd,xhh] = tai2utcSergio(yavg.rtime);        %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
    time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
    n2oppm = 315 + n2oslope*time_so_far; %% ppb
    n2oppm = n2oppm/1000;             %% change to ppm
  
    yavg.co2ppm = co2ppm;
    yavg.n2oppm = n2oppm;
    yavg.ch4ppm = ch4ppm;
  
  scatter_coast(yavg.rlon,yavg.rlat,50,yavg.co2ppm)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if exist('hERAI')
    havg = hERAI;
  else
    havg = load('/home/sergio/git/umbc_singlefootprint/h2645structure.mat');
    havg = havg.h;
    ha = [];
    pa = [];
  end    
  havg.ngas = 2;
  havg.glist = [1  3]';
  havg.gunit = [21 21]';
  havg.ptype = 0;
  havg.pfields = 1;
  yavg.plat = yavg.rlat;
  yavg.plon = yavg.rlon;
  yavg.nlevs = 37 * ones(size(yavg.stemp));
  rtpwrite(foutNyearaverageIP,havg,ha,yavg,pa);
  klayerser = ['!' topts.klayers ' fin=' foutNyearaverageIP ' fout=' foutNyearaverageOP]; eval(klayerser);

  [hnew,hax,pnew,pax] = rtpread(foutNyearaverageOP);
  addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
  ppmv2 = layers2ppmv(hnew,pnew,1:length(pnew.stemp),2);
  ppmv4 = layers2ppmv(hnew,pnew,1:length(pnew.stemp),4);
  ppmv6 = layers2ppmv(hnew,pnew,1:length(pnew.stemp),6);

  i500mb = find(pnew.plevs(:,3000) >= 500,1)
  i850mb = find(pnew.plevs(:,3000) >= 850,1)

  scatter_coast(yavg.rlon,yavg.rlat,50,ppmv2(i500mb,:))
  scatter_coast(yavg.rlon,yavg.rlat,50,ppmv2(i500mb,:)-co2ppm)

  pnew.gas_6 = pnew.gas_6 .* (ones(101,1) * ch4ppm./ppmv6(i500mb,:));
  pnew.gas_4 = pnew.gas_4 .* (ones(101,1) * n2oppm./ppmv4(i500mb,:));
  ppmv4 = layers2ppmv(hnew,pnew,1:length(pnew.stemp),4);
  ppmv6 = layers2ppmv(hnew,pnew,1:length(pnew.stemp),6);
  plot(ppmv6(i500mb,:),ch4ppm)
  plot(ppmv4(i500mb,:),n2oppm)

  scatter_coast(yavg.rlon,yavg.rlat,50,ppmv2(i500mb,:))
  scatter_coast(yavg.rlon,yavg.rlat,50,ppmv4(i500mb,:))
  scatter_coast(yavg.rlon,yavg.rlat,50,ppmv6(i500mb,:))
  rtpwrite(foutNyearaverageOP,hnew,hax,pnew,pax)
 
  sartaer = ['!' topts.sarta ' fin=' foutNyearaverageOP '     fout=' foutNyearaverageRP]; eval(sartaer);

  [hnew,hax,pnew,pax] = rtpread(foutNyearaverageRP);
  scatter_coast(yavg.rlon,yavg.rlat,50,rad2bt(1231,pnew.rcalc(1520,:)))  
  scatter_coast(yavg.rlon,yavg.rlat,50,pnew.stemp-rad2bt(1231,pnew.rcalc(1520,:)))

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% compare pnew profiles with pERAI
  addpath /home/sergio/MATLABCODE/COLORMAP/
  scatter_coast(yavg.rlon,yavg.rlat,50,pnew.stemp-pERAI.stemp); title('stemp : ERA5 - ERA-Interim'); colormap(usa2); caxis([-10 +10]/1)

  scatter_coast(yavg.rlon,yavg.rlat,50,pnew.ptemp(i500mb,:)-pERAI.ptemp(i500mb,:)); title('500 mb T(z) : ERA5 - ERA-Interim'); colormap(usa2); caxis([-10 +10]/10)
  scatter_coast(yavg.rlon,yavg.rlat,50,pnew.ptemp(i850mb,:)-pERAI.ptemp(i850mb,:)); title('850 mb T(z) : ERA5 - ERA-Interim'); colormap(usa2); caxis([-10 +10]/10)

  mmwI = mmwater_rtp(hERAI,pERAI);
  mmw5 = mmwater_rtp(hnew,pnew);
  scatter_coast(yavg.rlon,yavg.rlat,50,mmw5-mmwI); title('mmw : ERA5 - ERA-Interim'); colormap(usa2); caxis([-1 +1]/1)

  mmwI_300 = mmwater_rtp(hERAI,pERAI,300);
  mmw5_300 = mmwater_rtp(hnew,pnew,300);
  scatter_coast(yavg.rlon,yavg.rlat,50,mmw5_300-mmwI_300); title('mmw to 300 mb : ERA5 - ERA-Interim'); colormap(usa2); caxis([-1 +1]/100)

end

disp(' ')
lser = ['!ls -lt ' foutNyearaverageIP ' ' foutNyearaverageOP];
eval(lser)

disp(' ')
disp('to make kCARTA jacs : see ../AIRS_gridded_STM_May2021_trendsonlyCLR/Readme_make_avg_profs_N_jacs_4608')

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
erro('check bad emis')

%{
TURNS OUT KCARTA fails on a number of profiles, since the emissivity of first 5 freqs hinges is NaN over land (Africa, Asia, the Americas)
BUT LUCKILY ONLY first 5 hingepoimts, which are almost all high alt channels

see /home/sergio/git/kcarta_gen/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/check_files_N_sizes.m

%% fix bad NAN for 4608 profiles
[h,ha,p,pa] = rtpread('summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.op.rtp');

figure(3); plot(p.efreq(:,baddy),p.emis(:,baddy))
[badX,badY] = find(isnan(p.emis));
if ~exist('baddy')
  baddy = unique(badY);
  baddy = baddy';
end
whos baddy
unique(badY) - baddy'
figure(4); clf; scatter_coast(p.rlon(baddy),p.rlat(baddy),50,p.stemp(baddy)); colormap jet; %% all over land - Africa, Asia, N/S America
unique(badX)  %%% first five channels [1 2 3 4 5]
p.efreq(unique(badX),baddy(1))   %% 645   673   701   729   757    so these are almost  all high alt channels

%% plot(p.stemp(baddy(1)) - rad2bt(h.vchan,p.rcalc(:,baddy(1))))   BAH no sarta
p.efreq(6,1)    %%% 785

%%%%%%%%%%

pnew = p;
for iii = 1 : length(baddy)
  ii = baddy(iii);
  pnew.emis(1:5,ii) = p.emis(6,ii) * ones(5,1);
  pnew.rho(1:5,ii)  = p.rho(6,ii) * ones(5,1);
end
[badX2,badY2] = find(isnan(pnew.emis)); whos badX2 badY2  YAY

mver = ['!mv summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.op.rtp summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.BADEMIS.op.rtp']; eval(mver)
rtpwrite('summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.op.rtp',h,ha,pnew,pa);

%%%%%%%%%%

[h,ha,p,pa] = rtpread('summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.ip.rtp');
[badX,badY] = find(isnan(p.emis));
unique(badY) - baddy'
sum(unique(badY) - baddy')

pnew = p;
for iii = 1 : length(baddy)
  ii = baddy(iii);
  pnew.emis(1:5,ii) = p.emis(6,ii) * ones(5,1);
  pnew.rho(1:5,ii)  = p.rho(6,ii) * ones(5,1);
end
[badX2,badY2] = find(isnan(pnew.emis)); whos badX2 badY2  YAY

mver = ['!mv summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.ip.rtp summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.BADEMIS.ip.rtp']; eval(mver)
rtpwrite('summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.ip.rtp',h,ha,pnew,pa);

%}
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%


disp(' ')
disp('now make kCARTA jacs and look at eg /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly/clust_put_together_jacs_clrERA5.m')
disp('now make kCARTA jacs and look at eg /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly/clust_put_together_jacs_clrERA5.m')
disp('now make kCARTA jacs and look at eg /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly/clust_put_together_jacs_clrERA5.m')
disp('now make kCARTA jacs and look at eg /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Aug2022_startSept2002_endAug2014_trendsonly/clust_put_together_jacs_clrERA5.m')
disp('now make kCARTA jacs and look at eg /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Aug2022_startSept2002_endAug2014_trendsonly/clust_put_together_jacs_clrERA5.m')
disp('now make kCARTA jacs and look at eg /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Aug2022_startSept2002_endAug2014_trendsonly/clust_put_together_jacs_clrERA5.m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
[~,~,junk19] = rtpread('summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp');
[~,~,junk12] = rtpread('summary_12years_all_lat_all_lon_2002_2014_monthlyERA5.rp.rtp');

plevs = flipud(load('/home/sergio/MATLABCODE/airslevels.dat')); 

figure(1); pcolor(1:64,plevs,squeeze(nanmean(reshape(junk19.ptemp,101,72,64),2))); shading interp; caxis([200 300]); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap jet; colorbar; ylim([1 1000])
figure(2); pcolor(1:64,plevs,squeeze(nanmean(reshape(junk19.ptemp - junk12.ptemp,101,72,64),2))); shading interp; caxis([-1 +1]); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(usa2); colorbar; ylim([1 1000])
 title('\delta T(z) : 19 yrs - 12 years')

figure(3); pcolor(1:64,plevs,log10(squeeze(nanmean(reshape(junk19.gas_1,101,72,64),2)))); shading interp; caxis([15 22]); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap jet; colorbar; ylim([1 1000])
figure(4); pcolor(1:64,plevs,1-squeeze(nanmean(reshape(junk19.gas_1 ./ junk12.gas_1,101,72,64),2))); shading interp; caxis([-1 +1]/10); set(gca,'ydir','reverse'); set(gca,'yscale','log'); colormap(usa2); colorbar; ylim([1 1000])
  title('1 - WV19/WV12 : 19 yrs - 12 years')

figure(5); scatter_coast(junk12.rlon,junk12.rlat,50,junk19.stemp-junk12.stemp); caxis([-1 +1]); colormap(usa2); title('\delta stemp : 19 yrs - 12 years')
%}
