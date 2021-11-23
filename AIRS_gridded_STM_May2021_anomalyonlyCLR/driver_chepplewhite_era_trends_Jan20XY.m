addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /asl/matlib/maps

%% looks like the rtp files were made either by 
%%   /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/Code_For_AIRS_STM_Oct2020_allstarts_Jan20XY OR
%%   /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/Code_starts_Sept2002
%%
%% example rtp filename is fnamelast = /umbc/xfs2/strow/asl/s1/sergio/MakeAvgProfs2002_2020/pall_16daytimestep_001_2003_01_01_to_2003_01_16_1.rp.rtp
%% example rtp filename is fnamelast = /umbc/xfs2/strow/asl/s1/sergio/MakeAvgProfs2002_2020/pall_16daytimestep_015_2019_08_13_to_2019_08_28_15.rp.rtp
%%                                    where eg 015 is JOBID 
%% NOTE this pall_16daytimestep_XYZ goes from 001 to 023 .. suspicious .. this is Code_For_AIRS_STM_Oct2020_allstarts_Jan20XY!!!!
%% I am assuming the latter, so look at clust_make_profs_data_howard_bins_startSept2002.m in there
%% there are 18 years to work on, so JOBID = 001,002 ... 018

%{
/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/Readme
%%% USES HOWARDS 64 LATBINS, but starts Jan1 for every year
-rw-rw-r-- 1 sergio pi_strow  3374 Sep 16 14:58 job_done_already_howard.m
-rw-rw-r-- 1 sergio pi_strow  7251 Sep 18 13:53 loop_gather_pall_howardV3.m         %% even spiffier than v2 as it puts together obs cldcal clrcal ****** and makes rtp files to make 2645 chans ********
-rw-rw-r-- 1 sergio pi_strow  6057 Sep 16 18:22 loop_gather_pall_howardV2.m         %% this just looks for mising rtps and runs sarta to****** make the missing 2645 rcalc ******
-rw-rw-r-- 1 sergio pi_strow  5829 Sep 16 10:11 loop_gather_pall_howard.m           %% this assumes all pavg*.mat files done and ****** makes rtp files with 2645 rcalc in them ******
-rw-rw-r-- 1 sergio pi_strow 11799 Sep 15 09:26 clust_make_profs_data_howard_bins.m %% makes rtp files with 7 window chans (cind1)
-rw-rw-r-- 1 sergio pi_strow 11063 Sep 12 11:12 test_make_profs_data_howard_bins.m
-rw-rw-r-- 1 sergio pi_strow  1778 Sep 14 18:56 testplot_howard_bins.m
%}

d.home  = '/asl/s1/sergio/MakeAvgProfs2002_2020/';
fnpatt  = 'pall_16daytimestep_*.rtp';
d.dir   = dir([d.home fnpatt]);
% ------------------
%  Load in the data
% ------------------
for ifn = 1:length(d.dir)
  if(ifn == 1)
    [hd, ~, pd, ~] = rtpread([d.dir(ifn).folder '/' d.dir(ifn).name]);
    rlat = pd.rlat;
    rlon = pd.rlon;
  end
  junk = strsplit(d.dir(ifn).name,{'_','.'});
  dt1(ifn)  = datenum(cell2mat(junk(4:6)),'YYYYmmdd');
  dt2(ifn)  = datenum(cell2mat(junk(8:10)),'YYYYmmdd');
  iset(ifn) = str2num(junk{3});
end
[ixt iyt] = sort(dt1);
dnum = (dt1(iyt) + dt2(iyt))/2;

rtim = []; skt = []; wv = []; o3 = []; o3 = []; rh = []; rh1km = []; t = [];

k = 1;
for ifn = iyt
  [hd, ~, pd, ~] = rtpread([d.dir(ifn).folder '/' d.dir(ifn).name]);  
  disp(d.dir(ifn).name)
  
  rtim(k,:)   = pd.rtime;
  skt(k,:)    = pd.stemp;
  t(k,:,:)    = pd.ptemp;
  wv(k,:,:)   = pd.gas_1;
  co2(k,:,:)  = pd.gas_2;
  o3(k,:,:)   = pd.gas_3;
  % convert wv to relative humidity
  hd.ptype = 2;
  pd.palts(101,:) = pd.palts(100,:) - 240;
  [xh, xh_1km, colwater] = layeramt2RH(hd,pd);
  rh(k,:,:)    = xh;
  rh1km(k,:,:) = xh_1km;
  % convert o3 amt to ppmv
  np = size(pd.gas_3,2);
  [ppmvLAY,ppmvAVG,~,pavgLAY,tavgLAY,~,ppmv75,~] = layers2ppmv(hd,pd,[1:np],3);
  o3pp(k,:,:)  = ppmvLAY;
    
  k = k + 1;
  %fprintf(1,'.')
end

fnamelast = [d.dir(ifn).folder '/' d.dir(ifn).name];
fprintf(1,'fnamelast = %s \n',fnamelast);
error('kjgdkjd')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
load('llsmap5.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off
xtime = (1:386)*16;
for ii = 1 : 4608
  data = skt(:,ii);
  [B,stats] = Math_tsfit_lin_robust(xtime,data,4);
  skt_trends(ii)     = B(2);
  skt_trends_err(ii) = stats(2);
end
warning on

figure(1); aslmap(1,rlat65,rlon73,smoothn((reshape(skt_trends,64,72)),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt ERA SKT');  
caxis([-0.15 +0.15])

comment = 'see driver_chepplewhite_era_trends.m';
%save era_trends_Jan20XY.mat skt_* X Y comment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off
xtime = (1:386)*16;
jj0 = 79;
jj0 = 1;
for jj = jj0 : 100
  if mod(jj,10) == 0
    fprintf('+');
  else
    fprintf('.');
  end
  for ii = 1 : 4608
    data = rh(:,jj,ii);
    oo = find(isfinite(data));
    if length(oo) > 20
      [B,stats] = Math_tsfit_lin_robust(xtime(oo),data(oo),4);
      rh_trends(jj,ii)     = B(2);
      rh_trends_err(jj,ii) = stats(2);
    else
      rh_trends(jj,ii)     = NaN;
      rh_trends_err(jj,ii) = NaN;
    end
  end
end
fprintf(1,'\n');
warning on

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN = plevs(1:100)-plevs(2:101); playsD = log(plevs(1:100)./plevs(2:101)); plays = flipud(playsN./playsD); plays = plays(1:100);
 
comment = 'see driver_chepplewhite_era_trends.m';
%save era_trends_Jan20XY.mat rh_* skt_* X Y comment plevs plays

i800 = find(plays >= 800,1);
figure(1); aslmap(1,rlat65,rlon73,smoothn((reshape(rh_trends(i800,:),64,72)),1),[-90 +90],[-180 +180]); colormap(llsmap5);  
title('d/dt ERA RH(800 mb)');  caxis([-0.25 +0.25])

for ii = 1 : 64
  boo = find(abs(pd.rlat-rlat(ii)) <= 0.5);
  iaFound(ii) = length(boo);
  rh_trends_zonal(ii,:) = nanmean(rh_trends(:,boo),2);
end
clf; pcolor(rlat,plays,rh_trends_zonal'); colormap(llsmap5); caxis([-0.25 +0.25])
  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off
xtime = (1:386)*16;
jj0 = 79;
jj0 = 1;
for jj = jj0 : 100
  if mod(jj,10) == 0
    fprintf('+');
  else
    fprintf('.');
  end
  for ii = 1 : 4608
    data = t(:,jj,ii);
    oo = find(isfinite(data));
    if length(oo) > 20
      [B,stats] = Math_tsfit_lin_robust(xtime(oo),data(oo),4);
      t_trends(jj,ii)     = B(2);
      t_trends_err(jj,ii) = stats(2);
    else
      t_trends(jj,ii)     = NaN;
      t_trends_err(jj,ii) = NaN;
    end
  end
end
fprintf(1,'\n');
warning on

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN = plevs(1:100)-plevs(2:101); playsD = log(plevs(1:100)./plevs(2:101)); plays = flipud(playsN./playsD); plays = plays(1:100);
 
comment = 'see driver_chepplewhite_era_trends.m';
%save era_trends_Jan20XY.mat t_* rh_* skt_* X Y comment plevs plays

i800 = find(plays >= 800,1);
figure(1); aslmap(1,rlat65,rlon73,smoothn((reshape(t_trends(i800,:),64,72)),1),[-90 +90],[-180 +180]); colormap(llsmap5);  
title('d/dt ERA T(800 mb)');  caxis([-0.25 +0.25])

for ii = 1 : 64
  boo = find(abs(pd.rlat-rlat(ii)) <= 0.5);

  iaFound(ii) = length(boo);
  boo = find(abs(pd.rlat-rlat(ii)) <= 0.5 & pd.landfrac == 0);
  t_trends_zonal(ii,:) = nanmean(t_trends(:,boo),2);

  iaFoundOcean(ii) = length(boo);
  t_trends_ocean_zonal(ii,:) = nanmean(t_trends(:,boo),2);
end
plot(iaFoundOcean/72,rlat,'o-'); title('fraction of ocean tiles')

clf; pcolor(rlat,plays,t_trends_zonal'); colormap(llsmap5); caxis([-0.05 +0.05])
  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
clf; pcolor(rlat,plays,t_trends_ocean_zonal'); colormap(llsmap5); caxis([-0.05 +0.05])
  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
