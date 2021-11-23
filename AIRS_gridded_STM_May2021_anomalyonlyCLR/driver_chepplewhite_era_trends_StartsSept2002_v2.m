addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /asl/matlib/maps

%% looks like the rtp files were made either by 
%%   /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/Code_For_AIRS_STM_Oct2020_allstarts_Jan20XY OR
%%   /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/Code_starts_Sept2002
%%
%% example rtp filename is fnamelast = /umbc/xfs2/strow/asl/s1/sergio/MakeAvgProfsOneYear/RTP_PROFSV2/TrueERA/trueERA_timestep_lonbin_72_latbin_64_JOB_4608_cld.rtp
%%                                    where eg 4608 is TileIndex (1 .. 4608)

%{
/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/Readme
%%% USES HOWARDS 64 LATBINS, but starts Sept1, 2002
-rw-rw-r-- 1 sergio pi_strow 6824 Nov  1 21:59 clust_gather_pall_howardV3_startSept2002.m
  and then loop_gather_pall_howardV3_startSept2002.m to make the "not done" ones
-rw-rw-r-- 1 sergio pi_strow 11799 Sep 15 09:26 clust_make_profs_data_howard_bins_startSept2002.m %% makes rtp files with 7 window chans (cind1)
%}

d.home  = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/LookAtTimeSeries/RTP_PROFSV2/TrueERA/';
fnpatt  = 'trueERA_timestep_lonbin*.rtp';
d.dir   = dir([d.home fnpatt]);

%{
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
%}

homedir = pwd;

rtim = []; skt = []; wv = []; o3 = []; o3pp = zeros(4608,100,387); rh = []; rh1km = []; t = [];

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';

iyt = 1 : 4608;
k = 1;
cder = ['cd ' d.home]; eval(cder);
 
for ifn = iyt
  thename = d.dir(ifn).name;
  kind = str2num(thename(42:45)); fprintf(1,'%4i --> %4i \n',ifn,kind);

  fin = [d.dir(ifn).folder '/' d.dir(ifn).name];
  fin = [d.dir(ifn).name]; 

  klayerser = ['!' klayers ' fin=' fin ' fout=/asl/s1/sergio/JUNK/junk.op.rtp >& /asl/s1/sergio/JUNK/ugh'];
  eval(klayerser);
  [hd, ~, pd, ~] = rtpread('/asl/s1/sergio/JUNK/junk.op.rtp');
  
  allrlon(kind)  = pd.rlon(1);
  allrlat(kind)  = pd.rlat(1);

  rtim(kind,:)   = pd.rtime;
  skt(kind,:)    = pd.stemp;
  t(kind,:,:)    = pd.ptemp;
  wv(kind,:,:)   = pd.gas_1;
  co2(kind,:,:)  = pd.gas_2;
  o3(kind,:,:)   = pd.gas_3;
  % convert wv to relative humidity
  hd.ptype = 2;
  pd.palts(101,:) = pd.palts(100,:) - 240;
  [xh, xh_1km, colwater] = layeramt2RH(hd,pd);
  rh(kind,:,:)    = xh;
  rh1km(kind,:,:) = xh_1km;
  % convert o3 amt to ppmv
  np = size(pd.gas_3,2);
  [ppmvLAY,ppmvAVG,~,pavgLAY,tavgLAY,~,ppmv75,~] = layers2ppmv(hd,pd,[1:np],3);
  [mm0,nn0] = size(ppmvLAY);
  o3pp(kind,1:mm0,:)  = ppmvLAY;
    
  k = k + 1;
  %fprintf(1,'.')
end

cd /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/

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
ind = 1:386;
for ii = 1 : 4608
  data = skt(ii,:);
  data = data(ind);
  oo = find(isfinite(data));
  [B,stats] = Math_tsfit_lin_robust(xtime,data(ind),4);
  skt_trends(ii)     = B(2);
  skt_trends_err(ii) = stats(2);
end
warning on

figure(1); aslmap(1,rlat65,rlon73,smoothn((reshape(skt_trends,72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d/dt ERA SKT');  
caxis([-0.15 +0.15])
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/era_skt_trends_StartsSept2002_v2.pdf');

comment = 'see driver_chepplewhite_era_StartsSept2002_v2.m';
%save era_trends_StartSept2002_v2.mat skt_* X Y comment

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
    data = rh(ii,jj,:);
    data = data(ind);
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
%save era_trends_StartSept2002_v2.mat rh_* skt_* X Y comment plevs plays

i800 = find(plays >= 800,1);
i500 = find(plays >= 500,1);
i200 = find(plays >= 200,1);
i600 = find(plays >= 600,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); aslmap(1,rlat65,rlon73,smoothn((reshape(rh_trends(i600,:),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  
title('d/dt ERA RH(600 mb)');  caxis([-0.25 +0.25])
title('d/dt ERA RH(600 mb)');  caxis([-1 +1])
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/era_rh_600mb_global_trends_StartsSept2002_v2.pdf');

figure(1); aslmap(1,rlat65,rlon73,smoothn((reshape(rh_trends(i500,:),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  
title('d/dt ERA RH(500 mb)');  caxis([-0.25 +0.25])
title('d/dt ERA RH(500 mb)');  caxis([-1 +1])
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/era_rh_500mb_global_trends_StartsSept2002_v2.pdf');

for ii = 1 : 64
  boo = find(abs(pd.rlat-rlat(ii)) <= 0.5);
  boo = find(allrlat >= rlat65(ii) & allrlat < rlat65(ii));
  boo = (1:72) + (ii-1)*72;
  iaFound(ii) = length(boo);
  rh_trends_zonal(ii,:) = nanmean(rh_trends(:,boo),2);
end
figure(2); clf;  pcolor(rlat,plays,rh_trends_zonal'); colormap(llsmap5); caxis([-0.25 +0.25])
  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar('horizontal')
  title('Zonal d/dt RH ERA')
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/era_rh_zonal_trends_StartsSept2002_v2.pdf');

%{
clear data dataMap
data    = smoothn(rh_trends_zonal',1); 
dataMap = smoothn((reshape(rh_trends(i500,:),72,64)'),1); 
save era_rh_zonal_trends.mat rlat plays data dataMap rlat65 rlon73 
%}

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
    data = t(ii,jj,:);
    data = data(ind);
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
 
comment = 'see driver_chepplewhite_era_trends_StartsSept2002_v2.m';
%save era_trends_StartSept2002_v2.mat t_* rh_* skt_* X Y comment plevs plays

i800 = find(plays >= 800,1);
figure(1); aslmap(1,rlat65,rlon73,smoothn((reshape(t_trends(i800,:),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  
title('d/dt ERA T(800 mb)');  caxis([-0.25 +0.25])

i500 = find(plays >= 500,1);
figure(1); aslmap(1,rlat65,rlon73,smoothn((reshape(t_trends(i500,:),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  
title('d/dt ERA T(500 mb)');  caxis([-0.25 +0.25])

for ii = 1 : 64
  boo = find(abs(allrlat-rlat(ii)) <= 0.5);
  boo = find(allrlat >= rlat65(ii) & allrlat < rlat65(ii));
  boo = (1:72) + (ii-1)*72;

  iaFound(ii) = length(boo);
  t_trends_zonal(ii,:) = nanmean(t_trends(:,boo),2);

%  iaFoundOcean(ii) = length(boo);
%  t_trends_ocean_zonal(ii,:) = nanmean(t_trends(:,boo),2);
end
%plot(iaFoundOcean/72,rlat,'o-'); title('fraction of ocean tiles')

figure(2); clf; pcolor(rlat,plays,t_trends_zonal'); colormap(llsmap5); caxis([-0.05 +0.05])
  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar('horizontal')
  title('Zonal d/dt T ERA')
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/era_T_zonal_trends_StartsSept2002_v2.pdf');
%figure(2); clf; pcolor(rlat,plays,t_trends_ocean_zonal'); colormap(llsmap5); caxis([-0.05 +0.05])
%  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar

%{
clear data dataMap
data    = smoothn(t_trends_zonal',1); 
dataMap = smoothn((reshape(t_trends(i500,:),72,64)'),1); 
save era_T_zonal_trends.mat rlat plays data dataMap rlat65 rlon73 
%}

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
    data = o3pp(ii,jj,:);
    data = data(ind);
    oo = find(isfinite(data));
    if length(oo) > 20
      [B,stats] = Math_tsfit_lin_robust(xtime(oo),data(oo),4);
      o3pp_trends(jj,ii)     = B(2);
      o3pp_trends_err(jj,ii) = stats(2);
    else
      o3pp_trends(jj,ii)     = NaN;
      o3pp_trends_err(jj,ii) = NaN;
    end
  end
end
fprintf(1,'\n');
warning on

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN = plevs(1:100)-plevs(2:101); playsD = log(plevs(1:100)./plevs(2:101)); plays = flipud(playsN./playsD); plays = plays(1:100);
 
comment = 'see driver_chepplewhite_era_trends_StartsSept2002_v2.m';
%save era_trends_StartSept2002_v2.mat o3pp_* t_* rh_* skt_* X Y comment plevs plays

i025 = find(plays >= 025,1);
figure(1); aslmap(1,rlat65,rlon73,smoothn((reshape(o3pp_trends(i025,:),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  
title('d/dt ERA O3(025 mb)');  caxis([-0.01 +0.01])
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/era_o3_025mb_global_trends_StartsSept2002_v2.pdf');

for ii = 1 : 64
  boo = find(abs(allrlat-rlat(ii)) <= 0.5);
  boo = find(allrlat >= rlat65(ii) & allrlat < rlat65(ii));
  boo = (1:72) + (ii-1)*72;

  iaFound(ii) = length(boo);
  o3pp_trends_zonal(ii,:) = nanmean(o3pp_trends(:,boo),2);

end

figure(2); clf; pcolor(rlat,plays,smoothn(o3pp_trends_zonal',1)); colormap(llsmap5); caxis([-0.05 +0.05]*2)
  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([1 100]); colorbar('horizontal')
  title('Zonal d/dt O3 ERA')
%% aslprint('/home/sergio/PAPERS/AIRS/AIRS-STM-May-2021/tiletrends/Figs/era_o3_zonal_trends_StartsSept2002_v2.pdf');

%figure(2); clf; pcolor(rlat,plays,o3pp_trends_ocean_zonal'); colormap(llsmap5); caxis([-0.05 +0.05])
%  set(gca,'ydir','reverse');   set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar

%{
clear data dataMap
data    = smoothn(o3pp_trends_zonal',1); 
dataMap = smoothn((reshape(o3pp_trends(i025,:),72,64)'),1); 
save era_o3_zonal_trends.mat rlat plays data dataMap rlat65 rlon73 
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
