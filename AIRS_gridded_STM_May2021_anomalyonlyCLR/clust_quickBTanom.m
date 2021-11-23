addpath /asl/matlib/aslutil
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies
addpath /home/sergio/MATLABCODE/COLORMAP

JOBIN = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% tile number
thelist = load('do_these_tiles1.txt');
thelist = load('do_these_tiles2.txt');
thelist = load('do_these_tiles3.txt');
JOB = thelist(JOBIN);

fprintf(1,'clust slurm index JOBIN = %4i will be doing anomalies for tile %4i \n',JOBIN,JOB);

%JOB = 3000
%JOB = 4608
%JOB = 3000

if JOB > 4608
  JOB
  error('JOB cannot be larger than 4608')
end

YY = floor((JOB-1)/72) + 1;
XX = JOB-(YY-1)*72;

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');

%JOB = 32
%JOB = 47
%JOB  = 47

%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
load latB64.mat
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;
Xstretch = X(:); 
Ystretch = Y(:); 
%Xuse = Xstretch(JOB); Xstr = num2str(Xuse,'%02d');
%Yuse = Ystretch(JOB); Ystr = num2str(Yuse,'%02d');

Yuse = YY; Ystr = num2str(Yuse,'%02d'); fprintf(1,'LAT : Y(JOB) p.rlat(JOB) Yind = %8.4f %8.4f %2i \n',[Ystretch(JOB)  p.rlat(JOB) YY])
Xuse = XX; Xstr = num2str(Xuse,'%02d'); fprintf(1,'LON : X(JOB) p.rlon(JOB) Xind = %8.4f %8.4f %2i \n',[Xstretch(JOB)  p.rlon(JOB) XX])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdirIN = ['/umbc/xfs2/strow/asl/s1/strow/home/Work/Airs/Tiles/Data/Quantv1/'];
fdirIN = [fdirIN '/LatBin' Ystr '/LonBin' Xstr '/'];
filenameIN = [fdirIN '/summarystats_LatBin' Ystr '_LonBin' Xstr '_timesetps_001_412_V1.mat'];

iQuantile = 08; iQstr = num2str(iQuantile,'%02d');
iQuantile = 16; iQstr = num2str(iQuantile,'%02d');

fdirOUT     = ['ANOM/LatBin' Ystr '/LonBin' Xstr '/v1_Quantile' iQstr '/'];
if ~exist(fdirOUT)
  mker = ['!mkdir -p ' fdirOUT];
  eval(mker);
end
filenameOUT = [fdirOUT '/desc_anaomaly_LatBin' Ystr '_LonBin' Xstr '_Quantile' iQstr '_timesetps_001_412_V1.mat'];

if exist(filenameOUT)
  fprintf(1,'%s exists, exit \n',filenameOUT);
else
  tic
  load h2645structure.mat
  loader = ['load ' filenameIN];
  eval(loader);
  rads = squeeze(rad_quantile_desc(:,:,iQuantile))';
  bt = rad2bt(h.vchan,rads);

  [mm,nn] = size(rads);
  xin = (1:nn);
  xin = xin-1;
  xin = xin*16;
  k = ones(1,nn);
  k = true(size(k));
  warning off
  for ii = 1 : 2645
    if mod(ii,1000) == 0
      fprintf(1,'+');
    elseif mod(ii,100) == 0
      fprintf(1,'.');
    end
    [B, stats] = Math_tsfit_lin_robust_finitewrapper(xin,rads(ii,:),4);
    deriv = drdbt(h.vchan(ii),mean(rads(ii,:)));
    rad_trend(ii) = B(2);
    BT_trend(ii)  = B(2)/deriv;
    [bt_anom rad_anom] = compute_anomaly(k,xin,B,h.vchan(ii),rads(ii,:)');
    Ball(ii,:) = B;
    dBall(ii,:) = stats.se;
    bt_anom_all(ii,:)  = bt_anom; 
    rad_anom_all(ii,:) = rad_anom; 
  end
  fprintf(1,'\n');
  warning on
  toc
  comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_May2021_anomalyonlyCLR/clust_quickBTanom.m';
  imagesc(bt_anom_all); colorbar; caxis([-1 +1]); title('BT anomaly'); colormap(usa2);
  pcolor(h.vchan,1:nn,bt_anom_all'); shading flat; colorbar; caxis([-1 +1]); title('BT anomaly'); colormap(usa2);

  saver = ['save ' filenameOUT ' *all comment rad_trend rad_trend'];
  eval(saver)
end


