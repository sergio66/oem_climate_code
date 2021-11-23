addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /asl/matlib/aslutil/
addpath /asl/matlib/maps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = '/asl/models/gistemp4/gistemp1200_GHCNv4_ERSSTv5.nc';
giss = read_netcdf_lls(fname);

%% want to do doy since 01/01/1800
doyS = change2days(2002,09,01,1800);
doyE = change2days(2019,08,31,1800);

plot(giss.time)
  line([0 length(giss.time)],[doyS doyS]);
  line([0 length(giss.time)],[doyE doyE]);
oo = find(giss.time >= doyS & giss.time <= doyE);

[Y,X] = meshgrid(giss.lat,giss.lon);

warning off
for ii = 1 : 180
  if mod(ii,50) == 0
    fprintf(1,'+')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end
  for jj = 1 : 90
    data = squeeze(giss.tempanomaly(ii,jj,:));
    aha = find(isfinite(data(oo)));
    if length(aha) > 20
      [B stats] = Math_tsfit_lin_robust(double(giss.time(oo(aha))),double(data(oo(aha))),4);
      giss_trend(ii,jj) = B(2);  
      giss_trend_err(ii,jj) = stats(2);  
    else
      giss_trend(ii,jj) = NaN;  
      giss_trend_err(ii,jj) = NaN;
    end
  end
end
warning on

simplemap(Y(:),X(:),giss_trend(:))
colormap(usa2);
caxis([-0.1 +0.1])

aslmap(1,-90:2:+90,-180:2:+180,smoothn(giss_trend',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('d/dt GISS K/yr'); 
caxis([-0.1 +0.1]);

%% so now interp2 to the Howard tiles
%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2; 
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y4608,X4608] = meshgrid(rlat,rlon);
giss_trend4608     = interp2(Y,X,giss_trend,Y4608,X4608);
giss_trend_err4608 = interp2(Y,X,giss_trend_err,Y4608,X4608);

aslmap(1,rlat65,rlon73,smoothn(giss_trend4608',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('d/dt GISS K/yr'); 
caxis([-0.1 +0.1]);

comment = 'see find_GISS_AIRSL3_trends.m';
save ChrisHTrends/giss_trends.mat giss_* X Y comment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = load(‘/home/chepplew/data/rates_anomalies/tiled/airs_l3_skt_d_fit_results.mat’)
airs_time = (1:220)*30;
oo = 1:216;  %% 12months*18years = 216
warning off
for ii = 1 : 360
  if mod(ii,100) == 0
    fprintf(1,'+')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end
  for jj = 1 : 180
    data = squeeze(x.airs.skt_d(jj,ii,:));
    aha = find(isfinite(data(oo)));
    if length(aha) > 20
      [B stats] = Math_tsfit_lin_robust(airs_time(oo(aha)),double(data(oo(aha))),4);
      airsL3_trend(ii,jj) = B(2);  
      airsL3_trend_err(ii,jj) = stats(2);  
    else
      airsL3_trend(ii,jj) = NaN;  
      airsL3_trend_err(ii,jj) = NaN;
    end
  end
end
warning on

[YL3,XL3] = meshgrid(-89.5:+89.5,-179.5:+179.5);

figure(2);
simplemap(YL3(:),XL3(:),airsL3_trend(:))
colormap(usa2);
caxis([-0.1 +0.1])

aslmap(2,-90:+90,-180:+180,smoothn(airsL3_trend',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('d/dt AIRS L3 K/yr'); 
caxis([-0.1 +0.1]);

%% so now interp2 to the Howard tiles
%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2; 
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y4608,X4608] = meshgrid(rlat,rlon);
airsL3_trend4608     = interp2(YL3,XL3,airsL3_trend,Y4608,X4608);
airsL3_trend_err4608 = interp2(YL3,XL3,airsL3_trend_err,Y4608,X4608);

aslmap(2,rlat65,rlon73,smoothn(airsL3_trend4608',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('d/dt AIRSL3 K/yr'); 
caxis([-0.1 +0.1]);

comment = 'see find_GISS_AIRSL3_trends.m';
save ChrisHTrends/airsL3_trends.mat airsL3* XL3 YL3 comment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
addpath /asl/matlib/plotutils

load /home/sergio/MATLABCODE/COLORMAP/LLS/llsmap5.mat

% MAP trend
%addpath /asl/matlib/maps          % aslmap
%addpath /home/strow/Matlab/Extra/
%addpath /asl/matlib/plotutils
%load llsmap5

fign = 2;
mopts.color = 'k';
mopts.title = 'AIRS L3 v7 SKT.D rate 2003:20  K/yr';
mopts.caxis = [-Inf Inf];
%mopts.caxis = [-0.25 0.25];
mopts.cmap  = llsmap5;
mopts.titlesize = 14;
fh = aslmap(fign, rot90(rlat65), rlon73, junk2, [-90 90], [-180 180], mopts);
%}
