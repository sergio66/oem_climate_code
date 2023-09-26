addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /asl/matlib/aslutil/
addpath /asl/matlib/maps

iDoAnom = -1;
iDoAnom = +1;

iDoGiss   = -1;
iDoAIRSL3 = -1;
iDoERA5   = +1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% so now interp2 to the Howard tiles
    %load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
    load latB64.mat
    rlat65 = latB2; rlon73 = -180 : 5 : +180;
    rlon = -180 : 5 : +180;  rlat = latB2; 
    rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
    rlat = 0.5*(rlat(1:end-1)+rlat(2:end));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iDoGiss > 0

  iNumCyclesFit = 4; %% standard fit
  iNumCyclesFit = 1; %% giss data is already an anomaly

  fname = '/asl/models/gistemp4/gistemp1200_GHCNv4_ERSSTv5.nc';
  fname = '/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/GISTEMP/F77/gistemp1200_ERSST.nc';
  
  figure(1); clf; 
  
  giss = read_netcdf_lls(fname);
  giss.lat = double(giss.lat);
  giss.lon = double(giss.lon);
  giss.time = double(giss.time);
  giss.time_bnds = double(giss.time_bnds);

  %% want to do doy since 01/01/1800
  yS = 2002; doyS = change2days(yS,09,01,1800);
  yE = 2022; doyE = change2days(yE,08,31,1800);
  
  iaNumYears = [05 10 15 20];
  iaNumYears = [20];
  for tt = 1 : length(iaNumYears)
    iNumYears = iaNumYears(tt);
    yE = yS + iNumYears; doyE = change2days(yE,08,31,1800);
    
    plot(giss.time)
      line([0 length(giss.time)],[doyS doyS]);
      line([0 length(giss.time)],[doyE doyE]);
    oo = find(giss.time >= doyS & giss.time <= doyE);
    
    [Y,X] = meshgrid(giss.lat,giss.lon);
    
    warning off
    disp(' "+" = 50, "." = 10 .. till 180 lonbins')
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
          [B stats] = Math_tsfit_lin_robust(double(giss.time(oo(aha))),double(data(oo(aha))),iNumCyclesFit);
          if iDoAnom > 0
            [giss_anom(ii,jj,:) Bjunk statsjunk] = generic_compute_anomaly(giss.time,data,oo(aha),1,iNumCyclesFit);
          end
          giss_trend(ii,jj) = B(2);  
          giss_trend_err(ii,jj) = stats.se(2);  
        else
          if iDoAnom > 0
            giss_anom(ii,jj,:) = NaN;
          end
          giss_trend(ii,jj) = NaN;  
          giss_trend_err(ii,jj) = NaN;
        end
      end
    end
    fprintf(1,'\n');
    warning on
    
    simplemap(Y(:),X(:),giss_trend(:))
    colormap(usa2);
    caxis([-1 +1]*0.15)
    
    aslmap(1,-90:2:+90,-180:2:+180,smoothn(giss_trend',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('d/dt GISS K/yr'); 
    caxis([-1 +1]*0.15);
    
    addpath /home/sergio/MATLABCODE/COLORMAP/LLS
    load llsmap5
    colormap(llsmap5)
    
    [Y4608,X4608] = meshgrid(rlat,rlon);
    giss_trend4608     = interp2(Y,X,giss_trend,Y4608,X4608);
    giss_trend_err4608 = interp2(Y,X,giss_trend_err,Y4608,X4608);
    if iDoAnom > 0
      for ii = 1 : length(aha)
        boo = giss_anom(:,:,ii);
        moo = interp2(Y,X,boo,Y4608,X4608);
        giss_anom4608(:,:,ii) = moo;
      end
      giss_time = giss.time(oo(aha));
    end

    aslmap(1,rlat65,rlon73,smoothn(giss_trend4608',1), [-90 +90],[-180 +180]);  colormap(llsmap5);  title('d/dt GISS K/yr'); 
    caxis([-1 +1] * 0.15);
    
    comment = 'see find_GISS_AIRSL3_trends.m';
    if iDoAnom < 0
      saver = ['save ChrisHTrends/giss_trends.mat giss_* X Y comment'];                       %% this is 2002-2019
      saver = ['save ChrisHTrends/giss_trends_2002_' num2str(yE) '.mat giss_* X Y comment'];  %% this is 2002-2019
    else
      saver = ['save ChrisHTrends/giss_trends_plus_anom.mat giss_* X Y comment'];                       %% this is 2002-2019
      saver = ['save ChrisHTrends/giss_trends_plus_anom_2002_' num2str(yE) '.mat giss_* X Y comment'];  %% this is 2002-2019
    end

    eval(saver);

    %% larrabee says this tile is in TWP so it is typically cloudy
    wawooY = find(giss.lat >= rlat(35),1) - 1;
    wawooX = find(giss.lon >= rlon(67),1) - 1;
    wawooY = find(giss.lat >= rlat(35),1);
    wawooX = find(giss.lon >= rlon(67),1);
    figure(2); clf
    plot(giss.time,squeeze(giss.tempanomaly(wawooX,wawooY,:)),giss.time,squeeze(giss_anom(67,35,:)))
    plot(giss.time,squeeze(giss.tempanomaly(wawooX,wawooY,:)),giss.time,squeeze(giss_anom(67,35,:)),'g',giss.time,smooth(squeeze(giss_anom(67,35,:)),23*0.25),'r','linewidth',2)
    %keyboard_nowindow

  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iDoAIRSL3 > 0
  x = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2022_20yr_desc.mat','save64x72_stemp','days');
  iaNumYears = [05 10 15 20];
  iaNumYears = [20];

  for tt = 1 : length(iaNumYears)
    iNumYears = iaNumYears(tt);
    airs_time = x.days(1:iNumYears*12);

    warning off
    disp(' "." = 10 .. till 64 latbins')
    for jj = 1 : 64
      if mod(jj,10) == 0
        fprintf(1,'.');
      end
      for ii = 1 : 72
        oo = 1 : length(airs_time);
        data = squeeze(x.save64x72_stemp(jj,ii,1:iNumYears*12));
        aha = find(isfinite(data(oo)));
        if length(aha) > 20
          [B stats] = Math_tsfit_lin_robust(airs_time(oo(aha)),double(data(oo(aha))),4);
          if iDoAnom > 0
            [airsL3_anom4608(ii,jj,:) Bjunk statsjunk] = generic_compute_anomaly(airs_time,data,oo(aha),1,4);
          end
          airsL3_trend4608(ii,jj) = B(2);  
          airsL3_trend_err4608(ii,jj) = stats.se(2);  
        else
          if iDoAnom > 0
            airsL3_anom4608(ii,jj,:) = NaN;
          end
          airsL3_trend4608(ii,jj) = NaN;  
          airsL3_trend_err4608(ii,jj) = NaN;
        end
      end
    end
    warning on
    aslmap(3,rlat65,rlon73,smoothn(airsL3_trend4608',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('d/dt AIRS L3 K/yr'); 
    caxis([-1 +1]*0.15);
        
    comment = 'see find_GISS_AIRSL3_trends.m';
    if iDoAnom < 0
      saver = ['save ChrisHTrends/airsL3_newtrends_2002_' num2str(2002+iNumYears,'%02d') '.mat airsL3* comment'];   %% this was default for 18 years
    else
      saver = ['save ChrisHTrends/airsL3_newtrends_and_anom_2002_' num2str(2002+iNumYears,'%02d') '.mat airsL3* comment'];   %% this was default for 18 years
    end
    eval(saver)

    figure(4); clf
    plot(2002.75+(1:240)/12,squeeze(airsL3_anom4608(67,35,:)),'g',2002.75+(1:240)/12,smooth(squeeze(airsL3_anom4608(67,35,:)),23*0.25),'r','linewidth',2)
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iDoERA5
  era5 = load('ERA5_atm_N_cld_data_2002_09_to_2022_08_desc.mat');
  era5daysSince2002 = change2days(era5.all.yy,era5.all.mm,era5.all.dd,2002);
  alldata = era5.all.stemp;
  disp(' "+" = 100, "." = 100 .. till 4608 tiles')

  for ii = 1 : 4608
    if mod(ii,1000) == 0
      fprintf(1,'+');
    elseif mod(ii,100) == 0
      fprintf(1,'.');
    end

    iNumYears = 20;
    oo = 1 : 12*iNumYears;
    data = alldata(:,ii);
    aha = find(isfinite(data(oo)));
    if length(aha) > 20
      [B stats] = Math_tsfit_lin_robust(era5daysSince2002(oo(aha)),double(data(oo(aha))),4);
      if iDoAnom > 0
        [era5_anom4608(ii,:) Bjunk statsjunk] = generic_compute_anomaly(era5daysSince2002,data,oo(aha),1,4);
      end
      era5_trend4608(ii) = B(2);  
      era5_trend_err4608(ii) = stats.se(2);  
    else
      if iDoAnom > 0
        era5_anom4608(ii,:) = NaN;
      end
      era5_trend4608(ii) = NaN;  
      era5_trend_err4608(ii) = NaN;
    end
  end

  warning on
  aslmap(5,rlat65,rlon73,smoothn(reshape(era5_trend4608,72,64)',1), [-90 +90],[-180 +180]);  colormap(usa2);  title('d/dt ERA5 K/yr'); 
  caxis([-1 +1]*0.15);

  comment = 'see find_GISS_AIRSL3_trends.m';
  if iDoAnom < 0
    saver = ['save ChrisHTrends/era5_newtrends_2002_' num2str(2002+iNumYears,'%02d') '.mat era5_* comment'];   %% this was default for 18 years
  else
    saver = ['save ChrisHTrends/era5_newtrends_and_anom_2002_' num2str(2002+iNumYears,'%02d') '.mat era5_* comment'];   %% this was default for 18 years
  end
  eval(saver)

  figure(6); clf
  i67_35 = (35-1)*72 + 67;
  plot(2002.75+(1:240)/12,squeeze(era5_anom4608(i67_35,:)),'g',2002.75+(1:240)/12,smooth(squeeze(era5_anom4608(i67_35,:)),23*0.25),'r','linewidth',2)

end

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
