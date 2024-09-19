addpath /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/

disp('plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2 : did preliminaries, starting the spectral and geophysical trends ...')

%% see also /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_put_together_QuantileChoose_anomalies.m
%% https://www.arm.gov/publications/proceedings/conf05/extended_abs/mlawer_ej.pdf
RRTM_bands0 = [10 250 500 630 700 820 980 1080 1180 1390 1480 1800 2080 2250 2380 2600 3000];
wx = 10:1:3000; rx = ttorad(wx,300); 
fluxx = trapz(wx,rx)/1000;
for ijunk = 2 : length(RRTM_bands0)
  junk = find(wx > RRTM_bands0(ijunk-1) & wx <= RRTM_bands0(ijunk));
  fluxx(ijunk) = trapz(wx(junk),rx(junk))/1000;
end

RRTM_bands = RRTM_bands0(4:end);
thesave.RRTM_bands = RRTM_bands;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thesave.bad_lonbins = [];
if isfield(p72x,'verybad')
  bad_lonbins = find(p72x.verybad > 0);
  if length(bad_lonbins) > 0
    disp('the are the verybad lonbins, so you may not want to believe the spectral or geophysical trends from this set of routines')
    thesave.bad_lonbins = unique(p72x.lonbin(bad_lonbins))
  end
end

fchanx = h72x.vchan;
plevsx = p72x.plevs(1:97,20);
rlatx  = unique(p72.rlat);
if iType == 6 | iType == 61
  plevsnwp = squeeze(cmip6_64x72.all.nwp_plevs(1,:,3000));
  iNlev = 19;
elseif iType == 5 | iType == 51
  plevsnwp = squeeze(era5_64x72.all.nwp_plevs(1,:,3000));
  iNlev = 37;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' .... stemp trends ...')
warning off
t = p72x.stemp; t = reshape(t,72,numtimesteps);
for iii = 1 : 72
  data = real(t(iii,:));
%  junk = polyfit(dayOFtime/365,data,1);
%  thesave.st_trend(iii) = junk(1);
%  st_constr(iii) = junk(2);

  zoo = find(isfinite(data));
  if length(find(isfinite(data))) > 16
    %junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
    [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
    thesave.sxt_anom(iii,:) = junkanom;
    thesave.xst_trend(iii)  = junk(2);
    xst_constr(iii) = junk(1);
  else
    thesave.xst_anom(iii,:) = NaN;
    thesave.xst_trend(iii) = NaN;
    xst_constr(iii) = NaN;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%

  data = real(squeeze(tcalc(1520,iii,:)));
%  junk = polyfit(dayOFtime/365,data,1);
%  thesave.bt1231_trend(iii) = junk(1);
%  bt1231_constr(iii) = junk(2);

  zoo = find(isfinite(data));
  if length(find(isfinite(data))) > 16
    %junk = Math_tsfit_lin_robust(dayOFtime,data,4);
    [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
    thesave.xbt1231_anom(iii,:) = junkanom;
    thesave.xbt1231_trend(iii) = junk(2);
    xbt1231_constr(iii) = junk(1);
  else
    thesave.xbt1231_anom(iii,:) = NaN;
    thesave.xbt1231_trend(iii) = NaN;
    xbt1231_constr(iii) = NaN;
  end

end
warning on

figure(6); plot(rlon,thesave.xst_trend,rlon,thesave.xbt1231_trend); title('dST/dt and dBT1231/dt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' .... fluxes for each of the 72 lon bins ... ')
warning off
for ilon = 1 : 72
  % old
  % FAFA = squeeze(rcalc(:,ilon,:));  
  % data = real(squeeze(trapz(h.vchan,FAFA))/1000);  

  %% see cluster_do_the_fits_airsL3_ratesv7_tiles_radiances.m
  [Y,I] = sort(h.vchan);
  fchan = Y;
  FAFA = squeeze(rcalc(I,ilon,:));
  data = real(squeeze(trapz(fchan,FAFA))/1000);  

  zoo = find(isfinite(data));
  if length(zoo) > 20
    %[junk,err] = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
    [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
    thesave.xanomflux(ilon,1,:)    = junkanom;
    thesave.xtrendflux(ilon,1)     = junk(2);
    thesave.xtrendflux_unc(ilon,1) = err.se(2);
  else
    thesave.xanomflux(ilon,1,:)    = NaN;
    thesave.xtrendflux(ilon,1)     = NaN;
    thesave.xtrendflux_unc(ilon,1) = NaN;
  end

  for flfl = 1 : length(RRTM_bands)-1
    junk = find(h.vchan >= RRTM_bands(flfl) & h.vchan < RRTM_bands(flfl+1));
    %data = squeeze(trapz(h.vchan(junk),FAFA(junk,:)))/1000;      %% OLD
    data = real(squeeze(trapz(fchan(junk),FAFA(junk,:)))/1000);      
    zoo = find(isfinite(data));
    if length(zoo) > 20
      %[junk,err] = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
      [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
      thesave.xanomflux(ilon,flfl+1,:)    = junkanom;
      thesave.xtrendflux(ilon,flfl+1)     = junk(2);
      thesave.xtrendflux_unc(ilon,flfl+1) = err.se(2);
    else
      thesave.xanomflux(ilon,flfl+1,:)    = NaN;
      thesave.xtrendflux(ilon,flfl+1)     = NaN;
      thesave.xtrendflux_unc(ilon,flfl+1) = NaN;
    end
  end

end  
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' .... 2645 channel trends and anomalies for each of the 72 lon bins ... + = 1000 . = 100')
warning off
for iii = 1 : 2645
  if mod(iii,1000) == 0
    fprintf(1,'+');
  elseif mod(iii,100) == 0
    fprintf(1,'.');
  end
  for ilon = 1 : 72
    data = real(squeeze(tcalc(iii,ilon,:)));

%    junk = polyfit(dayOFtime/365,data,1);
%    thesave.trendSpectral(iii,ilon) = junk(1);
%    constr72(iii,ilon) = junk(2);

    zoo = find(isfinite(data));
    if length(zoo) > 20
      %[junk,err] = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
      [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
      thesave.xanomSpectral(iii,ilon,:)    = junkanom;
      thesave.xtrendSpectral(iii,ilon)     = junk(2);
      thesave.xtrendSpectral_unc(iii,ilon) = err.se(2);
%      xconstr72(iii,ilon) = junk(1);
    else
      thesave.xanomSpectral(iii,ilon,:)    = NaN;
      thesave.xtrendSpectral(iii,ilon)     = NaN;
      thesave.xtrendSpectral_unc(iii,ilon) = NaN;
%      xconstr72(iii,ilon) = NaN;
    end
  end
end
thesave.freq =  h72x.vchan;

fprintf(1,'\n');
figure(1); pcolor(h72x.vchan,1:72,thesave.xtrendSpectral'); shading interp; title('Spectral Rate K/yr')
  ylabel('72 LonBins'); xlabel('Wavenumber cm-1'); caxis([-1 +1]); colorbar; colormap(llsmap5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' .... 2645 channel trends for the average latbin (over the 72 lon bins) ... should be fast')
for iii = 1 : 2645
  data = real(tcalcavg(iii,:));
%  junk = polyfit(dayOFtime/365,data,1);
%  thesave.trend(iii) = junk(1);
%  constr(iii) = junk(2);

  zoo = find(isfinite(data));
  if length(zoo) > 20
    %[junk,err] = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
    [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
    thesave.xanom(iii,:)    = junkanom;
    thesave.xtrend(iii)     = junk(2);
    thesave.xtrend_unc(iii) = err.se(2);
    xconstr(iii) = junk(1);
  else
    thesave.xanom(iii,:)    = NaN;
    thesave.xtrend(iii)     = NaN;
    thesave.xtrend_unc(iii) = NaN;
    xconstr(iii) = NaN;
  end
end

figure(2); plot(fchanx,thesave.xtrend); grid;  axis([640 1640 -0.1 +0.05]); plotaxis2; title('dBT/dt from CMIP6 or ERA5')
%figure(2); plot(fchanx,thesave.trend,fchanx,thesave.xtrend); grid;  axis([640 1640 -0.1 +0.05]); plotaxis2; title('dBT/dt from CMIP6 or ERA5')

%figure(2); plot(fchanx,thesave.xtrend,fchanx,thesave.xtrend,fchanx,nanmean(thesave.xtrendSpectral')); grid;
%cmean = nanmean(cos(pi/180*rlatx));
%figure(2); plot(fchanx,thesave.trend,fchanx,thesave.xtrend,fchanx,1/cmean*nanmean((thesave.trendSpectral .* (ones(2645,1)*cos(pi/180*rlatx)))')); grid;
%figure(2); plot(fchanx,thesave.trend,fchanx,thesave.xtrend,fchanx,1/cmean*nanmean((thesave.xtrendSpectral .* (ones(2645,1)*cos(pi/180*rlatx)))')); grid;
%  hl = legend('trend','xtrend','weighted xtrend','location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' .... Tz trends ...')
t = p72.ptemp; t = reshape(t,iNlev,72,numtimesteps); 
for jj = 1 : iNlev
  for ll = 1 : 72
    data = real(squeeze(t(jj,ll,:)));
%    junk = polyfit(dayOFtime/365,data,1);
%    thesave.t2d_trendnwp(jj,ll) = junk(1);
%    t2d_constnwp(jj,ll) = junk(2);

    zoo = find(isfinite(data));
    if length(zoo) > 20    
      %junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
      [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
      thesave.t2d_xanomnwp(jj,ll,:) = junkanom;
      thesave.t2d_xtrendnwp(jj,ll)  = junk(2);
      t2d_xconstrnwp(jj,ll) = junk(1);
    else
      thesave.t2d_xanomnwp(jj,ll,:) = NaN;
      thesave.t2d_xtrendnwp(jj,ll)  = NaN;
      t2d_xconstrnwp(jj,ll) = NaN;
    end
  end
end

figure(7); pcolor(rlon,plevsnwp,thesave.t2d_xtrendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); 
  title('T rates straight from CMIP6 or ERA5 zonal levels')

t = p72x.ptemp; t = reshape(t,101,72,numtimesteps); 
for jj = 1 : 101
  for ll = 1 : 72
    data = real(squeeze(t(jj,ll,:)));
%    junk = polyfit(dayOFtime/365,data,1);
%    thesave.t2d_trend(jj,ll) = junk(1);
%    t2d_constr(jj,ll) = junk(2);

    zoo = find(isfinite(data));
    if length(zoo) > 20    
      %junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
      [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
      thesave.t2d_xanom(jj,ll,:) = junkanom;
      thesave.t2d_xtrend(jj,ll) = junk(2);
      t2d_xconstr(jj,ll) = junk(1);
    else
      thesave.t2d_xanom(jj,ll,:) = NaN;
      thesave.t2d_xtrend(jj,ll) = NaN;
      t2d_xconstr(jj,ll) = NaN;
    end
  end
end

figure(8); pcolor(rlon,plevsx,thesave.t2d_xtrend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('T rates from rtp after klayers')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' .... WV(z) trends ...')
t = p72.gas_1; t = reshape(t,iNlev,72,numtimesteps); tX = squeeze(nanmean(t,3));
[~,~,kk] = size(t); clear tXX
for kkk = 1 : kk
 tXX(:,:,kkk) = tX;
end
t = t./tXX;

for jj = 1 : iNlev
  for ll = 1 : 72
    data = real(squeeze(t(jj,ll,:)));
    zoo = find(isfinite(data));
    if length(zoo) > 20    
      %junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
      [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
      thesave.wv2d_xanomnwp(jj,ll,:) = junkanom;
      thesave.wv2d_xtrendnwp(jj,ll) = junk(2);
      wv2d_xconstrnwp(jj,ll) = junk(1);
    else
      thesave.wv2d_xanomnwp(jj,ll,:) = NaN;
      thesave.wv2d_xtrendnwp(jj,ll) = NaN;
      wv2d_xconstrnwp(jj,ll) = NaN;
    end
  end
end

figure(9); pcolor(rlon,plevsnwp,thesave.wv2d_xtrendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-1 +1]*0.015); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('WV frac rates straight from CMIP6 or ERA5 zonal levels')

t = p72x.gas_1; t = reshape(t,101,72,numtimesteps); tX = squeeze(nanmean(t,3));
[~,~,kk] = size(t); clear tXX
for kkk = 1 : kk
 tXX(:,:,kkk) = tX;
end
t = t./tXX;
for jj = 1 : 101
  for ll = 1 : 72
    data = real(squeeze(t(jj,ll,:)));
    zoo = find(isfinite(data));
    if length(zoo) > 20    
      %junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
      [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);  
      thesave.wv2d_xanom(jj,ll,:) = junkanom;
      thesave.wv2d_xtrend(jj,ll) = junk(2);
      wv2d_xconstr(jj,ll) = junk(1);
    else
      thesave.wv2d_xanom(jj,ll,:) = NaN;
      thesave.wv2d_xtrend(jj,ll) = NaN;
      wv2d_xconstr(jj,ll) = NaN;
    end
  end
end

figure(10); pcolor(rlon,plevsx,thesave.wv2d_xtrend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-1 +1]*0.015); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('T rates from rtp after klayers')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' .... other trends if needed ...')
if iType ~= 4
  t = p72.gas_3; t = reshape(t,iNlev,72,numtimesteps); tX = squeeze(nanmean(t,3));
  [~,~,kk] = size(t); clear tXX
  for kkk = 1 : kk
   tXX(:,:,kkk) = tX;
  end
  t = t./tXX;
  
  for jj = 1 : iNlev
    for ll = 1 : 72
      data = real(squeeze(t(jj,ll,:)));
      zoo = find(isfinite(data));
      if length(zoo) > 20    
        %junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
        [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
        thesave.oz2d_xanomnwp(jj,ll,:) = junkanom;
        thesave.oz2d_xtrendnwp(jj,ll) = junk(2);
        oz2d_xconstrnwp(jj,ll) = junk(1);
      else
        thesave.oz2d_xanomnwp(jj,ll,:) = NaN;
        thesave.oz2d_xtrendnwp(jj,ll) = NaN;
        oz2d_xconstrnwp(jj,ll) = NaN;
      end
    end
  end
  
  figure(3); pcolor(rlon,plevsnwp,thesave.oz2d_xtrendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-1 +1]*0.015); 
    xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); 
    title('OZ frac rates straight from CMIP6 or ERA5 zonal levels')
  
  t = p72x.gas_3; t = reshape(t,101,72,numtimesteps); tX = squeeze(nanmean(t,3));
  [~,~,kk] = size(t); clear tXX
  for kkk = 1 : kk
   tXX(:,:,kkk) = tX;
  end
  t = t./tXX;
  for jj = 1 : 101
    for ll = 1 : 72
      data = real(squeeze(t(jj,ll,:)));
      zoo = find(isfinite(data));
      if length(zoo) > 20    
        %junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
        [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
        thesave.oz2d_xanom(jj,ll,:) = junkanom;
        thesave.oz2d_xtrend(jj,ll) = junk(2);
        oz2d_xconstr(jj,ll) = junk(1);
      else
        thesave.oz2d_xanom(jj,ll,:) = NaN;
        thesave.oz2d_xtrend(jj,ll) = NaN;
        oz2d_xconstr(jj,ll) = NaN;
      end
    end
  end
  
  figure(4); pcolor(rlon,plevsx,thesave.oz2d_xtrend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-1 +1]*0.015); 
    xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('T rates from rtp after klayers')
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' .... RH trends ...')
rh = p72.rh; rh = reshape(rh,iNlev,72,numtimesteps); 
for jj = 1 : iNlev
  for ll = 1 : 72
    data = real(squeeze(rh(jj,ll,:)));
%    junk = polyfit(dayOFtime/365,data,1);
%    thesave.rh2d_trendnwp(jj,ll) = junk(1);
%    rh2d_constnwp(jj,ll) = junk(2);

    zoo = find(isfinite(data));
    if length(zoo) > 20    
      %junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
      [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
      thesave.rh2d_xanomnwp(jj,ll,:) = junkanom;      
      thesave.rh2d_xtrendnwp(jj,ll) = junk(2);      
      rh2d_xconstrnwp(jj,ll) = junk(1);
    else
      thesave.rh2d_xanomnwp(jj,ll,:) = NaN;      
      thesave.rh2d_xtrendnwp(jj,ll) = NaN;      
      rh2d_xconstrnwp(jj,ll) = NaN;
    end
  end
end

thesave.rh2d_xtrendnwp = real(thesave.rh2d_xtrendnwp);
figure(11); pcolor(rlon,plevsnwp,thesave.rh2d_xtrendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); 
  title('RH rates straight from CMIP6 or ERA5 zonal levels')

rh = p72x.rh; rh = reshape(rh,100,72,numtimesteps); 
for jj = 1 : 100
  for ll = 1 : 72
    data = real(squeeze(rh(jj,ll,:)));
%    junk = polyfit(dayOFtime/365,data,1);
%    thesave.rh2d_trend(jj,ll) = junk(1);
%    rh2d_constr(jj,ll) = junk(2);

    zoo = find(isfinite(data));
    if length(zoo) > 20    
      %junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
      [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
      thesave.rh2d_xanom(jj,ll,:) = real(junkanom);
      thesave.rh2d_xtrend(jj,ll) = real(junk(2));
      rh2d_xconstr(jj,ll) = junk(1);
    else
      thesave.rh2d_xanom(jj,ll,:) = NaN;
      thesave.rh2d_xtrend(jj,ll) = NaN;
      rh2d_xconstr(jj,ll) = NaN;
    end

  end
end

thesave.rh2d_xtrendnwp = real(thesave.rh2d_xtrendnwp);
figure(12); pcolor(rlon,plevsx,thesave.rh2d_xtrend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  xlabel('Longitude (deg)'); set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('RH rates from rtp after klayers')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' .... NOW 101 trends ...')
disp(' .... RH 101 trends ...')
rh = p72x.rh; rh = reshape(rh,100,72,numtimesteps); rh = squeeze(nanmean(rh,2)); 
for iii = 1 : 97
  data = real(rh(iii,:));
%  junk = polyfit(dayOFtime/365,data,1);
%  thesave.rh_trend(iii) = junk(1);
%  rh_constr(iii) = junk(2);

  zoo = find(isfinite(data));
  if length(zoo) > 20    
    %junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
    [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
    thesave.rh_xanom(iii,:) = junkanom;
    thesave.rh_xtrend(iii) = junk(2);
    rh_xconstr(iii) = junk(1);
  else
    thesave.rh_xanom(iii,:) = NaN;
    thesave.rh_xtrend(iii) = NaN;
    rh_xconstr(iii) = NaN;
  end

end

figure(5); subplot(122); semilogy(thesave.rh_xtrend(1:97),plevsx,'r','linewidth',2); ylim([100 1000]); set(gca,'ydir','reverse'); plotaxis2;
  title('dRH/dt'); xlim([-0.5 +0.5])
%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' .... TZ 101 trends ...')

t = p72x.ptemp; t = reshape(t,101,72,numtimesteps); t = squeeze(nanmean(t,2)); 
for iii = 1 : 97
  data = real(t(iii,:));
%  junk = polyfit(dayOFtime/365,data,1);
%  thesave.t_trend(iii) = junk(1);
%  t_constr(iii) = junk(2);

  zoo = find(isfinite(data));
  if length(zoo) > 20    
    %junk = Math_tsfit_lin_robust(dayOFtime(zoo),data(zoo),4);
    [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dayOFtime,data,4,-1,-1);
    thesave.t_xanom(iii,:) = junkanom;
    thesave.t_xtrend(iii) = junk(2);
    t_xconstr(iii) = junk(1);
  else
    thesave.t_xanom(iii,:) = NaN;
    thesave.t_xtrend(iii) = NaN;
    t_xconstr(iii) = NaN;
  end
end

figure(5); subplot(121); semilogy(thesave.t_xtrend(1:97),plevsx,'r','linewidth',2); ylim([10 1000]); set(gca,'ydir','reverse');  plotaxis2;
  title('dT/dt'); xlim([-0.5 +0.5])
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ... finally done!! saving ....')

diroutX = dirout;
%diroutX = ' ';

ii = JOB;

disp('here I am doing this save')
if iType == 1
  saver = ['save ' diroutX 'reconstruct_umbc_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx *xconstr*'];
  %saver = [saver ' zonalrlat zonalplays zonalRHUMBCrate zonalTUMBCrate *xconstr*'];

elseif iType == 2
  saver = ['save ' diroutX 'reconstruct_merra2_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHMERRA2rate zonalTMERRA2rate *xconstr*'];
elseif iType == 5
  if ~exist('idRH')
    saver = ['save ' diroutX 'reconstruct_era5_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
    saver = [saver ' zonalrlat zonalplays zonalRHERA5rate zonalTERA5rate *xconstr*'];
  else
    saver = ['save ' diroutX 'reconstruct_era5_spectra_geo_idRH_' num2str(idRH) '_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
    saver = [saver ' zonalrlat zonalplays zonalRHERA5rate zonalTERA5rate *xconstr*'];
  end
elseif iType == 51
  saver = ['save ' diroutX 'reconstruct_era5_const_tracegas_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHERA5rate zonalTERA5rate *xconstr*'];

elseif iType == 6
  saver = ['save ' diroutX 'reconstruct_cmip6_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHCMIP6rate zonalTCMIP6rate *xconstr*'];
elseif iType == 61
  saver = ['save ' diroutX 'reconstruct_cmip6_const_tracegas_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHCMIP6rate zonalTCMIP6rate *xconstr*'];
elseif iType == 7
  saver = ['save ' diroutX 'reconstruct_amip6_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHAMIP6rate zonalTAMIP6rate *xconstr*'];

elseif iType == 3
  saver = ['save ' diroutX 'reconstruct_airsL3_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHAIRSL3rate zonalTAIRSL3rate *xconstr*'];
elseif iType == 4
  saver = ['save ' diroutX 'reconstruct_climcapsL3_spectra_geo_rlat' num2str(ii,'%02i') fstr '.mat fchanx thesave rlon rlatx '];
  saver = [saver ' zonalrlat zonalplays zonalRHAIRSCLIMCAPSL3rate zonalTAIRSCLIMCAPSL3rate *xconstr*'];

end

foutname = findstr(saver,'.mat');
foutname = [saver(6:foutname) 'mat'];

if ~exist(foutname)
  saver = [saver ' plevsnwp plevsx dayOFtime'];
  fprintf(1,'%s \n',saver');
  eval(saver)
else
  fprintf(1,'%s already exists, not saving \n',foutname)
  eval(saver)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('now use driver_read64x72flux72_save64zonalflux_ERA5_AIRSL3.m to put the individual latbins together')
disp('now use driver_read64x72flux72_save64zonalflux_ERA5_AIRSL3.m to put the individual latbins together')
disp('now use driver_read64x72flux72_save64zonalflux_ERA5_AIRSL3.m to put the individual latbins together')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
make_all_replots_check_WV_T_RH_CMIP6_geo_and_spectral_rates2

