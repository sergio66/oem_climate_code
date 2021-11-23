addpath ../../FIND_TRENDS/

fchanx = h64x.vchan;
plevsx = p64x.plevs(1:97,20);
rlatx  = unique(p64.rlat);
plevsnwp = squeeze(cmip6_64x72.all.nwp_plevs(1,:,3000));

warning off
days = (1:numtimesteps)*30/365;
for ii = 1 : 2645
  if mod(ii,1000) == 0
    fprintf(1,'+');
  elseif mod(ii,100) == 0
    fprintf(1,'.');
  end
  for ilat = 1 : 64
    data = squeeze(tcalc(ii,ilat,:));
    junk = polyfit(days,data,1);
    trend64(ii,ilat) = junk(1);
    constr64(ii,ilat) = junk(2);

    junk = Math_tsfit_lin_robust(days*365,data,4);
    xtrend64(ii,ilat) = junk(2);
    xconstr64(ii,ilat) = junk(1);
  end
end
warning on

warning off
days = (1:numtimesteps)*30/365;
for ii = 1 : 2645
  data = tcalcavg(ii,:);
  junk = polyfit(days,tcalcavg(ii,:),1);
  trend(ii) = junk(1);
  constr(ii) = junk(2);

  junk = Math_tsfit_lin_robust(days*365,tcalcavg(ii,:),4);
  xtrend(ii) = junk(2);
  xconstr(ii) = junk(1);
end
warning on

figure(2); plot(fchanx,trend); grid;  axis([640 1640 -0.1 +0.05]); plotaxis2; title('dBT/dt from CMIP6')
figure(2); plot(fchanx,trend,fchanx,xtrend); grid;  axis([640 1640 -0.1 +0.05]); plotaxis2; title('dBT/dt from CMIP6')

figure(2); plot(fchanx,trend,fchanx,xtrend,fchanx,nanmean(trend64')); grid;
cmean = nanmean(cos(pi/180*rlatx));
figure(2); plot(fchanx,trend,fchanx,xtrend,fchanx,1/cmean*nanmean((trend64 .* (ones(2645,1)*cos(pi/180*rlatx)))')); grid;
figure(2); plot(fchanx,trend,fchanx,xtrend,fchanx,1/cmean*nanmean((xtrend64 .* (ones(2645,1)*cos(pi/180*rlatx)))')); grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = p64.ptemp; t = reshape(t,19,64,numtimesteps); 
for jj = 1 : 19
  for ll = 1 : 64
    data = squeeze(t(jj,ll,:));
    junk = polyfit(days,data,1);
    t2d_trendnwp(jj,ll) = junk(1);
    t2d_constnwp(jj,ll) = junk(2);
  end
end

figure(3); pcolor(rlatx,plevsnwp,t2d_trendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('T rates straight from CMIP6 zonal levels')

t = p64x.ptemp; t = reshape(t,101,64,numtimesteps); 
for jj = 1 : 101
  for ll = 1 : 64
    data = squeeze(t(jj,ll,:));
    junk = polyfit(days,data,1);
    t2d_trend(jj,ll) = junk(1);
    t2d_constr(jj,ll) = junk(2);
  end
end

figure(4); pcolor(rlatx,plevsx,t2d_trend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('T rates from rtp after klayers')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rh = p64.rh; rh = reshape(rh,19,64,numtimesteps); 
for jj = 1 : 19
  for ll = 1 : 64
    data = squeeze(rh(jj,ll,:));
    junk = polyfit(days,data,1);
    rh2d_trendnwp(jj,ll) = junk(1);
    rh2d_constnwp(jj,ll) = junk(2);
  end
end

figure(3); pcolor(rlatx,plevsnwp,rh2d_trendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('RH rates straight from CMIP6 zonal levels')

rh = p64x.rh; rh = reshape(rh,100,64,numtimesteps); 
for jj = 1 : 100
  for ll = 1 : 64
    data = squeeze(rh(jj,ll,:));
    junk = polyfit(days,data,1);
    rh2d_trend(jj,ll) = junk(1);
    rh2d_constr(jj,ll) = junk(2);
  end
end

figure(4); pcolor(rlatx,plevsx,rh2d_trend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('RH rates from rtp after klayers')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rh = p64x.rh; rh = reshape(rh,100,64,numtimesteps); rh = squeeze(nanmean(rh,2)); 
for ii = 1 : 97
  data = rh(ii,:);
  junk = polyfit(days,data,1);
  rh_trend(ii) = junk(1);
  rh_constr(ii) = junk(2);
end

figure(5); subplot(122); semilogy(rh_trend(1:97),plevsx,'r','linewidth',2); ylim([100 1000]); set(gca,'ydir','reverse'); plotaxis2;
  title('dRH/dt'); xlim([-0.25 +0.25])
%%%%%%%%%%%%%%%%%%%%%%%%%

t = p64x.ptemp; t = reshape(t,101,64,numtimesteps); t = squeeze(nanmean(t,2)); 
for ii = 1 : 97
  data = t(ii,:);
  junk = polyfit(days,data,1);
  t_trend(ii) = junk(1);
  t_constr(ii) = junk(2);
end

figure(5); subplot(121); semilogy(t_trend(1:97),plevsx,'r','linewidth',2); ylim([10 1000]); set(gca,'ydir','reverse');  plotaxis2;
  title('dT/dt'); xlim([-0.05 +0.05])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off
t = p64x.stemp; t = reshape(t,64,numtimesteps);
for ii = 1 : 64
  data = t(ii,:);
  junk = polyfit(days,data,1);
  st_trend(ii) = junk(1);
  st_constr(ii) = junk(2);
  junk = Math_tsfit_lin_robust(days*365,data,4);
  xst_trend(ii) = junk(2);
  xst_constr(ii) = junk(1);

  data = squeeze(tcalc(1520,ii,:));
  junk = polyfit(days,data,1);
  bt1231_trend(ii) = junk(1);
  bt1231_constr(ii) = junk(2);
  junk = Math_tsfit_lin_robust(days*365,data,4);
  xbt1231_trend(ii) = junk(2);
  xbt1231_constr(ii) = junk(1);
end
warning on
figure(6); plot(rlatx,st_trend,rlatx,bt1231_trend); title('dST/dt and dBT1231/dt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saver = ['save reconstruct_cmip6_spectra_geo.mat fchanx trend* xtrend* rlatx '];
saver = [saver ' plevsnwp t2d_trendnwp rh2d_trendnwp '];
saver = [saver ' plevsx t2d_trend rh2d_trend rh_trend t_trend '];
saver = [saver ' zonalrlat zonalplays zonalRHCMIP6rate zonalTCMIP6rate st_trend bt1231_trend '];
eval(saver)

