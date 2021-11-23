addpath ../../FIND_TRENDS/

fchanx = h40x.vchan;
plevsx = p40x.plevs(1:97,20);
rlatx  = unique(p40.rlat);
plevsnwp = airsL3zonal.Tlevs;

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

figure(2); plot(h40x.vchan,trend); grid;  axis([640 1640 -0.1 +0.05]); plotaxis2; title('dBT/dt from AIRS L3')
figure(2); plot(h40x.vchan,trend,h40x.vchan,xtrend); grid;  axis([640 1640 -0.1 +0.05]); plotaxis2; title('dBT/dt from L3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = p40.ptemp; t = reshape(t,24,40,numtimesteps); 
for jj = 1 : 24
  for ll = 1 : 40
    data = squeeze(t(jj,ll,:));
    junk = polyfit(days,data,1);
    t2d_trendnwp(jj,ll) = junk(1);
    t2d_constnwp(jj,ll) = junk(2);
  end
end

figure(3); pcolor(rlatx,plevsnwp,t2d_trendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('T rates straight from L3 zonal levels')

t = p40x.ptemp; t = reshape(t,101,40,numtimesteps); 
for jj = 1 : 101
  for ll = 1 : 40
    data = squeeze(t(jj,ll,:));
    junk = polyfit(days,data,1);
    t2d_trend(jj,ll) = junk(1);
    t2d_constr(jj,ll) = junk(2);
  end
end

figure(4); pcolor(rlatx,p40x.plevs(1:97,20),t2d_trend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('T rates from rtp after klayers')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rh = p40.rh; rh = reshape(rh,24,40,numtimesteps); 
for jj = 1 : 24
  for ll = 1 : 40
    data = squeeze(rh(jj,ll,:));
    junk = polyfit(days,data,1);
    rh2d_trendnwp(jj,ll) = junk(1);
    rh2d_constnwp(jj,ll) = junk(2);
  end
end

figure(3); pcolor(rlatx,plevsnwp,rh2d_trendnwp); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('RH rates straight from L3 zonal levels')

rh = p40x.rh; rh = reshape(rh,100,40,numtimesteps); 
for jj = 1 : 100
  for ll = 1 : 40
    data = squeeze(rh(jj,ll,:));
    junk = polyfit(days,data,1);
    rh2d_trend(jj,ll) = junk(1);
    rh2d_constr(jj,ll) = junk(2);
  end
end

figure(4); pcolor(rlatx,p40x.plevs(1:97,20),rh2d_trend(1:97,:)); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('RH rates from rtp after klayers')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rh = p40x.rh; rh = reshape(rh,100,40,numtimesteps); rh = squeeze(nanmean(rh,2)); 
for ii = 1 : 97
  data = rh(ii,:);
  junk = polyfit(days,data,1);
  rh_trend(ii) = junk(1);
  rh_constr(ii) = junk(2);
end

figure(5); subplot(122); semilogy(rh_trend(1:97),p40x.plevs(1:97,20),'r','linewidth',2); ylim([100 1000]); set(gca,'ydir','reverse'); plotaxis2;
  title('dRH/dt'); xlim([-0.25 +0.25])
%%%%%%%%%%%%%%%%%%%%%%%%%

t = p40x.ptemp; t = reshape(t,101,40,numtimesteps); t = squeeze(nanmean(t,2)); 
for ii = 1 : 97
  data = t(ii,:);
  junk = polyfit(days,data,1);
  t_trend(ii) = junk(1);
  t_constr(ii) = junk(2);
end

figure(5); subplot(121); semilogy(t_trend(1:97),p40x.plevs(1:97,20),'r','linewidth',2); ylim([10 1000]); set(gca,'ydir','reverse');  plotaxis2;
  title('dT/dt'); xlim([-0.05 +0.05])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = p40x.stemp; t = reshape(t,40,numtimesteps);
for ii = 1 : 40
  data = t(ii,:);
  junk = polyfit(days,data,1);
  st_trend(ii) = junk(1);
  st_constr(ii) = junk(2);

  data = squeeze(tcalc(1520,ii,:));
  junk = polyfit(days,data,1);
  bt1231_trend(ii) = junk(1);
  bt1231_constr(ii) = junk(2);
end
figure(6); plot(rlatx,st_trend,rlatx,bt1231_trend); title('dST/dt and dBT1231/dt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saver = ['save reconstruct_L3_spectra_geo.mat fchanx trend xtrend rlatx '];
saver = [saver ' plevsnwp t2d_trendnwp rh2d_trendnwp '];
saver = [saver ' plevsx t2d_trend rh2d_trend rh_trend t_trend '];
saver = [saver ' zonalrlat zonalplays zonalRHL3rate zonalTL3rate st_trend bt1231_trend '];
eval(saver)
