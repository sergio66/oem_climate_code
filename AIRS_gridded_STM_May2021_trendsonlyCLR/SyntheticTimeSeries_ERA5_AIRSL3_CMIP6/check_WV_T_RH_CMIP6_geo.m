addpath /home/sergio/MATLABCODE

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /asl/matlib/maps
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE
addpath ../../FIND_TRENDS/

dirout = '../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries';

load('llsmap5.mat');

%% see  FIND_NWP_MODEL_TRENDS/driver_computeCMIP6_monthly_trends.m  and do_the_AIRSL3_trends.m
cmip6_64x72 = load('../FIND_NWP_MODEL_TRENDS/CMIP6_atm_data_2002_09_to_2014_08.mat');

addpath ../FIND_NWP_MODEL_TRENDS/
all = cmip6_64x72.all;

dayOFtime = change2days(all.yy,all.mm,all.dd,2002);
numtimesteps = length(dayOFtime);

%computeERA5_surface_trends

%%%%%%%%%%%%%%%%%%%%%%%%%

stemp_all = [];
for ii = 1 : 64
  
  fip = [dirout '/simulate64binsCMIP6' num2str(ii) '.ip.rtp'];
  fop = [dirout '/simulate64binsCMIP6' num2str(ii) '.op.rtp'];
  frp = [dirout '/simulate64binsCMIP6' num2str(ii) '.rp.rtp'];

  ind = (1:72) + (ii-1)*72;
  [h72x,~,p72x,~] = rtpread(frp);
  junk = reshape(p72x.stemp,72,numtimesteps); stemp_all(ind,:) = junk;
end

meanST = nanmean(stemp_all,2); pcolor(reshape(meanST,72,64)'); shading interp; colorbar
boo = 1:4608; pcolor(reshape(boo,72,64)'); shading interp; colorbar

origdata = load('../FIND_NWP_MODEL_TRENDS/CMIP6_atm_data_2002_09_to_2014_08.mat');
origtrends = load('../FIND_NWP_MODEL_TRENDS/CMIP6_atm_data_2002_09_to_2014_08_trends.mat');

%plot(origtrends.trend_stemp - trend_stemp)
plot(stemp_all - cmip6_64x72.all.stemp')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obsrates = load('gather_tileCLRnight_rates_fits.mat');
zonalavg = load('reconstruct_cmip6_spectra_geo.mat');

[numtimesteps,~] = size(cmip6_64x72.all.mmw);
rlat = load('latB64.mat'); rlat65 = rlat.latB2; rlat = 0.5*(rlat.latB2(1:end-1)+rlat.latB2(2:end));
rlon73 = (1:73); rlon73 = -180 + (rlon73-1)*5;  rlon = (1:72); rlon = -177.5 + (rlon-1)*5;
[Y,X] = meshgrid(rlat,rlon); X = X(:); Y = Y(:);

yy = []; mm = []; dd = [];
for ii = 2002 : 2014
  clear yyx mmx ddx
  if ii == 2002
    inum = 4;
    yyx(1:inum) = ii;
    mmx = 9:12;
    ddx = ones(size(mmx)) * 15;
  elseif ii == 2014
    inum = 8;
    yyx(1:inum) = ii;
    mmx = 1 : 8;
    ddx = ones(size(mmx)) * 15;
  else
    inum = 12;
    yyx(1:inum) = ii;
    mmx = 1:12;
    ddx = ones(size(mmx)) * 15;
  end
  fprintf(1,'%4i %2i \n',[ii inum])
  yy = [yy yyx];
  mm = [mm mmx];
  dd = [dd ddx];
end
rtime = utc2taiSergio(yy,mm,dd,ones(size(yy))*12.0);
co2ppm = 370 + 2.2*((yy+mm/12)-2002);

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stemp = [];
spres = [];
ptemp = [];
gas_1 = [];
gas_3 = [];
plevs = [];
rcalc = [];

for ii = 1 : 64
  
  fip = [dirout '/simulate64binsCMIP6' num2str(ii) '.ip.rtp'];
  fop = [dirout '/simulate64binsCMIP6' num2str(ii) '.op.rtp'];
  frp = [dirout '/simulate64binsCMIP6' num2str(ii) '.rp.rtp'];

  ind = (1:72) + (ii-1)*72;
  [h72x,~,p72x,~] = rtpread(frp);
  junk = reshape(p72x.stemp,72,numtimesteps);     junk = nanmean(junk,2); stemp(ind) = junk;
  junk = reshape(p72x.spres,72,numtimesteps);     junk = nanmean(junk,2); spres(ind) = junk;
  junk = reshape(p72x.plevs,101,72,numtimesteps); junk = nanmean(junk,3); plevs(:,ind) = junk;
  junk = reshape(p72x.ptemp,101,72,numtimesteps); junk = nanmean(junk,3); ptemp(:,ind) = junk;
  junk = reshape(p72x.gas_1,101,72,numtimesteps); junk = nanmean(junk,3); gas_1(:,ind) = junk;
  junk = reshape(p72x.gas_3,101,72,numtimesteps); junk = nanmean(junk,3); gas_3(:,ind) = junk;
  junk = reshape(p72x.rcalc,2645,72,numtimesteps); junk = nanmean(junk,3); rcalc(:,ind) = junk;
end

whos stemp ptemp rcalc

for ii = 1 : 4; figure(ii); colormap jet; end
figure(1); pcolor(reshape(X,72,64)',reshape(Y,72,64)',reshape(spres,72,64)'); shading interp; colorbar
figure(2); pcolor(reshape(X,72,64)',reshape(Y,72,64)',reshape(stemp,72,64)'); shading interp; colorbar
figure(3); pcolor(reshape(X,72,64)',reshape(Y,72,64)',reshape(rad2bt(1231,rcalc(1520,:)),72,64)'); shading interp; colorbar
figure(4); pcolor(reshape(X,72,64)',reshape(Y,72,64)',reshape(stemp - rad2bt(1231,rcalc(1520,:)),72,64)'); shading interp; colorbar;  colormap(llsmap5); caxis([-6 +6])

figure(5); junk = reshape(ptemp,101,72,64); junk = squeeze(nanmean(junk,2)); pcolor(junk); caxis([200 300]); colormap jet; shading interp;  colorbar
figure(5); junk = reshape(log10(gas_1),101,72,64); junk = squeeze(nanmean(junk,2)); pcolor(junk); caxis([10 22]); colormap jet; shading interp;  colorbar
figure(5); junk = reshape(log10(gas_3),101,72,64); junk = squeeze(nanmean(junk,2)); pcolor(junk); caxis([15 18]); colormap jet; shading interp;  colorbar
pause(0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 64
  ind = (1:72) + (ii-1)*72;
  fmat = [dirout '/reconstruct_cmip6_spectra_geo_rlat' num2str(ii,'%02d') '.mat'];
  x = load(fmat);
  %xtrendSpectral(:,ind) = x.thesave.xtrend72;
  xtrendSpectral(:,ind) = x.thesave.xtrendSpectral;
  xtrend_bt1231(ind) = x.thesave.xbt1231_trend;
  xtrend_st(ind) = x.thesave.xst_trend;
  xtrend_Tz(:,ind) = x.thesave.t2d_xtrend;
  xtrend_RHz(:,ind) = x.thesave.t2d_xtrend;
end
figure(6); clf; colormap(llsmap5); pcolor(reshape(X,72,64)',reshape(Y,72,64)',reshape(xtrend_st,72,64)'); shading interp; colorbar; caxis([-0.15 +0.15]);
figure(7); clf; colormap(llsmap5); pcolor(reshape(X,72,64)',reshape(Y,72,64)',reshape(xtrend_bt1231,72,64)'); shading interp; colorbar

figure(6); clf; aslmap(6,rlat65,rlon73,smoothn(reshape(xtrend_st,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-0.15 +0.15]);
figure(7); clf; aslmap(7,rlat65,rlon73,smoothn(reshape(xtrend_bt1231,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-0.15 +0.15]);

%old_cmip6_trends = load('../FIND_NWP_MODEL_TRENDS/CMIP6_atm_data_2002_09_to_2021_07_trends_desc.mat');
%figure(8); clf; aslmap(8,rlat65,rlon73,smoothn(reshape(old_cmip6_trends.trend_stemp,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-0.15 +0.15]);

figure(9);  clf; pcolor(origdata.all.stemp'-stemp_all); shading interp; colorbar; title('ORIG STEMP - what is in the rtp files')
figure(10); clf; plot(origtrends.trend_stemp - xtrend_st);

figure(11); plot(h72x.vchan,mean(xtrendSpectral,2),h72x.vchan,mean(obsrates.rates,2)); grid; xlim([640 1640]); title('CMIP6'); plotaxis2; 
  hl = legend('CMIP6','AIRS obs','location','best')

pause(0.1);
error('lgksg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE/TIME
days = (1:numtimesteps)*30/365;
[sum(yy-origdata.all.yy) sum(mm-origdata.all.mm) sum(dd-origdata.all.dd)]
dayOFtime = change2days(yy,mm,dd,2002);
warning off
for ii = 1 : 4608
  if mod(ii,1000) == 0
    fprintf(1,'+')
  elseif mod(ii,100) == 0
    fprintf(1,'.')
  end
  addpath ../../FIND_TRENDS/
  junk = Math_tsfit_lin_robust(days*365,stemp_all(ii,:),4);; 
  trend1(ii) = junk(2);

  addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
  junk = Math_tsfit_lin_robust(days*365,stemp_all(ii,:),4);; 
  trend2(ii) = junk(2);  

  addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
  junk = Math_tsfit_lin_robust(dayOFtime,stemp_all(ii,:),4);; 
  trend3(ii) = junk(2);  
end
fprintf(1,'\n');
warning on
plot(1:length(trend1),trend1-trend2,'+-',1:length(trend1),trend1-trend3)
plot(origtrends.trend_stemp-trend3)

plot(trend1-origtrends.trend_stemp)
plot(trend3-origtrends.trend_stemp)
