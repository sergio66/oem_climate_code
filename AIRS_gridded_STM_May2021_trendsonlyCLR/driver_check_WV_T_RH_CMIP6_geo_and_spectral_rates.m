load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE

load('llsmap5.mat');

RH000 = layeramt2RH(h,p);

pCMIP6 = p;
pCMIP6.stemp          = pCMIP6.stemp          + nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.stemp;
pCMIP6.ptemp(1:100,:) = pCMIP6.ptemp(1:100,:) + nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp;
pCMIP6.gas_1(1:100,:) = pCMIP6.gas_1(1:100,:).*(1 + nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.gas_1);
pCMIP6.gas_3(1:100,:) = pCMIP6.gas_3(1:100,:).*(1 + nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.gas_3);
RHCMIP6 = layeramt2RH(h,pCMIP6);

RHCMIP6rate = RHCMIP6 - RH000;
zonalRHCMIP6rate = reshape(RHCMIP6rate,100,72,64);
zonalRHCMIP6rate = squeeze(nanmean(zonalRHCMIP6rate,2));

zonalrlat = rlat;
zonalplays = p.plays(1:100,3000);
figure(1); pcolor(zonalrlat,zonalplays,zonalRHCMIP6rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('reconstruct Trate,WVrate \newline -> RH rate')

TCMIP6rate = nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp;
zonalTCMIP6rate = reshape(TCMIP6rate,100,72,64);
zonalTCMIP6rate = squeeze(nanmean(zonalTCMIP6rate,2));
figure(2); pcolor(zonalrlat,zonalplays,zonalTCMIP6rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.15 +0.15]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('reconstruct Trate')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see  FIND_NWP_MODEL_TRENDS/driver_computeCMIP6_monthly_trends.m  and do_the_AIRSL3_trends.m
cmip6_64x72 = load('../FIND_NWP_MODEL_TRENDS/CMIP6_atm_data_2002_09_to_2014_08.mat');

[numtimesteps,~] = size(cmip6_64x72.all.mmw);
rlat = load('latB64.mat'); rlat = 0.5*(rlat.latB2(1:end-1)+rlat.latB2(2:end));

yy = []; mm = []; dd = [];
for ii = 2002 : 2014
  clear yyx mmx ddx
  if ii == 2002
    inum = 4;
    yyx(1:inum) = ii;
    mmx = 9:12;
    ddx = ones(size(mmx)) * 15;
  elseif ii == 2021
    inum = 7;
    yyx(1:inum) = ii;
    mmx = 1 : 7;
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

h64 = struct;
p64 = struct;

h64 = h;
h64.ptype = 0;
h64.pfields = 1;
h64.ngas = 2;
h64.gunit = [21]';  %% g/g
h64.glist = [ 1]';

p64.rtime = [];
p64.co2ppm = [];
for ii = 1 : numtimesteps
  p64.rtime  = [p64.rtime ones(1,64)*rtime(ii)];
  p64.co2ppm = [p64.co2ppm ones(1,64)*co2ppm(ii)];
end

p64.nlevs = ones(size(p64.rtime)) * 19;
p64.plevs = squeeze(cmip6_64x72.all.nwp_plevs(1,:,3000))' * ones(1,64*numtimesteps);

p64.rlat  = rlat * ones(1,1*numtimesteps); p64.rlat = p64.rlat(:)';
p64.rlon = ones(size(p64.rlat))*0;
p64.plat = p64.rlat;
p64.plon = p64.rlon;

junk = cmip6_64x72.all.stemp;     junk = reshape(junk,numtimesteps,72,64);    junk = squeeze(nanmean(junk,2)); junk = junk'; p64.stemp = reshape(junk,1,64*numtimesteps);
junk = cmip6_64x72.all.nwp_ptemp; junk = reshape(junk,numtimesteps,19,72,64); junk = squeeze(nanmean(junk,3)); junk = permute(junk,[2 3 1]); p64.ptemp = reshape(junk,19,64*numtimesteps);
junk = cmip6_64x72.all.nwp_gas_1; junk = reshape(junk,numtimesteps,19,72,64); junk = squeeze(nanmean(junk,3)); junk = permute(junk,[2 3 1]); p64.gas_1 = reshape(junk,19,64*numtimesteps);
%junk = cmip6_64x72.all.nwp_gas_3; junk = reshape(junk,numtimesteps,19,72,64); junk = squeeze(nanmean(junk,3)); junk = permute(junk,[2 3 1]); p64.gas_3 = reshape(junk,19,64*numtimesteps);
junk = cmip6_64x72.all.nwp_rh;    junk = reshape(junk,numtimesteps,19,72,64); junk = squeeze(nanmean(junk,3)); junk = permute(junk,[2 3 1]); p64.rh    = reshape(junk,19,64*numtimesteps);

p64.scanang = zeros(size(p64.stemp));
p64.satzen = zeros(size(p64.stemp));
p64.solzen = 150 * ones(size(p64.stemp));
%p64.spres = 1000 * ones(size(p64.stemp));
%p64.salti = 0 * ones(size(p64.stemp));
p64.spres = nanmean(reshape(p.spres,72,64),1)' * ones(1,1*numtimesteps); p64.spres = p64.spres(:)';
p64.salti = nanmean(reshape(p.salti,72,64),1)' * ones(1,1*numtimesteps); p64.salti = p64.salti(:)';
p64.zobs = 705000 * ones(size(p64.stemp));

pcolor(rlat,1:numtimesteps,reshape(p64.stemp,64,numtimesteps)'); colormap jet; colorbar
pcolor(rlat,1:numtimesteps,reshape(p64.spres,64,numtimesteps)'); colormap jet; colorbar
pcolor(rlat,1:numtimesteps,reshape(p64.salti,64,numtimesteps)'); colormap jet; colorbar

p64.cngwat = zeros(size(p64.stemp));
p64.cngwat2 = zeros(size(p64.stemp));
p64.cfrac = zeros(size(p64.stemp));
p64.cfrac2 = zeros(size(p64.stemp));
p64.cfrac12 = zeros(size(p64.stemp));

p64.nemis = ones(size(p64.stemp)) * 2;
p64.efreq(1,:) = ones(size(p64.stemp)) * 200;
p64.efreq(2,:) = ones(size(p64.stemp)) * 3200;
p64.emis(1,:) = ones(size(p64.stemp)) * 0.98;
p64.emis(2,:) = ones(size(p64.stemp)) * 0.98;
p64.rho = (1-p64.emis)/pi;

rtpwrite('simulate64binsCMIP6.ip.rtp',h64,[],p64,[]);

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';;

klayerser = ['!' klayers ' fin=simulate64binsCMIP6.ip.rtp fout=simulate64binsCMIP6.op.rtp'];
sartaer   = ['!' sarta '   fin=simulate64binsCMIP6.op.rtp fout=simulate64binsCMIP6.rp.rtp'];

eval(klayerser);
eval(sartaer);
[h64x,~,p64x,~] = rtpread('simulate64binsCMIP6.rp.rtp');
p64x.rh = layeramt2RH(h64x,p64x);
p64x.mmw = mmwater_rtp(h64x,p64x);
plot(p64x.rlat,p64x.mmw,'.')
plot(p64x.rlat,p64x.stemp,'.')
plot(p64x.rlat,p64x.ptemp(80,:),'.')
plot(p64x.rlat,p64x.gas_1(80,:),'.')

ppmvLAY_1 = layers2ppmv(h64x,p64x,1:length(p64x.ptemp),1);
ppmvLAY_2 = layers2ppmv(h64x,p64x,1:length(p64x.ptemp),2);
ppmvLAY_3 = layers2ppmv(h64x,p64x,1:length(p64x.ptemp),3);
ppmvLAY_4 = layers2ppmv(h64x,p64x,1:length(p64x.ptemp),4);
ppmvLAY_5 = layers2ppmv(h64x,p64x,1:length(p64x.ptemp),5);
ppmvLAY_6 = layers2ppmv(h64x,p64x,1:length(p64x.ptemp),6);

tcalc = reshape(rad2bt(h64x.vchan,p64x.rcalc),2645,64,numtimesteps);;
tcalcavg = squeeze(nanmean(tcalc,2));

plot(squeeze(nanmean(tcalc(1520,:,:),3)))
plot(1:numtimesteps,tcalcavg(1520,:),1:numtimesteps,squeeze(tcalc(1520,:,:)))
plot(1:numtimesteps,squeeze(tcalc(1520,:,:)),'b.-',1:numtimesteps,tcalcavg(1520,:),'r')
plot(1:numtimesteps,nanmean(squeeze(tcalc(1520,:,:))),'b.-',1:numtimesteps,tcalcavg(1520,:),'r')
  days = (1:numtimesteps)*30/365;
  polyfit(days,nanmean(squeeze(tcalc(1520,:,:))),1); ans(1)
  addpath ../../FIND_TRENDS/
  Math_tsfit_lin_robust(days*365,nanmean(squeeze(tcalc(1520,:,:))),4); ans(2)

plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates
