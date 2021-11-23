load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/

load('llsmap5.mat');

RH000 = layeramt2RH(h,p);

pL3 = p;
pL3.stemp          = pL3.stemp          + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
pL3.ptemp(1:100,:) = pL3.ptemp(1:100,:) + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp;
pL3.gas_1(1:100,:) = pL3.gas_1(1:100,:).*(1 + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1);
pL3.gas_3(1:100,:) = pL3.gas_3(1:100,:).*(1 + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_3);
RHL3 = layeramt2RH(h,pL3);

RHL3rate = RHL3 - RH000;
zonalRHL3rate = reshape(RHL3rate,100,72,64);
zonalRHL3rate = squeeze(nanmean(zonalRHL3rate,2));

zonalrlat = rlat;
zonalplays = p.plays(1:100,3000);
figure(1); pcolor(zonalrlat,zonalplays,zonalRHL3rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('reconstruct Trate,WVrate \newline -> RH rate')

TL3rate = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp;
zonalTL3rate = reshape(TL3rate,100,72,64);
zonalTL3rate = squeeze(nanmean(zonalTL3rate,2));
figure(2); pcolor(zonalrlat,zonalplays,zonalTL3rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.15 +0.15]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('reconstruct Trate')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see  FIND_NWP_MODEL_TRENDS/driver_compute_AIRSL3_trends_desc_or_asc.m and do_the_AIRSL3_trends.m
airsL3zonal = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_zonal_rates_Sept2002_Jul2021_19yr_desc.mat');

numtimesteps = length(airsL3zonal.days);
rlat = airsL3zonal.latbins; rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
yy = []; mm = []; dd = [];
for ii = 2002 : 2021
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

h40 = struct;
p40 = struct;

h40 = h;
h40.ptype = 0;
h40.pfields = 1;
h40.ngas = 2;
h40.gunit = [20 20]';  %% g/kg
h40.glist  = [ 1  3]';

p40.rtime = [];
p40.co2ppm = [];
for ii = 1 : numtimesteps
  p40.rtime  = [p40.rtime ones(1,40)*rtime(ii)];
  p40.co2ppm = [p40.co2ppm ones(1,40)*co2ppm(ii)];
end

p40.nlevs = ones(size(p40.rtime)) * length(airsL3zonal.Tlevs);
p40.plevs = airsL3zonal.Tlevs' * ones(1,40*numtimesteps);

p40.rlat  = rlat' * ones(1,1*numtimesteps); p40.rlat = p40.rlat(:)';
p40.rlon = ones(size(p40.rlat))*0;
p40.plat = p40.rlat;
p40.plon = p40.rlon;

p40.stemp = reshape(airsL3zonal.save_stemp,1,40*numtimesteps);
p40.ptemp = reshape(permute(airsL3zonal.save_T,[2 1 3]),24,40*numtimesteps);
p40.rh    = reshape(permute(airsL3zonal.save_RH,[2 1 3]),12,40*numtimesteps);
  p40.rh(13:24,:) = ones(12,1) * p40.rh(12,:);
p40.gas_1 = reshape(permute(airsL3zonal.save_Q,[2 1 3]),12,40*numtimesteps);  %% g/kg
  p40.gas_1(13:24,:) = ones(12,1)*p40.gas_1(12,:);
p40.gas_3 = reshape(permute(airsL3zonal.save_O3,[2 1 3]),24,40*numtimesteps); %% g/kg
p40.txover = ones(size(p40.stemp)) * 1;
p40.gxover(1,:) = ones(size(p40.stemp)) * 100;  %% WV
p40.gxover(2,:) = ones(size(p40.stemp)) * 1;    %% O3

p40.scanang = zeros(size(p40.stemp));
p40.satzen = zeros(size(p40.stemp));
p40.solzen = 150 * ones(size(p40.stemp));
p40.spres = 1000 * ones(size(p40.stemp));
p40.salti = 0 * ones(size(p40.stemp));
p40.zobs = 705000 * ones(size(p40.stemp));

p40.cngwat = zeros(size(p40.stemp));
p40.cngwat2 = zeros(size(p40.stemp));
p40.cfrac = zeros(size(p40.stemp));
p40.cfrac2 = zeros(size(p40.stemp));
p40.cfrac12 = zeros(size(p40.stemp));

p40.nemis = ones(size(p40.stemp)) * 2;
p40.efreq(1,:) = ones(size(p40.stemp)) * 200;
p40.efreq(2,:) = ones(size(p40.stemp)) * 3200;
p40.emis(1,:) = ones(size(p40.stemp)) * 0.98;
p40.emis(2,:) = ones(size(p40.stemp)) * 0.98;
p40.rho = (1-p40.emis)/pi;

rtpwrite('simulate40binsAIRSL3.ip.rtp',h40,[],p40,[]);

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';;

klayerser = ['!' klayers ' fin=simulate40binsAIRSL3.ip.rtp fout=simulate40binsAIRSL3.op.rtp'];
sartaer   = ['!' sarta '   fin=simulate40binsAIRSL3.op.rtp fout=simulate40binsAIRSL3.rp.rtp'];

eval(klayerser);
eval(sartaer);
[h40x,~,p40x,~] = rtpread('simulate40binsAIRSL3.rp.rtp');
p40x.rh = layeramt2RH(h40x,p40x);
p40x.mmw = mmwater_rtp(h40x,p40x);
plot(p40x.rlat,p40x.mmw,'.')

ppmvLAY_1 = layers2ppmv(h40x,p40x,1:length(p40x.ptemp),1);
ppmvLAY_2 = layers2ppmv(h40x,p40x,1:length(p40x.ptemp),2);
ppmvLAY_3 = layers2ppmv(h40x,p40x,1:length(p40x.ptemp),3);
ppmvLAY_4 = layers2ppmv(h40x,p40x,1:length(p40x.ptemp),4);
ppmvLAY_5 = layers2ppmv(h40x,p40x,1:length(p40x.ptemp),5);
ppmvLAY_6 = layers2ppmv(h40x,p40x,1:length(p40x.ptemp),6);

tcalc = reshape(rad2bt(h40x.vchan,p40x.rcalc),2645,40,numtimesteps);;
tcalcavg = squeeze(nanmean(tcalc,2));

plot(squeeze(nanmean(tcalc(1520,:,:),3)))
plot(1:numtimesteps,tcalcavg(1520,:),1:numtimesteps,squeeze(tcalc(1520,:,:)))
plot(1:numtimesteps,squeeze(tcalc(1520,:,:)),'b.-',1:numtimesteps,tcalcavg(1520,:),'r')
plot(1:numtimesteps,nanmean(squeeze(tcalc(1520,:,:))),'b.-',1:numtimesteps,tcalcavg(1520,:),'r')
  days = (1:numtimesteps)*30/365;
  polyfit(days,nanmean(squeeze(tcalc(1520,:,:))),1); ans(1)
  addpath ../../FIND_TRENDS/
  Math_tsfit_lin_robust(days*365,nanmean(squeeze(tcalc(1520,:,:))),4); ans(2)

plot_check_WV_T_RH_L3_geo_and_spectral_rates
