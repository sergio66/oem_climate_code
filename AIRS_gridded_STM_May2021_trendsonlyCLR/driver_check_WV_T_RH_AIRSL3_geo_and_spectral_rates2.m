addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCOD
%addpath ../../../FIND_TRENDS/
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME

system_slurm_stats

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 64 for the 64 latbins
%JOB = 1

load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat

load('llsmap5.mat');

RH000 = layeramt2RH(h,p);

pAIRSL3 = p;
pAIRSL3.stemp          = pAIRSL3.stemp          + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
pAIRSL3.ptemp(1:100,:) = pAIRSL3.ptemp(1:100,:) + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp;
pAIRSL3.gas_1(1:100,:) = pAIRSL3.gas_1(1:100,:).*(1 + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1);
pAIRSL3.gas_3(1:100,:) = pAIRSL3.gas_3(1:100,:).*(1 + nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_3);
RHAIRSL3 = layeramt2RH(h,pAIRSL3);

RHAIRSL3rate = RHAIRSL3 - RH000;
zonalRHAIRSL3rate = reshape(RHAIRSL3rate,100,72,64);
zonalRHAIRSL3rate = squeeze(nanmean(zonalRHAIRSL3rate,2));

zonalrlat = rlat;
zonalplays = p.plays(1:100,3000);
figure(1); pcolor(zonalrlat,zonalplays,zonalRHAIRSL3rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('reconstruct Trate,WVrate \newline -> RH rate')

TAIRSL3rate = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp;
zonalTAIRSL3rate = reshape(TAIRSL3rate,100,72,64);
zonalTAIRSL3rate = squeeze(nanmean(zonalTAIRSL3rate,2));
figure(2); pcolor(zonalrlat,zonalplays,zonalTAIRSL3rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.15 +0.15]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('reconstruct Trate')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see  FIND_NWP_MODEL_TRENDS/driver_computeAIRSL3_monthly_trends.m  and do_the_AIRSL3_trends.m
%airsl3_64x72 = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_Sept2002_Jul2021_19yr_desc.mat');
airsl3_64x72 = load('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_Sept2002_Aug2021_19yr_desc.mat');

[numtimesteps] = length(airsl3_64x72.days);
fprintf(1,'driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m : numtimesteps = %3i \n',numtimesteps)
rlat = load('latB64.mat'); rlat = 0.5*(rlat.latB2(1:end-1)+rlat.latB2(2:end));
rlon = (1:72); rlon = -177.5 + (rlon-1)*5;

yy = []; mm = []; dd = [];
for ii = 2002 : 2021
  clear yyx mmx ddx
  if ii == 2002
    inum = 4;
    yyx(1:inum) = ii;
    mmx = 9:12;
    ddx = ones(size(mmx)) * 15;
  elseif ii == 2021
    %inum = 7;
    %yyx(1:inum) = ii;
    %mmx = 1 : 7;
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

%% see ~/MATLABCODE/CRODGERS_FAST_CLOUD/driver_stage2_ESRL_set_CO2_CH4_N2O.m
time_so_far = (yy-2000) + ((mm-1)+1)/12;
co2ppm = 368 + 2.1*time_so_far;
n2oppm = 315  + (332-315)/(2020-2000)*time_so_far; n2oppm = n2oppm/1000;
ch4ppm = 1.75 + (1.875-1.750)/(2020-2000)*time_so_far;

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';;

dirout = '../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries';

co2ppm_t = [];
n2oppm_t = [];
ch4ppm_t = [];

for ii = JOB

  h72 = struct;
  p72 = struct;
  
  h72 = h;
  h72.ptype = 0;
  h72.pfields = 1;

  h72.ngas = 2;
  h72.gunit = [20 12]';  %% g/kg and VMR
  h72.glist = [ 1 3 ]';

  h72.ngas = 4;
  h72.gunit = [20 12 12 12]';  %% g/kg and VMR
  h72.glist = [ 1  3  5  6]';
  
  p72.rtime = [];
  p72.co2ppm = [];
  for iii = 1 : numtimesteps
    co2ppm_t   = [co2ppm_t ones(1,72)*co2ppm(iii)];
    n2oppm_t   = [n2oppm_t ones(1,72)*n2oppm(iii)];
    ch4ppm_t   = [ch4ppm_t ones(1,72)*ch4ppm(iii)];

    p72.rtime  = [p72.rtime ones(1,72)*rtime(iii)];
    p72.co2ppm = [p72.co2ppm ones(1,72)*co2ppm(iii)];
  end
  
  iNlev = 24;
  plevsnwp = airsl3_64x72.Tlevs;
  p72.nlevs = ones(size(p72.rtime)) * iNlev;
  p72.plevs = airsl3_64x72.Tlevs' * ones(1,72*numtimesteps);
  
  p72.rlat = rlat(ii) * ones(1,72*numtimesteps);  p72.rlat = p72.rlat(:)';
  p72.rlon = rlon' * ones(1,numtimesteps);        p72.rlon = p72.rlon(:)';
  p72.plat = p72.rlat;
  p72.plon = p72.rlon;
  
  junk = airsl3_64x72.save64x72_stemp;     junk = permute(junk,[3 2 1]);   junk = squeeze(junk(:,:,ii));    junk = junk'; p72.stemp = reshape(junk,1,72*numtimesteps);
  junk = airsl3_64x72.save64x72_T;         junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.ptemp = reshape(junk,iNlev,72*numtimesteps);
  junk = airsl3_64x72.save64x72_Q;         junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.gas_1 = reshape(junk,iNlev/2,72*numtimesteps);
  junk = airsl3_64x72.save64x72_O3;        junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.gas_3 = reshape(junk,iNlev,72*numtimesteps);
  junk = airsl3_64x72.save64x72_CO;        junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.gas_5 = reshape(junk,iNlev,72*numtimesteps);
  junk = airsl3_64x72.save64x72_CH4;       junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.gas_6 = reshape(junk,iNlev,72*numtimesteps);
  junk = airsl3_64x72.save64x72_RH;        junk = permute(junk,[4 3 2 1]); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.rh    = reshape(junk,iNlev/2,72*numtimesteps);

  junk = p72.gas_1; p72.gas_1 = zeros(size(p72.ptemp)); p72.gas_1(1:iNlev/2,:) = junk; for jj = iNlev/2+1 : iNlev; frac = (iNlev-jj+1)/(iNlev/2+1); p72.gas_1(jj,:) = junk(iNlev/2,:) * frac; end;
  junk = p72.rh;    p72.rh = zeros(size(p72.ptemp));    p72.rh(1:iNlev/2,:) = junk;    for jj = iNlev/2+1 : iNlev; frac = (iNlev-jj+1)/(iNlev/2+1); p72.rh(jj,:) = junk(iNlev/2,:) * frac; end;

%  p72.scanang = zeros(size(p72.stemp));
%  p72.satzen = zeros(size(p72.stemp));
  p72.zobs = 705000 * ones(size(p72.stemp));
  p72.scanang = ones(size(p72.stemp)) * 22;
  p72.satzen = vaconv(p72.scanang, p72.zobs, zeros(size(p72.zobs)));
  p72.solzen = 150 * ones(size(p72.stemp));
  %p72.spres = 1000 * ones(size(p72.stemp));
  %p72.salti = 0 * ones(size(p72.stemp));
  p72.spres = reshape(p.spres,72,64); p72.spres = p72.spres(:,ii) * ones(1,1*numtimesteps); p72.spres = p72.spres(:)';
  p72.salti = reshape(p.salti,72,64); p72.salti = p72.salti(:,ii) * ones(1,1*numtimesteps); p72.salti = p72.salti(:)';
  
  pcolor(rlon,1:numtimesteps,reshape(p72.stemp,72,numtimesteps)'); colormap jet; colorbar
  pcolor(rlon,1:numtimesteps,reshape(p72.spres,72,numtimesteps)'); colormap jet; colorbar
  pcolor(rlon,1:numtimesteps,reshape(p72.salti,72,numtimesteps)'); colormap jet; colorbar
  
  p72.cngwat = zeros(size(p72.stemp));
  p72.cngwat2 = zeros(size(p72.stemp));
  p72.cfrac = zeros(size(p72.stemp));
  p72.cfrac2 = zeros(size(p72.stemp));
  p72.cfrac12 = zeros(size(p72.stemp));
  
  p72.nemis = ones(size(p72.stemp)) * 2;
  p72.efreq(1,:) = ones(size(p72.stemp)) * 200;
  p72.efreq(2,:) = ones(size(p72.stemp)) * 3200;
  p72.emis(1,:) = ones(size(p72.stemp)) * 0.98;
  p72.emis(2,:) = ones(size(p72.stemp)) * 0.98;
  p72.rho = (1-p72.emis)/pi;
  
  fip = [dirout '/simulate64binsAIRSL3_' num2str(ii) '.ip.rtp'];
  fop = [dirout '/simulate64binsAIRSL3_' num2str(ii) '.op.rtp'];
  frp = [dirout '/simulate64binsAIRSL3_' num2str(ii) '.rp.rtp'];

  rtpwrite(fip,h72,[],p72,[]);
    
  klayerser = ['!' klayers ' fin=' fip ' fout=' fop];
  sartaer   = ['!' sarta '   fin=' fop ' fout=' frp];

  %%%%%%%%%%%%%%%%%%%%%%%%%
  eval(klayerser);
  [h72I,ha72I,p72I,pa72I] = rtpread(fop);
  ppmvLAY_2 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),2);
  ppmvLAY_4 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),4);
  ppmvLAY_6 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),6);
  
  i500 = find(p72I.plevs(:,1) >= 500,1);
  p72I.gas_4 = p72I.gas_4 .* (ones(101,1)*(n2oppm_t./ppmvLAY_4(i500,:)));
  %p72I.gas_6 = p72I.gas_6 .* (ones(101,1)*(ch4ppm_t./ppmvLAY_6(i500,:)));  %% use AIRS L3 CH4

  ppmvLAY_2 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),2);
  ppmvLAY_4 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),4);
  ppmvLAY_6 = layers2ppmv(h72I,p72I,1:length(p72I.stemp),6);
  
  rtpwrite(fop,h72I,ha72I,p72I,pa72I);
  eval(sartaer);
  %%%%%%%%%%%%%%%%%%%%%%%%%

  % numtimesteps = 144
  % [h72,~,p72,~] = rtpread(fip);
  [h72x,~,p72x,~] = rtpread(frp);
  p72x.rh = layeramt2RH(h72x,p72x);
  p72x.mmw = mmwater_rtp(h72x,p72x);
  plot(p72x.rlat,p72x.mmw,'.')
  plot(p72x.rlat,p72x.stemp,'.')
  plot(p72x.rlat,p72x.ptemp(80,:),'.')
  plot(p72x.rlat,p72x.gas_1(80,:),'.')
  
  ppmvLAY_1 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),1);
  ppmvLAY_2 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),2);
  ppmvLAY_3 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),3);
  ppmvLAY_4 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),4);
  ppmvLAY_5 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),5);
  ppmvLAY_6 = layers2ppmv(h72x,p72x,1:length(p72x.ptemp),6);
  
  tcalc = reshape(rad2bt(h72x.vchan,p72x.rcalc),2645,72,numtimesteps);;
  tcalcavg = squeeze(nanmean(tcalc,2));
  
  plot(squeeze(nanmean(tcalc(1520,:,:),3)))
  plot(1:numtimesteps,tcalcavg(1520,:),1:numtimesteps,squeeze(tcalc(1520,:,:)))
  plot(1:numtimesteps,squeeze(tcalc(1520,:,:)),'b.-',1:numtimesteps,tcalcavg(1520,:),'r')
  plot(1:numtimesteps,nanmean(squeeze(tcalc(1520,:,:))),'b.-',1:numtimesteps,tcalcavg(1520,:),'r')

  %%%%%%%%%%%%%%%%%%%%%%%%%
    days = (1:numtimesteps)*30/365;

    polyfit(days,nanmean(squeeze(tcalc(1520,:,:))),1); ans(1)
    %Math_tsfit_lin_robust(days*365,nanmean(squeeze(tcalc(1520,:,:))),4); ans(2)
  
    stempjunk = reshape(p72x.stemp,72,numtimesteps);
    polyfit(days,nanmean(stempjunk,1),1); ans(1)
    %Math_tsfit_lin_robust(days*365,nanmean(stempjunk),4); ans(2)

  %%%%%%%%%%%%%%%%%%%%%%%%%
    dayOFtime = change2days(yy,mm,dd,2002);
    disp('dude I just computed dayOFtime')

    polyfit(dayOFtime/365.25,nanmean(squeeze(tcalc(1520,:,:))),1); ans(1)
    %Math_tsfit_lin_robust(dayOFtime,nanmean(squeeze(tcalc(1520,:,:))),4); ans(2)
  
    stempjunk = reshape(p72x.stemp,72,numtimesteps);
    polyfit(dayOFtime/365.25,nanmean(stempjunk,1),1); ans(1)
    %Math_tsfit_lin_robust(dayOFtime,nanmean(stempjunk),4); ans(2)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%

  iType = 3;
  plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2

end
