addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE
addpath ../../../FIND_TRENDS/
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/

system_slurm_stats

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 64 for the 64 latbins
%JOB = 1

load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat

load('llsmap5.mat');

RH000 = layeramt2RH(h,p);

pAMIP6 = p;
xmiptrend = getdata_XMIP6(1);
pAMIP6.stemp          = pAMIP6.stemp          + xmiptrend.trend_stemp;
pAMIP6.ptemp(1:100,:) = pAMIP6.ptemp(1:100,:) + xmiptrend.trend_ptemp;
pAMIP6.gas_1(1:100,:) = pAMIP6.gas_1(1:100,:).*(1 + xmiptrend.trend_gas_1);
pAMIP6.gas_3(1:100,:) = pAMIP6.gas_3(1:100,:).*(1 + xmiptrend.trend_gas_3);

%{
pAMIP6.stemp          = pAMIP6.stemp          + nwp_spectral_trends_amip6_era5_airsL3_umbc.amip6_100_layertrends.stemp;
pAMIP6.ptemp(1:100,:) = pAMIP6.ptemp(1:100,:) + nwp_spectral_trends_amip6_era5_airsL3_umbc.amip6_100_layertrends.ptemp;
pAMIP6.gas_1(1:100,:) = pAMIP6.gas_1(1:100,:).*(1 + nwp_spectral_trends_amip6_era5_airsL3_umbc.amip6_100_layertrends.gas_1);
pAMIP6.gas_3(1:100,:) = pAMIP6.gas_3(1:100,:).*(1 + nwp_spectral_trends_amip6_era5_airsL3_umbc.amip6_100_layertrends.gas_3);
%}
RHAMIP6 = layeramt2RH(h,pAMIP6);

RHAMIP6rate = RHAMIP6 - RH000;
zonalRHAMIP6rate = reshape(RHAMIP6rate,100,72,64);
zonalRHAMIP6rate = squeeze(nanmean(zonalRHAMIP6rate,2));

kaboom = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat','rlat');
rlat = kaboom.rlat; 

zonalrlat = rlat;
zonalplays = p.plays(1:100,3000);
figure(1); pcolor(zonalrlat,zonalplays,zonalRHAMIP6rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.25 +0.25]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000]); title('reconstruct Trate,WVrate \newline -> RH rate')

%TAMIP6rate = nwp_spectral_trends_amip6_era5_airsL3_umbc.amip6_100_layertrends.ptemp;
TAMIP6rate = xmiptrend.trend_ptemp;
zonalTAMIP6rate = reshape(TAMIP6rate,100,72,64);
zonalTAMIP6rate = squeeze(nanmean(zonalTAMIP6rate,2));
figure(2); pcolor(zonalrlat,zonalplays,zonalTAMIP6rate); shading interp; colorbar; colormap(llsmap5); caxis([-0.15 +0.15]); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000]); title('reconstruct Trate')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hkcarta_emis,~,kcarta_emis,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp');
ind = (1:72) + (JOB-1)*72;
[hkcarta_emis,kcarta_emis] = subset_rtp(hkcarta_emis,kcarta_emis,[],[],ind);

%% see  FIND_NWP_MODEL_TRENDS/driver_computeAMIP6_monthly_trends.m  and do_the_AIRSL3_trends.m
amip6_64x72 = load('../../FIND_NWP_MODEL_TRENDS/AMIP6_atm_data_2002_09_to_2014_08.mat');

[numtimesteps,~] = size(amip6_64x72.all.mmw);
rlat = load('latB64.mat'); rlat = 0.5*(rlat.latB2(1:end-1)+rlat.latB2(2:end));
rlon = (1:72); rlon = -177.5 + (rlon-1)*5;

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
time_so_far = (yy-2000) + ((mm-1)+1)/12;
co2ppm = 370 + 2.2*((yy+mm/12)-2002);
%% see ~/MATLABCODE/CRODGERS_FAST_CLOUD/driver_stage2_ESRL_set_CO2_CH4_N2O.m
co2ppm = 368 + 2.1*time_so_far;
n2oppm = 315  + (332-315)/(2020-2000)*time_so_far; n2oppm = n2oppm/1000;
ch4ppm = 1.75 + (1.875-1.750)/(2020-2000)*time_so_far;

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';;

dirout = '../../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries';
dirout = 'SimulateTimeSeries/AMIP6/';

co2ppm_t = [];
n2oppm_t = [];
ch4ppm_t = [];

nemis_t = [];
emis_t  = [];
efreq_t = [];
rho_t   = [];

for ii = JOB

  h72 = struct;
  p72 = struct;
  
  h72 = h;
  h72.ptype = 0;
  h72.pfields = 1;
  h72.ngas = 1;
  h72.gunit = [21]';  %% g/g
  h72.glist = [ 1]';
  
  p72.rtime = [];
  p72.co2ppm = [];
  for iii = 1 : numtimesteps
    nemis_t  = [nemis_t  kcarta_emis.nemis];
    efreq_t  = [efreq_t  kcarta_emis.efreq];
    emis_t   = [emis_t   kcarta_emis.emis];
    rho_t    = [rho_t    kcarta_emis.rho];

    co2ppm_t   = [co2ppm_t ones(1,72)*co2ppm(iii)];
    n2oppm_t   = [n2oppm_t ones(1,72)*n2oppm(iii)];
    ch4ppm_t   = [ch4ppm_t ones(1,72)*ch4ppm(iii)];

    p72.rtime  = [p72.rtime ones(1,72)*rtime(iii)];
    p72.co2ppm = [p72.co2ppm ones(1,72)*co2ppm(iii)];
  end

  iNlev = 19;
  plevsnwp = squeeze(amip6_64x72.all.nwp_plevs(1,:,3000))';  
  p72.nlevs = ones(size(p72.rtime)) * iNlev;
  p72.plevs = squeeze(amip6_64x72.all.nwp_plevs(1,:,3000))' * ones(1,72*numtimesteps);
  
  p72.rlat  = rlat(ii) * ones(1,72*numtimesteps); p72.rlat = p72.rlat(:)';
  p72.rlon = rlon' * ones(1,numtimesteps);        p72.rlon = p72.rlon(:)';
  p72.plat = p72.rlat;
  p72.plon = p72.rlon;
  
  junk = amip6_64x72.all.stemp;     junk = reshape(junk,numtimesteps,72,64);    junk = squeeze(junk(:,:,ii)); junk = junk'; p72.stemp = reshape(junk,1,72*numtimesteps);
  junk = amip6_64x72.all.nwp_ptemp; junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.ptemp = reshape(junk,iNlev,72*numtimesteps);
  junk = amip6_64x72.all.nwp_gas_1; junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.gas_1 = reshape(junk,iNlev,72*numtimesteps);
  %junk = amip6_64x72.all.nwp_gas_3; junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.gas_3 = reshape(junk,iNlev,72*numtimesteps);
  junk = amip6_64x72.all.nwp_rh;    junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); p72.rh    = reshape(junk,iNlev,72*numtimesteps);

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
  
  p72.nemis = nemis_t;
  p72.efreq = efreq_t;
  p72.emis  = emis_t;
  p72.rho   = rho_t;

  fip = [dirout '/simulate64binsAMIP6_' num2str(ii) '.ip.rtp'];
  fop = [dirout '/simulate64binsAMIP6_' num2str(ii) '.op.rtp'];
  frp = [dirout '/simulate64binsAMIP6_' num2str(ii) '.rp.rtp'];

  klayers_sarta_check_WV_T_RH_geo_and_spectral_rates2

  iType = 7;  
  plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2

end
