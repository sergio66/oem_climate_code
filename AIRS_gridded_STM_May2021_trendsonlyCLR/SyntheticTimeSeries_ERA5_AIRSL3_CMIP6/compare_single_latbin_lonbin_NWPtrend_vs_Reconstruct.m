iLonBin = find(chisqr(5,:) == max(chisqr(5,:)) );
fprintf(1,'worst WV chisqr for Latbin %2i is at LonBin %2i \n',iLatBin,iLonBin);

%iLonBin = find(chisqr(5,:) == min(chisqr(5,:)) );
%fprintf(1,'best WV chisqr for Latbin %2i is at LonBin %2i \n',iLatBin,iLonBin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% lot of this is from driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE
%addpath ../../../FIND_TRENDS/
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/

JOB = iLatBin;
[hkcarta_emis,~,kcarta_emis,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp');
ind = (1:72)  + (JOB-1)*72;
ind = iLonBin + (JOB-1)*72;
[hkcarta_emis,kcarta_emis] = subset_rtp(hkcarta_emis,kcarta_emis,[],[],ind);

%load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat
wah = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat','h');
h = wah.h;
wah = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_v2_unc.mat','p');
p = wah.p;

%% see  FIND_NWP_MODEL_TRENDS/driver_computeERA5_monthly_trends.m  and do_the_AIRSL3_trends.m
if ~exist('era5_64x72')
  %era5_64x72 = load('../../FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_07_desc.mat');
  era5_64x72 = load('../../FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_08_desc.mat');
end

[numtimesteps,~] = size(era5_64x72.all.mmw);
fprintf(1,'driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m : numtimesteps = %3i \n',numtimesteps)

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

yy2002 = yy + (mm)/12;

rtime = utc2taiSergio(yy,mm,dd,ones(size(yy))*12.0);
time_so_far = (yy-2000) + ((mm-1)+1)/12;
co2ppm = 370 + 2.2*((yy+mm/12)-2002);
%% see ~/MATLABCODE/CRODGERS_FAST_CLOUD/driver_stage2_ESRL_set_CO2_CH4_N2O.m
co2ppm = 368 + 2.1*time_so_far;
n2oppm = 315  + (332-315)/(2020-2000)*time_so_far; n2oppm = n2oppm/1000;
ch4ppm = 1.75 + (1.875-1.750)/(2020-2000)*time_so_far;

iRampCO2_CH4_N2O = input('ramp up CO2/CH4/N2O (-1/+1) : ');
if length(iRampCO2_CH4_N2O) == 0
  iRampCO2_CH4_N2O = +1;
end
if iRampCO2_CH4_N2O < 0
  co2ppm = ones(size(co2ppm)) * mean(co2ppm);
  n2oppm = ones(size(n2oppm)) * mean(n2oppm);
  ch4ppm = ones(size(ch4ppm)) * mean(ch4ppm);
end

klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta   = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';;

dirout = 'JUNK/';

co2ppm_t = [];
n2oppm_t = [];
ch4ppm_t = [];

nemis_t = [];
emis_t  = [];
efreq_t = [];
rho_t   = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the profile timeseries
for ii = JOB

  h72 = struct;
  p72 = struct;
  
  h72 = h;
  h72.ptype = 0;
  h72.pfields = 1;
  h72.ngas = 2;
  h72.gunit = [21 21]';  %% g/g
  h72.glist = [ 1 3]';
  
  p72.rtime = [];
  p72.co2ppm = [];
  i72 = 1;

  for iii = 1 : numtimesteps
    nemis_t  = [nemis_t  kcarta_emis.nemis];
    efreq_t  = [efreq_t  kcarta_emis.efreq];
    emis_t   = [emis_t   kcarta_emis.emis];
    rho_t    = [rho_t    kcarta_emis.rho];

    co2ppm_t = [co2ppm_t ones(1,i72)*co2ppm(iii)];
    n2oppm_t = [n2oppm_t ones(1,i72)*n2oppm(iii)];
    ch4ppm_t = [ch4ppm_t ones(1,i72)*ch4ppm(iii)];

    p72.rtime  = [p72.rtime ones(1,i72)*rtime(iii)];
    p72.co2ppm = [p72.co2ppm ones(1,i72)*co2ppm(iii)];
  end
  
  iNlev = 37;
  p72.nlevs = ones(size(p72.rtime)) * iNlev;
  p72.plevs = squeeze(era5_64x72.all.nwp_plevs(1,:,3000))' * ones(1,i72*numtimesteps);
  
  p72.rlat = rlat(ii) * ones(1,i72*numtimesteps);   p72.rlat = p72.rlat(:)';
  p72.rlon = rlon(iLonBin) * ones(1,numtimesteps);  p72.rlon = p72.rlon(:)';
  p72.plat = p72.rlat;
  p72.plon = p72.rlon;
  
  junk = era5_64x72.all.stemp;     junk = reshape(junk,numtimesteps,72,64);       junk = squeeze(junk(:,:,ii));   junk = junk';                 stemp = junk(iLonBin,:);              p72.stemp = stemp;
  junk = era5_64x72.all.nwp_ptemp; junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); ptemp = squeeze(junk(:,iLonBin,:));   p72.ptemp = ptemp;
  junk = era5_64x72.all.nwp_gas_1; junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); gas_1 = squeeze(junk(:,iLonBin,:));   p72.gas_1 = gas_1;
  junk = era5_64x72.all.nwp_gas_3; junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); gas_3 = squeeze(junk(:,iLonBin,:));   p72.gas_3 = gas_3;
  junk = era5_64x72.all.nwp_rh;    junk = reshape(junk,numtimesteps,iNlev,72,64); junk = squeeze(junk(:,:,:,ii)); junk = permute(junk,[2 3 1]); rh    = squeeze(junk(:,iLonBin,:));   p72.rh = rh;

%  p72.scanang = zeros(size(p72.stemp));
%  p72.satzen = zeros(size(p72.stemp));
  p72.zobs = 705000 * ones(size(p72.stemp));
  p72.scanang = ones(size(p72.stemp)) * 22;
  p72.satzen = vaconv(p72.scanang, p72.zobs, zeros(size(p72.zobs)));
  p72.satzen = ones(size(p72.stemp)) * 24; %%% to match what is in home/sergio/KCARTA/WORK/RUN/RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp which was used to make jacs
  p72.solzen = 150 * ones(size(p72.stemp));
  %p72.spres = 1000 * ones(size(p72.stemp));
  %p72.salti = 0 * ones(size(p72.stemp));
  spres = reshape(p.spres,72,64); p72.spres = spres(iLonBin,iLatBin) * ones(1,1*numtimesteps);
  salti = reshape(p.salti,72,64); p72.salti = salti(iLonBin,iLatBin) * ones(1,1*numtimesteps);

  plot(yy2002,p72.stemp);
  plot(yy2002,p72.spres);
  plot(yy2002,p72.salti);
  
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
    
  size_struct_fields(p72)

  fip = [dirout '/simulate64binsERA5_' num2str(ii) '.ip.rtp'];
  fop = [dirout '/simulate64binsERA5_' num2str(ii) '.op.rtp'];
  frp = [dirout '/simulate64binsERA5_' num2str(ii) '.rp.rtp'];

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
  p72I.gas_6 = p72I.gas_6 .* (ones(101,1)*(ch4ppm_t./ppmvLAY_6(i500,:)));

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
  plot(1:length(p72x.rlat),p72x.mmw)
  plot(1:length(p72x.rlat),p72x.stemp)
  plot(1:length(p72x.rlat),p72x.ptemp(80,:))
  plot(1:length(p72x.rlat),p72x.gas_1(80,:))
  
  ppmvLAY_1 = layers2ppmv(h72x,p72x,1:length(p72x.stemp),1);
  ppmvLAY_2 = layers2ppmv(h72x,p72x,1:length(p72x.stemp),2);
  ppmvLAY_3 = layers2ppmv(h72x,p72x,1:length(p72x.stemp),3);
  ppmvLAY_4 = layers2ppmv(h72x,p72x,1:length(p72x.stemp),4);
  ppmvLAY_5 = layers2ppmv(h72x,p72x,1:length(p72x.stemp),5);
  ppmvLAY_6 = layers2ppmv(h72x,p72x,1:length(p72x.stemp),6);
  
  tcalc = rad2bt(h72x.vchan,p72x.rcalc);
  tcalcavg = squeeze(nanmean(tcalc,2));

  plot(yy2002,tcalc(1520,:),yy2002,tcalcavg(1520)*ones(size(p72.stemp)))  

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %{
    days = (1:numtimesteps)*30/365;

    polyfit(days,tcalc(1520,:),1); ans(1)
    Math_tsfit_lin_robust(days*365,tcalc(1520,:),4); ans(2)
  
    stempjunk = p72x.stemp;
    polyfit(days,stempjunk,1); ans(1)
    Math_tsfit_lin_robust(days*365,stempjunk,4); ans(2)

  %%%%%%%%%%%%%%%%%%%%%%%%%
    dayOFtime = change2days(yy,mm,dd,2002);
    disp('dude I just computed dayOFtime')

    polyfit(dayOFtime/365.25,tcalc(1520,:),1); ans(1)
    Math_tsfit_lin_robust(dayOFtime,tcalc(1520,:),4); ans(2)
  
    stempjunk = p72x.stemp;
    polyfit(dayOFtime/365.25,stempjunk,1); ans(1)
    Math_tsfit_lin_robust(dayOFtime,stempjunk,4); ans(2)
  %}
  %%%%%%%%%%%%%%%%%%%%%%%%%
end

error('lklklk')
%{
stempjunk = p72x.stemp;
if iRampCO2_CH4_N2O > 0
  saver = ['save test' num2str(iLatBin,'%02i') '.mat'];
else
  saver = ['save test' num2str(iLatBin,'%02i') '_tracegasconst.mat'];
end
saver = [saver ' yy mm dd tcalc h72x p72x chisqr stempjunk ind nwp_trends fKc sartatrend iLonBin raaReconstruct fop frp  ha72I pa72I jac iaNlays iLonBin iRampCO2_CH4_N2O'];
eval(saver);
%}

%%% now do geophysical trends and jacs for all 12x19 years
%%% make_single_tile_timeseries_jacs_reconstruct

%%% now do geophysical trends and jacs for all 12x19x4 years/seasons (spring summer fall winter)
%%% this is more general purpose, and uses log jacs for gases d(BT)/log(1+dQ)
make_single_tile_timeseries_jacs_reconstruct_SSFW_ln

%%% now do geophysical trends and jacs for all 12x19x4 years/seasons (spring summer fall winter)
%%% this is more general purpose, and uses abs jacs d(BT)/dQ)
make_single_tile_timeseries_jacs_reconstruct_SSFW_noln

