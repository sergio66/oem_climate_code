%% monthly, 18 years x 12 months/year = 216
%% monthly, 19 years x 12 months/year = 228
%% monthly, 12 years x 12 months/year = 144

addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/Strow_humidity/convert_humidity/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/IDL_WV_ROUTINES/atmos_phys/MATLAB/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/create_ecrad_inputSergio/
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil

load('llsmap5.mat');

%% note this code only handles 2002_09 to YYYY_08 sp
%%      this code does not currently handle eg OCO2 ERA5_atm_N_cld_data_2012_05_to_2019_04_trends_desc.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('assumes you have run driver_compute_AIRSL3_trends_desc_or_asc.m');

[h0,ha,p0,pa] = rtpread('summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.ip.rtp');

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 240 for months
if length(JOB) == 0
  JOB = 97;
  JOB = 78;
end

epss = eps;
epss = eps *1e10;
epss = 1e-10;

iNumYears = 20;
timespan  = iNumYears;
iDorA = +1;

yy0 = 2002; mm0 = 8;
yy = yy0; mm = mm0;
for JOBB = 1 : JOB
  yy = yy; mm = mm + 1;
  if mod(JOBB-5,12) == 0
    yy = yy + 1;
    mm = 1;
  end
  saveYY(JOBB) = yy; saveMM(JOBB) = mm;
  fprintf(1,'JOBB = %3i %4i/%2i \n',JOBB,yy,mm)
end

%for JOB = 1 : iNumYears*12

  fnameIN = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/AIRSL3/Tile_Center/DESC/airsL3_tile_center_monthly_' num2str(JOB,'%03d') '.mat'];

  fnameOUT = ['OLR_ecRad/AIRSL3v0/airsL3_olr_' num2str(JOB,'%03d') '.mat'];
  if ~exist(fnameOUT)    

    if timespan <= 20
      savestr_version_big = 'Sept2002_Aug2021_19yr_';
      savestr_version_big = 'Sept2002_Aug2022_20yr_';
    end
  
    if iDorA > 0
      fnameIN = ['/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_' savestr_version_big 'desc.mat'];
      loader = ['load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_' savestr_version_big 'desc.mat'];
    else
      fnameIN = ['/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_' savestr_version_big 'asc.mat'];
      loader = ['load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_' savestr_version_big 'asc.mat'];
    end

    a = load(fnameIN);

    clear fnameIN

    xpnew_ip = struct;
    xpnew_ip.olrclr = squeeze(a.save64x72_clrolr(:,:,JOB));
    xpnew_ip.ptemp = squeeze(a.save64x72_T(:,:,:,JOB));
    xpnew_ip.stemp = squeeze(a.save64x72_stemp(:,:,JOB));
    xpnew_ip.gas_1 = squeeze(a.save64x72_Q(:,:,:,JOB));
    for ii = 13 : 24
      xpnew_ip.gas_1(:,:,ii) = xpnew_ip.gas_1(:,:,12);
    end
    xpnew_ip.gas_3 = squeeze(a.save64x72_O3(:,:,:,JOB));
    xpnew_ip.gas_5 = squeeze(a.save64x72_CO(:,:,:,JOB));
    xpnew_ip.gas_6 = squeeze(a.save64x72_CH4(:,:,:,JOB));
    xpnew_ip.plevs = (a.Tlevs)' * ones(1,4608);
    [xpnew_ip.rlat,xpnew_ip.rlon]  = meshgrid(a.save_lat64x72',a.save_lon64x72);
    xpnew_ip.rlat = xpnew_ip.rlat'; xpnew_ip.rlon = xpnew_ip.rlon';

%    [badX,badY]  = find(isnan(xpnew_ip.stemp));
%    [goodX,goodY] = find(isfinite(xpnew_ip.stemp));
%    %Vq = interp2(xpnew_ip.rlon(goodX,goodY),xpnew_ip.rlat(goodX,goodY),xpnew_ip.stemp(goodX,goodY),xpnew_ip.rlon(badX,badY),xpnew_ip.rlat(badX,badY));
%    [bad]  = find(isnan(xpnew_ip.stemp));
%    xpnew_ip.stemp(bad) = xpnew_ip.ptemp(1,bad);

    pnew_ip.plevs = xpnew_ip.plevs;
    pnew_ip.rlat = reshape(xpnew_ip.rlat,1,4608);
    pnew_ip.rlon = reshape(xpnew_ip.rlon,1,4608);
    pnew_ip.olrclr = reshape(xpnew_ip.olrclr,1,4608);
    pnew_ip.stemp = reshape(xpnew_ip.stemp,1,4608);
    for ii = 1 : 24
      junk = squeeze(xpnew_ip.ptemp(:,:,ii)); pnew_ip.ptemp(ii,:) = reshape(junk,1,4608);
      junk = squeeze(xpnew_ip.gas_1(:,:,ii)); pnew_ip.gas_1(ii,:) = reshape(junk,1,4608);
      junk = squeeze(xpnew_ip.gas_3(:,:,ii)); pnew_ip.gas_3(ii,:) = reshape(junk,1,4608);
      junk = squeeze(xpnew_ip.gas_6(:,:,ii)); pnew_ip.gas_6(ii,:) = reshape(junk,1,4608);
    end

    [badX,badY]  = find(isnan(pnew_ip.stemp));
    [goodX,goodY] = find(isfinite(pnew_ip.stemp));
    [bad]  = find(isnan(pnew_ip.stemp));
    pnew_ip.stemp(bad) = pnew_ip.ptemp(1,bad);
    
    pnew_ip00 = pnew_ip;

    bad = find(isnan(pnew_ip.ptemp)); pnew_ip.ptemp(bad) = 0;
    bad = find(isnan(pnew_ip.gas_1)); pnew_ip.gas_1(bad) = 0;
    bad = find(isnan(pnew_ip.gas_3)); pnew_ip.gas_3(bad) = 0;
    bad = find(isnan(pnew_ip.gas_6)); pnew_ip.gas_6(bad) = 0;

    scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.stemp); colormap jet
    scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.olrclr); colormap jet

    pnew_ip.rtime = utc2taiSergio(saveYY(JOB),saveMM(JOB),15,01.50) * ones(size(pnew_ip.rlat));
    pnew_ip.nemis = p0.nemis;
    pnew_ip.emis  = p0.emis;
    pnew_ip.rho   = p0.rho;
    pnew_ip.efreq = p0.efreq;
    pnew_ip.spres = p0.spres;
    pnew_ip.zobs  = p0.zobs;
    pnew_ip.salti = p0.salti;
    %pnew_ip.landfrac = p0.landfrac;
    pnew_ip.scanang = p0.scanang;
    pnew_ip.satzen = p0.scanang;
    pnew_ip.solzen = p0.solzen;
    pnew_ip.nlevs = 24 * ones(size(pnew_ip.stemp));

    pnew_ip.plevs = flipud(pnew_ip.plevs);
    pnew_ip.ptemp = flipud(pnew_ip.ptemp);
    pnew_ip.gas_1 = flipud(pnew_ip.gas_1);
    pnew_ip.gas_3 = flipud(pnew_ip.gas_3);
    pnew_ip.gas_6 = flipud(pnew_ip.gas_6);

    pnew_ip0 = pnew_ip;

    for ii = 1 : 4608
      plevs = pnew_ip.plevs(:,ii);
      spres = pnew_ip.spres(ii);
      boo = find(plevs <= spres);
      nlevs = min(max(boo) + 1,24);
      pnew_ip.nlevs(ii) = nlevs;
      tz = pnew_ip.ptemp(:,ii);
      wz = pnew_ip.gas_1(:,ii);
      oz = pnew_ip.gas_3(:,ii);

      bad = find(tz(1:nlevs) < 180); good = find(tz(1:nlevs) >= 180);
      if length(bad) > 0
        tz(bad) = interp1(log(plevs(good)),tz(good),log(plevs(bad)),[],'extrap');
      end
      bad = find(wz <= epss); good = find(wz > epss);
      if length(bad) > 0
        wz(bad) = interp1(log(plevs(good)),wz(good),log(plevs(bad)),[],'extrap');
      end
      bad = find(wz <= epss); good = find(wz > epss);
      if length(bad) > 0
        wz(bad) = wz(good(end))*ones(size(bad));
      end
      bad = find(oz <= epss); good = find(oz > epss);
      if length(good) == 0
        oz = pnew_ip.gas_3(:,ii-1);
        bad = find(oz <= epss); good = find(oz > epss);  
      end
      if length(bad) > 0
        oz(bad) = interp1(log(plevs(good)),oz(good),log(plevs(bad)),[],'extrap');
      end
      pnew_ip.ptemp(:,ii) = tz;
      pnew_ip.gas_1(:,ii) = wz;
      pnew_ip.gas_3(:,ii) = oz;
    end

    %% https://docserver.gesdisc.eosdis.nasa.gov/repository/Mission/AIRS/3.3_ScienceDataProductDocumentation/3.3.4_ProductGenerationAlgorithms/V6_Released_Processing_Files_Description.pdf, pg 222
    hnew_ip.ngas = 4;
    hnew_ip.glist = [1 3 5 6]';
    hnew_ip.gunit = [20 10 10 10]';  %% WV = mass mix ratio in g/kg, O3,CO,CH4 = volume mix ratio in ppv

    hnew_ip.ngas = 2;
    hnew_ip.glist = [1 3]';
    hnew_ip.gunit = [20 10 ]';  %% WV = mass mix ratio in g/kg, O3,CO,CH4 = volume mix ratio in ppv

    hnew_ip.ptype = 0;
    hnew_ip.pfields = 1;
    hnew_ip.nchan = h0.nchan;
    hnew_ip.ichan = h0.ichan;
    hnew_ip.vchan = h0.vchan;

    pnew_ip.plat = pnew_ip.rlat;
    pnew_ip.plon = pnew_ip.rlon;

   scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,pnew_ip.stemp-pnew_ip.ptemp(1,:)); title('Stemp-AIrTemp')
   colormap(usa2); caxis([-1 +1]*5)

    fip = ['junk.ip.rtp_' num2str(JOB)];
    fop = ['junk.op.rtp_' num2str(JOB)];
    frp = ['junk.rp.rtp_' num2str(JOB)];
    rtpwrite(fip,hnew_ip,[],pnew_ip,[]);

 %  run through klayers to get hnew_op,pnew_op
 %  run through sarta   to get hnew_op,pnew_op with rads

    klayerser = ['!/asl/packages/klayersV205/BinV201/klayers_airs fin=' fip ' fout=' fop ' >& ugh'];
    eval(klayerser);
    sarta = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';
    sarta = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';
    sartaer = ['!' sarta ' fin=' fop ' fout=' frp];
    eval(sartaer);

    [hnew_op,~,pnew_op,~] = rtpread(frp);
    scatter_coast(pnew_ip.rlon,pnew_ip.rlat,50,rad2bt(1231,pnew_op.rcalc(1520,:))); colormap jet
    
    h = hnew_op;
    px = pnew_op;
    bad = find(isnan(px.ptemp) | isnan(px.gas_1) | isnan(px.gas_2) | isnan(px.gas_3));
    fprintf(1,'found %8i of %8i NAN profile entries \n',length(bad),101*length(px.stemp));
    px.ptemp(bad) = 0;
    px.gas_1(bad) = 0;
    px.gas_2(bad) = 0;
    px.gas_3(bad) = 0;
    bad = find(isnan(px.efreq) | isnan(px.emis) | isnan(px.rho));
    px.efreq(bad) = 0;
    px.emis(bad) = 0;
    px.rho(bad) = 0;

    olr = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);
    
    stemp = px.stemp;
    olr_AIRSL3 = pnew_ip.olrclr;
    saver = ['save ' fnameOUT ' stemp olr  olr_AIRSL3'];
    eval(saver)
    rmer = ['!/bin/rm ' fip ' ' fop ' ' frp];
    eval(rmer);
    fprintf(1,'saved %s \n',fnameOUT)

  end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computeAIRSL3_OLR_trend
