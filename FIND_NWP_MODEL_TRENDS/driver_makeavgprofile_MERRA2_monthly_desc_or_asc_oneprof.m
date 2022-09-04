addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /asl/matlib/aslutil/

%% see eg ~/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_AIRS_STM_Oct2020_allstarts_Jan20XY/call_save_split_apart_rtp_howard_bins.m
[hERAI,ha,pERAI,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/summary_17years_all_lat_all_lon_2002_2019.rtp');

topts.klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
topts.sarta = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';

if ~exist('all')
  %load ERA5_atm_data_2002_09_to_2021_08_desc.mat
  load MERRA2_atm_data_2002_09_to_2021_08_desc.mat
end

%% these are actually DESC

ii0 = 2276;
%4133        4134        4204        4205        4274        4276        4277        4346        4349        4416        4417        4418        4421        4490        4491        4492
ii0 = 4492
%          32         261         266        4134        4419        4420
ii0 = 4420
ii0 = 4378

foutNyearaverageIP = ['summary_19years_all_lat_all_lon_2002_2021_monthlyMERRA2_' num2str(ii0,'%4i') '.ip.rtp'];
foutNyearaverageOP = ['summary_19years_all_lat_all_lon_2002_2021_monthlyMERRA2_' num2str(ii0,'%4i') '.op.rtp'];
foutNyearaverageRP = ['summary_19years_all_lat_all_lon_2002_2021_monthlyMERRA2_' num2str(ii0,'%4i') '.rp.rtp'];

if ~exist(foutNyearaverageIP)
  pNyearaverage = [];
  for ii = ii0
    yavg.stemp = nanmean(all.stemp(:,ii));
    yavg.ptemp(:,1) = nanmean(squeeze(all.nwp_ptemp(:,:,ii)),1);
    yavg.gas_1(:,1) = nanmean(squeeze(all.nwp_gas_1(:,:,ii)),1);
    yavg.gas_3(:,1) = nanmean(squeeze(all.nwp_gas_3(:,:,ii)),1);
    yavg.plevs(:,1) = nanmean(squeeze(all.nwp_plevs(:,:,ii)),1);

    yavg.rlon  = all.rlon(ii);
    yavg.rlat  = all.rlat(ii);
  end

  yavg.spres = pERAI.spres(ii);
  yavg.solzen = pERAI.solzen(ii);
  yavg.satzen = pERAI.satzen(ii);
  yavg.rtime = pERAI.rtime(ii);
  yavg.zobs = pERAI.zobs(ii);
  yavg.salti = pERAI.salti(ii);
  yavg.scanang = pERAI.scanang(ii);
  yavg.nemis = pERAI.nemis(ii);
  yavg.emis = pERAI.emis(:,ii);
  yavg.efreq = pERAI.efreq(:,ii);
  yavg.rho = pERAI.rho(:,ii)

  iFix = +1;
  if iFix > 0
    iN = find(yavg.plevs <= yavg.spres);
    iN = [iN; iN(end)+1];
    bad = find(yavg.ptemp(iN) < 0 | yavg.gas_1(iN) < 0 | yavg.gas_3(iN) < 0);
    bad = find(yavg.ptemp(iN) < 150);
    if length(bad) > 0
      good = setdiff(iN,bad);
      yavg.ptemp(bad) = interp1(log(yavg.plevs(good)),yavg.ptemp(good),log(yavg.plevs(bad)),[],'extrap');
    end
    bad = find(yavg.gas_1(iN) < 0);
    if length(bad) > 0
      good = setdiff(iN,bad);
      yavg.gas_1(bad) = interp1(log(yavg.plevs(good)),log(yavg.gas_1(good)),log(yavg.plevs(bad)),[],'extrap');
      yavg.gas_1(bad) = exp(yavg.gas_1(bad));
    end
    bad = find(yavg.gas_3(iN) < 0);
    if length(bad) > 0
      good = setdiff(iN,bad);
      yavg.gas_3(bad) = interp1(log(yavg.plevs(good)),log(yavg.gas_3(good)),log(yavg.plevs(bad)),[],'extrap');
      yavg.gas_3(bad) = exp(yavg.gas_3(bad));
    end
  end

  [yy,mm,dd,hh] = tai2utcSergio(yavg.rtime);
  tS = utc2taiSergio(2002,09,01,00);
  tE = utc2taiSergio(2021,08,31,23.99);
  [yyX,mmX,ddX,hhX] = tai2utcSergio((tS+tE)/2);
  yavg.rtime = ones(size(yavg.rtime)) * (tS+tE)/2;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  yavg.co2ppm = yyX + (mmX-1)/12; yavg.co2ppm = ones(size(yavg.rtime)) * (yavg.co2ppm-2002)*2.2 + 370;
  
    addpath /home/sergio/MATLABCODE/ESRL_TRACE_GAS
    disp('putting in CO2,N2O.CH4 ESRL')
    co2ppm = read_trace_gas(double(yavg.rlat),double(yavg.rtime),2);
    ch4ppm = read_trace_gas(double(yavg.rlat),double(yavg.rtime),6)/1000;
  
    n2oslope = (332-315)/(2020-2000); %% https://www.esrl.noaa.gov/gmd/hats/combined/N2O.html
    [xyy,xmm,xdd,xhh] = tai2utcSergio(yavg.rtime);        %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<< for SdSM old time
    time_so_far = (xyy-2000) + ((xmm-1)+1)/12;
    n2oppm = 315 + n2oslope*time_so_far; %% ppb
    n2oppm = n2oppm/1000;             %% change to ppm
  
    yavg.co2ppm = co2ppm;
    yavg.n2oppm = n2oppm;
    yavg.ch4ppm = ch4ppm;
  
  scatter_coast(yavg.rlon,yavg.rlat,50,yavg.co2ppm)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  havg = hERAI;
  havg.ngas = 2;
  havg.glist = [1  3]';
  havg.gunit = [21 21]';
  havg.ptype = 0;
  havg.pfields = 1;
  yavg.plat = yavg.rlat;
  yavg.plon = yavg.rlon;
  %yavg.nlevs = 37 * ones(size(yavg.stemp));
  boo = find(yavg.plevs <= yavg.spres);
  yavg.nlevs = max(boo);
  rtpwrite(foutNyearaverageIP,havg,ha,yavg,pa);
  klayerser = ['!' topts.klayers ' fin=' foutNyearaverageIP ' fout=' foutNyearaverageOP]; eval(klayerser);

  [hnew,hax,pnew,pax] = rtpread(foutNyearaverageOP);
  addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
  ppmv2 = layers2ppmv(hnew,pnew,1:length(pnew.stemp),2);
  ppmv4 = layers2ppmv(hnew,pnew,1:length(pnew.stemp),4);
  ppmv6 = layers2ppmv(hnew,pnew,1:length(pnew.stemp),6);

  i500mb = find(pnew.plevs >= 500,1)
  i850mb = find(pnew.plevs >= 850,1)

  ppmv2(i500mb)
  ppmv2(i500mb)-co2ppm

  pnew.gas_6 = pnew.gas_6 .* (ones(101,1) * ch4ppm./ppmv6(i500mb));
  pnew.gas_4 = pnew.gas_4 .* (ones(101,1) * n2oppm./ppmv4(i500mb));
  ppmv4 = layers2ppmv(hnew,pnew,1:length(pnew.stemp),4);
  ppmv6 = layers2ppmv(hnew,pnew,1:length(pnew.stemp),6);

 rtpwrite(foutNyearaverageOP,hnew,hax,pnew,pax)
 
  sartaer = ['!' topts.sarta ' fin=' foutNyearaverageOP '     fout=' foutNyearaverageRP]; eval(sartaer);

  [hnew,hax,pnew,pax] = rtpread(foutNyearaverageRP);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% compare pnew profiles with pERAI
  addpath /home/sergio/MATLABCODE/COLORMAP/

end

disp('now make kCARTA jacs and look at eg /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly/clust_put_together_jacs_clrMERRA2.m')
disp('now make kCARTA jacs and look at eg /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly/clust_put_together_jacs_clrMERRA2.m')
disp('now make kCARTA jacs and look at eg /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly/clust_put_together_jacs_clrMERRA2.m')

figure(1); plot(yavg.ptemp,yavg.plevs); set(gca,'yscale','log'); set(gca,'ydir','reverse'); xlim([130 300])
figure(2); plot(pnew.ptemp,pnew.plevs); set(gca,'yscale','log'); set(gca,'ydir','reverse'); xlim([130 300])
figure(3); plot(yavg.gas_1,yavg.plevs); set(gca,'yscale','log'); set(gca,'ydir','reverse'); xlim([0 1e-3])
figure(4); plot(pnew.gas_1,pnew.plevs); set(gca,'yscale','log'); set(gca,'ydir','reverse'); xlim([0 1e20])
figure(5); plot(hnew.vchan,rad2bt(hnew.vchan,pnew.rcalc))
