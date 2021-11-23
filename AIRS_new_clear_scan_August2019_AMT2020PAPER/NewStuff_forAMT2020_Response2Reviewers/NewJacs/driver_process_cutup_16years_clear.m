addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /home/sergio/MATLABCODE/ESRL_TRACE_GAS

iCount = 0;

results.rlon = [];
results.rlat = [];
results.rtime = [];
results.QA = [];

results.btobs = [];
results.btcal0 = [];
results.btcalF = [];

results.stemp_umbc = [];
results.colwv_umbc = [];
results.o3_umbc    = [];
results.co2_umbc   = [];
results.ch4_umbc   = [];
results.tz_umbc    = [];
results.wv_umbc    = [];

results.stemp_era = [];
results.colwv_era = [];
results.o3_era    = [];
results.co2_esrl  = [];
results.ch4_esrl  = [];
results.tz_era    = [];
results.wv_era    = [];

for ii = 1 : 200
  fin = ['/asl/s1/sergio/rtp/singlefootprintretrievals_airicrad_v6/2002/08/31/cutup_16years_clear_' num2str(ii,'%03d') '_5_iDET_4_iStemp_ColWV_5_iCenterFov_-1_retrclr_no_clouds.mat'];
  if exist(fin)    
    iCount = iCount + 1;
    iaFound(iCount) = ii;
    fprintf(1,'%3i : %3i %s \n',ii,iCount,fin)

    clear poemNew poem
    a = load(fin);

    poemNew = a.poemNew;
    hoemNew = a.hoemNew;
    g       = a.g;
    [junk2] = layers2ppmv(hoemNew,poemNew,1:length(poemNew.stemp),2);
    [junk6] = layers2ppmv(hoemNew,poemNew,1:length(poemNew.stemp),6);

    poem = poemNew;
    poem.gas_1 = poem.gas_1_orig;
    poem.gas_2 = poem.gas_2_orig;
    poem.gas_3 = poem.gas_3_orig;
    %poem.gas_4 = poem.gas_4_orig;
    poem.gas_5 = poem.gas_5_orig;
    poem.gas_6 = poem.gas_6_orig;
    poem.ptemp = poem.ptemp_orig;
    poem.stemp = poem.stemp_orig;

    co2ppm = read_trace_gas(double(poemNew.rlat),double(poemNew.rtime),2);
    ch4ppm = read_trace_gas(double(poemNew.rlat),double(poemNew.rtime),6)/1000;

    mmU = mmwater_rtp(hoemNew,poemNew);
    mmE = mmwater_rtp(hoemNew,poem);
    o3U = dobson_rtp(hoemNew,poemNew);
    o3E = dobson_rtp(hoemNew,poem);

    if ~isfield(hoemNew,'vchan')
      hoemNew.vchan = instr_chans;
    end
    i735 = find(hoemNew.vchan >= 735,1);
    i900 = find(hoemNew.vchan >= 900,1);
    i1231 = find(hoemNew.vchan >= 1231,1);
    i1419 = find(hoemNew.vchan >= 1419,1);
    ix = [i735 i900 i1231 i1419];

    j850 = find(poemNew.plevs(:,1) >= 850,1);
    j500 = find(poemNew.plevs(:,1) >= 500,1);
    j250 = find(poemNew.plevs(:,1) >= 250,1);
    jx = [j850 j500 j250];

    results.rlon = [results.rlon poemNew.rlon];
    results.rlat = [results.rlat poemNew.rlat];
    results.rtime = [results.rtime poemNew.rtime];
    results.QA    = [results.QA poemNew.QA];
    results.btobs  = [results.btobs  rad2bt(hoemNew.vchan(ix),poemNew.robs1(ix,:))];
    results.btcal0 = [results.btcal0 rad2bt(hoemNew.vchan(ix),poemNew.rcalc_prof0(ix,:))];
    results.btcalF = [results.btcalF rad2bt(hoemNew.vchan(ix),poemNew.rcalc(ix,:))];
    
    results.tz_umbc    = [results.tz_umbc     poemNew.ptemp(jx,:)];
    results.wv_umbc    = [results.wv_umbc     poemNew.gas_1(jx,:)];
    results.stemp_umbc = [results.stemp_umbc  poemNew.stemp];
    results.colwv_umbc = [results.colwv_umbc  mmU];
    results.o3_umbc    = [results.o3_umbc     o3U];
    results.co2_umbc   = [results.co2_umbc    junk2(80,:)];
    results.ch4_umbc   = [results.ch4_umbc    junk6(80,:)];

    results.tz_era    = [results.tz_era     poem.ptemp(jx,:)];
    results.wv_era    = [results.wv_era     poem.gas_1(jx,:)];
    results.stemp_era = [results.stemp_era  poem.stemp];
    results.colwv_era = [results.colwv_era  mmE];
    results.o3_era    = [results.o3_era     o3E];
    results.co2_esrl  = [results.co2_esrl  co2ppm];
    results.ch4_esrl  = [results.ch4_esrl  ch4ppm];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[yy,mm,dd,hh] = tai2utcSergio(results.rtime);
doy2002 = change2days(yy,mm,dd,2002);
tdoy = 2002 + doy2002/365;

figure(1); plot(tdoy,smooth(results.stemp_era-results.stemp_umbc,1000)); title('ERA-UMBC stemp (K)');
figure(2); plot(tdoy,smooth(results.colwv_era-results.colwv_umbc,1000)); title('ERA-UMBC colwv (mmWV)');
figure(3); plot(tdoy,smooth(results.o3_era-results.o3_umbc,1000)); title('ERA-UMBC o3 (du) ');
figure(4); plot(tdoy,smooth(results.co2_esrl-results.co2_umbc,1000)); title('ESRL-UMBC co2 (ppm)');
figure(5); plot(tdoy,smooth(results.ch4_esrl-results.ch4_umbc,1000)); title('ESRL-UMBC ch4 (ppm)');

latbins = equal_area_spherical_bands(20);

latbins = -70 : 10 : +70;
clear latbins_results
for ii = 1 : length(latbins)-1
  boo = find(results.rlat >= latbins(ii) & results.rlat < latbins(ii+1) & results.QA > 100);
  boo = find(results.rlat >= latbins(ii) & results.rlat < latbins(ii+1) & results.QA == 127);
  if length(boo) > 1

    [xyy,xmm,xdd,xhh] = tai2utcSergio(results.rtime(boo));
    xdoy2002 = change2days(xyy,xmm,xdd,2002);
    xtdoy = 2002 + xdoy2002/365;

    latbins_results(ii).doy   = xtdoy;
    latbins_results(ii).rtime = results.rtime(boo);
    latbins_results(ii).rlat = results.rlat(boo);
    latbins_results(ii).rlon = results.rlon(boo);
    
    latbins_results(ii).btobs  = results.btobs(:,boo);
    latbins_results(ii).btcal0 = results.btcal0(:,boo);
    latbins_results(ii).btcalF = results.btcalF(:,boo);

    latbins_results(ii).tz_umbc  = results.tz_umbc(:,boo);
    latbins_results(ii).wv_umbc  = results.wv_umbc(:,boo);
    latbins_results(ii).stemp_umbc = results.stemp_umbc(boo);
    latbins_results(ii).colwv_umbc = results.colwv_umbc(boo);
    latbins_results(ii).o3_umbc = results.o3_umbc(boo);
    latbins_results(ii).co2_umbc = results.co2_umbc(boo);
    latbins_results(ii).ch4_umbc = results.ch4_umbc(boo);

    latbins_results(ii).tz_era  = results.tz_era(:,boo);
    latbins_results(ii).wv_era  = results.wv_era(:,boo);
    latbins_results(ii).stemp_era = results.stemp_era(boo);
    latbins_results(ii).colwv_era = results.colwv_era(boo);
    latbins_results(ii).o3_era = results.o3_era(boo);
    latbins_results(ii).co2_era = results.co2_esrl(boo);
    latbins_results(ii).ch4_era = results.ch4_esrl(boo);
  end
end

figure(1); plot(tdoy,smooth(results.stemp_era-results.stemp_umbc,1000)); title('ERA-UMBC stemp (K)');

%%%%%%%%%%%%%%%%%%%%%%%%%

plot(latbins_results(02).doy,latbins_results(02).co2_era,'b',latbins_results(02).doy,smooth(latbins_results(02).co2_umbc,100),'c',....
     latbins_results(07).doy,latbins_results(07).co2_era,'k',latbins_results(07).doy,smooth(latbins_results(07).co2_umbc,100),'g',....
     latbins_results(12).doy,latbins_results(12).co2_era,'r',latbins_results(12).doy,smooth(latbins_results(12).co2_umbc,100),'m','linewidth',2)

[BCO2_era12, stats]=Math_tsfit_lin_robust((latbins_results(12).doy-2002)*365,latbins_results(12).co2_era,4); 
[BCO2_umbc12, stats]=Math_tsfit_lin_robust((latbins_results(12).doy-2002)*365,latbins_results(12).co2_umbc,4); 

%%%%%%%%%%%%%%%%%%%%%%%%%

plot(latbins_results(02).doy,latbins_results(02).ch4_era,'b',latbins_results(02).doy,smooth(latbins_results(02).ch4_umbc,100),'c',....
     latbins_results(07).doy,latbins_results(07).ch4_era,'k',latbins_results(07).doy,smooth(latbins_results(07).ch4_umbc,100),'g',....
     latbins_results(12).doy,latbins_results(12).ch4_era,'r',latbins_results(12).doy,smooth(latbins_results(12).ch4_umbc,100),'m','linewidth',2)

[BCH4_era12, stats]=Math_tsfit_lin_robust((latbins_results(12).doy-2002)*365,latbins_results(12).ch4_era,4); 
[BCH4_umbc12, stats]=Math_tsfit_lin_robust((latbins_results(12).doy-2002)*365,latbins_results(12).ch4_umbc,4); 

%%%%%%%%%%%%%%%%%%%%%%%%%

plot(latbins_results(02).doy,latbins_results(02).o3_era,'b',latbins_results(02).doy,smooth(latbins_results(02).o3_umbc,100),'c',....
     latbins_results(07).doy,latbins_results(07).o3_era,'k',latbins_results(07).doy,smooth(latbins_results(07).o3_umbc,100),'g',....
     latbins_results(12).doy,latbins_results(12).o3_era,'r',latbins_results(12).doy,smooth(latbins_results(12).o3_umbc,100),'m','linewidth',2)

[BO3_era12, stats]=Math_tsfit_lin_robust((latbins_results(12).doy-2002)*365,latbins_results(12).o3_era,4); 
[BO3_umbc12, stats]=Math_tsfit_lin_robust((latbins_results(12).doy-2002)*365,latbins_results(12).o3_umbc,4); 

%%%%%%%%%%%%%%%%%%%%%%%%%

plot(latbins_results(02).doy,latbins_results(02).stemp_era,'b',latbins_results(02).doy,smooth(latbins_results(02).stemp_umbc,100),'c',....
     latbins_results(07).doy,latbins_results(07).stemp_era,'k',latbins_results(07).doy,smooth(latbins_results(07).stemp_umbc,100),'g',....
     latbins_results(12).doy,latbins_results(12).stemp_era,'r',latbins_results(12).doy,smooth(latbins_results(12).stemp_umbc,100),'m','linewidth',2)

[BSTEMP_era12, stats]=Math_tsfit_lin_robust((latbins_results(12).doy-2002)*365,latbins_results(12).stemp_era,4); 
[BSTEMP_umbc12, stats]=Math_tsfit_lin_robust((latbins_results(12).doy-2002)*365,latbins_results(12).stemp_umbc,4); 

%%%%%%%%%%%%%%%%%%%%%%%%%

plot(latbins_results(02).doy,latbins_results(02).colwv_era,'b',latbins_results(02).doy,smooth(latbins_results(02).colwv_umbc,100),'c',....
     latbins_results(07).doy,latbins_results(07).colwv_era,'k',latbins_results(07).doy,smooth(latbins_results(07).colwv_umbc,100),'g',....
     latbins_results(12).doy,latbins_results(12).colwv_era,'r',latbins_results(12).doy,smooth(latbins_results(12).colwv_umbc,100),'m','linewidth',2)

[BCOLwv_era12, stats]=Math_tsfit_lin_robust((latbins_results(12).doy-2002)*365,latbins_results(12).colwv_era,4); 
[BCOLwv_umbc12, stats]=Math_tsfit_lin_robust((latbins_results(12).doy-2002)*365,latbins_results(12).colwv_umbc,4); 
