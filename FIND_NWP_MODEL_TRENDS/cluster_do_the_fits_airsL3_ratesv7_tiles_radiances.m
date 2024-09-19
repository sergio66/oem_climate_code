addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/emis/
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/

%% iaFound = check_all_jobs_done('/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_latbin_',64,'.mat');
JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 64 for the 64 latbins

if length(JOB) == 0
  JOB = 32;
  JOB = 39;
  JOB = 29;
  JOB = 41;
end

filename = '/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc.mat';
loader = ['load ' filename];
eval(loader)

fnameout = [filename(1:end-4) '_btanom_latbin_' num2str(JOB,'%02i') '.mat'];
%if exist(fnameout)
%  fprintf(1,'%s already exists, skipping \n',fnameout)
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hraw,ha,praw,pa] = rtpread('summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.op.rtp');

%% see also /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_put_together_QuantileChoose_anomalies.m
%% https://www.arm.gov/publications/proceedings/conf05/extended_abs/mlawer_ej.pdf
RRTM_bands0 = [10 250 500 630 700 820 980 1080 1180 1390 1480 1800 2080 2250 2380 2600 3000];
wx = 10:1:3000; rx = ttorad(wx,300); 
fluxx = trapz(wx,rx)/1000;
for ijunk = 2 : length(RRTM_bands0)
  junk = find(wx > RRTM_bands0(ijunk-1) & wx <= RRTM_bands0(ijunk));
  fluxx(ijunk) = trapz(wx(junk),rx(junk))/1000;
end

RRTM_bands = RRTM_bands0(4:end);
thestatsradtrend64x72.RRTM_bands = RRTM_bands;

do_XX_YY_from_X_Y
klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
sarta = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin//jac_airs_l1c_2834_cloudy_may19_prod';

Xprime = X';
Yprime = Y';

[~,~,lenT] = size(save64x72_stemp);
iNumYears = ceil(lenT/12);
iCntJunk = 0;
yy = 2002;
mm = 9;
for iJunk = 1 : lenT  
  yyx(iJunk) = yy;
  mmx(iJunk) = mm;
  mm = mm +1;
  if mm == 13
    mm = 1;
    yy = yy + 1;
  end
end

%iioo = 19;
for jj = JOB
  for ii = 1 : 72
  %for ii = iioo
    p = struct;
    p.spres = praw.spres((jj-1)*72 + ii) * ones(1,lenT);;
    p.gas_6 = squeeze(save64x72_CH4(jj,ii,:,:));
    p.gas_5 = squeeze(save64x72_CO(jj,ii,:,:));
    p.gas_3 = squeeze(save64x72_O3(jj,ii,:,:));

    p.gas_1 = 0 * p.gas_5;
    p.gas_1(01:12,:) = squeeze(save64x72_Q(jj,ii,:,:)); 
    iNlev = 24;
    for ll = 13 : 24
      frac = (iNlev-ll+1)/(iNlev/2+1);
      p.gas_1(ll,:) = p.gas_1(12,:)*frac;
    end
    p.ptemp = squeeze(save64x72_T(jj,ii,:,:));
    p.stemp = squeeze(save64x72_stemp(jj,ii,:,:))';
    p.wspeed = ones(size(p.stemp))*5;
    p.rlon = Xprime(jj,ii) * ones(1,lenT);
    p.rlat = Yprime(jj,ii) * ones(1,lenT);

    p.zobs = 705000 * ones(size(p.stemp));
    p.solzen = 150 * ones(1,lenT);
    p.scanang = 22  * ones(1,lenT);
    p.satzen = vaconv(p.scanang, p.zobs, zeros(size(p.zobs)));

    p.rtime = utc2taiSergio(yyx,mmx,15,12);
    p.co2ppm = 370 + (1:lenT)/12*2.2;
    p.rcalc = zeros(2645,lenT);
    [p.salti,p.landfrac] =  usgs_deg10_dem(p.rlat,p.rlon);
    p.plat = p.rlat;
    p.plon = p.rlon;
    p.plon = p.rlon;
    p.plon = p.rlon;

    p.nlevs = 24 * ones(size(p.stemp));
    p.plevs = Tlevs' * ones(1,lenT);

    %%%%%%%%%%%%%%%%%%%%%%%%%

    avgwv = nanmean(p.gas_1,2);
      bad = find(isnan(avgwv) | isinf(avgwv));
      if length(bad) > 0
        good = find(isfinite(avgwv));
        avgwv(bad) = interp1(log(Tlevs(good)),avgwv(good),log(Tlevs(bad)),[],'extrap');
        avgwv(bad) = abs(avgwv(bad));
      end
    avgoz = nanmean(p.gas_3,2);
      bad = find(isnan(avgoz) | isinf(avgoz));
      if length(bad) > 0
        good = find(isfinite(avgoz));
        avgoz(bad) = interp1(log(Tlevs(good)),avgoz(good),log(Tlevs(bad)),[],'extrap');
        avgoz(bad) = abs(avgoz(bad));
      end
    avgtz = nanmean(p.ptemp,2);
      bad = find(isnan(avgtz) | isinf(avgtz));
      if length(bad) > 0
        good = find(isfinite(avgtz));
        avgtz(bad) = interp1(log(Tlevs(good)),avgtz(good),log(Tlevs(bad)),[],'extrap');
        avgtz(bad) = abs(avgtz(bad));
      end

    for iii = 1 : length(p.stemp)
      plevs = p.plevs(:,iii);
      nlevs = p.nlevs(iii);
      spres = p.spres(iii);
      moo = find(plevs > spres);
      %moo(length(moo)+1) = moo(end)+1;
      moo = moo(1:end-1);

      wv = p.gas_1(:,iii);
      oz = p.gas_3(:,iii);
      tz = p.ptemp(:,iii);

      bad = find(isnan(wv)); wv(bad) = 0; 
        badbad = setdiff(bad,moo); 
        good = max(bad)+1 : 24;
        if length(badbad) > 0  & length(good) > 15
          wv(badbad) = interp1(log(plevs(good)),wv(good),log(plevs(badbad)),[],'extrap');
        elseif length(badbad) > 0 & length(good) <= 15
          wv(badbad) = p.gas_1(badbad,iii-1);
        end
      bad = find(isnan(oz)); oz(bad) = 0;
        badbad = setdiff(bad,moo); 
        good = max(bad)+1 : 24;
        if length(badbad) > 0 & length(good) > 15
          oz(badbad) = interp1(log(plevs(good)),oz(good),log(plevs(badbad)),[],'extrap');
        elseif length(badbad) > 0 & length(good) <= 15          
          oz(badbad) = avgoz(badbad);
        end
      bad = find(isnan(tz)); tz(bad) = 0; 
        badbad = setdiff(bad,moo); 
        good = max(bad)+1 : 24;
        if length(badbad) > 0 & length(good) > 15
          tz(badbad) = interp1(log(plevs(good)),tz(good),log(plevs(badbad)),[],'extrap');
        elseif length(badbad) > 0 & length(good) <= 15           
          tz(badbad) = p.ptemp(badbad,iii-1);
        end

      p.gas_1(:,iii) = wv;
      p.gas_3(:,iii) = oz;
      p.ptemp(:,iii) = tz;

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%

%{
    p.verybad = zeros(size(p.salti));
    iCheck = 12; %% AIRS L3 has 24 levels for T, 12 for WV
    verybad1 = find(isnan(p.ptemp(iCheck,:)) | isnan(p.gas_1(iCheck,:)) | isnan(p.gas_3(iCheck,:)))
    iCheck = 18; %% AIRS L3 has 24 levels for T, 12 for WV, but we have ""augmented" the latter
    verybad2 = find(isnan(p.ptemp(iCheck,:)) | isnan(p.gas_1(iCheck,:)) | isnan(p.gas_3(iCheck,:)))
    verybad = union(verybad1,verybad2);
    p.verybad(verybad) = 1;
    if length(verybad) > 0  
      for jjj = 1 : length(verybad)
        bah = verybad(jjj);
        badah = (1:12*iNumYears);
        [Y,I1,I2] = intersect(bah,badah);
        if I2 > 1
          p.stemp(bah)   = p.stemp(badah(I2-1));
          p.ptemp(:,bah) = p.ptemp(:,badah(I2-1));
          p.gas_1(:,bah) = p.gas_1(:,badah(I2-1));
          p.gas_3(:,bah) = p.gas_3(:,badah(I2-1));
  %        p.rh(:,bah)    = p.rh(:,badah(I2-1));
        else
          p.stemp(bah)   = p.stemp(badah(I2+1));
          p.ptemp(:,bah) = p.ptemp(:,badah(I2+1));
          p.gas_1(:,bah) = p.gas_1(:,badah(I2+1));
          p.gas_3(:,bah) = p.gas_3(:,badah(I2+1));
  %        p.rh(:,bah)    = p.rh(:,badah(I2+1));
        end
      end    
    end
%}

    %%%%%%%%%%%%%%%%%%%%%%%%%
    h.ptype = 0;
    h.pfields = 5; % (1=prof + 4=IRobs);
    h.pfields = 1; %% 1 = profile
    h.nchan = 2645;

    h.ichan = (1:2645)';
    h.vchan = instr_chans2645;

    junk = load('/home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/h2645structure.mat');
    h.ichan = junk.h.ichan;
    h.vchan = junk.h.vchan;

    %% see driver_check_WV_T_RH_AIRSCLIMCAPSL3_geo_and_spectral_rates2.m
    h.ngas = 4;    
    h.glist = [1  3  5  6]';
    h.gunit = [21 12 12 12]'; %% kg/kg and VMR
    %% see driver_check_WV_T_RH_AIRSL3_geo_and_spectral_rates2.m
    h.ngas = 2;
    h.gunit = [20 12]';  %% g/kg and VMR
    h.glist = [ 1 3 ]';

    pa = {{'profiles','rtime','seconds since 1993'}};
    ha = {{'header','hdf file',filename}};
    [p,pa] = rtp_add_emis(p,pa);

    fip = ['junk_' num2str(ii,'%02i') '_' num2str(jj,'%02i')  '.ip.rtp'];
    fop = ['junk_' num2str(ii,'%02i') '_' num2str(jj,'%02i')  '.op.rtp'];
    frp = ['junk_' num2str(ii,'%02i') '_' num2str(jj,'%02i')  '.rp.rtp'];

    rtpwrite(fip,h,ha,p,pa);
    klayerser = ['!' klayers ' fin=' fip ' fout=' fop ' >& ughklayers'];
    sartaer   = ['!' sarta '   fin=' fop ' fout=' frp ' >& ughsarta'];
    fprintf(1,'lonbin %2i of 72 : did sarta and klayers \n',ii)
    eval(klayerser);
    eval(sartaer);
    [hcalc,ha,pcalc,pa] = rtpread(frp);
    rmer = ['!rm ' fip ' ' fop ' ' frp]; eval(rmer)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rcalc = pcalc.rcalc;
    tcalc = rad2bt(hcalc.vchan,pcalc.rcalc);

    iNumSineCosCycles = 4;
    fairs = h.vchan;
    k = 1 : lenT;
    [yyx,mmx,ddx,hhx] = tai2utcSergio(pcalc.rtime);
    dtime = change2days(yyx,mmx,ddx,2002);        

    clear moo woo
    [Y,I] = sort(h.vchan);
    fchan = Y;
    FAFA = pcalc.rcalc(I,:);
    bonk = find(diff(fchan(I)) > 10); bonk = [bonk bonk+1];
    data = trapz(fchan,FAFA)/1000 - trapz(fchan(bonk),FAFA(bonk,:))/1000;  moo(1,:) = data;
    chuse =[];
    for flfl = 1 : length(RRTM_bands)-1
      junk = find(fchan >= RRTM_bands(flfl) & fchan < RRTM_bands(flfl+1));
      chuse = [chuse; junk];
      data = trapz(fchan(junk),FAFA(junk,:))/1000; moo(flfl+1,:) = data;
    end
    woo = sum(moo(2:14,:),1); 
    figure(1); clf; plot(1:262,woo,1:262,moo(1,:)); title('Compare flux(all) vs sum(bandFlux)')
    figure(2); clf; plot(nanmean(moo,2)); title('Mean flux per RRTM band')
    figure(3); clf; plot(fchan,nanmean(FAFA,2),'b.-',h.vchan,nanmean(pcalc.rcalc,2))
    figure(4); clf; plot(fchan,rad2bt(fchan,nanmean(FAFA,2)),'b.-',h.vchan,rad2bt(hcalc.vchan,nanmean(pcalc.rcalc,2)))
    printarray([sum(woo-moo(1,:)) mean(woo-moo(1,:)) std(woo-moo(1,:))],'difference between trapz(all) - sum(trapz(RRTM bands)) : total, mean, stdev')

    fprintf(1,'doing fluxes for latbin jj = %2i lonbin ii = %2i \n',jj,ii);
    data = trapz(fchan,FAFA)/1000;         
    zoo = find(isfinite(data));
    if length(zoo) > 20
      %[junk,err] = Math_tsfit_lin_robust(dtime(zoo),data(zoo),4);
      [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dtime,data,4,-1,-1);
      % thestatsradtrend64x72.anomflux(ii,jj,1,:)    = junkanom;
      % thestatsradtrend64x72.trendflux(ii,jj,1)     = junk(2);
      % thestatsradtrend64x72.trendflux_unc(ii,jj,1) = err.se(2);
      thestatsradtrend64x72.anomflux(ii,1,:)    = junkanom;
      thestatsradtrend64x72.trendflux(ii,1)     = junk(2);
      thestatsradtrend64x72.trendflux_unc(ii,1) = err.se(2);
    else
      % thestatsradtrend64x72.anomflux(ii,jj,1,:)    = NaN;
      % thestatsradtrend64x72.trendflux(ii,jj,1)     = NaN;
      % thestatsradtrend64x72.trendflux_unc(ii,jj,1) = NaN;
      thestatsradtrend64x72.anomflux(ii,1,:)    = NaN;
      thestatsradtrend64x72.trendflux(ii,1)     = NaN;
      thestatsradtrend64x72.trendflux_unc(ii,1) = NaN;
    end
    for flfl = 1 : length(RRTM_bands)-1
      junk = find(fchan >= RRTM_bands(flfl) & fchan < RRTM_bands(flfl+1));
      data = squeeze(trapz(fchan(junk),FAFA(junk,:)))/1000;      
      zoo = find(isfinite(data));
      if length(zoo) > 20
        %[junk,err] = Math_tsfit_lin_robust(dtime(zoo),data(zoo),4);
        [junk,err,junkanom] = compute_anomaly_wrapper(zoo,dtime,data,4,-1,-1);
        % thestatsradtrend64x72.anomflux(ii,jj,flfl+1,:)    = junkanom;
        % thestatsradtrend64x72.trendflux(ii,jj,flfl+1)     = junk(2);
        % thestatsradtrend64x72.trendflux_unc(ii,jj,flfl+1) = err.se(2);
        thestatsradtrend64x72.anomflux(ii,flfl+1,:)    = junkanom;
        thestatsradtrend64x72.trendflux(ii,flfl+1)     = junk(2);
        thestatsradtrend64x72.trendflux_unc(ii,flfl+1) = err.se(2);
      else
        % thestatsradtrend64x72.anomflux(ii,jj,flfl+1,:)    = NaN;
        % thestatsradtrend64x72.trendflux(ii,jj,flfl+1)     = NaN;
        % thestatsradtrend64x72.trendflux_unc(ii,jj,flfl+1) = NaN;
        thestatsradtrend64x72.anomflux(ii,flfl+1,:)    = NaN;
        thestatsradtrend64x72.trendflux(ii,flfl+1)     = NaN;
        thestatsradtrend64x72.trendflux_unc(ii,flfl+1) = NaN;
      end
    end

    moo = squeeze(thestatsradtrend64x72.anomflux(ii,:,:));
    woo1 = moo(1,:); woo2 = sum(moo(2:14,:),1);
    hist(100*(woo1 - woo2)./woo1,100)


    fprintf(1,'doing 2645 channels trend and anom for jj = %2i ii = %2i \n',jj,ii);
    for ch = 1 : 2645
      if mod(ch,1000) == 0 
        fprintf(1,'+')
      elseif mod(ch,100) == 0 
        fprintf(1,'.')
      end
      [junkB junkstats junkbtanom junkradanom] = compute_anomaly_wrapper(k,dtime,pcalc.rcalc(ch,:),iNumSineCosCycles,fairs(ch),+1,-1);
      % thestatsradtrend64x72.BTtrend(ii,jj,ch)    = junkB(2);
      % thestatsradtrend64x72.BTtrenderr(ii,jj,ch) = junkstats.se(2);
      % thestatsradtrend64x72.radanom(ii,jj,ch,:)  = junkradanom;
      % thestatsradtrend64x72.BTanom(ii,jj,ch,:)   = junkbtanom;
      thestatsradtrend64x72.BTtrend(ii,ch)    = junkB(2);
      thestatsradtrend64x72.BTtrenderr(ii,ch) = junkstats.se(2);
      thestatsradtrend64x72.radanom(ii,ch,:)  = junkradanom;
      thestatsradtrend64x72.BTanom(ii,ch,:)   = junkbtanom;
    end

    fprintf(1,'\n');

  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% can put together the files using eg 
for ibah = 1 : 64; 
  if mod(ii,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end
  filein = ['/asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_fastgrib_Sept2002_Aug2024_22yr_desc_btanom_latbin_' num2str(ibah,'%02i') '.mat'];
  a = load(filein);
  flux64x72ta.anomflux(:,ibah,:,:) = a.thestatsradtrend64x72.anomflux;  
  flux64x72ta.trendflux(:,ibah,:)  = a.thestatsradtrend64x72.trendflux;  
  flux64x72ta.BTtrend(:,ibah,:)    = a.thestatsradtrend64x72.BTtrend;  
  flux64x72ta.radanom(:,ibah,:,:)  = a.thestatsradtrend64x72.radanom;
  flux64x72ta.BTanom(:,ibah,:,:)   = a.thestatsradtrend64x72.BTanom;
end
fprintf(1,'\n');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
comment = 'see cluster_do_the_fits_airsL3_ratesv7_tiles_radiances.m <-----> do_the_AIRSL3_trends_OQuestioN.m (which is called by eg driver_compute_AIRSL3_trends_desc_or_ascNOQuestioN.m)';
saver = ['save ' fnameout ' thestatsradtrend64x72 filename'];
if ~exist(fnameout)
  eval(saver);
  fprintf(1,'saving to %s \n',fnameout)
else
  fprintf(1,'already have saved %s so not again \n',fnameout)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('now use driver_read64x72flux72_save64zonalflux_ERA5_AIRSL3.m to put the individual latbins together')
disp('now use driver_read64x72flux72_save64zonalflux_ERA5_AIRSL3.m to put the individual latbins together')
disp('now use driver_read64x72flux72_save64zonalflux_ERA5_AIRSL3.m to put the individual latbins together')

