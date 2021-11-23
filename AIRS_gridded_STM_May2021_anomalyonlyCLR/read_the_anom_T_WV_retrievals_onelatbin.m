hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
f = vchan2834;
load sarta_chans_for_l1c.mat
f = f(ichan);

%%%%%%%%%%%%%%%%%%%%%%%%%

disp('checking number of running ANOM jobs')
sqr = ['!squeue | grep -in '' ANOMALY '' | grep -in ''batch'' | grep -in '' R '' | wc -l']; eval(sqr);
disp('checking number of delayed ANOM jobs')
sqr = ['!squeue | grep -in '' ANOMALY '' | grep -in ''batch'' | grep -in '' PD '' | wc -l']; eval(sqr);

fprintf(1,'sum(aaFound) at start = %6i \n',sum(iaFound));
sum0 = sum(iaFound);

if ~exist('co2')
  co2 = zeros(1,365);
  n2o = zeros(1,365);
  ch4 = zeros(1,365);
  cfc11 = zeros(1,365);
  cfc12 = zeros(1,365);
  stemp = zeros(1,365);

  co2_sig = zeros(1,365);
  n2o_sig = zeros(1,365);
  ch4_sig = zeros(1,365);
  cfc11_sig = zeros(1,365);
  cfc12_sig = zeros(1,365);
  stemp_sig = zeros(1,365);

  dofs  = zeros(1,365);
  %cdofs = zeros(1,365,66);
  %ak = zeros(1,365,66,66);
end

for ii = iReadLatbin
  iix = 1;
  fprintf(1,'reading in data for latbin %2i ',ii);
  mapp = save_days_map;
  mapp = mapp(mapp > 0);
  for iTime = 1 : 365
    if mod(iTime,23) == 0    %% there are 365 days per year, 16 day steps ==> 23 timesteps per year
      fprintf(1,'.');
    end

    if iOBSorCAL == 0
      fname = ['OutputAnomaly_OBS/' num2str(ii,'%02d') '/anomtest_timestep' num2str(iTime) '.mat'];
    elseif iOBSorCAL == 1
      fname = ['OutputAnomaly_CAL/' num2str(ii,'%02d') '/anomtest_timestep' num2str(iTime) '.mat'];
    end

    if exist(fname) & iaFound(mapp(iTime)) == 0
      loader = ['a = load(''' fname ''');'];

      eval(loader)
      iaFound(mapp(iTime))  = +1;

      %% a.jacobian.water_i a.jacobian.temp_i
      nlays = length(a.oem.wunc);
      if co2lays == 1
        %%% a.jacobian.wvjaclays_offset = 6!!!
        indoffset = a.jacobian.wvjaclays_offset;
        wvind = (1:nlays) + 0*nlays + indoffset;
        tzind = (1:nlays) + 1*nlays + indoffset;
        o3ind = (1:nlays) + 2*nlays + indoffset;
      elseif co2lays == 3
        wvind = (1:nlays) + 0*nlays + indoffset;
        tzind = (1:nlays) + 1*nlays + indoffset;
        o3ind = (1:nlays) + 2*nlays + indoffset;
      end

      if co2lays == 1
        dofs(mapp(iTime))      = a.oem.dofs;
        %cdofs(mapp(iTime),:)   = a.oem.cdofs;
        %ak(mapp(iTime),:,:)    = a.oem.ak;

        co2(mapp(iTime))       = a.oem.finalrates(1);
        co2_sig(mapp(iTime))   = a.oem.finalsigs(1);
        n2o(mapp(iTime))       = a.oem.finalrates(2);
        n2o_sig(mapp(iTime))   = a.oem.finalsigs(2);
        ch4(mapp(iTime))       = a.oem.finalrates(3);
        ch4_sig(mapp(iTime))   = a.oem.finalsigs(3);
        cfc11(mapp(iTime))     = a.oem.finalrates(4);
        cfc11_sig(mapp(iTime)) = a.oem.finalsigs(4);
        cfc12(mapp(iTime))     = a.oem.finalrates(5);
        cfc12_sig(mapp(iTime)) = a.oem.finalsigs(5);
        stemp(mapp(iTime))     = a.oem.finalrates(6);
        stemp_sig(mapp(iTime)) = a.oem.finalsigs(6);
        bestloop(mapp(iTime))  = a.oem.bestloop;
        wv(mapp(iTime),:)      = a.oem.finalrates(wvind);
        wv_sig(mapp(iTime),:)  = a.oem.finalsigs(wvind);
        tz(mapp(iTime),:)      = a.oem.finalrates(tzind);
        tz_sig(mapp(iTime),:)  = a.oem.finalsigs(tzind);
        o3(mapp(iTime),:)      = a.oem.finalrates(o3ind);
        o3_sig(mapp(iTime),:)  = a.oem.finalsigs(o3ind);

      elseif co2lays == 3
        dofs(mapp(iTime))      = a.oem.dofs;
        aco2(mapp(iTime))       = a.oem.finalrates(1);
        aco2_sig(mapp(iTime))   = a.oem.finalsigs(1);
        bco2(mapp(iTime))       = a.oem.finalrates(2);
        bco2_sig(mapp(iTime))   = a.oem.finalsigs(2);
        cco2(mapp(iTime))       = a.oem.finalrates(3);
        cco2_sig(mapp(iTime))   = a.oem.finalsigs(3);
        n2o(mapp(iTime))       = a.oem.finalrates(4);
        n2o_sig(mapp(iTime))   = a.oem.finalsigs(4);
        ch4(mapp(iTime))       = a.oem.finalrates(5);
        ch4_sig(mapp(iTime))   = a.oem.finalsigs(5);
        cfc11(mapp(iTime))     = a.oem.finalrates(6);
        cfc11_sig(mapp(iTime)) = a.oem.finalsigs(6);
        cfc12(mapp(iTime))     = a.oem.finalrates(7);
        cfc12_sig(mapp(iTime)) = a.oem.finalsigs(7);
        stemp(mapp(iTime))     = a.oem.finalrates(8);
        stemp_sig(mapp(iTime)) = a.oem.finalsigs(8);
        bestloop(mapp(iTime))  = a.oem.bestloop;
        wv(mapp(iTime),:)      = a.oem.finalrates(wvind);
        wv_sig(mapp(iTime),:)  = a.oem.finalsigs(wvind);
        tz(mapp(iTime),:)      = a.oem.finalrates(tzind);
        tz_sig(mapp(iTime),:)  = a.oem.finalsigs(tzind);
        o3(mapp(iTime),:)      = a.oem.finalrates(o3ind);
        o3_sig(mapp(iTime),:)  = a.oem.finalsigs(o3ind);
      end
    end
  end
  fprintf(1,' done \n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%
okrtime = save_rtime(1:length(co2));
okdates = 2002+save_days/365;
okdates = okdates(1:length(co2));

latbinsx = equal_area_spherical_bands(20);
latbins = (latbinsx(1:end-1) + latbinsx(2:end)) * 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7); plot(iaFound,'o-'); xlabel('Latbins'); ylabel('Made so far')
figure(8); plot(iaFound,'o-');  xlabel('TimeStep'); ylabel('Made so far')
[bbm,bbn] = size(bestloop);
if bbm == 40 & bbn == length(okdates)
  figure(9); pcolor(okdates,latbins,bestloop); title('bestloop'); colorbar; colormap jet; shading flat
end

fprintf(1,'sum(iaFound) at start = %6i \n',sum0);
fprintf(1,'sum(iaFound) at end   = %6i \n',sum(iaFound));
fprintf(1,'expect about 365 files \n')

if co2lays == 3
  co2 = aco2;
  co2 = cco2;
  co2 = bco2;
end

figure(2); clf; plot(okdates,co2);   caxis([-10 +40]);  colorbar; title('ANOM CO2 ppm'); shading flat
figure(3); clf; plot(okdates,n2o);   caxis([-2 +20]);   colorbar; title('ANOM N2O ppm'); shading flat
figure(4); clf; plot(okdates,ch4);   caxis([-10 +100]); colorbar; title('ANOM CH4 ppb'); shading flat
figure(5); clf; plot(okdates,cfc11); caxis([-20 +2]);   colorbar; title('ANOM CFC11 ppb'); shading flat
figure(6); clf; plot(okdates,cfc12); caxis([-20 +2]);   colorbar; title('ANOM CFC12 ppb'); shading flat
figure(7); clf; plot(okdates,stemp); caxis([-1 +1]);    colorbar; title('ANOM STEMP K'); shading flat

for ii = 1 : 6
  figure(ii); colormap jet; shading interp
end

[yyok,mmok,ddok,hhok] = tai2utcSergio(okrtime);
topts = a.topts;
topts.nloop = a.oem.nloop;

if exist('ameanco2')
  clear *mean*
end

if iOBSorCAL == 0
  saver = ['save anomaly_' num2str(iAvgNumDays) 'dayavg_results.mat okdates okrtime latbins *co2 n2o ch4 cfc11 cfc12 stemp topts dofs bestloop'];
elseif iOBSorCAL == 1
  saver = ['save anomaly_' num2str(iAvgNumDays) 'dayavg_cal_results.mat okdates okrtime latbins *co2 n2o ch4 cfc11 cfc12 stemp topts dofs bestloop'];
end
%eval(saver)

disp('still need to get the cluster to do these timesteps')
[badrow,badcol] = find(iaFound == 0);
baddy = unique(badcol)'
if length(baddy) > 0
  %for bb = 1 : length(baddy)
  %  str = ['sbatch -exclude=cnode203,cnode204,cnode260,cnode267 --array=' num2str(baddy(bb)) ' sergio_matlab_jobB.sbatch'];
  %  fprintf(fid,'%s \n',str);
  %end

  fid = fopen('badanom.sc','w');
  fprintf(1,'found that %3i of 365 timesteps did not finish : see badanom.sc \n',length(baddy))
  str = ['sbatch --account=pi_strow --exclude=cnode[204,225,267] --array='];
  fprintf(fid,'%s',str);
  iX = nice_output(fid,baddy);   %% now put in continuous strips
  fprintf(1,'length(badanom) = %4i Num Continuous Strips = %4i \n',length(baddy),iX)
  str = [' sergio_matlab_jobB.sbatch '];
  fprintf(fid,'%s \n',str);
  fclose(fid);
end

disp('now keep running read_the_anom_T_WV_retrievals  (and/or read_the_anom_retrievals_spectra)')

if co2lays == 1
  %disp('ret to see all co2'); pause;
  pause(1)
  figure(2); clf; plot(okdates,co2);   caxis([-10 +40]);  colorbar; title('ANOM CO2 ppm'); shading flat

  latbinsx = equal_area_spherical_bands(20);
  latbins = (latbinsx(1:end-1) + latbinsx(2:end)) * 0.5;
  iaTropics = find(abs(latbins) <= 30);

  bad2 = find(abs(co2) > 45);
  bad1 = find(abs(stemp) > 1);
  badco2 = union(bad1,bad2);
  co2(badco2) = NaN;
  n2o(badco2) = NaN;
  ch4(badco2) = NaN;
  cfc11(badco2) = NaN;
  cfc12(badco2) = NaN;
  stemp(badco2) = NaN;

  meanco2 = co2;
  meann2o = n2o;
  meanch4 = ch4;
  meancfc11 = cfc11;
  meancfc12 = cfc12;
  meanstemp = stemp;

  figure(4); plot(okdates,smooth(meanco2,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical CO2 smoothed over 10 years')

  figure(5); plot(okdates,smooth(meann2o,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical N2O smoothed over 10 years')

  figure(6); plot(okdates,smooth(meanch4,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical CH4 smoothed over 10 years')

  figure(7); plot(okdates,smooth(meancfc11,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical CFC11 smoothed over 10 years')

  figure(8); plot(okdates,smooth(meancfc12,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical CFC12 smoothed over 10 years')

  figure(9); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical stemp smoothed over 10 years')

  figure(10); clf; plot(okdates,dofs);   caxis([15 20]);  colorbar; title('DOFS'); shading flat  

  for ii=4:9; figure(ii); disp('showing individual trace gas/stemp; ret to continue'); pause; end

  figure(11); plot(a.jacobian.chanset,f(a.jacobian.chanset),'o-'); title('chans used')
  %disp('ret to continue to T/WV/O3 trends'); pause

elseif co2lays == 3
  %disp('ret to see all co2'); pause;
  pause(1)
  figure(1); clf; plot(okdates,aco2);   caxis([-10 +40]);  colorbar; title('ANOM ACO2 ppm'); shading flat
  figure(2); clf; plot(okdates,bco2);   caxis([-10 +40]);  colorbar; title('ANOM BCO2 ppm'); shading flat
  figure(3); clf; plot(okdates,cco2);   caxis([-10 +40]);  colorbar; title('ANOM CCO2 ppm'); shading flat

  latbinsx = equal_area_spherical_bands(20);
  latbins = (latbinsx(1:end-1) + latbinsx(2:end)) * 0.5;
  iaTropics = find(abs(latbins) <= 30);

  ameanco2 = nanmean(aco2(iaTropics,:));
  bmeanco2 = nanmean(bco2(iaTropics,:));
  cmeanco2 = nanmean(cco2(iaTropics,:));
  meann2o = nanmean(n2o(iaTropics,:));
  meanch4 = nanmean(ch4(iaTropics,:));
  meancfc11 = nanmean(cfc11(iaTropics,:));
  meancfc12 = nanmean(cfc12(iaTropics,:));
  meanstemp = nanmean(stemp(iaTropics,:));

  figure(4); plot(okdates,smooth(ameanco2,2*5),'b',okdates,smooth(bmeanco2,2*5),'r',okdates,smooth(cmeanco2,2*5),'g','linewidth',2); 
  hl = legend('a mean tropical co2','b mean tropical co2','c mean tropical co2','location','best'); set(hl,'fontsize',10);
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical CO2 smoothed over 10 years')

  figure(5); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical stemp smoothed over 10 years')

  figure(6); clf; plot(okdates,dofs);   caxis([15 20]);  colorbar; title('DOFS'); shading flat  

  figure(7); plot(a.jacobian.chanset,f(a.jacobian.chanset),'o-'); title('chans used')

  disp('ret to continue to T/WV/O3 trends'); pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
anom_T_WV_plots_onelatbin
