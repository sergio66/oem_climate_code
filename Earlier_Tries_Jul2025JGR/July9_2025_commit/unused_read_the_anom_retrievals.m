fprintf(1,'sum(iaaFound) at start = %6i \n',sum(sum(iaaFound)));
sum0 = sum(sum(iaaFound));

if ~exist('co2')
  co2 = zeros(40,365);
  n2o = zeros(40,365);
  ch4 = zeros(40,365);
  cfc11 = zeros(40,365);
  cfc12 = zeros(40,365);
  stemp = zeros(40,365);

  co2_sig = zeros(40,365);
  n2o_sig = zeros(40,365);
  ch4_sig = zeros(40,365);
  cfc11_sig = zeros(40,365);
  cfc12_sig = zeros(40,365);
  stemp_sig = zeros(40,365);
end

for ii = 1 : 40
  fprintf(1,'reading in data for latbin %2i ',ii);
  mapp = save_days_map(:,ii);
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
    if exist(fname) & iaaFound(ii,mapp(iTime)) == 0
      loader = ['a = load(''' fname ''');'];

      eval(loader)
      iaaFound(ii,mapp(iTime))  = +1;

      if co2lays == 1
        dofs(ii,mapp(iTime))      = a.oem.dofs;
        co2(ii,mapp(iTime))       = a.oem.finalrates(1);
        co2_sig(ii,mapp(iTime))   = a.oem.finalsigs(1);
        n2o(ii,mapp(iTime))       = a.oem.finalrates(2);
        n2o_sig(ii,mapp(iTime))   = a.oem.finalsigs(2);
        ch4(ii,mapp(iTime))       = a.oem.finalrates(3);
        ch4_sig(ii,mapp(iTime))   = a.oem.finalsigs(3);
        cfc11(ii,mapp(iTime))     = a.oem.finalrates(4);
        cfc11_sig(ii,mapp(iTime)) = a.oem.finalsigs(4);
        cfc12(ii,mapp(iTime))     = a.oem.finalrates(5);
        cfc12_sig(ii,mapp(iTime)) = a.oem.finalsigs(5);
        stemp(ii,mapp(iTime))     = a.oem.finalrates(6);
        stemp_sig(ii,mapp(iTime)) = a.oem.finalsigs(6);
        bestloop(ii,mapp(iTime))  = a.oem.bestloop;
      elseif co2lays == 3
        dofs(ii,mapp(iTime))      = a.oem.dofs;
        aco2(ii,mapp(iTime))       = a.oem.finalrates(1);
        aco2_sig(ii,mapp(iTime))   = a.oem.finalsigs(1);
        bco2(ii,mapp(iTime))       = a.oem.finalrates(2);
        bco2_sig(ii,mapp(iTime))   = a.oem.finalsigs(2);
        cco2(ii,mapp(iTime))       = a.oem.finalrates(3);
        cco2_sig(ii,mapp(iTime))   = a.oem.finalsigs(3);
        n2o(ii,mapp(iTime))       = a.oem.finalrates(4);
        n2o_sig(ii,mapp(iTime))   = a.oem.finalsigs(4);
        ch4(ii,mapp(iTime))       = a.oem.finalrates(5);
        ch4_sig(ii,mapp(iTime))   = a.oem.finalsigs(5);
        cfc11(ii,mapp(iTime))     = a.oem.finalrates(6);
        cfc11_sig(ii,mapp(iTime)) = a.oem.finalsigs(6);
        cfc12(ii,mapp(iTime))     = a.oem.finalrates(7);
        cfc12_sig(ii,mapp(iTime)) = a.oem.finalsigs(7);
        stemp(ii,mapp(iTime))     = a.oem.finalrates(8);
        stemp_sig(ii,mapp(iTime)) = a.oem.finalsigs(8);
        bestloop(ii,mapp(iTime))  = a.oem.bestloop;
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
figure(7); plot(sum(iaaFound'),'o-'); xlabel('Latbins'); ylabel('Made so far')
figure(8); plot(sum(iaaFound),'o-');  xlabel('TimeStep'); ylabel('Made so far')
[bbm,bbn] = size(bestloop);
if bbm == 40 & bbn == length(okdates)
  figure(9); pcolor(okdates,latbins,bestloop); title('bestloop'); colorbar; colormap jet; shading flat
end

fprintf(1,'sum(iaaFound) at start = %6i \n',sum0);
fprintf(1,'sum(iaaFound) at end   = %6i \n',sum(sum(iaaFound)));
fprintf(1,'expect about 365*40 = 14600 files (shoot for >= 14579) \n');

if co2lays == 3
  co2 = aco2;
  co2 = cco2;
  co2 = bco2;
end


figure(2); clf; pcolor(okdates,latbins,co2);   caxis([-10 +40]);  colorbar; title('ANOM CO2 ppm'); shading flat
figure(3); clf; pcolor(okdates,latbins,n2o);   caxis([-2 +20]);   colorbar; title('ANOM N2O ppm'); shading flat
figure(4); clf; pcolor(okdates,latbins,ch4);   caxis([-10 +100]); colorbar; title('ANOM CH4 ppb'); shading flat
figure(5); clf; pcolor(okdates,latbins,cfc11); caxis([-20 +2]);   colorbar; title('ANOM CFC11 ppb'); shading flat
figure(6); clf; pcolor(okdates,latbins,cfc11); caxis([-20 +2]);   colorbar; title('ANOM CFC12 ppb'); shading flat
figure(7); clf; pcolor(okdates,latbins,stemp); caxis([-1 +1]);    colorbar; title('ANOM STEMP K'); shading flat

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
[badrow,badcol] = find(iaaFound == 0);
baddy = unique(badcol)'
if length(baddy) > 0
  fid = fopen('badanom.sc','w');
  for bb = 1 : length(baddy)
    str = ['sbatch -exclude=cnode260,cnode267 --array=' num2str(baddy(bb)) ' sergio_matlab_jobB.sbatch'];
    fprintf(fid,'%s \n',str);
  end
  fclose(fid);
end

disp('now keep running read_the_anom_retrievals  (and/or read_the_anom_retrievals_spectra)')
if co2lays == 1
  disp('ret to see all co2'); pause;
  figure(2); clf; pcolor(okdates,latbins,co2);   caxis([-10 +40]);  colorbar; title('ANOM CO2 ppm'); shading flat

  latbinsx = equal_area_spherical_bands(20);
  latbins = (latbinsx(1:end-1) + latbinsx(2:end)) * 0.5;
  iaTropics = find(abs(latbins) <= 30);

  meanco2 = nanmean(co2(iaTropics,:));
  meann2o = nanmean(n2o(iaTropics,:));
  meanch4 = nanmean(ch4(iaTropics,:));
  meancfc11 = nanmean(cfc11(iaTropics,:));
  meancfc12 = nanmean(cfc12(iaTropics,:));
  meanstemp = nanmean(stemp(iaTropics,:));

  figure(4); plot(okdates,smooth(meanco2,2*5),'b','linewidth',2); 
		    hl = legend('mean tropical co2','location','best'); set(hl,'fontsize',10);
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical CO2 smoothed over 10 years')

  figure(5); plot(okdates,smooth(meanstemp,2*5),'b','linewidth',2); 
  ax = axis; ax(1) = min(okdates); ax(2) = max(okdates); axis(ax); grid; title('tropical stemp smoothed over 10 years')

  figure(6); clf; pcolor(okdates,latbins,dofs);   caxis([15 20]);  colorbar; title('DOFS'); shading flat  

elseif co2lays == 3
  disp('ret to see all co2'); pause;
  figure(1); clf; pcolor(okdates,latbins,aco2);   caxis([-10 +40]);  colorbar; title('ANOM ACO2 ppm'); shading flat
  figure(2); clf; pcolor(okdates,latbins,bco2);   caxis([-10 +40]);  colorbar; title('ANOM BCO2 ppm'); shading flat
  figure(3); clf; pcolor(okdates,latbins,cco2);   caxis([-10 +40]);  colorbar; title('ANOM CCO2 ppm'); shading flat

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

  figure(6); clf; pcolor(okdates,latbins,dofs);   caxis([15 20]);  colorbar; title('DOFS'); shading flat  
end
