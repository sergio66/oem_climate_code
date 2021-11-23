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

fprintf(1,'sum(iaaFound) at start = %6i \n',sum(sum(iaaFound)));
sum0 = sum(sum(iaaFound));

iDoDOFS = -1;
iDoDOFS = +1;

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

  dofs  = zeros(40,365);
  if iDoDOFS > 0
    cdofs = zeros(40,365,66);
    ak = zeros(40,365,66,66);
  end
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
        dofs(ii,mapp(iTime))      = a.oem.dofs;
        if iDoDOFS > 0
          cdofs(ii,mapp(iTime),:)   = a.oem.cdofs;
          ak(ii,mapp(iTime),:,:)    = a.oem.ak;
        end

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
        wv(ii,mapp(iTime),:)      = a.oem.finalrates(wvind);
        wv_sig(ii,mapp(iTime),:)  = a.oem.finalsigs(wvind);
        tz(ii,mapp(iTime),:)      = a.oem.finalrates(tzind);
        tz_sig(ii,mapp(iTime),:)  = a.oem.finalsigs(tzind);
        o3(ii,mapp(iTime),:)      = a.oem.finalrates(o3ind);
        o3_sig(ii,mapp(iTime),:)  = a.oem.finalsigs(o3ind);

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
        wv(ii,mapp(iTime),:)      = a.oem.finalrates(wvind);
        wv_sig(ii,mapp(iTime),:)  = a.oem.finalsigs(wvind);
        tz(ii,mapp(iTime),:)      = a.oem.finalrates(tzind);
        tz_sig(ii,mapp(iTime),:)  = a.oem.finalsigs(tzind);
        o3(ii,mapp(iTime),:)      = a.oem.finalrates(o3ind);
        o3_sig(ii,mapp(iTime),:)  = a.oem.finalsigs(o3ind);
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
figure(6); clf; pcolor(okdates,latbins,cfc12); caxis([-20 +2]);   colorbar; title('ANOM CFC12 ppb'); shading flat
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
  figure(2); clf; pcolor(okdates,latbins,co2);   caxis([-10 +40]);  colorbar; title('ANOM CO2 ppm'); shading flat

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

  meanco2 = nanmean(co2(iaTropics,:));
  meann2o = nanmean(n2o(iaTropics,:));
  meanch4 = nanmean(ch4(iaTropics,:));
  meancfc11 = nanmean(cfc11(iaTropics,:));
  meancfc12 = nanmean(cfc12(iaTropics,:));
  meanstemp = nanmean(stemp(iaTropics,:));

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

  figure(10); clf; pcolor(okdates,latbins,dofs);   caxis([15 20]);  colorbar; title('DOFS'); shading flat  

  for ii=4:9; figure(ii); disp('showing individual trace gas/stemp; ret to continue'); pause; end

  figure(11); plot(a.jacobian.chanset,f(a.jacobian.chanset),'o-'); title('chans used')
  %disp('ret to continue to T/WV/O3 trends'); pause

elseif co2lays == 3
  %disp('ret to see all co2'); pause;
  pause(1)
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

  figure(7); plot(a.jacobian.chanset,f(a.jacobian.chanset),'o-'); title('chans used')

  disp('ret to continue to T/WV/O3 trends'); pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iDoDOFS > 0
  cxdofs = squeeze(sum(cdofs,3));
  figure(1); clf
    plot(1:40,dofs-cxdofs,'r'); axis([2 40 -1 +1])
    plot(1:40,dofs-cxdofs,'r')

  tdof  = squeeze(sum(cdofs(:,:,a.jacobian.temp_i),3));
  wvdof = squeeze(sum(cdofs(:,:,a.jacobian.water_i),3));
  o3dof = squeeze(sum(cdofs(:,:,a.jacobian.ozone_i),3));
  scalardof = cdofs(:,:,a.jacobian.scalar_i);

  figure(2)
  plot(1:40,nanmean(tdof,2),'r+-',1:40,nanmean(wvdof,2),'b+-',1:40,nanmean(o3dof,2),'g+-'); hl = legend('T','WV','O3'); ylabel('DOFs'); xlabel('Latitude bin');
  grid

  junkdof = squeeze(nanmean(scalardof,2));
  figure(3)
  plot(1:40,junkdof,'+-')
  hl = legend('CO2','N2O','CH4','CFC11','CFC12','Stemp',;ocation','best');
  grid

  figure(3)
  plot(1:40,junkdof(:,[1 2 3 5 6]),'+-','linewidth',2)
  hl = legend('CO2','N2O','CH4','CFC12','Stemp','location','best');
  grid

  disp('ret to continue'); pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
anom_T_WV_plots
