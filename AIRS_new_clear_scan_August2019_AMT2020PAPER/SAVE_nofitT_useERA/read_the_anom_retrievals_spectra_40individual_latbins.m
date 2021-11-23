raaObs = nan * ones(2645,365); 
raaCal = nan * ones(2645,365); 

latbinsx = equal_area_spherical_bands(20);
latbins = (latbinsx(1:end-1) + latbinsx(2:end)) * 0.5;

iaTropics = find(abs(latbins) <= 30);
iaTropics = 1 : 40;

iOBSorCAL = input('Enter (0) OBS (1) CAL fits : ');

load save_365datemaps.mat

iY = 0;
for iTime = 1 : 365
  if mod(iTime,23) == 0    %% there are 365 days per year, 16 day steps ==> 23 timesteps per year
    iY = iY + 1;
    if iY ~= 10 
      fprintf(1,'.');
    else
      fprintf(1,'+');
    end
  end
  raa1 = nan * ones(2645,length(iaTropics));
  raa2 = nan * ones(2645,length(iaTropics));
  for iii = 1 : length(iaTropics)
    ii = iaTropics(iii);
    mapp = save_days_map(:,ii);
    mapp = mapp(mapp > 0);
    if iOBSorCAL == 0
      fname = ['OutputAnomaly_OBS/' num2str(ii,'%02d') '/anomtest_timestep' num2str(iTime) '.mat'];
    elseif iOBSorCAL == 1
      fname = ['OutputAnomaly_CAL/' num2str(ii,'%02d') '/anomtest_timestep' num2str(iTime) '.mat'];
    end
    %if exist(fname) & iaaFound(ii,mapp(iTime)) == 1
    if exist(fname)
      loader = ['a = load(''' fname ''');'];
      eval(loader)
      raa1ALL(:,ii,iTime) = a.rateset.rates;
      raa2ALL(:,ii,iTime) = a.oem.fit';
      raa1(:,iii) = a.rateset.rates;
      raa2(:,iii) = a.oem.fit';
    end
  end
  raaObs(:,iTime) = nanmean(raa1');
  raaCal(:,iTime) = nanmean(raa2');
end
fprintf(1,' done \n');

error('ooooo')

%%%%%%%%%%%%%%%%%%%%%%%%%
hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
f = vchan2834;

load sarta_chans_for_l1c.mat
f = f(ichan);
%%%%%%%%%%%%%%%%%%%%%%%%%

%tObs = rad2bt(f,raaObs);
%tCal = rad2bt(f,raaCal);
%% the anoms are already BT
tObs = raaObs;
tCal = raaCal;
chanset = a.jacobian.chanset;

figure(6)
plot(f(chanset),nanmean(tObs(chanset,:)'-tCal(chanset,:)'),'b',f(chanset),nanstd(tObs(chanset,:)'-tCal(chanset,:)'),'r',...
     f(chanset),nanmean(tObs(chanset,:)'),'k',f(chanset),nanstd(tObs(chanset,:)'),'k--',f(chanset),-nanstd(tObs(chanset,:)'),'k--','linewidth',2)
  ylabel('Obs-Cal'); hl = legend('mean','std','mean signal','+std signal','-std signal','location','best'); grid
  set(hl,'fontsize',8); 

topts = a.topts;
%if iOBSorCAL == 0
%  saver = ['save anomaly_' num2str(iAvgNumDays) 'dayavg_results_spectra.mat okdates okrtime iaTropics raaObs raaCal chanset topts'];
%elseif iOBSorCAL == 1
%  saver = ['save anomaly_' num2str(iAvgNumDays) 'dayavg_results_spectra_cal.mat okdates okrtime iaTropics raaObs raaCal chanset topts'];
%end

%eval(saver)

