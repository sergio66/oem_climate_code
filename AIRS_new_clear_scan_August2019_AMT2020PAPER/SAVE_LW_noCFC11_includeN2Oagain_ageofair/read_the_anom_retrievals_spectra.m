raaObs = nan * ones(2645,365); 
raaCal = nan * ones(2645,365); 

latbinsx = equal_area_spherical_bands(20);
latbins = (latbinsx(1:end-1) + latbinsx(2:end)) * 0.5;

iaTropics = find(abs(latbins) <= 30);

for iTime = 1 : 365
  if mod(iTime,23) == 0    %% there are 365 days per year, 16 day steps ==> 23 timesteps per year
    fprintf(1,'.');
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
    if exist(fname) & iaaFound(ii,mapp(iTime)) == 1
      loader = ['a = load(''' fname ''');'];
      eval(loader)
      raa1(:,iii) = a.rateset.rates;
      raa2(:,iii) = a.oem.fit';
    end
  end
  allraa1(:,:,iTime) = raa1;
  allraa2(:,:,iTime) = raa2;
  raaObs(:,iTime) = nanmean(raa1');
  raaCal(:,iTime) = nanmean(raa2');
end
fprintf(1,' done \n');

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

ii961 = find(f(chanset) >= 961,1)
figure(7); plot(tObs(chanset(ii961),:)-tCal(chanset(ii961),:))
figure(7); hist(tObs(chanset(ii961),:)-tCal(chanset(ii961),:),100)

topts = a.topts;
if iOBSorCAL == 0
  saver = ['save anomaly_' num2str(iAvgNumDays) 'dayavg_results_spectra.mat okdates okrtime iaTropics raaObs raaCal chanset topts'];
elseif iOBSorCAL == 1
  saver = ['save anomaly_' num2str(iAvgNumDays) 'dayavg_results_spectra_cal.mat okdates okrtime iaTropics raaObs raaCal chanset topts'];
end

%eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is to see bad points

%figure(8)
i1231 = find(f >= 1231,1);
i18 = find(iaTropics == 18);
t18Obs = squeeze(allraa1(i1231,i18,:));
plot(okdates,tObs(i1231,:),'r',okdates,tCal(i1231,:),'b',okdates,stemp(18,:),'kx-',okdates,t18Obs,'g'); 
  hl=legend('mean tropic data','mean tropic fit','STEMP(18,:)','Raw DATA(18,:)','location','best');

plot(okdates,stemp(18,:),'kx-',okdates,t18Obs,'g'); 
  hl=legend('STEMP(18,:)','Raw DATA(18,:)','location','best');

oo2016 = find(okdates >= 2016,1);
plot(f,squeeze(allraa1(:,18,oo2016-5:oo2016+5)),'.-')

%{
whos okdates tObs tCal allraa* stemp f iaTropics co2
save junkspectra_vs_stemp.mat okdates tObs tCal allraa* stemp f iaTropics co2
%}

