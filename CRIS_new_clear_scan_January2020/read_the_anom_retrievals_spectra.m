raaObs = nan * ones(1305,157); 
raaCal = nan * ones(1305,157); 

latbinsx = equal_area_spherical_bands(20);
latbins = (latbinsx(1:end-1) + latbinsx(2:end)) * 0.5;

iaTropics = find(abs(latbins) <= 30);

for iTime = 1 : 157
  if mod(iTime,23) == 0    %% there are 365 days per year, 16 day steps ==> 23 timesteps per year
    fprintf(1,'.');
  end
  raa1 = nan * ones(1305,length(iaTropics));
  raa2 = nan * ones(1305,length(iaTropics));
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
xyz = load('f1305.mat');
fcris = xyz.f1305;
f = fcris;
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
figure(7); clf; plot(tObs(chanset(ii961),:)-tCal(chanset(ii961),:))
figure(7); clf; hist(tObs(chanset(ii961),:)-tCal(chanset(ii961),:),100)

topts = a.topts;
if iOBSorCAL == 0
  saver = ['save anomaly_' num2str(iAvgNumDays) 'dayavg_results_spectra.mat okdates okrtime iaTropics raaObs raaCal chanset topts'];
elseif iOBSorCAL == 1
  saver = ['save anomaly_' num2str(iAvgNumDays) 'dayavg_results_spectra_cal.mat okdates okrtime iaTropics raaObs raaCal chanset topts'];
end

iSave = input('save the output (-1 default/+1) : ');
if length(iSave == 0)
  iSave = -1;
end
if iSave > 0
  eval(saver)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is to see bad points

figure(8)
i1231 = find(f >= 1231,1);
i18 = find(iaTropics == 18);
t18Obs = squeeze(allraa1(i1231,i18,:));
plot(okdates,tObs(i1231,:),'r',okdates,tCal(i1231,:),'b',okdates,stemp(18,:),'kx-',okdates,t18Obs,'g','linewidth',2); 
  hl=legend('mean BT1231 tropic data','mean BT1231 tropic fit','STEMP fit (18,:)','Raw DATA(18,:)','location','best');
title('Latbin 18')

plot(okdates,stemp(18,:),'kx-',okdates,t18Obs,'g'); 
  hl=legend('STEMP(18,:)','Raw DATA(18,:)','location','best');

figure(9)
plot(f,squeeze(allraa1(:,18,80:90)),'.-')
oo2016 = find(okdates >= 2016,1);
plot(f,squeeze(allraa1(:,18,oo2016-5:oo2016+5)),'.-')
title('Latbin 18 : 2016 +/- 5*16 days')

figure(10);
pcolor(okdates,f,squeeze(allraa1(:,18,:))); shading flat; title('Signal, latbin 18')
pcolor(okdates,f,squeeze(allraa2(:,18,:))); shading flat; title('Fit, latbin18')
pcolor(okdates,f,squeeze(allraa1(:,18,:))-squeeze(allraa2(:,18,:))); shading flat; title('Signal-Fit, latbin 18')
pcolor(okdates,f(chanset),squeeze(allraa1(chanset,18,:))-squeeze(allraa2(chanset,18,:))); shading flat; title('Signal-Fit, latbin 18')
   caxis([-1 +1]); colorbar

pcolor(okdates,f(chanset),raaObs(chanset,:)-raaCal(chanset,:)); shading flat; title('Signal-Fit, <tropical>')
   caxis([-0.1 +0.1]); colorbar
pcolor(okdates,f(chanset),raaObs(chanset,:)); shading flat; title('Signal-Fit, <tropical>')
   caxis([-1 +1]); colorbar
  hold on
    plot(okdates,790+50*smooth(meanstemp,2*5),'k',okdates,790+50*smooth(mean(save_dat_1231(:,iaTropics)'),10),'rx-',okdates,790+50*nanmean(smoothtz,2),'bs-',okdates,790*ones(size(okdates)),'k--','linewidth',4);
  hold off
axis([min(okdates) max(okdates) 640 940])
hl = legend('mean spectra','mean Stemp','mean BT1231','mean T800 mb');

figure(11)
clf; [mmbad,nnbad] = find(stemporig > 1.75); plot(mmbad,nnbad,'o'); whos mmbad nnbad; lala = find(mmbad > 3 & mmbad < 35); [mmbad(lala) nnbad(lala)]
for xyz = 1 : length(lala)
  wah1(:,xyz) = squeeze(allraa1(:,mmbad(lala(xyz)),nnbad(lala(xyz)))); 
  wah2(:,xyz) = squeeze(allraa2(:,mmbad(lala(xyz)),nnbad(lala(xyz)))); 
end
if length(lala) > 0
  plot(f,wah1,'r.-',f,wah2,'b.-'); grid
end
title('Bad points')

%{
whos okdates tObs tCal allraa* stemp f iaTropics co2
save junkspectra_vs_stemp.mat okdates tObs tCal allraa* stemp f iaTropics co2
%}

% [mmbad(lala) nnbad(lala)]
%    19    84
%    20    84
%    20    85

