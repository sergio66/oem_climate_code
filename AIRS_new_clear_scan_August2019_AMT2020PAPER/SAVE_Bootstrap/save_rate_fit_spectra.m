%% copied from read_the_anom_retrievals_spectra_40individual_latbins.m
raaObs = nan * ones(2645,365); 
raaCal = nan * ones(2645,365); 

latbinsx = equal_area_spherical_bands(20);
latbins = (latbinsx(1:end-1) + latbinsx(2:end)) * 0.5;

iaTropics = find(abs(latbins) <= 30);
iaTropics = 1 : 40;

iOBSorCAL = input('Enter (0) OBS (1) CAL fits : ');

load save_365datemaps.mat

fdir0 = 'SAVE_nofitT_useERA/';
fdir0 = 'SAVE_BESTRUNv1/';

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
      fname = [fdir0 'OutputAnomaly_OBS/' num2str(ii,'%02d') '/anomtest_timestep' num2str(iTime) '.mat'];
    elseif iOBSorCAL == 1
      fname = [ddir0 'OutputAnomaly_CAL/' num2str(ii,'%02d') '/anomtest_timestep' num2str(iTime) '.mat'];
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

clear a fname iTime iY latbinsx iaTropics raa1 raa2 raaCal raaObs save_dat_day ii iii iOBSorCAL

iSave = -1;
if iSave > 0
  for ii = 1 : 40
    rates = squeeze(raa1ALL(:,ii,:));
    fits  = squeeze(raa1ALL(:,ii,:));
    saver = ['save  ' fdir0 '/SAVESPECTRA/save_spectra_latbin_' num2str(ii) '.mat rates fits save_days save_rtime'];
    eval(saver);
  end
end
