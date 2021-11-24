function allrates = get_all_387_rates(driver,iLen)

%{
if nargin == 1
  iLen = 387;
end

if iLen == -1
  iLen = length(dtime);
else
  iLen = 387;
end
%}

addpath StrowCode
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP

if nargin == 1
  iLen = 387;
end

airs = instr_chans2645;

latstr = num2str(driver.latlon.latbin,'%02d');
lonstr = num2str(driver.latlon.lonbin,'%02d');

%% sergio quantiles
d = load(['/home/strow/Work/Airs/Tiles/Data/Quantv1/LatBin' latstr '/LonBin' lonstr '/summarystats_LatBin' latstr '_LonBin' lonstr '_timesetps_001_412_V1.mat']);

%% strow anomalies
loader = ['load /home/strow/Work/Airs/Tiles/Data/Quantv1_fits/LatBin' latstr  '/LonBin' lonstr '/fits_LonBin' lonstr '_LatBin' latstr '_V1.mat'];
eval(loader);

if iLen == -1
  iLen = length(dtime);
end

dtime = dtime(1:iLen); %% fixing this to iLen
k_desc = k_desc(1:iLen);

ch = 1520;
qi = 16;
r = squeeze(d.rad_quantile_desc(:,ch,qi)); r = r(1:iLen,1);
[xbt_anom xr_anom] = compute_anomaly(k_desc,dtime,b_desc(ch,qi,:),fairs(ch),r);

bt_anom = nan(2645,2,iLen);
r_anom  = nan(2645,2,iLen);
for iich = 1:2645
  if mod(iich,1000) == 0
    fprintf(1,'++')
  elseif mod(iich,500) == 0
    fprintf(1,'+')
  elseif mod(iich,100) == 0
    fprintf(1,'.')
  end
  for jjqi = 15 : 16
    r = squeeze(d.rad_quantile_desc(:,iich,jjqi)); r = r(1:iLen,1);
    [ bt_anom(iich, jjqi-14,:)  r_anom(iich,jjqi-14,:)] = compute_anomaly(k_desc,dtime,b_desc(iich,jjqi,:),fairs(iich),r);
    [dbt_anom(iich, jjqi-14,:) dr_anom(iich,jjqi-14,:)] = compute_anomaly(k_desc,dtime,berr_desc(iich,jjqi,:),fairs(iich),r);
  end
end

fprintf(1,'\n');

%whos bt_anom dbt_anom

%allrates.rates      = squeeze(nanmean(bt_anom(:,15:16,:),2));   %% err check the size, dummy, before you do that
%allrates.xunc_rates = squeeze(nanmean(dbt_anom(:,15:16,:),2));

allrates.rates      = squeeze(nanmean(bt_anom,2));
allrates.xunc_rates = squeeze(nanmean(dbt_anom,2));
allrates.xlagcorr   = nanmean(lag_desc(:,15:16),2);

allrates.unc_rates = ones(size(allrates.rates)) * 0.001;
allrates.unc_rates = ones(size(allrates.rates)) * 0.01;
allrates.unc_rates = ones(size(allrates.rates)) * 0.1;

moo = find(abs(allrates.rates(1520,:)) < eps);
fprintf(1,'lon/lat = %s %s has zero anomaly for ch1520 %3i out of %3i .. so expect only %3i output retrievals \n',lonstr,latstr,length(moo),iLen,iLen-length(moo));
pause(0.1);

%allrates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf; pcolor(fairs,1:iLen,allrates.rates'); shading flat
caxis([-10 +10]); colormap(usa2); colorbar; axis([640 1640 0 iLen])

figure(2); clf; pcolor(fairs,1:iLen,allrates.unc_rates'); shading flat
caxis([0 +1]); colormap(usa2); colorbar; axis([640 1640 0 iLen])

%% good try but looks awful
figure(3); clf; pcolor(fairs,1:iLen,allrates.xunc_rates'); shading flat
caxis([0 100]); colormap(usa2); colorbar; axis([640 1640 0 iLen])

