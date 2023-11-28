addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/

%% https://acp.copernicus.org/articles/22/10603/2022/ Christian Borger et al      Jan 2005 - Dec2021  OMI total column water
omi = read_netcdf_lls('/home/sergio/PAPERS/SUBMITPAPERS/trends//MPIC_OMI_TCWV_v1.2.nc');

days = (1:12*16)*30;
for jj = 1 : 180
  if mod(jj,100) == 0
    fprintf(1,'+');
  elseif mod(jj,10) == 0
    fprintf(1,'.');
  end
  for ii = 1 : 360
    data = squeeze(omi.tcwv(ii,jj,:));
    zoo = find(isfinite(data));
    if length(zoo) > 20
      [junk,err] = Math_tsfit_lin_robust(days(zoo),data(zoo),4);
      omi_colwv.trend(ii,jj)     = junk(2);
      omi_colwv.trend_err(ii,jj) = err.se(2);
    else
      omi_colwv.trend(ii,jj) = NaN;
      omi_colwv.trend_err(ii,jj) = NaN;
    end
  end
end
fprintf(1,'\n');

omi_colwv.lat = omi.latitude;
omi_colwv.lon = omi.longitude;
save omi_tcwv_trends.mat omi_colwv 

figure(1); clf; simplemap(omi_colwv.trend'); caxis([-1 +1]*0.15)
figure(2); clf; plot(omi.latitude,nanmean(omi_colwv.trend'))
