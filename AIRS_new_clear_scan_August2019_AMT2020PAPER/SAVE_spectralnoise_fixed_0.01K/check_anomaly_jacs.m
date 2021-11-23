addpath /home/sergio/MATLABCODE/PLOTTER
latbinsx = equal_area_spherical_bands(20);
latbins = (latbinsx(1:end-1) + latbinsx(2:end)) * 0.5;
iaTropics = find(abs(latbins) < 30);

for ii = 1 : 365
  if mod(ii,25) == 0
    fprintf(1,'%3i of 365 \n',ii);
  end
 driver.jacobian.filename = ['../../oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/Anomaly365_16/RESULTS/kcarta_' num2str(ii,'%03d') '_M_TS_jac_all_5_97_97_97_2645.mat'];
  a = load(driver.jacobian.filename);
  f = a.f;
  co2_20(ii,:) = squeeze(a.M_TS_jac_all(20,:,1));

  junk = squeeze(a.M_TS_jac_all(iaTropics,:,1));
  co2_mean(ii,:) = nanmean(junk);
  co2_std(ii,:)  = nanstd(junk);
end

iaChan = find(f <= 820);
iaDate = 2002.8 + (1:365)*(16/365);
figure(1); pcolor(iaDate,f(iaChan),co2_20(:,iaChan)'); shading flat; colorbar; colormap jet; title('latbin 20');
figure(2); pcolor(iaDate,f(iaChan),co2_mean(:,iaChan)'); shading flat; colorbar; colormap jet; title('latbin trop mean');
figure(3); pcolor(iaDate,f(iaChan),co2_std(:,iaChan)'); shading flat; colorbar; colormap jet; title('latbin trop std');
for ii = 1 : 3
  figure(ii); xlabel('time'); ylabel('Jac dBT/dCO2')
end
