iTorO3 = -1;  %% O3
iTorO3 = 0;   %% WV
iTorO3 = +1;  %% T

for ii = 1 : 157
  x = ['OutputAnomaly_CAL/27/anomtest_timestep' num2str(ii) '.mat'];
  if exist(x)
    xx = load(x);
    if iTorO3 > 0
      fit_anom(:,ii) = xx.oem.finalrates(xx.jacobian.temp_i);
    elseif iTorO3 < 0 
      fit_anom(:,ii) = xx.oem.finalrates(xx.jacobian.ozone_i);
    elseif iTorO3 == 0
      fit_anom(:,ii) = xx.oem.finalrates(xx.jacobian.water_i);
    end
  end
end
figure(2); pcolor(fit_anom); shading flat; colorbar; colormap jet
caxis([-5 5])
title('oem fits to SARTA cals')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
era_model_file = ['/asl/s1/sergio/home/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Aug20_2019_Clr/'];
era_model_file = [era_model_file '/Desc/16DayAvgNoS_withanom/latbin27_16day_avg.rp.mat'];

era_anom = load(era_model_file);
if iTorO3 > 0
  figure(1); pcolor(era_anom.p16anomaly.ptempanomaly); shading flat; colorbar; colormap jet; caxis([-5 5])
  junk = era_anom.p16anomaly.ptempanomaly;  
elseif iTorO3 < 0
  figure(1); pcolor(era_anom.p16anomaly.gas_3anomalyfrac); shading flat; colorbar; colormap jet; caxis([-5 5])
  junk = era_anom.p16anomaly.gas_3anomalyfrac;  
elseif iTorO3 == 0
  figure(1); pcolor(era_anom.p16anomaly.gas_1anomalyfrac); shading flat; colorbar; colormap jet; caxis([-5 5])
  junk = era_anom.p16anomaly.gas_1anomalyfrac;  
end
for ii = 1 : 20
  if iTorO3 > 0
    ix = xx.jacobian.wvjaclays_used{ii}-6;  %% 6 = xx.jacobian.scalar_i below
  elseif iTorO3 < 0
    ix = xx.jacobian.wvjaclays_used{ii}-6;  %% 6 = xx.jacobian.scalar_i below, and then add T also
  elseif iTorO3 == 0
    ix = xx.jacobian.wvjaclays_used{ii}-6;  %% 6 = xx.jacobian.scalar_i below, and then add T+WV also
  end
  if max(ix) <= 98
    era_ptempanom(ii,:) = mean(junk(ix,:));
  end
end
figure(1); pcolor(era_ptempanom); shading flat; colorbar; colormap jet; caxis([-5 5])
title('ERA data anomaly')

whos era_ptempanom fit_anom
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iTorO3 > 0
  era_ptempanom = zeros(40,20,157);
elseif iTorO3 < 0
  era_gas_3anom = zeros(40,20,157);
elseif iTorO3 == 0
  era_gas_1anom = zeros(40,20,157);
end

for iLat = 1 : 40
  era_model_file = ['/asl/s1/sergio/home/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/CLO_LAT40_avg_made_Dec2019_Clr/'];
  era_model_file = [era_model_file '/Desc/16DayAvgNoS/latbin' num2str(iLat) '_16day_avg.rp.mat'];
  
  if exist(era_model_file)
    era_anom = load(era_model_file);
    if iTorO3 > 0
      figure(1); pcolor(era_anom.p16anomaly.ptempanomaly); shading flat; colorbar; colormap jet; caxis([-5 5])
      junk = era_anom.p16anomaly.ptempanomaly;  
    elseif iTorO3 < 0
      figure(1); pcolor(era_anom.p16anomaly.gas_3anomalyfrac); shading flat; colorbar; colormap jet; caxis([-1 1])
      junk = era_anom.p16anomaly.gas_3anomalyfrac;  
    elseif iTorO3 == 0
      figure(1); pcolor(era_anom.p16anomaly.gas_1anomalyfrac); shading flat; colorbar; colormap jet; caxis([-1 1])
      junk = era_anom.p16anomaly.gas_1anomalyfrac;  
    end
    for ii = 1 : 20
      if iTorO3 > 0
        ix = xx.jacobian.wvjaclays_used{ii}-6;  %% 6 = xx.jacobian.scalar_i below
      elseif iTorO3 < 0
        ix = xx.jacobian.wvjaclays_used{ii}-6;  %% 6 = xx.jacobian.scalar_i below
      elseif iTorO3 == 0
        ix = xx.jacobian.wvjaclays_used{ii}-6;  %% 6 = xx.jacobian.scalar_i below
      end
      if max(ix) <= 98
        wah = mean(junk(ix,:));
        if iTorO3 > 0
          era_ptempanom(iLat,ii,1:length(wah)) = wah;
        elseif iTorO3 < 0
          era_gas_3anom(iLat,ii,1:length(wah)) = wah;
        elseif iTorO3 == 0
          era_gas_1anom(iLat,ii,1:length(wah)) = wah;
        end
      end
    end
    if iTorO3 > 0
      figure(2); pcolor(squeeze(era_ptempanom(iLat,:,:))); shading flat; colorbar; colormap jet; caxis([-5 5])
    elseif iTorO3 < 0
      figure(2); pcolor(squeeze(era_gas_3anom(iLat,:,:))); shading flat; colorbar; colormap jet; caxis([-1 1])
    elseif iTorO3 == 0
      figure(2); pcolor(squeeze(era_gas_1anom(iLat,:,:))); shading flat; colorbar; colormap jet; caxis([-1 1])
    end
    title(['ERA data anomaly iLat = ' num2str(iLat)])
    pause(0.1)
  else
    fprintf(1,'%s DNE \n',era_model_file)
  end
end

%{
comment = 'see compare_era_anomaly_from_fit_and_model.m';
if iTorO3 > 0
  save era_ptempanom.mat era_ptempanom comment
elseif iTorO3 < 0
  save era_gas_3anom.mat era_gas_3anom comment
elseif iTorO3 == 0
  save era_gas_1anom.mat era_gas_1anom comment
end
%}
