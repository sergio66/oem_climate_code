for ii = 1 : 365
  x = ['OutputAnomaly_CAL/27/anomtest_timestep' num2str(ii) '.mat'];
  if exist(x)
    xx = load(x);
    fit_anom(:,ii) = xx.oem.finalrates(xx.jacobian.temp_i);
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
figure(1); pcolor(era_anom.p16anomaly.ptempanomaly); shading flat; colorbar; colormap jet; caxis([-5 5])
junk = era_anom.p16anomaly.ptempanomaly;
for ii = 1 : 20
  ix = xx.jacobian.wvjaclays_used{ii}-6;  %% 6 = xx.jacobian.scalar_i below
  if max(ix) <= 98
    era_ptempanom(ii,:) = mean(junk(ix,:));
  end
end
figure(1); pcolor(era_ptempanom); shading flat; colorbar; colormap jet; caxis([-5 5])
title('ERA data anomaly')

whos era_ptempanom fit_anom
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
era_ptempanom = zeros(40,20,365);
for iLat = 1 : 40
  era_model_file = ['/asl/s1/sergio/home/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Aug20_2019_Clr/'];
  era_model_file = [era_model_file '/Desc/16DayAvgNoS_withanom/latbin' num2str(iLat) '_16day_avg.rp.mat'];
  
  if exist(era_model_file)
    era_anom = load(era_model_file);
    figure(1); pcolor(era_anom.p16anomaly.ptempanomaly); shading flat; colorbar; colormap jet; caxis([-5 5])
    junk = era_anom.p16anomaly.ptempanomaly;
    for ii = 1 : 20
      ix = xx.jacobian.wvjaclays_used{ii}-6;  %% 6 = xx.jacobian.scalar_i below
      if max(ix) <= 98
        wah = mean(junk(ix,:));
        era_ptempanom(iLat,ii,1:length(wah)) = wah;
      end
    end
    figure(2); pcolor(squeeze(era_ptempanom(iLat,:,:))); shading flat; colorbar; colormap jet; caxis([-5 5])
    title(['ERA data anomaly iLat = ' num2str(iLat)])
    pause(0.1)
  else
    fprintf(1,'%s DNE \n',era_model_file)
  end
end

comment = 'see compare_era_anomaly_from_fit_and_model.m';
%save era_ptempanom.mat era_ptempanom comment
