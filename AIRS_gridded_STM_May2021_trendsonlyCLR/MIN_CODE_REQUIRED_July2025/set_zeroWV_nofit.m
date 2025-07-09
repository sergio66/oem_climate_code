if exist('iFixWV_NoFit','var')
  %% from strow_override_defaults_latbins_AIRS_fewlays.m
  if iFixWV_NoFit >= 0 & (strcmp(driver.rateset.ocb_set,'obs') | strcmp(driver.rateset.ocb_set,'cal'))
    if iZeroWVVers == 0
      izname = ['SAVE_BESTRUNv1/OutputAnomaly_CAL/' num2str(driver.iibin,'%02d') '/anomtest_timestep' num2str(driver.i16daytimestep) '.mat'];
      if exist(izname)
        fprintf(1,'WARNING setting dWV(z,lat,t)/dt using CAL anom in %s \n',izname);       
        izt = load(izname);
        cal_WV_rates = izt.oem.finalrates(driver.jacobian.ozone_i);

        spectra_due_to_WV_jac = zeros(size(driver.rateset.rates));
        renorm_cal_WV_rates = cal_WV_rates./driver.qrenorm(driver.jacobian.ozone_i)';
        for iiWV = 1 : length(driver.jacobian.ozone_i)
          spectra_due_to_WV_jac = spectra_due_to_WV_jac + renorm_cal_WV_rates(iiT)*m_ts_jac(:,driver.jacobian.ozone_i(iiWV));
        end
        load f2645.mat
        plot(f2645,spectra_due_to_WV_jac,'b',f2645,driver.rateset.rates,'m.-',...
             f2645,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2645,driver.rateset.rates-spectra_due_to_WV_jac,'r.-')
        plot(f2645,driver.rateset.rates,'m.-',f2645,driver.rateset.rates-spectra_due_to_WV_jac,'r.-',f2645,m_ts_jac(:,1)*10+0.85,'k.-')
      else
        fprintf(1,'WARNING trying to setting dWV(z,lat,t)/dt using CAL anom but %s DNE so set xb(WV) = 0 \n',izname);
        cal_WV_rates = zeros(size(driver.jacobian.ozone_i));
        spectra_due_to_WV_jac = zeros(size(driver.rateset.rates));
        xb(driver.jacobian.ozone_i) = 0.0;
      end

    elseif iZeroWVVers == 1
      era_model_file = ['/asl/s1/sergio/home/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Aug20_2019_Clr/'];
      era_model_file = [era_model_file '/Desc/16DayAvgNoS_withanom/latbin' num2str(driver.iibin,'%02d') '_16day_avg.rp.mat'];
      izname = era_model_file;
      x = ['SAVE_BESTRUNv1/OutputAnomaly_CAL/27/anomtest_timestep' num2str(123) '.mat'];
      xx = load(x);

      if exist(izname)
        fprintf(1,'WARNING setting dWV(z,lat,t)/dt using raw ERA snom in %s \n',izname);       
        era_anom = load(era_model_file);
        figure(1); pcolor(era_anom.p16anomaly.gas_3anomaly); shading flat; colorbar; colormap jet; caxis([-5 5])
        junk = era_anom.p16anomaly.ozoneanomaly;
        for iiii = 1 : iNlays_retrieve   %% before we had 1 : 20
          ixix = xx.jacobian.wvjaclays_used{iiii}-6;  %% 6 = xx.jacobian.scalar_i below
          if max(ixix) <= 98
            era_ozoneanom(iiii,:) = mean(junk(ixix,:));
          end
        end        
        cal_ozone_rates = era_ozoneanom(:,driver.i16daytimestep);

        spectra_due_to_ozone_jac = zeros(size(driver.rateset.rates));
        renorm_cal_WV_rates = cal_WV_rates./driver.qrenorm(driver.jacobian.ozone_i)';
        for iiWV = 1 : length(driver.jacobian.ozone_i)
          spectra_due_to_WV_jac = spectra_due_to_WV_jac + renorm_cal_WV_rates(iiT)*m_ts_jac(:,driver.jacobian.ozone_i(iiWV));
        end
        load f2645.mat
        plot(f2645,spectra_due_to_WV_jac,'b',f2645,driver.rateset.rates,'m.-',...
             f2645,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2645,driver.rateset.rates-spectra_due_to_WV_jac,'r.-')
        plot(f2645,driver.rateset.rates,'m.-',f2645,driver.rateset.rates-spectra_due_to_WV_jac,'r.-',f2645,m_ts_jac(:,1)*10+0.85,'k.-')

      else
        fprintf(1,'WARNING trying to setting dWV(z,lat,t)/dt using CAL anom but %s DNE so set xb(WV) = 0 \n',izname);
        cal_WV_rates = zeros(size(driver.jacobian.ozone_i));
        spectra_due_to_WV_jac = zeros(size(driver.rateset.rates));
        xb(driver.jacobian.ozone_i) = 0.0;
      end

    elseif iZeroWVVers == 2
      era_model_file = 'era_gas_3anom.mat';
      izname = era_model_file;
      if exist(izname)
        fprintf(1,'WARNING setting dWV(z,lat,t)/dt using raw ERA anom in %s \n',izname);       
        era_anom = load(era_model_file);
        junk = era_anom.era_gas_3anom;
        era_ozoneanom = squeeze(junk(driver.iibin,:,:));
        cal_WV_rates = era_ozoneanom(:,driver.i16daytimestep);

        spectra_due_to_WV_jac = zeros(size(driver.rateset.rates));
        renorm_cal_WV_rates = cal_WV_rates./driver.qrenorm(driver.jacobian.ozone_i)';
        for iiT = 1 : length(driver.jacobian.ozone_i)
          spectra_due_to_WV_jac = spectra_due_to_WV_jac + renorm_cal_WV_rates(iiT)*m_ts_jac(:,driver.jacobian.ozone_i(iiT));
        end
        load f2645.mat
        plot(f2645,spectra_due_to_WV_jac,'b',f2645,driver.rateset.rates,'m.-',...
             f2645,sum(m_ts_jac(:,driver.jacobian.ozone_i)')*100,'k.-',f2645,driver.rateset.rates-spectra_due_to_WV_jac,'r.-')
        plot(f2645,driver.rateset.rates,'m.-',f2645,driver.rateset.rates-spectra_due_to_WV_jac,'r.-',f2645,m_ts_jac(:,1)*10+0.85,'k.-')

      else
        fprintf(1,'WARNING trying to setting dT(z,lat,t)/dt using CAL anom but %s DNE so set xb(T) = 0 \n',izname);
        cal_WV_rates = zeros(size(driver.jacobian.ozone_i));
        spectra_due_to_WV_jac = zeros(size(driver.rateset.rates));
        xb(driver.jacobian.ozone_i) = 0.0;
      end

    end  %% if iZeroWVVers == 0 ! if iZeroWVVers == 1 | if iZeroWVVers == 2
    if iFixWV_NoFit == 0
      disp('hmm in SW no need to worry about WV, zero all this WV stuff ...')
      spectra_due_to_WV_jac = zeros(size(spectra_due_to_WV_jac));
      cal_WV_rates = zeros(size(cal_WV_rates));
    end 
  end    %% if iFixWV_NoFit > 0 & strcmp(driver.rateset.ocb_set,'obs')
end      %% if exist('iFixWV_NoFit','var')
