%% copied from simple_get_the_model_trends_do_feedbacks.m
%% copied from simple_get_the_model_trends_do_feedbacks.m
%% copied from simple_get_the_model_trends_do_feedbacks.m

    junk2 = +1;
    if exist(feedbacknameUMBC)
      lser = ['!ls -lth ' feedbacknameUMBC];
      eval(lser);
      junk2 = input('file already exists, overwrite (-1 default/+1) : ');
      if length(junk2) == 0
        junk2 = -1;
      end
    end       
    saver = ['save ' feedbacknameUMBC ' umbc_spectral_olr results resultsWV resultsT resultsO3 pavg plays'];  %% if you only want to save UMBC
    if junk2 > 0
      fprintf(1,'saving to %s \n',feedbacknameUMBC);
      eval(saver);
    else
      fprintf(1,'this already exists %s not saving \n',feedbacknameUMBC);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isfield('nwp_spectral_trends_cmip6_era5_airsL3_umbc','airsL3_100_layertrends')
      %% this test is more because of deriver_quick_redo_regressions_umbc_feedbacks_timeseries.m
      aL3trend.stemp = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp;
      c6trend.stemp = cmip6.trend_stemp;
    end

      junk2 = +1;
      if exist(feedbacknameNWP_ERA5)
        lser = ['!ls -lth ' feedbacknameNWP_ERA5];
        eval(lser);
        junk2 = input('file already exists, overwrite (-1 default/+1) : ');
        if length(junk2) == 0
          junk2 = -1;
        end
      end       
      stemptrend.era5   = era5.trend_stemp;
      stemptrend.airsL3 = aL3trend.stemp;
      stemptrend.cmip6  = c6trend.stemp;
      saver = ['save ' feedbacknameNWP_ERA5 ' era5_spectral_olr cmip6_spectral_olr airsL3_spectral_olr pavg plays stemptrend'];  %% if you want to save models/NWP only
      if junk2 > 0
        fprintf(1,'saving to %s \n',feedbacknameNWP_ERA5);
        eval(saver);
      else
        fprintf(1,'this already exists %s not saving \n',feedbacknameNWP_ERA5);
      end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('might fail here if you do not have MERRA2/CLIMCAPSL3/AMIP6')
    if isfield('nwp_spectral_trends_amip6_merra2_climcapsL3_umbc','airsL3_100_layertrends')
      %% this test is more because of deriver_quick_redo_regressions_umbc_feedbacks_timeseries.m
      cL3trend.stemp = nwp_spectral_trends_amip6_merra2_climcapsL3_umbc.airsL3_100_layertrends.stemp;
      a6trend.stemp = amip6.trend_stemp;
    end
      junk2 = +1;
      if exist(feedbacknameNWP_MERRA2)
        lser = ['!ls -lth ' feedbacknameNWP_MERRA2];
        eval(lser);
        junk2 = input('file already exists, overwrite (-1 default/+1) : ');
        if length(junk2) == 0
          junk2 = -1;
        end
      end       
      stemptrend.era5   = merra2.trend_stemp;
      stemptrend.airsL3 = cL3trend.stemp;
      stemptrend.cmip6  = a6trend.stemp;
      saver = ['save ' feedbacknameNWP_MERRA2 ' merra2_spectral_olr amip6_spectral_olr climcapsL3_spectral_olr pavg plays stemptrend'];  %% if you want to save models/NWP only
      if junk2 > 0
        fprintf(1,'saving to %s \n',feedbacknameNWP_MERRA2);
        eval(saver);
      else
        fprintf(1,'this already exists %s not saving \n',feedbacknameNWP_MERRA2);
      end
