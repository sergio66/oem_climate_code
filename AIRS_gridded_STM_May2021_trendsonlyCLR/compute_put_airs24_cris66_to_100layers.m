for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end


  ind = (ii-1)*72 + (1:72);

  %%%%%%%%%%%%%%%%%%%%%%%%%

  %% note we use 1/co2ppmv,1/n2oppb,1/ch4ppb instead of 2.2/co2ppmv,0.8/n2oppb,5/ch4ppb -- because we are really doing a jacobian normalized by 2.2 ppmv change etc etc
  %% see get_jac_fast.m

  %% remember first arg is either ERA or CMIP6 rates
  junkrate = era.trend_stemp(ind);  
    era_100_layertrends.stemp(ind) = junkrate;
  junkrate = era.trend_ptemp(:,ind); junkrate(isnan(junkrate)) = 0; 
    era_100_layertrends.ptemp(:,ind) = junkrate;
  junkrate = era.trend_gas_1(:,ind); junkrate(isnan(junkrate)) = 0; 
    era_100_layertrends.gas_1(:,ind) = junkrate;
  junkrate = era.trend_gas_3(:,ind); junkrate(isnan(junkrate)) = 0; 
    era_100_layertrends.gas_3(:,ind) = junkrate;

  if iERAorCMIP6 > 0
    % this basically works
    % era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCO2 * 1.0./m_ts_jac.subjac.ppmv2;
    % era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacN2O * 1.0./(m_ts_jac.subjac.ppmv4*1000);
    % era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCH4 * 1.0./(m_ts_jac.subjac.ppmv6*1000) *5.0/4.5;
    % era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC11 * (-1.0./300);
    % era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC12 * (-1.0/600);
  
    % so does this
    junkrate = era5.trend_stemp(ind);  
      era5_100_layertrends.stemp(ind) = junkrate;
    junkrate = era5.trend_ptemp(:,ind); junkrate(isnan(junkrate)) = 0; 
      era5_100_layertrends.ptemp(:,ind) = junkrate;
    junkrate = era5.trend_gas_1(:,ind); junkrate(isnan(junkrate)) = 0; 
      era5_100_layertrends.gas_1(:,ind) = junkrate;
    junkrate = era5.trend_gas_3(:,ind); junkrate(isnan(junkrate)) = 0; 
      era5_100_layertrends.gas_3(:,ind) = junkrate;
  
    %%%%%%%%%%%%%%%%%%%%%%%%%
    Qlevs = airsL3.Qlevs;
    Tlevs = airsL3.Tlevs;
    [~,~,OZlevs] = size(airsL3.thestats64x72.ozonerate);
    if OZlevs == length(Qlevs)
      OZlevs = Qlevs;
    elseif OZlevs == length(Tlevs)
      OZlevs = Tlevs;
    end
    junkrate = airsL3.thestats64x72.stemprate(:,ii)';   
      airsL3_100_layertrends.stemp(ind) = junkrate;
    xjunkrate = airsL3.thestats64x72.ptemprate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      airsL3_100_layertrends.ptemp(:,ind) = junkrate;
    xjunkrate = airsL3.thestats64x72.waterrate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Qlevs) | plays >= max(Qlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      airsL3_100_layertrends.gas_1(:,ind) = junkrate;
    xjunkrate = airsL3.thestats64x72.ozonerate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(OZlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(OZlevs) | plays >= max(OZlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      airsL3_100_layertrends.gas_3(:,ind) = junkrate;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    Qlevs = pavg;
    Tlevs = pavg;
      umbc_20_layertrends.stemp(ind) = results(ind,6);
    xjunkrate = resultsT(ind,:); clear junkrate
      junkrate = xjunkrate; 
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      umbc_20_layertrends.ptemp(:,ind) = junkrate;
    xjunkrate = resultsWV(ind,:); clear junkrate
      junkrate = xjunkrate; 
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      umbc_20_layertrends.gas_1(:,ind) = junkrate;
    xjunkrate = resultsO3(ind,:); clear junkrate
      junkrate = xjunkrate;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; 
      umbc_20_layertrends.gas_3(:,ind) = junkrate;
  end
end
