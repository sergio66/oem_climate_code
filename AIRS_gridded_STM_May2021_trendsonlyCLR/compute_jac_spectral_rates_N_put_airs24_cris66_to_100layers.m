for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end

% really not needed
%  jacfile = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//kcarta_clr_subjacLatBin_newSARTA_' num2str(ii,'%02d') '.mat'];
%  jac = load(jacfile);

  for jjj = 1 : 72
    if iERA5orERAI == 2019
      jacfilex = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//usethisjac_clear_reconstructcode_latbin_' num2str(ii,'%02d') '_lonbin_' num2str(jjj,'%02d') '.mat'];
    elseif iERA5orERAI == 2021
      jacfilex = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//usethisjac_clear_reconstructcode_ERA5_2021_latbin_' num2str(ii,'%02d') '_lonbin_' num2str(jjj,'%02d') '.mat'];
    end

    jacx = load(jacfilex);    
    %% THESE ARE NORMALIZED JACS   JTRUE*RENORM
    m_ts_jac.subjac.coljacCO2(:,jjj)   = jacx.m_ts_jac(:,1);
    m_ts_jac.subjac.coljacN2O(:,jjj)   = jacx.m_ts_jac(:,2);
    m_ts_jac.subjac.coljacCH4(:,jjj)   = jacx.m_ts_jac(:,3);
    m_ts_jac.subjac.coljacCFC11(:,jjj) = jacx.m_ts_jac(:,4);
    m_ts_jac.subjac.coljacCFC12(:,jjj) = jacx.m_ts_jac(:,5);
    m_ts_jac.subjac.jacST(:,jjj)       = jacx.m_ts_jac(:,6);
    wah = ((1:jacx.nlays)+6+jacx.nlays*0); boo = zeros(100,2645); boo(1:jacx.nlays,:) =  jacx.m_ts_jac(:,wah)';  x100m_ts_jac.subjac.jacWV(:,:,jjj) = boo;
    wah = ((1:jacx.nlays)+6+jacx.nlays*1); boo = zeros(100,2645); boo(1:jacx.nlays,:) =  jacx.m_ts_jac(:,wah)';  x100m_ts_jac.subjac.jacT(:,:,jjj)  = boo;
    wah = ((1:jacx.nlays)+6+jacx.nlays*2); boo = zeros(100,2645); boo(1:jacx.nlays,:) =  jacx.m_ts_jac(:,wah)';  x100m_ts_jac.subjac.jacO3(:,:,jjj) = boo;
  end
%  m_ts_jac.subjac.ppmv2 = jac.subjac.ppmv2;
%  m_ts_jac.subjac.ppmv4 = jac.subjac.ppmv4;
%  m_ts_jac.subjac.ppmv6 = jac.subjac.ppmv6;

  % forget the renorm
  m_ts_jac.subjac.jacST = m_ts_jac.subjac.jacST * 1;
  m_ts_jac.subjac.jacWV = quick_combinejaclays_make_profile_spectral_trends(x100m_ts_jac.subjac.jacWV * 1.00,1:100,nnUMBC);
  m_ts_jac.subjac.jacT  = quick_combinejaclays_make_profile_spectral_trends(x100m_ts_jac.subjac.jacT  * 1.00,1:100,nnUMBC);
  m_ts_jac.subjac.jacO3 = quick_combinejaclays_make_profile_spectral_trends(x100m_ts_jac.subjac.jacO3 * 1.00,1:100,nnUMBC);

  ind = (ii-1)*72 + (1:72);

  %%%%%%%%%%%%%%%%%%%%%%%%%

  %% note we use 1/co2ppmv,1/n2oppb,1/ch4ppb instead of 2.2/co2ppmv,0.8/n2oppb,5/ch4ppb -- because we are really doing a jacobian normalized by 2.2 ppmv change etc etc
  %% see get_jac_fast.m

  %% remember first arg is either ERA or CMIP6 rates
  era_spectral_rates(:,ind) = era_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCO2;
  era_spectral_rates(:,ind) = era_spectral_rates(:,ind) + m_ts_jac.subjac.coljacN2O;
  era_spectral_rates(:,ind) = era_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCH4;
  era_spectral_rates(:,ind) = era_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC11*0;
  era_spectral_rates(:,ind) = era_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC12*0;
  junkrate = era.trend_stemp(ind);  era_spectral_rates(:,ind) = era_spectral_rates(:,ind) + m_ts_jac.subjac.jacST .* (ones(2645,1)*junkrate/0.1);
    era_100_layertrends.stemp(ind) = junkrate;
  junkrate = era.trend_ptemp(:,ind); junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100; era_spectral_rates(:,ind) = era_spectral_rates(:,ind) + squeeze(x100m_ts_jac.subjac.jacT(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
    era_100_layertrends.ptemp(:,ind) = junkrate;
  junkrate = era.trend_gas_1(:,ind); junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100; era_spectral_rates(:,ind) = era_spectral_rates(:,ind) + squeeze(x100m_ts_jac.subjac.jacWV(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
    era_100_layertrends.gas_1(:,ind) = junkrate;
  junkrate = era.trend_gas_3(:,ind); junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100; era_spectral_rates(:,ind) = era_spectral_rates(:,ind) + squeeze(x100m_ts_jac.subjac.jacO3(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
    era_100_layertrends.gas_3(:,ind) = junkrate;

  if iERAorCMIP6 > 0
    % this basically works
    % era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCO2 * 1.0./m_ts_jac.subjac.ppmv2;
    % era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacN2O * 1.0./(m_ts_jac.subjac.ppmv4*1000);
    % era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCH4 * 1.0./(m_ts_jac.subjac.ppmv6*1000) *5.0/4.5;
    % era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC11 * (-1.0./300);
    % era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC12 * (-1.0/600);
  
    % so does this
    era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCO2;
    era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacN2O;
    era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCH4;
    era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC11*0;
    era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC12*0;
    junkrate = era5.trend_stemp(ind);  era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + m_ts_jac.subjac.jacST .* (ones(2645,1)*junkrate/0.1);
      era5_100_layertrends.stemp(ind) = junkrate;
    junkrate = era5.trend_ptemp(:,ind); junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100; era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + squeeze(x100m_ts_jac.subjac.jacT(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
      era5_100_layertrends.ptemp(:,ind) = junkrate;
    junkrate = era5.trend_gas_1(:,ind); junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100; era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + squeeze(x100m_ts_jac.subjac.jacWV(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
      era5_100_layertrends.gas_1(:,ind) = junkrate;
    junkrate = era5.trend_gas_3(:,ind); junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100; era5_spectral_rates(:,ind) = era5_spectral_rates(:,ind) + squeeze(x100m_ts_jac.subjac.jacO3(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
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
    airsL3_spectral_rates(:,ind) = airsL3_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCO2;
    airsL3_spectral_rates(:,ind) = airsL3_spectral_rates(:,ind) + m_ts_jac.subjac.coljacN2O;
    airsL3_spectral_rates(:,ind) = airsL3_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCH4;
    airsL3_spectral_rates(:,ind) = airsL3_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC11*0;
    airsL3_spectral_rates(:,ind) = airsL3_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC12*0;
    junkrate = airsL3.thestats64x72.stemprate(:,ii)';   airsL3_spectral_rates(:,ind) = airsL3_spectral_rates(:,ind) + m_ts_jac.subjac.jacST .* (ones(2645,1)*junkrate/0.1);
      airsL3_100_layertrends.stemp(ind) = junkrate;
    xjunkrate = airsL3.thestats64x72.ptemprate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Tlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Tlevs) | plays >= max(Tlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100; airsL3_spectral_rates(:,ind) = airsL3_spectral_rates(:,ind) + squeeze(x100m_ts_jac.subjac.jacT(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
      airsL3_100_layertrends.ptemp(:,ind) = junkrate;
    xjunkrate = airsL3.thestats64x72.waterrate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(Qlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(Qlevs) | plays >= max(Qlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100; airsL3_spectral_rates(:,ind) = airsL3_spectral_rates(:,ind) + squeeze(x100m_ts_jac.subjac.jacWV(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
      airsL3_100_layertrends.gas_1(:,ind) = junkrate;
    xjunkrate = airsL3.thestats64x72.ozonerate(:,ii,:); xjunkrate = squeeze(xjunkrate); clear junkrate
      for jjj = 1 : 72; junkrate(jjj,:) = interp1(log(OZlevs),xjunkrate(jjj,:),log(plays),[],'extrap'); end; jjj = find(plays <= min(OZlevs) | plays >= max(OZlevs)); junkrate(:,jjj) = 0;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100; airsL3_spectral_rates(:,ind) = airsL3_spectral_rates(:,ind) + squeeze(x100m_ts_jac.subjac.jacO3(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
      airsL3_100_layertrends.gas_3(:,ind) = junkrate;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%  
    Qlevs = pavg;
    Tlevs = pavg;
    umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCO2 .* (ones(2645,1)*results(ind,1)'/2.2);
    umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.coljacN2O .* (ones(2645,1)*results(ind,2)'/1.0);
    umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCH4 .* (ones(2645,1)*results(ind,3)'/5.0);
    umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC11 * (0);
    umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC12 * (0);
    umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.jacST .* (ones(2645,1)*results(ind,6)'/0.1);
      umbc_20_layertrends.stemp(ind) = results(ind,6);
    xjunkrate = resultsT(ind,:); clear junkrate
      junkrate = xjunkrate; 
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : nnUMBC; umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + squeeze(m_ts_jac.subjac.jacT(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
      umbc_20_layertrends.ptemp(:,ind) = junkrate;
    xjunkrate = resultsWV(ind,:); clear junkrate
      junkrate = xjunkrate; 
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : nnUMBC; umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + squeeze(m_ts_jac.subjac.jacWV(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
      umbc_20_layertrends.gas_1(:,ind) = junkrate;
    xjunkrate = resultsO3(ind,:); clear junkrate
      junkrate = xjunkrate;
      junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : nnUMBC; umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + squeeze(m_ts_jac.subjac.jacO3(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
      umbc_20_layertrends.gas_3(:,ind) = junkrate;
  end
end
