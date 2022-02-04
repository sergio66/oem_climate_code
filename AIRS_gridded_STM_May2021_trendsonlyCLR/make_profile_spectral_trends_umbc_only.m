function umbc_spectral_trends = make_profile_spectral_trends_umbc_only(results,resultsWV,resultsT,resultsO3,pavg,iERA5orERAI);

%% inputs 
%%   results*     = 6 scalars, 20 or 25 or 49 or 97 layers
%%
%% note "mask" is not part of argument list, this routine computes spectral rates for all 4608 profiles
%% you apply the mask to the results externally after this routine has been called

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

umbc_spectral_rates     = zeros(2645,4608);

[h1_4608,~,p1_4608,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp'); %% line 23 of see_clust_put_together_jacs_clr.m
[~,kcarta.subjac.ppmv2] = layers2ppmv(h1_4608,p1_4608,1:4608,2);
[~,kcarta.subjac.ppmv4] = layers2ppmv(h1_4608,p1_4608,1:4608,4);
[~,kcarta.subjac.ppmv6] = layers2ppmv(h1_4608,p1_4608,1:4608,6);

[mmUMBC,nnUMBC] = size(resultsWV);    
if nnUMBC ~= 20
  fprintf(1,' <<<<<<<<<<< WARNING make_profile_spectral_trends_umbc_only.m has length of WV,T,O3 retrievals as %3i and not 20 \n',nnUMBC);
end

for ii = 1 : 64
  if mod(ii,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end

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

  Qlevs = pavg;
  Tlevs = pavg;
  umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCO2 .* (ones(2645,1)*results(ind,1)'/2.2);
  umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.coljacN2O .* (ones(2645,1)*results(ind,2)'/1.0);
  umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCH4 .* (ones(2645,1)*results(ind,3)'/5.0);
  umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC11 * (0);
  umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC12 * (0);
  umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + m_ts_jac.subjac.jacST .* (ones(2645,1)*results(ind,6)'/0.1);
    umbc_20_layertrends.stemp(ind) = results(ind,6);
  xjunkrate = resultsWV(ind,:); clear junkrate
    junkrate = xjunkrate; 
    junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : nnUMBC; umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + squeeze(m_ts_jac.subjac.jacWV(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
    umbc_20_layertrends.gas_3(:,ind) = junkrate;
  xjunkrate = resultsT(ind,:); clear junkrate
    junkrate = xjunkrate; 
    junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : nnUMBC; umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + squeeze(m_ts_jac.subjac.jacT(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
    umbc_20_layertrends.gas_3(:,ind) = junkrate;
  xjunkrate = resultsO3(ind,:); clear junkrate
    junkrate = xjunkrate;
    junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : nnUMBC; umbc_spectral_rates(:,ind) = umbc_spectral_rates(:,ind) + squeeze(m_ts_jac.subjac.jacO3(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
    umbc_20_layertrends.gas_3(:,ind) = junkrate;
end

fprintf(1,'\n')

umbc_spectral_trends.umbc_spectral_rates = umbc_spectral_rates;
umbc_spectral_trends.umbc_20_layertrends = umbc_20_layertrends;


