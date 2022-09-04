iNumYears = 19
iNumYears = 12

if iNumYears == 12
  load iType_5_convert_sergio_clearskygrid_obsonly_Q16.mat
elseif iNumYears == 19
  load iType_4_convert_sergio_clearskygrid_obsonly_Q16.mat
end

for jj = 1 : 64
  fprintf(1,'latbin %2i of 64 \n',jj);
  if iNumYears == 12
    jacname = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_12yr_' num2str(jj) '.mat'];
  elseif iNumYears == 19
    jacname = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_Dec2021_' num2str(jj) '.mat'];
  end
  for ii = 1 : 72
    iibin = (jj-1)*72 + ii;
    [m_ts_jac0,nlays,qrenorm,freq2645]  = get_jac_fast(jacname,iibin,ii,jj,2021);
    co2 = m_ts_jac0(:,1);
    n20 = m_ts_jac0(:,2);
    ch4 = m_ts_jac0(:,3);
    cfc11 = m_ts_jac0(:,4);
    cfc12 = m_ts_jac0(:,5);
    newrad(:,iibin) = squeeze(b_desc(ii,jj,:)) - co2 - ch4 - n20 - cfc11 - cfc12;
    nedt_b(:,iibin) = squeeze(b_err_desc(ii,jj,:));
    nedt_true(:,iibin) = squeeze(airs_noiseTtrue(ii,jj,:));
    %plot(1:2645,newrad(:,iibin),1:2645,squeeze(b_desc(ii,jj,:)),1:2645,co2)
  end
end

oldrad = reshape(permute(b_desc,[3 1 2]),2645,72*64);

comment = 'see driver_remove_CO2_CH4_N20.m';
if iNumYears == 12
  save no_tracegas_spectral_rate_12years.mat newrad nedt_b nedt_true oldrad comment freq2645
elseif iNumYears == 19
  save no_tracegas_spectral_rate_19years.mat newrad nedt_b nedt_true oldrad comment freq2645
end

figure(1);
pcolor(reshape(oldrad(50,:)-newrad(50,:),72,64)'); colormap(usa2); colorbar; shading interp; caxis([-1 +1]*0.15)

figure(2)
plot(freq2645,nanmean(oldrad,2),freq2645,nanmean(newrad,2))
plotaxis2; xlim([645 1640])

