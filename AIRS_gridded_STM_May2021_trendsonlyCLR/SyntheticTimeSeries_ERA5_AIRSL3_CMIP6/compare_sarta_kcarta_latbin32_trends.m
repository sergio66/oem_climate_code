for JOB = 1 : 72
  x = load(['KCARTA_latbin32_spectral_trends/kcarta_spectraltrends_latbin32_lonbin' num2str(JOB,'%02d') '.mat']);
  kctrendALL(:,JOB) = x.kctrend;
  fKc = x.fKc;
end

%% see eg plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
JOB = 32;
iLatBin = 32;
dirout = '../../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries';
fsarta = [dirout '/reconstruct_era5_spectra_geo_rlat' num2str(JOB,'%02i') '.mat'];
x = load(fsarta);
sartatrend = x.thesave.xtrendSpectral;

plot(fKc,nanmean(kctrendALL'),fKc,nanmean(kctrendALL'-sartatrend'),fKc,nanstd(kctrendALL'-sartatrend')); xlim([640 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/
for iLonBin = 1 : 72
  jacname = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/usethisjac_clear_reconstructcode_ERA5_2021_latbin_' num2str(iLatBin,'%02d') '_lonbin_' num2str(iLonBin,'%02d') '.mat'];
  boo = load(jacname);
  [mm,nn] = size(boo.m_ts_jac);
  jac(iLonBin,:,1:nn) = boo.m_ts_jac;
  iaNlays(iLonBin) = boo.nlays;
end

nwp_trends = load('../../AIRS_gridded_STM_May2021_trendsonlyCLR/nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat');
raaReconstruct = zeros(2645,72);
for iLonBin = 1 : 72
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,1))'*2.2/2.2;
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,2))'*0.8/1;
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,3))'*4.75/5;
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,6))'*nwp_trends.era5_100_layertrends.stemp((iLatBin-1)*72+iLonBin)/0.1;

  %% look at eg /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/make_profile_spectral_trends.m
  ix = (1:iaNlays(iLonBin));  
  ixx1 = 6 + ix + 0*iaNlays(iLonBin); wah = squeeze(jac(iLonBin,:,ixx1)); bah = ones(2645,1) * nwp_trends.era5_100_layertrends.gas_1(ix,(iLatBin-1)*72+iLonBin)'/0.01; 
    raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + sum(wah.*bah,2);
  ixxT = 6 + ix + 1*iaNlays(iLonBin); wah = squeeze(jac(iLonBin,:,ixxT)); bah = ones(2645,1) * nwp_trends.era5_100_layertrends.ptemp(ix,(iLatBin-1)*72+iLonBin)'/0.01; 
    raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + sum(wah.*bah,2);
  ixx3 = 6 + ix + 2*iaNlays(iLonBin); wah = squeeze(jac(iLonBin,:,ixx3)); bah = ones(2645,1) * nwp_trends.era5_100_layertrends.gas_3(ix,(iLatBin-1)*72+iLonBin)'/0.01; 
    raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + sum(wah.*bah,2);
end

ind = (1:72) + (iLatBin-1)*72;

plot(nwp_trends.era5_100_layertrends.stemp((57*72)+(1:72))); plotaxis2;
plot(nwp_trends.era5_spectral_rates(1520,(57*72)+(1:72))); plotaxis2;
figure(1); pcolor(reshape(nwp_trends.era5_100_layertrends.stemp,72,64)'); colormap(usa2); colorbar; shading flat; caxis([-0.15 +0.15]); title('Stemp rates');
figure(2); pcolor(reshape(nwp_trends.era5_spectral_rates(1520,:),72,64)'); colormap(usa2); colorbar; shading flat; caxis([-0.15 +0.15]); title('Orig BT1231 rates')

figure(3);
ijunk = (iLatBin-1)*72 + (1:72);
sum(ind - ijunk)
plot(1:72,nwp_trends.era5_spectral_rates(1520,ind),1:72,nwp_trends.era5_100_layertrends.stemp(ind),1:72,raaReconstruct(1520,:)); plotaxis2; title(num2str(iLatBin))
  hl = legend('Orig reconstruct 1231 cm-1','ERA5 stemp','New reconstruct 1231 cm-1','location','best');

previousReconstruct = nwp_trends.era5_spectral_rates(:,ind);
plot(fKc,nanmean(raaReconstruct'),fKc,nanmean(kctrendALL'),fKc,nanmean(previousReconstruct')); plotaxis2;
  hl = legend('New reconstruct','kCARTA trend','Orig reconstruct','location','best');
plot(fKc,nanmean(raaReconstruct'),'x-',fKc,nanmean(kctrendALL'),'.-',fKc,nanmean(sartatrend'),fKc,nanmean(previousReconstruct'),fKc,ones(size(fKc))*mean(nwp_trends.era5_100_layertrends.stemp(ind))); plotaxis2;; title(['Latbin ' num2str(iLatBin)])
  xlim([min(fKc) max(fKc)]);  hl = legend('New reconstruct','kCARTA trend','SARTA trend','Orig reconstruct','Mean ERA5 stemp','location','best','fontsize',10); 
  xlim([640 1640]);  hl = legend('New reconstruct','kCARTA trend','SARTA trend','Orig reconstruct','Mean ERA5 stemp','location','best','fontsize',10); 
