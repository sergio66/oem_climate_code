%% see eg plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
iLatBin = 32;
iLatBin = input('Enter latbin (1:64), default = 32 : ');
if length(iLatBin) == 0
  iLatBin = 32;
end

JOB = iLatBin;
dirout = '../../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries';
fsarta = [dirout '/reconstruct_era5_spectra_geo_rlat' num2str(JOB,'%02i') '.mat'];
x = load(fsarta);
sartatrend = x.thesave.xtrendSpectral;

for JOB = 1 : 72
  x = load(['KCARTA_latbin32_spectral_trends/kcarta_spectraltrends_latbin32_lonbin' num2str(JOB,'%02d') '.mat']);
  kctrendALL(:,JOB) = x.kctrend;
  fKc = x.fKc;
end

if iLatBin == 32
  plot(fKc,nanmean(kctrendALL'),fKc,nanmean(kctrendALL'-sartatrend'),fKc,nanstd(kctrendALL'-sartatrend')); xlim([640 1640])
  hl = legend('kCARTA trends','mean(kcata-sarta) trends','std(kcarta-sarta) trends','location','best','fontsize',8);
  disp('ret to continue'); pause
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/
for iLonBin = 1 : 72
  jacname = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/usethisjac_clear_reconstructcode_ERA5_2021_latbin_' num2str(iLatBin,'%02d') '_lonbin_' num2str(iLonBin,'%02d') '.mat'];
  boo = load(jacname);
  [mm,nn] = size(boo.m_ts_jac);
  jac(iLonBin,:,1:nn) = boo.m_ts_jac;
  iaNlays(iLonBin) = boo.nlays;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nwp_trends0 = load('../../AIRS_gridded_STM_May2021_trendsonlyCLR/nwp_spectral_trends_cmip6_era5_airsL3_umbc.mat');          %% wrong jacs! used log10
nwp_trends1 = load('../../AIRS_gridded_STM_May2021_trendsonlyCLR/nwp_spectral_trends_cmip6_era5_airsL3_umbc_fixedjac.mat'); %% correct jacs! used ln
nwp_trends = nwp_trends1;

for JOB = 1 : 64
  dirout = '../../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries';
  fsarta = [dirout '/reconstruct_era5_spectra_geo_rlat' num2str(JOB,'%02i') '.mat'];
  x = load(fsarta);
  sartatrend_all(:,(JOB-1)*72+(1:72)) = x.thesave.xtrendSpectral;
end 

%% ~/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/plot_ERA_ERA5_AIRSL3_CMIP6_trends.m  -- gives same rates as nwp_trends.era5_100_layertrends.stemp
era5   = load('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ERA5_atm_data_2002_09_to_2021_08_trends_desc.mat'); 

obs_rates_all = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_startwithERA5trends.mat','rates');
rlon_all = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_startwithERA5trends.mat','rlon');
rlat_all = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_startwithERA5trends.mat','rlat');
junk = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_startwithERA5trends.mat','results'); stemp_rate = junk.results(:,6);
[Rlat ,Rlon] = meshgrid(rlat_all.rlat,rlon_all.rlon);
plot(fKc,nanmean(obs_rates_all.rates'-sartatrend_all'),fKc,nanstd(obs_rates_all.rates'-sartatrend_all')); xlim([640 1640])
pcolor(Rlon',Rlat',reshape(obs_rates_all.rates(1520,:),72,64)'); caxis([-0.15 +0.15]); colormap(usa2); shading interp
pcolor(Rlon',Rlat',reshape(sartatrend_all(1520,:),72,64)'); caxis([-0.15 +0.15]); colormap(usa2); shading interp
figure(1); scatter_coast(Rlon(:),Rlat(:),50,obs_rates_all.rates(1520,:));            caxis([-0.15 +0.15]); colormap(usa2); shading interp; title('dBT1231/dt obs');
figure(2); scatter_coast(Rlon(:),Rlat(:),50,stemp_rate);                             caxis([-0.15 +0.15]); colormap(usa2); shading interp; title('retr dST/dt')
figure(3); scatter_coast(Rlon(:),Rlat(:),50,sartatrend_all(1520,:));                 caxis([-0.15 +0.15]); colormap(usa2); shading interp; title('ERA5-->SARTA dBT1231/dt')
figure(4); scatter_coast(Rlon(:),Rlat(:),50,nwp_trends.era5_100_layertrends.stemp);  caxis([-0.15 +0.15]); colormap(usa2); shading interp; title('ERA5 dST/dt')
disp('pause to continue');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raaReconstruct = zeros(2645,72);

%% I used these in driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m
%co2ppm = 368 + 2.1*time_so_far;
%n2oppm = 315  + (332-315)/(2020-2000)*time_so_far; n2oppm = n2oppm/1000;
%ch4ppm = 1.75 + (1.875-1.750)/(2020-2000)*time_so_far;

for iLonBin = 1 : 72
  %{
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,1))'*2.2/2.2;
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,2))'*0.8/1;
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,3))'*4.75/5;
  %}

  %{
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,1))'*2.2/2.2;
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,2))'*1/1;
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,3))'*5/5;
  %}

  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,1))'*2.1/2.2;
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,2))'*(332-315)/(2020-2000)/1;
  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,3))'*(1.875-1.750)/(2020-2000)*1000/5;

  raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + squeeze(jac(iLonBin,:,6))'*nwp_trends.era5_100_layertrends.stemp((iLatBin-1)*72+iLonBin)/0.1;

  %% look at eg /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/make_profile_spectral_trends.m
  ix = (1:iaNlays(iLonBin));  
  gas_fudge = 0.5;
  gas_fudge = 1/log(10);  %% this is before, when I had rescaled the jacs to be log10
  gas_fudge = 1.0;        %% this is after,  when I found that silly mistake!
  ixx1 = 6 + ix + 0*iaNlays(iLonBin); wah = squeeze(jac(iLonBin,:,ixx1)); bah = ones(2645,1) * nwp_trends.era5_100_layertrends.gas_1(ix,(iLatBin-1)*72+iLonBin)'/0.01*gas_fudge; 
    raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + sum(wah.*bah,2);
  ixxT = 6 + ix + 1*iaNlays(iLonBin); wah = squeeze(jac(iLonBin,:,ixxT)); bah = ones(2645,1) * nwp_trends.era5_100_layertrends.ptemp(ix,(iLatBin-1)*72+iLonBin)'/0.01; 
    raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + sum(wah.*bah,2);
  ixx3 = 6 + ix + 2*iaNlays(iLonBin); wah = squeeze(jac(iLonBin,:,ixx3)); bah = ones(2645,1) * nwp_trends.era5_100_layertrends.gas_3(ix,(iLatBin-1)*72+iLonBin)'/0.01*gas_fudge; 
    raaReconstruct(:,iLonBin) = raaReconstruct(:,iLonBin) + sum(wah.*bah,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = (1:72) + (iLatBin-1)*72;

plot(nwp_trends.era5_100_layertrends.stemp((57*72)+(1:72))); plotaxis2;
plot(nwp_trends.era5_spectral_rates(1520,(57*72)+(1:72))); plotaxis2;
figure(1); pcolor(reshape(nwp_trends.era5_100_layertrends.stemp,72,64)'); colormap(usa2); colorbar; shading flat; caxis([-0.15 +0.15]); title('Stemp rates');
figure(2); pcolor(reshape(nwp_trends.era5_spectral_rates(1520,:),72,64)'); colormap(usa2); colorbar; shading flat; caxis([-0.15 +0.15]); title('BT1231 rates')

figure(3);
ijunk = (iLatBin-1)*72 + (1:72);
%sum(ind - ijunk)
plot(1:72,nwp_trends.era5_spectral_rates(1520,ind),'.-',1:72,nwp_trends.era5_100_layertrends.stemp(ind),'o-',1:72,raaReconstruct(1520,:),1:72,sartatrend(1520,:),'x-','linewidth',2); plotaxis2; title(num2str(iLatBin))
  hl = legend('reconstruct 1231 cm-1 rate','ERA5 stemp rate','in place new reconstruct 1231 cm-1','sarta 1231 trend','location','best','fontsize',10);
plot(1:72,raaReconstruct(1520,:),'o--',1:72,nwp_trends.era5_spectral_rates(1520,ind),'x--',1:72,kctrendALL(1520,:),'o-',1:72,sartatrend(1520,:),'x-',1:72,nwp_trends.era5_100_layertrends.stemp(ind),'d-','linewidth',2); plotaxis2; title(num2str(iLatBin))
  hl = legend('In-situ jacobian reconstruct 1231 cm-1 rate','Previous jacobian (driver\_gather)  reconstruct 1231 cm-1','NWP->kcarta 1231 cm-1 rate','NWP->sarta 1231 cm-1 trend','ERA5 stemp rate','location','best','fontsize',10);

figure(4)
previousReconstruct = nwp_trends.era5_spectral_rates(:,ind);
plot(fKc,nanmean(raaReconstruct'),fKc,nanmean(kctrendALL'),fKc,nanmean(previousReconstruct')); plotaxis2;
  hl = legend('In-situ jacobian reconstruct','kCARTA trend','Orig jacobian reconstruct','location','best');
plot(fKc,nanmean(raaReconstruct'),'x-',fKc,nanmean(previousReconstruct'),fKc,nanmean(kctrendALL'),'.-',fKc,nanmean(sartatrend'),fKc,ones(size(fKc))*mean(nwp_trends.era5_100_layertrends.stemp(ind))); plotaxis2;; title(['Latbin ' num2str(iLatBin)])
  xlim([min(fKc) max(fKc)]);  hl = legend('In-situ jacobian reconstruct','Orig jacobian reconstruct','kCARTA trend','SARTA trend','Mean ERA5 stemp','location','best','fontsize',10);
  xlim([640 1640]);  hl = legend('In-situ jacobian reconstruct','Previous jacobian (driver\_gather)  reconstruct','NWP->kcarta 1231 cm-1 rate','NWP->sarta 1231 cm-1 trend','Mean ERA5 stemp','location','best','fontsize',10);

figure(5); subplot(311); plot(fKc,nanmean(raaReconstruct'-previousReconstruct'),fKc,nanstd(raaReconstruct'-previousReconstruct'),fKc,nanmean(raaReconstruct')); ylabel('jacobian');
figure(5); subplot(312); plot(fKc,nanmean(kctrendALL'-sartatrend'),fKc,nanstd(kctrendALL'-sartatrend'),fKc,nanmean(kctrendALL')); ylabel('NWP');
figure(5); subplot(313); plot(fKc,nanmean(kctrendALL'-raaReconstruct'),fKc,nanstd(kctrendALL'-raaReconstruct'),fKc,nanmean(kctrendALL')); ylabel('jacobians vs NWP');
figure(5); subplot(313); plot(fKc,nanmean(sartatrend'-raaReconstruct'),fKc,nanstd(sartatrend'-raaReconstruct'),fKc,nanmean(sartatrend')); ylabel('jacobians vs NWP');

figure(6);plot(fKc,nanmean(sartatrend'-raaReconstruct'),fKc,nanstd(sartatrend'-raaReconstruct'),fKc,nanmean(sartatrend')); ylabel('jacobians vs NWP'); plotaxis2;
xlim([640 1640]);
hl = legend('mean(NWP-jacobians)','std(NWP-jacobians)','mean(NWP)','location','best');

figure(7);plot(fKc,sartatrend-raaReconstruct,'c',fKc,nanmean(sartatrend'-raaReconstruct'),'b',fKc,nanstd(sartatrend'-raaReconstruct'),'r',fKc,nanmean(sartatrend'),'k'); ylabel('jacobians vs NWP'); plotaxis2;
xlim([640 1640]);

figure(8); plot(fKc,sartatrend(:,1)-raaReconstruct(:,1),fKc,sartatrend(:,1));
xlim([640 1640]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% now divide into regions to see diffs
boo = find(fKc <= 800);               junk = sartatrend(boo,:)-raaReconstruct(boo,:); junk = junk.*junk; chisqr(1,:) = sum(junk,1);
boo = find(fKc > 800  & fKc <= 960);  junk = sartatrend(boo,:)-raaReconstruct(boo,:); junk = junk.*junk; chisqr(2,:) = sum(junk,1);
boo = find(fKc > 960  & fKc <= 1140); junk = sartatrend(boo,:)-raaReconstruct(boo,:); junk = junk.*junk; chisqr(3,:) = sum(junk,1);
boo = find(fKc > 1140 & fKc <= 1240); junk = sartatrend(boo,:)-raaReconstruct(boo,:); junk = junk.*junk; chisqr(4,:) = sum(junk,1);
boo = find(fKc > 1240 & fKc <= 1640); junk = sartatrend(boo,:)-raaReconstruct(boo,:); junk = junk.*junk; chisqr(5,:) = sum(junk,1);
figure(9); plot(chisqr','linewidth',2); hl = legend('15 um','12um window','O3','8 um window','WV','location','best','fontsize',10); xlabel('LonBin'); ylabel('\chi^2 = sartaNWP-JacReconstruct'); title(['LatBin ' num2str(iLatBin)])

compare_single_latbin_lonbin_NWPtrend_vs_Reconstruct
