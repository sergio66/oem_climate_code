disp('making d/dtRH = 0 calcs')

%deltaTT = deltaT;
%for ii = 1 : 4608
%  nlevs = p.nlevs(ii);
%  deltaTT(nlevs:101,ii) = nan;
%end
if ~exist('compute_deltaRH')
  disp(' ')
  disp('computing UMBC    spectral closure'); 
    compute_deltaRH.umbc = const_rh_4608(results(:,6)',deltaT,fracWV,rates,mask); 
  disp(' ')
  disp('computing CMIP6   spectral closure'); 
    compute_deltaRH.cmip6 = const_rh_4608(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.stemp,nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp,nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.gas_1,rates,mask); 
  disp(' ')
  disp('computing ERA5    spectral closure'); 
    compute_deltaRH.era5 = const_rh_4608(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.stemp,nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp,nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.gas_1,rates,mask); 
  disp(' ')
  disp('computing AIRS L3 spectral closure'); 
    compute_deltaRH.airsL3 = const_rh_4608(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.stemp,nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp,nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.gas_1,rates,mask); 
  disp(' ')
  [nansum(compute_deltaRH.cmip6.BT_orig(:)-compute_deltaRH.era5.BT_orig(:)) nansum(compute_deltaRH.cmip6.BT_orig(:)-compute_deltaRH.airsL3.BT_orig(:)) nansum(compute_deltaRH.cmip6.BT_orig(:)-compute_deltaRH.umbc.BT_orig(:))]
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
comment = 'see plot_profile_trends2.m';
the_nwp_trends.airsL3 = nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates;
the_nwp_trends.era5 = nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates;
the_nwp_trends.cmip6 = nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates;
if iNorD > 0
  save reconstructed_spectral_trends_nwp_night.mat the_nwp_trends
elseif iNorD < 0
  save reconstructed_spectral_trends_nwp_day.mat the_nwp_trends
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk = compute_deltaRH.cmip6.BT_orig;
figure(1); clf
  plot(f,nanstd(rates(:,mask)'-(compute_deltaRH.umbc.BT_final(:,mask)'-junk(:,mask)')),...
       f,nanstd(rates(:,mask)'-(compute_deltaRH.airsL3.BT_final(:,mask)'-junk(:,mask)')),...
       f,nanstd(rates(:,mask)'-(compute_deltaRH.cmip6.BT_final(:,mask)'-junk(:,mask)')),...
       f,nanstd(rates(:,mask)'-(compute_deltaRH.era5.BT_final(:,mask)'-junk(:,mask)')))
  plotaxis2;
  hl = legend('UMBC reconstruct','AIRS L3 reconstruct','CMIP6 reconstruct','ERA5 reconstruct','location','best','fontsize',8);
  xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr'); axis([640 1640 0.0 0.075]); title('std(AIRS obs rates - reconstruct rates)')
figure(1);
  plot(f,nanmean(rates(:,mask)'-(compute_deltaRH.umbc.BT_final(:,mask)'-junk(:,mask)')),...
       f,nanmean(rates(:,mask)'-(compute_deltaRH.airsL3.BT_final(:,mask)'-junk(:,mask)')),...
       f,nanmean(rates(:,mask)'-(compute_deltaRH.cmip6.BT_final(:,mask)'-junk(:,mask)')),...
       f,nanmean(rates(:,mask)'-(compute_deltaRH.era5.BT_final(:,mask)'-junk(:,mask)')))
  plot(f,nanmean(rates(:,mask)'-nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates(:,mask)'),...
       f,nanmean(rates(:,mask)'-nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask)'),...
       f,nanmean(rates(:,mask)'-nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,mask)'),...
       f,nanmean(rates(:,mask)'-nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask)'))
  plotaxis2;
  xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr'); axis([640 1640 -0.05 0.05]); title('mean(AIRS obs rates - reconstruct rates)')
  hl = legend('UMBC reconstruct','AIRS L3 reconstruct','CMIP6 reconstruct','ERA5 reconstruct','location','best','fontsize',8);

%figure(2); clf
%  plot(f,nanmean(compute_deltaRH.umbc.BT_final(:,mask)'-compute_deltaRH.umbc.BT_orig(:,mask)'),'x-',...
%       f,nanmean(compute_deltaRH.airsL3.BT_final(:,mask)'-compute_deltaRH.airsL3.BT_orig(:,mask)'),...
%       f,nanmean(compute_deltaRH.cmip6.BT_final(:,mask)'-compute_deltaRH.cmip6.BT_orig(:,mask)'),...
%       f,nanmean(compute_deltaRH.era5.BT_final(:,mask)'-compute_deltaRH.era5.BT_orig(:,mask)'),...
%       f,nanmean(rates(:,mask)'),...
%          'linewidth',2);
%  plotaxis2; axis([640 1640 -0.1 0.05])
%  hl = legend('UMBC reconstruct','AIRS L3 reconstruct','CMIP6 reconstruct','ERA5 reconstruct','AIRS obs','location','best','fontsize',8);
%  xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr'); title('Spectral Rates');
%disp('OH OH looks like this is not as good as figure(9) on line 132 above ....')
figure(2); clf
  plot(f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates(:,mask),2),'x-',...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,mask),2),...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),f,nanmean(rates(:,mask),2),'linewidth',2)
  plotaxis2; hl = legend('UMBC','AIRS L3','CMIP6','ERA5','AIRS Obs','location','best','fontsize',8); axis([640 1640 -0.1 0.05]); title('Spectral Rates');

figure(3); clf
figure(3); subplot(121);
           semilogy(nanmean(compute_deltaRH.umbc.final(1:97,mask)'-compute_deltaRH.umbc.orig(1:97,mask)'),playsjunk,...
                    nanmean(compute_deltaRH.airsL3.final(1:97,mask)'-compute_deltaRH.airsL3.orig(1:97,mask)'),playsjunk,...
                    nanmean(compute_deltaRH.cmip6.final(1:97,mask)'-compute_deltaRH.cmip6.orig(1:97,mask)'),playsjunk,...
                    nanmean(compute_deltaRH.era5.final(1:97,mask)'-compute_deltaRH.era5.orig(1:97,mask)'),playsjunk,...
                    'linewidth',2);
  plotaxis2; set(gca,'ydir','reverse'); ylim([10 1000])
  hl = legend('UMBC reconstruct','AIRS L3 reconstruct','CMIP6 reconstruct','ERA5 reconstruct','location','best','fontsize',8); xlabel('d(RH)/dt'); ylabel('p(mb)'); xlim([-0.5 +0.5]);
figure(3); subplot(122);
           semilogy(nanmean(deltaT(1:97,mask)'),playsjunk,...
                    nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp(1:97,mask)'),playsjunk,...
                    nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp(1:97,mask)'),playsjunk,...
                    nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp(1:97,mask)'),playsjunk,...
                    'linewidth',2);
  plotaxis2; set(gca,'ydir','reverse'); ylim([10 1000])
  hl = legend('UMBC reconstruct','AIRS L3 reconstruct','CMIP6 reconstruct','ERA5 reconstruct','location','best','fontsize',8); xlabel('d(T)/dt'); ylabel('p(mb)'); xlim([-0.05 +0.05])

%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); clf
columnSST.umbc   = ones(101,1)*results(:,6)';
columnSST.airsL3 = ones(101,1)*reshape(airsL3.thestats64x72.stemprate,1,72*64);
columnSST.cmip6  = ones(101,1)*cmip6.trend_stemp;
columnSST.era5   = ones(101,1)*era5.trend_stemp;;

figure(4); subplot(121);
           semilogy(nanmedian((compute_deltaRH.umbc.final(1:97,mask)'-compute_deltaRH.umbc.orig(1:97,mask)')./(columnSST.umbc(1:97,mask)')),playsjunk,...
                    nanmedian((compute_deltaRH.airsL3.final(1:97,mask)'-compute_deltaRH.airsL3.orig(1:97,mask)')./(columnSST.airsL3(1:97,mask)')),playsjunk,...
                    nanmedian((compute_deltaRH.cmip6.final(1:97,mask)'-compute_deltaRH.cmip6.orig(1:97,mask)')./(columnSST.cmip6(1:97,mask)')),playsjunk,...
                    nanmedian((compute_deltaRH.era5.final(1:97,mask)'-compute_deltaRH.era5.orig(1:97,mask)')./(columnSST.era5(1:97,mask)')),playsjunk,...
                    'linewidth',2);
  plotaxis2; set(gca,'ydir','reverse'); ylim([10 1000])
  hl = legend('UMBC reconstruct','AIRS L3 reconstruct','CMIP6 reconstruct','ERA5 reconstruct','location','best','fontsize',8); xlabel('d(RH)/dSKT'); ylabel('p(mb)'); xlim([-10 +10]);
figure(4); subplot(122);
           semilogy(nanmedian((deltaT(1:97,mask)')./(columnSST.umbc(1:97,mask)')),playsjunk,...
                    nanmedian((nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp(1:97,mask)')./(columnSST.airsL3(1:97,mask)')),playsjunk,...
                    nanmedian((nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp(1:97,mask)')./(columnSST.cmip6(1:97,mask)')),playsjunk,...
                    nanmedian((nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp(1:97,mask)')./(columnSST.era5(1:97,mask)')),playsjunk,...
                    'linewidth',2);
  plotaxis2; set(gca,'ydir','reverse'); ylim([10 1000])
  hl = legend('UMBC reconstruct','AIRS L3 reconstruct','CMIP6 reconstruct','ERA5 reconstruct','location','best','fontsize',8); xlabel('d(T)/dSKT'); ylabel('p(mb)'); xlim([-2 +2])

epsx = 1e-4;
epsx = 1e-3;
oo1 = find(abs(results(mask,6)) > epsx);
oo2 = reshape(airsL3.thestats64x72.stemprate,1,72*64); oo2 = find(abs(oo2(mask)) > epsx);
oo3 = find(abs(cmip6.trend_stemp(mask)) > epsx);
oo4 = find(abs(era5.trend_stemp(mask)) > epsx);

figure(4); subplot(121);
           semilogy(nanmean((compute_deltaRH.umbc.final(1:97,mask(oo1))'-compute_deltaRH.umbc.orig(1:97,mask(oo1))')./(columnSST.umbc(1:97,mask(oo1))')),playsjunk,...
                    nanmean((compute_deltaRH.airsL3.final(1:97,mask(oo2))'-compute_deltaRH.airsL3.orig(1:97,mask(oo2))')./(columnSST.airsL3(1:97,mask(oo2))')),playsjunk,...
                    nanmean((compute_deltaRH.cmip6.final(1:97,mask(oo3))'-compute_deltaRH.cmip6.orig(1:97,mask(oo3))')./(columnSST.cmip6(1:97,mask(oo3))')),playsjunk,...
                    nanmean((compute_deltaRH.era5.final(1:97,mask(oo4))'-compute_deltaRH.era5.orig(1:97,mask(oo4))')./(columnSST.era5(1:97,mask(oo4))')),playsjunk,...
                    'linewidth',2);
  plotaxis2; set(gca,'ydir','reverse'); ylim([10 1000])
  hl = legend('UMBC reconstruct','AIRS L3 reconstruct','CMIP6 reconstruct','ERA5 reconstruct','location','best','fontsize',8); xlabel('d(RH)/dSKT'); ylabel('p(mb)'); xlim([-10 +10]);
figure(4); subplot(122);
           semilogy(nanmean((deltaT(1:97,mask(oo1))')./(columnSST.umbc(1:97,mask(oo1))')),playsjunk,...
                    nanmean((nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_100_layertrends.ptemp(1:97,mask(oo2))')./(columnSST.airsL3(1:97,mask(oo2))')),playsjunk,...
                    nanmean((nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_100_layertrends.ptemp(1:97,mask(oo3))')./(columnSST.cmip6(1:97,mask(oo3))')),playsjunk,...
                    nanmean((nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_100_layertrends.ptemp(1:97,mask(oo4))')./(columnSST.era5(1:97,mask(oo4))')),playsjunk,...
                    'linewidth',2);
  plotaxis2; set(gca,'ydir','reverse'); ylim([10 1000])
  hl = legend('UMBC reconstruct','AIRS L3 reconstruct','CMIP6 reconstruct','ERA5 reconstruct','location','best','fontsize',8); xlabel('d(T)/dSKT'); ylabel('p(mb)'); xlim([-2 +2])

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5); clf;
figure(5); subplot(121);
           semilogy(nanmean(RH0(1:97,mask)'),playsjunk,'linewidth',2);
           plotaxis2; set(gca,'ydir','reverse'); ylim([10 1000])
           xlabel('<RH>'); ylabel('p(mb)'); xlim([0 100]);
figure(5); subplot(122);
           semilogy(nanmean(p.ptemp(1:97,mask)'),playsjunk,'linewidth',2);
           plotaxis2; set(gca,'ydir','reverse'); ylim([10 1000])
           xlabel('<T>'); ylabel('p(mb)'); xlim([200 280]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
