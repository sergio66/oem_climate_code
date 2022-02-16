
iFig = 50;
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(nanmean(fracWV(1:100,mask),2),plays,'linewidth',2);          ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
                                     subplot(122); semilogy(1+nanmean(fracWVunc(1:100,mask),2),plays,'linewidth',2);     ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(nanmean(deltaRH(1:100,mask),2),plays,'linewidth',2);         ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
                                     subplot(122); semilogy(nanmean(abs(deltaRHunc(1:100,mask)),2),plays,'linewidth',2);      ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(nanmean(deltaT(1:100,mask),2),plays,'linewidth',2);          ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
                                     subplot(122); semilogy(nanmean(deltaTunc(1:100,mask),2),plays,'linewidth',2);       ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');
iFig = iFig + 1; figure(iFig); clf;  subplot(121); semilogy(nanmean(deltaO3(1:97,mask),2),plays(1:97),'linewidth',2);    ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
                                     subplot(122); semilogy(nanmean(deltaO3unc(1:97,mask),2),plays(1:97),'linewidth',2); ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');

if exist('airsL3')
  iFig = 50;
  Tlevs = airsL3.Tlevs;
  Qlevs = airsL3.Qlevs;

  boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = maskLFmatr'; end
  junk = airsL3.thestats64x72.waterrate.*boo; junk = reshape(junk,72*64,12);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,1)',Qlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
  junk = airsL3.thestats64x72.waterratestd.*boo; junk = reshape(junk,72*64,12);  
                                   subplot(122); hold on; semilogy(1+nanmean(junk,1)',Qlevs,'linewidth',2);     hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');

  boo = zeros(72,64,12); for ijunk = 1 : 12; boo(:,:,ijunk) = maskLFmatr'; end
  junk = airsL3.thestats64x72.RHrate.*boo; junk = reshape(junk,72*64,12);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,1)',Qlevs,'linewidth',2);        hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
  junk = airsL3.thestats64x72.RHratestd.*boo; junk = reshape(junk,72*64,12);  
                                   subplot(122); hold on; semilogy(nanmean(junk,1)',Qlevs,'linewidth',2);      hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');

  boo = zeros(72,64,12); for ijunk = 1 : 24; boo(:,:,ijunk) = maskLFmatr'; end
  junk = airsL3.thestats64x72.ptemprate.*boo; junk = reshape(junk,72*64,24);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,1)',Tlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
  junk = airsL3.thestats64x72.ptempratestd.*boo; junk = reshape(junk,72*64,24);  
                                   subplot(122); hold on; semilogy(nanmean(junk,1)',Tlevs,'linewidth',2);       hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');

  boo = zeros(72,64,12); for ijunk = 1 : 24; boo(:,:,ijunk) = maskLFmatr'; end
  junk = airsL3.thestats64x72.ozonerate.*boo; junk = reshape(junk,72*64,24);  
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,1)',Tlevs,'linewidth',2);   hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
  junk = airsL3.thestats64x72.ozoneratestd.*boo; junk = reshape(junk,72*64,24);  
                                   subplot(122); hold on; semilogy(nanmean(junk,1)',Tlevs,'linewidth',2); hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');
end

if exist('cmip6')
  iFig = 50;
  Tlevs = plays;
  Qlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end

  junk = cmip6.trend_gas_1; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Qlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
  junk = cmip6.trend_gas_1_err; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
                                   subplot(122); hold on; semilogy(1+nanmean(junk,2)',Qlevs,'linewidth',2);     hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');

  junk = cmip6.trend_RH; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Qlevs,'linewidth',2);        hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
  junk = cmip6.trend_RH_err; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
                                   subplot(122); hold on; semilogy(nanmean(junk,2)',Qlevs,'linewidth',2);      hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');

  junk = cmip6.trend_ptemp; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
  junk = cmip6.trend_ptemp_err; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
                                   subplot(122); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2);       hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');

  junk = cmip6.trend_gas_3; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2);   hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
  junk = cmip6.trend_gas_3_err; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
                                   subplot(122); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2); hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');
end

if exist('era5')
  iFig = 50;
  Tlevs = plays;
  Qlevs = plays;
  boo = zeros(100,72,64); for ijunk = 1 : 100; boo(ijunk,:,:) = maskLFmatr'; end

  junk = era5.trend_gas_1; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Qlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\mu} frac pert');
  junk = era5.trend_gas_1_err; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
                                   subplot(122); hold on; semilogy(1+nanmean(junk,2)',Qlevs,'linewidth',2);     hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('WV_{\sigma} frac pert');

  junk = era5.trend_RH; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Qlevs,'linewidth',2);        hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\mu} pert (%)');
  junk = era5.trend_RH_err; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
                                   subplot(122); hold on; semilogy(nanmean(junk,2)',Qlevs,'linewidth',2);      hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('RH_{\sigma} pert (%)');

  junk = era5.trend_ptemp; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2);         hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\mu} pert');
  junk = era5.trend_ptemp_err; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
                                   subplot(122); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2);       hold off; ylim([1 1000]); set(gca,'ydir','reverse'); title('T_{\sigma} pert');

  junk = era5.trend_gas_3; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
  iFig = iFig + 1; figure(iFig);   subplot(121); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2);   hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\mu} ppm pert');
  junk = era5.trend_gas_3_err; junk = reshape(junk,100,72,64); junk = junk.*boo; junk = reshape(junk,100,72*64);
                                   subplot(122); hold on; semilogy(nanmean(junk,2)',Tlevs,'linewidth',2); hold off; ylim([0.01 1000]); set(gca,'ydir','reverse'); title('O3_{\sigma} ppm pert');
end

iFig = 50;
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend('UMBC','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*2e-2)
                             subplot(122); hold on; plotaxis2; hl = legend('UMBC','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 0.01]+1)
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend('UMBC','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*0.25)
                             subplot(122); hold on; plotaxis2; hl = legend('UMBC','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 1]/2)
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend('UMBC','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*0.05) 
                             subplot(122); hold on; plotaxis2; hl = legend('UMBC','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 1]*0.1)
iFig = iFig+1; figure(iFig); subplot(121); hold on; plotaxis2; hl = legend('UMBC','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([-1 +1]*0.04)
                             subplot(122); hold on; plotaxis2; hl = legend('UMBC','AIRSL3','CMIP6','ERA5','location','best','fontsize',8); hold off; xlim([0 1]*0.04)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
junk = input('make_profile_spectral_trends, can take a while!!!! (-1/+1) : [default = +1] ');
if length(junk) == 0
  junk = 1;
end
if junk > 0
  iVersJac = input('Enter jac version (2019) = 2002/2019 ERAI or (2021) = 2002/2021 ERA5 default : ');
  if length(iVersJac) == 0
    iVersJac = 2021;
  end
  if ~exist('nwp_spectral_trends_cmip6_era5_airsL3_umbc')
    if exist('cmip6')

      %% see SyntheticTimeSeries_ERA5_AIRSL3_CMIP6/compare_sarta_kcarta_latbin32_trends.m
      for lll = 1 : 64
        %% see eg plot_check_WV_T_RH_CMIP6_geo_and_spectral_rates2.m
        dirout = '../FIND_NWP_MODEL_TRENDS/SimulateTimeSeries';
        fsarta = [dirout '/reconstruct_era5_spectra_geo_rlat' num2str(lll,'%02i') '.mat'];
        x = load(fsarta);
        ind = (1:72) + (lll-1)*72;
        sartaERA5trend(:,ind) = x.thesave.xtrendSpectral;
      end

      nwp_spectral_trends_cmip6_era5_airsL3_umbc = make_profile_spectral_trends(cmip6,era5,airsL3,results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,2,iVersJac);

      figure(1); plot(f,nanmean(sartaERA5trend'),f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates')); xlim([640 1640]); hl = legend('from SARTA trends','from jac x dX/dt','location','best');
      figure(1); lll = 1 : 4608; plot(f,sartaERA5trend(:,lll)-nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,lll));
      figure(1); plot(f,nanmean(sartaERA5trend'-nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'),f,nanstd(sartaERA5trend'-nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'),...
                      f,nanmean(sartaERA5trend'),'.-',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates'))

%{
%% this is DA PLOT
plot(f,nanmean(rates(:,mask),2),'b',...
      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),'c',...
      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era_spectral_rates(:,mask),2),'g',...
      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),'r',...
      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.fits(:,mask),2),'m','linewidth',2)
plot(f,nanmean(rates(:,mask),2),'b',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),'c',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),'r',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,mask),2),'g',...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates(:,mask),2),'m',...
     'linewidth',2)
  plotaxis2; hl = legend('AIRS Obs','AIRS L3','ERA5','CMIP6','UMBC','location','best'); axis([640 1640 -0.1 0.05]); title('Spectral Rates');
plot(f,nanmean(rates(:,mask),2),...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),...
     f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,mask),2),...
     'linewidth',0.5)
  plotaxis2; hl = legend('AIRS Obs','AIRS L3','ERA5','CMIP6','location','best','fontsize',10); ylabel('dBT/dt K/yr'); xlabel('Wavenumber cm^{-1}'); axis([640 1640 -0.1 0.05]); %title('Spectral Rates');
  % aslprint('Figs_STMOct2021/spectral_all_N_ratesV2.pdf')
%}

    end

    addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
      [mmw0,lps0] = mmwater_rtp(h,p);
      figure(1); aslmap(1,rlat65,rlon73,maskLFmatr.*smoothn((reshape(mmw0,72,64)'),1),[-90 +90],[-180 +180]); colormap(jet);  title('mmw');  caxis([0 70])
      figure(1); aslmap(1,rlat65,rlon73,maskLFmatr.*smoothn((reshape(lps0.lapse_othermethod,72,64)'),1),[-90 +90],[-180 +180]); colormap(jet);  title('est. lapse rate K/km');   caxis([2 8])
      figure(1); aslmap(1,rlat65,rlon73,maskLFmatr.*smoothn((reshape(lps0.trp_zHI/1000,72,64)'),1),[-90 +90],[-180 +180]);      colormap(jet);  title('Tropopause Hgt HighRes'); caxis([5 20])
      figure(1); aslmap(1,rlat65,rlon73,maskLFmatr.*smoothn((reshape(lps0.trp_z/1000,72,64)'),1),[-90 +90],[-180 +180]);        colormap(jet);  title('Tropopause Hgt 101Res');  caxis([5 20])
      figure(1); aslmap(1,rlat65,rlon73,maskLFmatr.*smoothn((reshape(lps0.scott_trp_z/1000,72,64)'),1),[-90 +90],[-180 +180]);  colormap(jet);  title('Tropopause Hgt Scott');   caxis([5 20])
      booHI    = maskLFmatr.*smoothn((reshape(lps0.trp_zHI/1000,72,64)'),1); 
      boo101   = maskLFmatr.*smoothn((reshape(lps0.trp_z/1000,72,64)'),1); 
      booScott = maskLFmatr.*smoothn((reshape(lps0.scott_trp_z/1000,72,64)'),1); 
        plot(rlat,nanmean(booHI,2),rlat,nanmean(boo101,2),rlat,nanmean(booScott,2),'linewidth',2); title('Tropopause Height'); hl = legend('201 levels','101 levels','Scott','location','best');
      booT     = reshape(p.ptemp(1:97,:),97,72,64); booT = squeeze(nanmean(booT,2)); pcolor(ones(97,1)*rlat',p.palts(1:97,2000)*ones(1,64)/1000,booT); shading interp; colormap jet; colorbar; caxis([200 300])
        hold on; plot(rlat,nanmean(booHI,2),'kx-'); hold off; title('T(z,lat)'); ylim([0 30]); ylabel('Hgt(km)'); %set(gca,'yscale','log'); set(gca,'ydir','reverse'); axis([0 90 0 25])
      %% compare to Fig 1 of https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2011RG000355 by 
      %%  THE EXTRATROPICAL UPPER TROPOSPHERE AND LOWER STRATOSPHERE
      %%  A. Gettelman, P. Hoor, L. L. Pan, W. J. Randel, M. I. Hegglin, T. Birner
      %%  First published: 09 August 2011 https://doi.org/10.1029/2011RG000355Citations: 224

      %% this is relatively fast
      [xall,x200,xtrop,xUTLS] = find_umbc_nwp_profile_trends_mmw_strat(h,p,pavg,nwp_spectral_trends_cmip6_era5_airsL3_umbc,resultsT,resultsWV,resultsO3,results(:,6),maskLFmatr,2);

  end
  figure(9); clf
  plot(f,nanmean(rates(:,mask),2),'k',f,nanmean(fits(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),'linewidth',2)
  plot(f,nanmean(rates(:,mask),2),'k',f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,mask),2),...
                                      f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),'linewidth',2)
  plotaxis2; hl = legend('AIRS Obs','UMBC','AIRS L3','CMIP6','ERA5','location','best'); axis([640 1640 -0.1 0.05]); title('Spectral Rates');

  plot(f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.umbc_spectral_rates(:,mask),2),'x-',...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.airsL3_spectral_rates(:,mask),2),...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.cmip6_spectral_rates(:,mask),2),...
       f,nanmean(nwp_spectral_trends_cmip6_era5_airsL3_umbc.era5_spectral_rates(:,mask),2),f,nanmean(rates(:,mask),2),'linewidth',2)
  plotaxis2; hl = legend('UMBC','AIRS L3','CMIP6','ERA5','AIRS Obs','location','best','fontsize',8); axis([640 1640 -0.1 0.05]); title('Spectral Rates');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

