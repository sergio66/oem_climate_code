iFig = 55;

waha = RHSurfpert - RHSurf0;
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,smoothn(reshape(waha,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*1); title('dRH Surf')

iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,smoothn((reshape(rates(i1231,:),72,64)'),1), [-90 +90],[-180 +180]); title('dBT1231/dt'); caxis([-1 +1]*0.15); colormap(llsmap5)

%iFig = iFig + 1; aslmap(iFig,rlat65,rlon73,reshape(nlays_straight_from_results,72,64)', [-90 +90],[-180 +180]); title('nlays'); caxis([40 50]); colormap(jet)

umbc_wv_frac_lowestlayer = save_cov_set.xf_wvz(2,:); %% retrieved
iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,smoothn(reshape(umbc_wv_frac_lowestlayer,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.015); title('UMBC RETRIEVED xb(WVfrac(gnd))')

%% see eg isaac_held_dRH_dST_pgorman_dcolwater_dST.m or set_CO2_CH4_N2O_ESRL.m
waha = rates(i1231,:)./pMean17years.stemp./pMean17years.stemp;
  Lo = 2.5e6;  %%% J/kg
  Rv = 461.52; %%% J/kg/K
  guess_wv_frac_lowestlayer_RH0 = exp(Lo/Rv * waha)-1;
  guess_wv_frac_lowestlayer_RH0 = (Lo/Rv * waha);
  iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,smoothn(reshape(guess_wv_frac_lowestlayer_RH0,72,64)',1), [-90 +90],[-180 +180]); 
  colormap(llsmap5); caxis([-1 +1]*0.015); title('GUESS xb(WVfrac(gnd)) if dRH = 0')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if max(save_cov_set.xb_wvz(2,:)) >= 1e-6
  if topts.set_era5_cmip6_airsL3 == 0
    strX = 'RH = 0 or 0.01 Guess!!';
  elseif topts.set_era5_cmip6_airsL3 == 5
    strX = 'ERA5';
  elseif topts.set_era5_cmip6_airsL3 == 6
    strX = 'CMIP6';
  elseif topts.set_era5_cmip6_airsL3 == -6
    strX = 'AMIP6';
  elseif topts.set_era5_cmip6_airsL3 == 3
    strX = 'AIRS L3';
  elseif topts.set_era5_cmip6_airsL3 == -3
    strX = 'CLIMCAPS L3';
  elseif topts.set_era5_cmip6_airsL3 == 2
    strX = 'MERRA2';
  elseif topts.set_era5_cmip6_airsL3 == 8
    strX = 'UA MLS';
  end

  strX2 = ['warning : save_cov_set.xb_wvz(2,:) is nonzero, set by ' strX ' (topts.set_era5_cmip6_airsL3 = ' num2str(topts.set_era5_cmip6_airsL3) ')'];
  disp(strX2);
  if topts.set_era5_cmip6_airsL3 == 0
    disp(['        : it was set by the first guess based on dBT1231/dt in set_CO2_CH4_N2O_ESRL.m (topts.iAdjLowerAtmWVfrac = ' num2str(topts.iAdjLowerAtmWVfrac) ')'])
  end

  guess2_wv_frac_lowestlayer = save_cov_set.xb_wvz(2,:); %% a priori
  iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,smoothn(reshape(guess2_wv_frac_lowestlayer,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*0.015); title([strX ' xb(WVfrac(gnd))'])

  waha = guess_wv_frac_lowestlayer_RH0./guess2_wv_frac_lowestlayer;
  waha = abs(waha);
  waha = 1./waha;
  iFig = iFig + 1; figure(iFig); clf;  aslmap(iFig,rlat65,rlon73,smoothn(reshape(waha,72,64)',1), [-90 +90],[-180 +180]); colormap(llsmap5); caxis([-1 +1]*5); 
  caxis([0 1]*2);    title(['GUESS xb : RH=used/RH=0']); colormap(jett)
  caxis([1 1.25]);   title(['GUESS xb : RH=used/RH=0']); colormap(jett)

  iFig = iFig + 1; figure(iFig); clf;  
  dn = 0:0.01:010; plot(dn,hist(waha,dn));     title(['GUESS xb : RH=used/RH=0'])
  dn = 0:0.01:100; semilogy(dn,hist(waha,dn)); title(['GUESS xb : RH=used/RH=0']); xlim([0 10])

  dn = 0 : 0.0001 : 1; semilogy(dn,hist(umbc_wv_frac_lowestlayer,dn),'b',dn,hist(guess2_wv_frac_lowestlayer,dn),'r',dn,hist(guess_wv_frac_lowestlayer_RH0,dn),'g','linewidth',2)
    hl = legend('UMBC RETRIEVED',strX,'guess RH=0'); title('dWVfrac/dt'); xlim([0 0.01])

else
  iFig = iFig + 1; figure(iFig); clf;  
  dn = 0 : 0.0001 : 1; semilogy(dn,hist(umbc_wv_frac_lowestlayer,dn),'b',dn,hist(guess_wv_frac_lowestlayer_RH0,dn),'g','linewidth',2)
    hl = legend('UMBC RETRIEVED','this guess'); title('dWVfrac/dt'); xlim([0 0.01])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iFig = iFig + 1; 
guess_delta_ILR

