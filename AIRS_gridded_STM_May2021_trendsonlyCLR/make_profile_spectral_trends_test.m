function [] = make_profile_spectral_trends_test(results,resultsWV,resultsT,resultsO3,fits,rates,pavg,plays,f,mask);

[h1_4608,~,p1_4608,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp'); %% line 23 of see_clust_put_together_jacs_clr.m
[~,kcarta.subjac.ppmv2] = layers2ppmv(h1_4608,p1_4608,1:4608,2);
[~,kcarta.subjac.ppmv4] = layers2ppmv(h1_4608,p1_4608,1:4608,4);
[~,kcarta.subjac.ppmv6] = layers2ppmv(h1_4608,p1_4608,1:4608,6);

xtest_fit_spectral_rates = zeros(2645,4608);

iS = 31; iE = 31;
for ii = iS : iE
  if mod(ii,10) == 0
    fprintf(1,'+')
  else
    fprintf(1,'.')
  end

% really not needed
%  jacfile = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//kcarta_clr_subjacLatBin_newSARTA_' num2str(ii,'%02d') '.mat'];
%  jac = load(jacfile);

  for jjj = 1 : 72
    jacfilex = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin//usethisjac_clear_reconstructcode_latbin_' num2str(ii,'%02d') '_lonbin_' num2str(jjj,'%02d') '.mat'];
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

  %{
  plot(1:2645,jac.subjac.coljacCO2,'b',1:2645,m_ts_jac.subjac.coljacCO2*400,'r')
  plot(1:2645,jac.subjac.coljacN2O,'b',1:2645,m_ts_jac.subjac.coljacN2O*300,'r')
  plot(1:2645,jac.subjac.coljacCH4,'b',1:2645,m_ts_jac.subjac.coljacCH4*1860,'r')
  plot(1:2645,jac.subjac.jacST,'b',1:2645,m_ts_jac.subjac.jacST./0.1,'r')
  plot(1:2645,jac.subjac.jacWV(:,:,1),'b',1:2645,m_ts_jac.subjac.jacWV(:,:,1)*100,'r')
    plot(1:2645,sum(jac.subjac.jacWV(:,:,1),1),'b.-',1:2645,sum(m_ts_jac.subjac.jacWV(:,:,1),1)*100,'r')
  plot(1:2645,jac.subjac.jacT(:,:,1),'b',1:2645,m_ts_jac.subjac.jacT(:,:,1)*100,'r')
    plot(1:2645,sum(jac.subjac.jacT(:,:,1),1),'b.-',1:2645,sum(m_ts_jac.subjac.jacT(:,:,1),1)*100,'r')
  plot(1:2645,jac.subjac.jacO3(:,:,1),'b',1:2645,m_ts_jac.subjac.jacO3(:,:,1)*100,'r')
    plot(1:2645,sum(jac.subjac.jacO3(:,:,1),1),'b.-',1:2645,sum(m_ts_jac.subjac.jacO3(:,:,1),1)*100,'r')
  error('kl;jsg;kgs')
  %}

  % forget the renorm
  m_ts_jac.subjac.jacST = m_ts_jac.subjac.jacST * 1;
  m_ts_jac.subjac.jacWV = quick_combinejaclays_make_profile_spectral_trends(x100m_ts_jac.subjac.jacWV * 1.00,1:100,20);
  m_ts_jac.subjac.jacT  = quick_combinejaclays_make_profile_spectral_trends(x100m_ts_jac.subjac.jacT  * 1.00,1:100,20);
  m_ts_jac.subjac.jacO3 = quick_combinejaclays_make_profile_spectral_trends(x100m_ts_jac.subjac.jacO3 * 1.00,1:100,20);

  ind = (ii-1)*72 + (1:72);

  %% note we use 1/co2ppmv,1/n2oppb,1/ch4ppb instead of 2.2/co2ppmv,0.8/n2oppb,5/ch4ppb -- because we are really doing a jacobian normalized by 2.2 ppmv change etc etc
  %% see get_jac_fast.m

  Qlevs = pavg;
  Tlevs = pavg;

  xtest_fit_spectral_rates(:,ind) = xtest_fit_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCO2 .* (ones(2645,1)*results(ind,1)'/2.2);
  xtest_fit_spectral_rates(:,ind) = xtest_fit_spectral_rates(:,ind) + m_ts_jac.subjac.coljacN2O .* (ones(2645,1)*results(ind,2)'/1.0);
  xtest_fit_spectral_rates(:,ind) = xtest_fit_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCH4 .* (ones(2645,1)*results(ind,3)'/5.0);

  xtest_fit_spectral_rates(:,ind) = xtest_fit_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC11 * (0);
  xtest_fit_spectral_rates(:,ind) = xtest_fit_spectral_rates(:,ind) + m_ts_jac.subjac.coljacCFC12 * (0);
  xtest_fit_spectral_rates(:,ind) = xtest_fit_spectral_rates(:,ind) + m_ts_jac.subjac.jacST .* (ones(2645,1)*results(ind,6)'/0.1);

  xjunkrate = resultsWV(ind,:); clear junkrate
    junkrate = xjunkrate; 
    junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100/5; xtest_fit_spectral_rates(:,ind) = xtest_fit_spectral_rates(:,ind) + squeeze(m_ts_jac.subjac.jacWV(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
  xjunkrate = resultsT(ind,:); clear junkrate
    junkrate = xjunkrate; 
    junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100/5; xtest_fit_spectral_rates(:,ind) = xtest_fit_spectral_rates(:,ind) + squeeze(m_ts_jac.subjac.jacT(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end
  xjunkrate = resultsO3(ind,:); clear junkrate
    junkrate = xjunkrate;
    junkrate = junkrate'; junkrate(isnan(junkrate)) = 0; for jjj = 1 : 100/5; xtest_fit_spectral_rates(:,ind) = xtest_fit_spectral_rates(:,ind) + squeeze(m_ts_jac.subjac.jacO3(jjj,:,:)) .* (ones(2645,1)*junkrate(jjj,:)/0.01); end

end

fprintf(1,'\n')

%% see gather_gridded_retrieval_results_plots
figure(9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
junk = [m_ts_jac.subjac.coljacCO2(1520,72) m_ts_jac.subjac.coljacN2O(1520,72) m_ts_jac.subjac.coljacCH4(1520,72) m_ts_jac.subjac.coljacCFC11(1520,72) m_ts_jac.subjac.coljacCFC12(1520,72) m_ts_jac.subjac.jacST(1520,72)];
punk = [results(2232,:)./[2.2 1 5 1 1 0.1]; junk]';
fprintf(1,'TraceGases  = rates/qrenorm jac*qrenorm %10.6f %10.6e \n',punk');

disp(' ')
junk = [m_ts_jac.subjac.jacWV(:,1520,72)];
punk = [resultsWV(2232,:)'/0.01 junk]';
fprintf(1,'Water Row = rates/qrenorm jac*qrenorm %10.6f %10.6e \n',punk);

disp(' ')
junk = [m_ts_jac.subjac.jacT(:,1520,72)];
punk = [resultsT(2232,:)'/0.01 junk]';
fprintf(1,'Tempr Row = rates/qrenorm jac*qrenorm %10.6f %10.6e \n',punk);

disp(' ')
junk = [m_ts_jac.subjac.jacO3(:,1520,72)];
punk = [resultsO3(2232,:)'/0.01 junk]';
fprintf(1,'Ozone Row = rates/qrenorm jac*qrenorm %10.6f %10.6e \n',punk);

%% look at oem_fit/oem_lls.m
%kapoo = load('/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/usethisjac_clear_reconstructcode_latbin_31_lonbin_72.mat');
%junk = kapoo.m_ts_jac(1520,1:6)./([2.2 1 5 1/300 1/600 0.1]./[kcarta.subjac.ppmv2(2232) kcarta.subjac.ppmv4(2232) kcarta.subjac.ppmv6(2232) 1 1 1]);
%junk = kapoo.m_ts_jac(1520,1:6);
%fprintf(1,'TraceGases Row 2X = jac*qrenorm   %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e \n',junk);
%junk = junk.*[2.2 1 5 1 1 0.1];
%fprintf(1,'C %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e \n',junk);
%junk = kapoo.m_ts_jac(1520,1:6)./([2.2 1 5 1/300 1/600 0.1]./[kcarta.subjac.ppmv2(2232) kcarta.subjac.ppmv4(2232) kcarta.subjac.ppmv6(2232) 1 1 1]);
%junk = kapoo.m_ts_jac(1520,1:6);
%fprintf(1,'Tempr Row 2X = jac*qrenorm   %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e \n',junk);
%junk = junk.*[2.2 1 5 1 1 0.1];
%fprintf(1,'C %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e \n',junk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
maskx = mask;
maskx = 2232;
if length(maskx) == 1
  plot(f,nanmean(fits(:,mask)-test_fit_spectral_rates(:,mask),2),'k',f,nanstd(fits(:,mask)-test_fit_spectral_rates(:,mask),[],2),'g--')
  plot(f,rates(:,maskx),'k',f,fits(:,maskx),'b',f,test_fit_spectral_rates(:,maskx),'g',f,xtest_fit_spectral_rates(:,maskx),'r')
  plot(f,fits(:,maskx)-test_fit_spectral_rates(:,maskx),'b',f,fits(:,maskx)-xtest_fit_spectral_rates(:,maskx),'r')
  plot(f,fits(:,maskx)-xtest_fit_spectral_rates(:,maskx),'r')
end

plot(f,nanmean(rates(:,mask),2),'b',f,nanmean(fits(:,mask),2),'r',f,nanmean(xtest_fit_spectral_rates(:,mask),2),'k')
plot(f,nanmean(rates(:,mask)-fits(:,mask),2),'b',f,nanstd(rates(:,mask)-fits(:,mask),[],2),'c',f,nanmean(rates(:,mask)-xtest_fit_spectral_rates(:,mask),2),'r',f,nanstd(rates(:,mask)-xtest_fit_spectral_rates(:,mask),[],2),'m')

plot(f,rates(:,mask),'b',f,airsL3_spectral_rates(:,mask),'c',f,era_spectral_rates(:,mask),'g',f,era5_spectral_rates(:,mask),'r',f,fits(:,mask),'m')
plot(f,nanmean(rates(:,mask),2),'b',f,nanmean(fits(:,mask),2),'m',f,nanmean(xtest_fit_spectral_rates(:,mask),2),'k')
plot(f,nanmean(fits(:,mask)-test_fit_spectral_rates(:,mask),2),'b',f,nanstd(fits(:,mask)-test_fit_spectral_rates(:,mask),[],2),'c',f,nanmean(fits(:,mask)-xtest_fit_spectral_rates(:,mask),2),'r',f,nanstd(fits(:,mask)-xtest_fit_spectral_rates(:,mask),[],2),'m')

%% this is DA PLOT
plot(f,nanmean(rates(:,mask),2),'b',f,nanmean(airsL3_spectral_rates(:,mask),2),'c',f,nanmean(era_spectral_rates(:,mask),2),'g',f,nanmean(era5_spectral_rates(:,mask),2),'r',f,nanmean(fits(:,mask),2),'m','linewidth',2)
plot(f,nanmean(rates(:,mask),2),'b',f,nanmean(airsL3_spectral_rates(:,mask),2),'c',f,nanmean(era_spectral_rates(:,mask),2),'g',f,nanmean(era5_spectral_rates(:,mask),2),'r',f,nanmean(fits(:,mask),2),'m','linewidth',2)
plot(f,nanmean(rates(:,mask),2),'b',f,nanmean(airsL3_spectral_rates(:,mask),2),'c',f,nanmean(era_spectral_rates(:,mask),2),'g',f,nanmean(era5_spectral_rates(:,mask),2),'r',f,nanmean(xtest_fit_spectral_rates(:,mask),2),'m','linewidth',2)
  plotaxis2; hl = legend('AIRS Obs','AIRS L3','ERA','ERA5','quick test','location','best'); axis([640 1640 -0.1 0.05]); title('Spectral Rates');

