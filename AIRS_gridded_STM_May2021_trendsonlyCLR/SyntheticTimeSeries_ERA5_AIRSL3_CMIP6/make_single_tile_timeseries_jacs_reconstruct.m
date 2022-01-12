iNtermsFit = 4; %% default
iNtermsFit = input('enter how many terms for math-tsfit? (4 = default) : ');
if length(iNtermsFit) == 0
  iNtermsFit = 4;
end

dayOFtime = change2days(yy,mm,dd,2002);
for ic = 1 : 2645
  junk = Math_tsfit_lin_robust(dayOFtime,tcalc(ic,:),iNtermsFit); 
  quickBTtrend(ic) = junk(2);
end

clear quick*rate

warning off
junk = Math_tsfit_lin_robust(dayOFtime,stempjunk,iNtermsFit); quickstemprate = junk(2);
fprintf(1,'stemp rates : old vs new %8.6f %8.6f \n',[nwp_trends.era5_100_layertrends.stemp(ind) quickstemprate])

for ic = 1 : p72avg.nlevs-1
  junkjunk = p72x.ptemp(ic,:);
  junk = Math_tsfit_lin_robust(dayOFtime,junkjunk,iNtermsFit); quickptemprate(ic) = junk(2);
end
quickptemprate = quickptemprate';
subplot(121); plot(nwp_trends.era5_100_layertrends.ptemp(1 : p72avg.nlevs-1,ind),1 : p72avg.nlevs-1,'o-',quickptemprate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('T(z) rates : old vs new');
subplot(122); plot(nwp_trends.era5_100_layertrends.ptemp(1 : p72avg.nlevs-1,ind) - quickptemprate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('T(z) rates : old vs new');
disp('ret to continue'); pause;

for ic = 1 : p72avg.nlevs-1
  junkjunk = p72x.gas_1(ic,:); junkjunk = junkjunk/mean(junkjunk);
  junk = Math_tsfit_lin_robust(dayOFtime,junkjunk,iNtermsFit); quickgas_1rate(ic) = junk(2);
end
quickgas_1rate = quickgas_1rate';
subplot(121); plot(nwp_trends.era5_100_layertrends.gas_1(1 : p72avg.nlevs-1,ind),1 : p72avg.nlevs-1,'o-',quickgas_1rate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('frac WV(z) rates : old vs new');
subplot(122); plot(nwp_trends.era5_100_layertrends.gas_1(1 : p72avg.nlevs-1,ind) ./ quickgas_1rate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('frac WV(z) rates : old vs new');
  xlim([0.75 1.25]);
disp('ret to continue'); pause;

for ic = 1 : p72avg.nlevs-1
  junkjunk = p72x.gas_3(ic,:); junkjunk = junkjunk/mean(junkjunk);
  junk = Math_tsfit_lin_robust(dayOFtime,junkjunk,iNtermsFit); quickgas_3rate(ic) = junk(2);
end
quickgas_3rate = quickgas_3rate';
subplot(121); plot(nwp_trends.era5_100_layertrends.gas_3(1 : p72avg.nlevs-1,ind),1 : p72avg.nlevs-1,'o-',quickgas_3rate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('frac O3(z) rates : old vs new');
subplot(122); plot(nwp_trends.era5_100_layertrends.gas_3(1 : p72avg.nlevs-1,ind) ./ quickgas_3rate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('frac O3(z) rates : old vs new');
  xlim([0.75 1.25]);
disp('ret to continue'); pause;
warning on

obs = load(['../Output/Quantile16/test' num2str(ind) '.mat']);
obsx   = obs.rateset.rates;
calcsx = obs.oem.fit;

figure(1); clf
plot(fKc,sartatrend(:,iLonBin),'b.-',fKc,raaReconstruct(:,iLonBin),'r',fKc,quickBTtrend,'k',fKc,obsx,fKc,calcsx); 
  hl = legend('SARTA NWP trend','Reconstructed jac x dX/dt','quicktrend','ACTUAL AIRS OBS','location','best');
  xlim([640 1640])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h72avg,p72avg] = subset_rtp(h72x,p72x,[],[],1);
p72avg.stemp = nanmean(p72x.stemp);
p72avg.ptemp = nanmean(p72x.ptemp,2);
p72avg.gas_1 = nanmean(p72x.gas_1,2);
p72avg.gas_2 = nanmean(p72x.gas_2,2);
p72avg.gas_3 = nanmean(p72x.gas_3,2);
p72avg.gas_4 = nanmean(p72x.gas_4,2);
p72avg.gas_5 = nanmean(p72x.gas_5,2);
p72avg.gas_6 = nanmean(p72x.gas_6,2);
p72avg.gas_9 = nanmean(p72x.gas_9,2);
p72avg.gas_12 = nanmean(p72x.gas_12,2);

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf;  plot(yy2002,p72x.stemp); title('stemp timeseries'); xlim([min(yy2002) max(yy2002)]); plotaxis2;

icc = input('find weighted avg profile??? (-1 = N, +1 = Y, +2 = Y^2, 2002-2021 = specific year) : ');
if icc > 2000
  boo = find(yy == icc);
  p72avg.stemp = nanmean(p72x.stemp(boo));
  p72avg.ptemp = nanmean(p72x.ptemp(:,boo),2);
  p72avg.gas_1 = nanmean(p72x.gas_1(:,boo),2);
  p72avg.gas_2 = nanmean(p72x.gas_2(:,boo),2);
  p72avg.gas_3 = nanmean(p72x.gas_3(:,boo),2);
  p72avg.gas_4 = nanmean(p72x.gas_4(:,boo),2);
  p72avg.gas_5 = nanmean(p72x.gas_5(:,boo),2);
  p72avg.gas_6 = nanmean(p72x.gas_6(:,boo),2);
  p72avg.gas_9 = nanmean(p72x.gas_9(:,boo),2);
  p72avg.gas_12 = nanmean(p72x.gas_12(:,boo),2);
  
elseif icc > 0 & icc < 2000
  dt = 295:0.5:305; y = histc(p72x.stemp,dt); meany = sum(y.^icc.*dt)/sum(y.^icc);
  plot(dt,histc(p72x.stemp,dt)); line([mean(p72x.stemp) mean(p72x.stemp)],[0 50],'color','k'); line([meany meany],[0 50],'color','r');
  [mean(p72x.stemp) meany]
  newavgST = meany;
  
  mmw = mmwater_rtp(h72x,p72x);
  dt = 4:2:30; y = histc(mmw,dt); meany = sum(y.^icc.*dt)/sum(y.^icc);
  plot(dt,histc(mmw,dt)); line([mean(mmw) mean(mmw)],[0 60],'color','k'); line([meany meany],[0 60],'color','r');
  [mean(mmw) meany]
  newavgMMW = meany;
  
  newavgT = zeros(101,1);
  for ic = 1 : p72avg.nlevs-1
    boo = p72x.ptemp(ic,:);
    mn = floor(min(boo)); mx = ceil(max(boo)); dx = (mx-mn)/51;
    dt = mn : dx : mx;  y = histc(boo,dt); meany = sum(y.^icc.*dt)/sum(y.^icc); 
    newavgT(ic) = meany;
  end
  semilogx(p72avg.ptemp(1:p72avg.nlevs-1),1:p72avg.nlevs-1,newavgT(1:p72avg.nlevs-1),1:p72avg.nlevs-1); set(gca,'ydir','reverse');
  plot(p72avg.ptemp(1:p72avg.nlevs-1)-newavgT(1:p72avg.nlevs-1),1:p72avg.nlevs-1); set(gca,'ydir','reverse');
  
  newavgWV = zeros(101,1);
  for ic = 1 : p72avg.nlevs-1
    boo = p72x.gas_1(ic,:);
    boo = log10(boo);
    mn = floor(min(boo)); mx = ceil(max(boo)); dx = (mx-mn)/51;
    dt = mn : dx : mx;  y = histc(boo,dt); meany = sum(y.^icc.*dt)/sum(y.^icc); meany = exp10(meany);
    newavgWV(ic) = meany;
  end
  semilogx(p72avg.gas_1(1:p72avg.nlevs-1),1:p72avg.nlevs-1,newavgWV(1:p72avg.nlevs-1),1:p72avg.nlevs-1); set(gca,'ydir','reverse');
  plot(p72avg.gas_1(1:p72avg.nlevs-1)./newavgWV(1:p72avg.nlevs-1),1:p72avg.nlevs-1); set(gca,'ydir','reverse');
  
  newavgO3 = zeros(101,1);
  for ic = 1 : p72avg.nlevs-1
    boo = p72x.gas_3(ic,:);
    boo = log10(boo);
    mn = floor(min(boo)); mx = ceil(max(boo)); dx = (mx-mn)/51;
    dt = mn : dx : mx;  y = histc(boo,dt); meany = sum(y.^icc.*dt)/sum(y.^icc); meany = exp10(meany);
    newavgO3(ic) = meany;
  end
  semilogx(p72avg.gas_3(1:p72avg.nlevs-1),1:p72avg.nlevs-1,newavgO3(1:p72avg.nlevs-1),1:p72avg.nlevs-1); set(gca,'ydir','reverse');
  plot(p72avg.gas_3(1:p72avg.nlevs-1)./newavgO3(1:p72avg.nlevs-1),1:p72avg.nlevs-1); set(gca,'ydir','reverse');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  p72avg.stemp = newavgST;
  p72avg.ptemp = newavgT;
  p72avg.gas_1 = newavgWV;
  p72avg.gas_3 = newavgO3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

rtpwrite(fop,h72avg,ha72I,p72avg,pa72I);
eval(sartaer);
[~,~,p0,~] = rtpread(frp); t0 = rad2bt(h72x.vchan,p0.rcalc);

dT = 0.1; dQ = 0.01;

p72pert = p72avg; p72pert.gas_2 = p72pert.gas_2 * (1 + dQ);
rtpwrite(fop,h72avg,ha72I,p72pert,pa72I);
eval(sartaer);
[~,~,pCO2,~] = rtpread(frp); tCO2 = rad2bt(h72x.vchan,pCO2.rcalc);

p72pert = p72avg; p72pert.gas_4 = p72pert.gas_4 * (1 + dQ);
rtpwrite(fop,h72avg,ha72I,p72pert,pa72I);
eval(sartaer);
[~,~,pN2O,~] = rtpread(frp); tN2O = rad2bt(h72x.vchan,pN2O.rcalc);

p72pert = p72avg; p72pert.gas_6 = p72pert.gas_6 * (1 + dQ);
rtpwrite(fop,h72avg,ha72I,p72pert,pa72I);
eval(sartaer);
[~,~,pCH4,~] = rtpread(frp); tCH4 = rad2bt(h72x.vchan,pCH4.rcalc);

p72pert = p72avg; p72pert.stemp = p72pert.stemp + dT;
rtpwrite(fop,h72avg,ha72I,p72pert,pa72I);
eval(sartaer);
[~,~,pST,~] = rtpread(frp); tST = rad2bt(h72x.vchan,pST.rcalc);

[~,p72pert] = replicate_rtp_headprof(h72avg,p72avg,1,100); for kk = 1 : 100; p72pert.ptemp(kk,kk) = p72pert.ptemp(kk,kk) + dT; end 
rtpwrite(fop,h72avg,ha72I,p72pert,pa72I);
eval(sartaer);
[~,~,pT,~] = rtpread(frp); tT = rad2bt(h72x.vchan,pT.rcalc);

[~,p72pert] = replicate_rtp_headprof(h72avg,p72avg,1,100); for kk = 1 : 100; p72pert.gas_1(kk,kk) = p72pert.gas_1(kk,kk) * (1+dQ); end 
rtpwrite(fop,h72avg,ha72I,p72pert,pa72I);
eval(sartaer);
[~,~,pWV,~] = rtpread(frp); tWV = rad2bt(h72x.vchan,pWV.rcalc);

[~,p72pert] = replicate_rtp_headprof(h72avg,p72avg,1,100); for kk = 1 : 100; p72pert.gas_3(kk,kk) = p72pert.gas_3(kk,kk) * (1+dQ); end 
rtpwrite(fop,h72avg,ha72I,p72pert,pa72I);
eval(sartaer);
[~,~,pO3,~] = rtpread(frp); tO3 = rad2bt(h72x.vchan,pO3.rcalc);

jST = (tST-t0)/dT;
jT  = (tT-t0*ones(1,100))/dT;
jWV = (tWV-t0*ones(1,100))/log(1+dQ);
jO3 = (tO3-t0*ones(1,100))/log(1+dQ);
jCO2 = (tCO2-t0)/log(1+dQ) * 2.2/385;
jN2O = (tN2O-t0)/log(1+dQ) * 1.0/320;
jCH4 = (tCH4-t0)/log(1+dQ) * 5/1860;
plot(fKc,jST,fKc,sum(jT'),fKc,sum(jWV'),fKc,sum(jO3'))

plot(fKc,jCO2,fKc,squeeze(jac(iLonBin,:,1)));    title('jacCO2 compare (b) new (r) kcarta');
plot(fKc,jN2O,fKc,squeeze(jac(iLonBin,:,2)));    title('jacN2O compare (b) new (r) kcarta');
plot(fKc,jCH4,fKc,squeeze(jac(iLonBin,:,3)));    title('jacCH4 compare (b) new (r) kcarta');
plot(fKc,jST,fKc,squeeze(jac(iLonBin,:,6))/0.1); title('jacST compare (b) new (r) kcarta');
  ix = (1:iaNlays(iLonBin)); 
  ixN = ix;
  ixN = 1 : p72avg.nlevs-1;
  ixx1 = 6 + ixN + 0*iaNlays(iLonBin); plot(fKc,sum(jWV'),fKc,sum(squeeze(jac(iLonBin,:,ixx1))/0.01,2)); title('jacWV compare (b) new (r) kcarta');
  ixxT = 6 + ixN + 1*iaNlays(iLonBin); plot(fKc,sum(jT'),fKc,sum(squeeze(jac(iLonBin,:,ixxT))/0.01,2)); title('jacT compare (b) new (r) kcarta');
  ixx3 = 6 + ixN + 2*iaNlays(iLonBin); plot(fKc,sum(jO3'),fKc,sum(squeeze(jac(iLonBin,:,ixx3))/0.01,2)); title('jacO3 compare (b) new (r) kcarta');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i200 = find(p72avg.plevs >= 200,1);
i500 = find(p72avg.plevs >= 500,1);
i800 = find(p72avg.plevs >= 800,1);

i200 = find(p72avg.plevs >= 200,1);
i500 = find(p72avg.plevs >= 500,1);
i800 = find(p72avg.plevs >= 800,1);

%%%%%%%%%%%%%%%%%%%%%%%%%

p72pert = p72avg;
  p72pert.ptemp(i200:i500-1) = p72pert.ptemp(i200:i500-1) + 2*dT;
  p72pert.ptemp(i500:i800-1) = p72pert.ptemp(i500:i800-1) - 1*dT;
  p72pert.ptemp(i800:101)  = p72pert.ptemp(i800:101)  + 4*dT;
rtpwrite(fop,h72avg,ha72I,p72pert,pa72I);
eval(sartaer);
[~,~,pTestT,~] = rtpread(frp); tTestT = rad2bt(h72x.vchan,pTestT.rcalc);

xraReconstruct = zeros(1,2645);
wah = zeros(1,101); wah(i200:i500-1) = 2*dT; wah(i500:i800-1) = -1*dT; wah(i800:101) = 4*dT; bah = ones(2645,1) * wah; xraReconstruct = xraReconstruct + sum(bah(:,ixN).*jT(:,ixN),2)';
plot(fKc,tTestT-t0,'.-',fKc,xraReconstruct)

%%%%%%%%%%%%%%%%%%%%%%%%%
p72pert = p72avg;
  p72pert.gas_1(i200:i500-1) = p72pert.gas_1(i200:i500-1)*(1 + 2*dQ);
  p72pert.gas_1(i500:i800-1) = p72pert.gas_1(i500:i800-1)*(1 - 1*dQ);
  p72pert.gas_1(i800:101)  = p72pert.gas_1(i800:101) *(1 + 4*dQ);
rtpwrite(fop,h72avg,ha72I,p72pert,pa72I);
eval(sartaer);
[~,~,pTestWV,~] = rtpread(frp); tTestWV = rad2bt(h72x.vchan,pTestWV.rcalc);

xraReconstruct = zeros(1,2645);
wah = zeros(1,101); wah(i200:i500-1) = 2*dQ; wah(i500:i800-1) = -1*dQ; wah(i800:101) = 4*dQ; bah = ones(2645,1) * wah; xraReconstruct = xraReconstruct + sum(bah(:,ixN).*jWV(:,ixN),2)';
plot(wah,1:101,'o-'); set(gca,'ydir','reverse'); plotaxis2;
plot(fKc,tTestWV-t0,'.-',fKc,xraReconstruct)

xraReconstruct = zeros(1,2645);
wah = zeros(1,101); wah(i200:i500-1) = 2*dQ; wah(i500:i800-1) = -1*dQ; wah(i800:101) = 4*dQ; bah = ones(2645,1) * (exp(wah)-1); xraReconstruct = xraReconstruct + sum(bah(:,ixN).*jWV(:,ixN),2)';
plot(wah,1:101,'o-'); set(gca,'ydir','reverse'); plotaxis2;
plot(fKc,tTestWV-t0,'.-',fKc,xraReconstruct)

%%%%%%%%%%%%%%%%%%%%%%%%%
p72pert = p72avg;
  p72pert.ptemp(i200:i500-1) = p72pert.ptemp(i200:i500-1) + 2*dT;
  p72pert.ptemp(i500:i800-1) = p72pert.ptemp(i500:i800-1) - 1*dT;
  p72pert.ptemp(i800:101)    = p72pert.ptemp(i800:101)  + 4*dT;
  p72pert.gas_1(i200:i500-1) = p72pert.gas_1(i200:i500-1)*(1 + 2*dQ);
  p72pert.gas_1(i500:i800-1) = p72pert.gas_1(i500:i800-1)*(1 - 1*dQ);
  p72pert.gas_1(i800:101)  = p72pert.gas_1(i800:101) *(1 + 4*dQ);
rtpwrite(fop,h72avg,ha72I,p72pert,pa72I);
eval(sartaer);
[~,~,pTestWVT,~] = rtpread(frp); tTestWVT = rad2bt(h72x.vchan,pTestWVT.rcalc);

xraReconstruct = zeros(1,2645);
wah = zeros(1,101); wah(i200:i500-1) = 2*dT; wah(i500:i800-1) = -1*dT; wah(i800:101) = 4*dT; bah = ones(2645,1) * wah; xraReconstruct = xraReconstruct + sum(bah(:,ixN).*jT(:,ixN),2)';
wah = zeros(1,101); wah(i200:i500-1) = 2*dQ; wah(i500:i800-1) = -1*dQ; wah(i800:101) = 4*dQ; bah = ones(2645,1) * wah; xraReconstruct = xraReconstruct + sum(bah(:,ixN).*jWV(:,ixN),2)';
plot(wah,1:101,'o-'); set(gca,'ydir','reverse'); plotaxis2;
plot(fKc,tTestWVT-t0,fKc,xraReconstruct)

xraReconstruct = zeros(1,2645);
wah = zeros(1,101); wah(i200:i500-1) = 2*dT; wah(i500:i800-1) = -1*dT; wah(i800:101) = 4*dT; bah = ones(2645,1) * wah; xraReconstruct = xraReconstruct + sum(bah(:,ixN).*jT(:,ixN),2)';
wah = zeros(1,101); wah(i200:i500-1) = 2*dQ; wah(i500:i800-1) = -1*dQ; wah(i800:101) = 4*dQ; bah = ones(2645,1) * (exp(wah)-1); xraReconstruct = xraReconstruct + sum(bah(:,ixN).*jWV(:,ixN),2)';
plot(wah,1:101,'o-'); set(gca,'ydir','reverse'); plotaxis2;
plot(fKc,tTestWVT-t0,fKc,xraReconstruct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ixN = 1 : p72avg.nlevs-1;
raReconstruct = zeros(1,2645);
raReconstruct = raReconstruct + jCO2';
raReconstruct = raReconstruct + jN2O';
raReconstruct = raReconstruct + jCH4';
%{
%% orig
raReconstruct = raReconstruct + jST'*nwp_trends.era5_100_layertrends.stemp(ind);
bah = ones(2645,1) * nwp_trends.era5_100_layertrends.gas_1(ixN,ind)'); raReconstruct = raReconstruct + sum(bah.*jWV(:,ixN),2)';
bah = ones(2645,1) * nwp_trends.era5_100_layertrends.gas_3(ixN,ind)'); raReconstruct = raReconstruct + sum(bah.*jO3(:,ixN),2)';
bah = ones(2645,1) * nwp_trends.era5_100_layertrends.ptemp(ixN,ind)'; raReconstruct = raReconstruct + sum(bah.*jT(:,ixN),2)';
%}

%{
%% very similar if rates are tiny
raReconstruct = raReconstruct + jST'*nwp_trends.era5_100_layertrends.stemp(ind);
bah = ones(2645,1) * (exp(nwp_trends.era5_100_layertrends.gas_1(ixN,ind)')-1); raReconstruct = raReconstruct + sum(bah.*jWV(:,ixN),2)';
bah = ones(2645,1) * (exp(nwp_trends.era5_100_layertrends.gas_3(ixN,ind)')-1); raReconstruct = raReconstruct + sum(bah.*jO3(:,ixN),2)';
bah = ones(2645,1) * nwp_trends.era5_100_layertrends.ptemp(ixN,ind)'; raReconstruct = raReconstruct + sum(bah.*jT(:,ixN),2)';
%}

raReconstruct = raReconstruct + jST'*quickstemprate;
bah = ones(2645,1) * (exp(quickgas_1rate(ixN)')-1); raReconstruct = raReconstruct + sum(bah.*jWV(:,ixN),2)';
bah = ones(2645,1) * (exp(quickgas_3rate(ixN)')-1); raReconstruct = raReconstruct + sum(bah.*jO3(:,ixN),2)';
bah = ones(2645,1) * quickptemprate(ixN)'; raReconstruct = raReconstruct + sum(bah.*jT(:,ixN),2)';

%%%%%%%%%%%%%%%%%%%%%%%%%

plot(fKc,sartatrend(:,iLonBin),'b.-',fKc,quickBTtrend,'c'); 
   hl = legend('SARTA NWP trend','quicktrend','location','best','fontsize',10);
   xlim([640 1640])

plot(fKc,raaReconstruct(:,iLonBin),'r.-',fKc,raReconstruct,'m-'); 
   hl = legend('Reconstructed jac x dX/dt','NEW RECONSTRUCT','location','best','fontsize',10);
   xlim([640 1640])

plot(fKc,sartatrend(:,iLonBin),'b.-',fKc,quickBTtrend,'c',fKc,raaReconstruct(:,iLonBin),'r.-',fKc,raReconstruct,'m-',fKc,obsx,'k.-'); 
   hl = legend('SARTA NWP trend','quicktrend','Reconstructed jac x dX/dt','NEW RECONSTRUCT','ACTUAL AIRS OBS','location','best','fontsize',10);
   xlim([640 1640])

plot(fKc,sartatrend(:,iLonBin),'b.-',fKc,quickBTtrend,'c',fKc,raaReconstruct(:,iLonBin),'r.-',fKc,raReconstruct,'m-'); 
   hl = legend('SARTA NWP trend','quicktrend','Reconstructed jac x dX/dt','NEW RECONSTRUCT','location','best','fontsize',10);
   xlim([640 1640])

plot(fKc,sartatrend(:,iLonBin),'b.-',fKc,raReconstruct,'r-'); 
   plotaxis2;
   hl = legend('SARTA NWP trend','RECONSTRUCT','location','best','fontsize',10);
   xlim([640 1640])

