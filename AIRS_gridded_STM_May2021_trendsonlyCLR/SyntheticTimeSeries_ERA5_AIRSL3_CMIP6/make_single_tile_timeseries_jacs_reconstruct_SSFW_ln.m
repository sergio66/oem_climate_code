%% load test32.mat
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
%addpath /asl/matlib/rtptoolsV201/
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil/

[h72avg,p72avg] = subset_rtp(h72x,p72x,[],[],iLonBin); %% aslo see below

iSeason = input('Enter season 0=ALL, 1=MAM Spring 2=JJA Summer 3=SON Fall 4=DJF Winter : ');
 
iNtermsFit = 1; %% default, only need linear!!!
if iSeason == 0
  iNtermsFit = input('enter how many terms for math-tsfit? (0 = linear only! +N for sin(Nwt)    -N for sin(1/N wt)      4 = default) : ');
  if length(iNtermsFit) == 0
    iNtermsFit = 4;
  end
elseif iSeason > 0
  iNtermsFit = input('enter how many terms for math-tsfit? (0 = linear only! +N for sin(Nwt)    -N for sin(1/N wt)      0 = default) : ');
  if length(iNtermsFit) == 0
    iNtermsFit = 0;
  end
end

iLinOrQuad = -1;
if iNtermsFit >= 0
  iLinOrQuad = input('Enter (1) linear DEFAULT or (2) quad fit : ');
  if length(iLinOrQuad) == 0
    iLinOrQuad = 1;
  end  
end

dayOFtime = change2days(yy,mm,dd,2002);
if iSeason == 0
  usethesemonths = 1 : length(mm);
elseif iSeason == 1
  usethesemonths = find(mm >= 3 & mm <= 5);
elseif iSeason == 2
  usethesemonths = find(mm >= 6 & mm <= 8);
elseif iSeason == 3
  usethesemonths = find(mm >= 9 & mm <= 11);
elseif iSeason == 4
  usethesemonths = find(mm == 1 | mm == 2 | mm == 12);
end
fprintf(1,'of the %4i yy/mm here, use %4i of them \n',length(yy),length(usethesemonths))

%{
fakesigH = 300 + dayOFtime/365*0.1; plot(dayOFtime,fakesigH);
fakesigH = fakesigH - 2*sin(2*pi*1*dayOFtime/365) + 0.5*sin(2*pi*2*dayOFtime/365) + 0.25*cos(2*pi*3*dayOFtime/365) - 0.125*cos(2*pi*4*dayOFtime/365); 
fakesigH = fakesigH + randn(size(fakesigH));
coeff = fft(fakesigH); coeff(228/2-50*1.5:228/2+50*1.5) = 0; newfakesigH = real(ifft(coeff)); 
junkH = Math_tsfit_lin_robust(dayOFtime,fakesigH,4)
xjunkH = Math_tsfit_lin_robust(dayOFtime(usethesemonths),fakesigH(usethesemonths),1); donk = Math_timeseries_2(dayOFtime(usethesemonths),xjunkH);
%  [x_anom,b] = generic_compute_anomaly(dayOFtime,fakesigH); plot(dayOFtime,fakesigH,dayOFtime,x_anom+b(1));
plot(dayOFtime,fakesigH,dayOFtime,newfakesigH,dayOFtime(usethesemonths),donk);  
plot(dayOFtime(usethesemonths),fakesigH(usethesemonths),dayOFtime(usethesemonths),donk);  
newjunkH = Math_tsfit_lin_robust(dayOFtime,newfakesigH,4)

fakesigL = 300 + dayOFtime/365*0.1; plot(dayOFtime,fakesigL);
fakesigL = fakesigL + 2*sin(2*pi*1*dayOFtime/365) - 0.5*sin(2*pi/2*dayOFtime/365) - 0.25*cos(2*pi/3*dayOFtime/365) + 0.125*cos(2*pi/4*dayOFtime/365); 
fakesigL = fakesigL + randn(size(fakesigL));
plot(dayOFtime,fakesigL);
junkL = Math_tsfit_lin_robust(dayOFtime,fakesigL,4)
xjunkL = Math_tsfit_lin_robust(dayOFtime(usethesemonths),fakesigL(usethesemonths),1)
junkL = Math_tsfit_lin_robust_lowfreq(dayOFtime,fakesigL,4)
%}

for ic = 1 : 2645
  if iNtermsFit >= 0
    if iLinOrQuad == 1
      [junk err] = Math_tsfit_lin_robust(dayOFtime(usethesemonths),tcalc(ic,usethesemonths),iNtermsFit); 
    else
      [junk err] = Math_tsfit_quad_robust(dayOFtime(usethesemonths),tcalc(ic,usethesemonths),iNtermsFit); 
    end
  else
    [junk err] = Math_tsfit_lin_robust_lowfreq(dayOFtime(usethesemonths),tcalc(ic,usethesemonths),iNtermsFit); 
  end
  quickBTtrend(ic) = junk(2);
  errBTtrend(ic) = err.se(2);
end

clear quick*rate err*rate

warning off
if iNtermsFit >= 0
  if iLinOrQuad == 1
    [junk err] = Math_tsfit_lin_robust(dayOFtime(usethesemonths),stempjunk(usethesemonths),iNtermsFit); 
  else
    [junk err] = Math_tsfit_quad_robust(dayOFtime(usethesemonths),stempjunk(usethesemonths),iNtermsFit); 
  end
else
  [junk err] = Math_tsfit_lin_robust_lowfreq(dayOFtime(usethesemonths),stempjunk(usethesemonths),iNtermsFit); 
end
quickstemprate = junk(2);
errstemprate = err.se(2);
clf; plot(dayOFtime(usethesemonths),stempjunk(usethesemonths),'.'); pause(0.1)
fprintf(1,'stemp rates : old vs new %8.6f %8.6f \n',[nwp_trends.era5_100_layertrends.stemp(ind) quickstemprate])

for ic = 1 : p72avg.nlevs-1
  junkjunk = p72x.ptemp(ic,usethesemonths);
  if iNtermsFit >= 0
    if iLinOrQuad == 1
      [junk err] = Math_tsfit_lin_robust(dayOFtime(usethesemonths),junkjunk,iNtermsFit); 
    else
      [junk err] = Math_tsfit_quad_robust(dayOFtime(usethesemonths),junkjunk,iNtermsFit); 
    end
  else
    [junk err] = Math_tsfit_lin_robust_lowfreq(dayOFtime(usethesemonths),junkjunk,iNtermsFit); 
  end
  quickptemprate(ic) = junk(2);
  errptemprate(ic) = err.se(2);
end

quickptemprate = quickptemprate';
errptemprate = errptemprate';
subplot(121); plot(nwp_trends.era5_100_layertrends.ptemp(1 : p72avg.nlevs-1,ind),1 : p72avg.nlevs-1,'o-',quickptemprate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('T(z) rates : old vs new');
subplot(122); plot(nwp_trends.era5_100_layertrends.ptemp(1 : p72avg.nlevs-1,ind) - quickptemprate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('T(z) rates : old vs new');
%disp('ret to continue'); pause;
pause(0.1);

for ic = 1 : p72avg.nlevs-1
  junkjunk = p72x.gas_1(ic,usethesemonths); junkjunk = junkjunk/mean(junkjunk);
  if iNtermsFit >= 0
    if iLinOrQuad == 1
      [junk err] = Math_tsfit_lin_robust(dayOFtime(usethesemonths),junkjunk,iNtermsFit); 
      donk = Math_timeseries_2(dayOFtime(usethesemonths),junk);
    else
      [junk err] = Math_tsfit_quad_robust(dayOFtime(usethesemonths),junkjunk,iNtermsFit); 
      donk = Math_timeseries_quad(dayOFtime(usethesemonths),junk);
    end
  else
    [junk err] = Math_tsfit_lin_robust_lowfreq(dayOFtime(usethesemonths),junkjunk,iNtermsFit); 
    donk = Math_timeseries_2_LOWnHIGH(dayOFtime(usethesemonths),junk);
  end
  quickgas_1rate(ic) = junk(2);
  errgas_1rate(ic) = err.se(2);
  clf; plot(dayOFtime(usethesemonths),junkjunk,'.-',dayOFtime(usethesemonths),donk,'x-'); 
  title(['WV ' num2str(ic)]); xlim([min(dayOFtime(usethesemonths)) max(dayOFtime(usethesemonths))]); pause(0.1)
  %pause
end
quickgas_1rate = quickgas_1rate';
errgas_1rate = errgas_1rate';
subplot(121); plot(nwp_trends.era5_100_layertrends.gas_1(1 : p72avg.nlevs-1,ind),1 : p72avg.nlevs-1,'o-',quickgas_1rate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('frac WV(z) rates : old vs new');
subplot(122); plot(nwp_trends.era5_100_layertrends.gas_1(1 : p72avg.nlevs-1,ind) ./ quickgas_1rate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('frac WV(z) rates : old vs new');
  xlim([0.75 1.25]);
%disp('ret to continue'); pause;
pause(0.1)

for ic = 1 : p72avg.nlevs-1
  junkjunk = p72x.gas_3(ic,usethesemonths); junkjunk = junkjunk/mean(junkjunk);
  if iNtermsFit >= 0
    if iLinOrQuad == 1
      [junk err] = Math_tsfit_lin_robust(dayOFtime(usethesemonths),junkjunk,iNtermsFit);
    else
      [junk err] = Math_tsfit_quad_robust(dayOFtime(usethesemonths),junkjunk,iNtermsFit);
    end
  else
    [junk err] = Math_tsfit_lin_robust_lowfreq(dayOFtime(usethesemonths),junkjunk,iNtermsFit); 
  end
  quickgas_3rate(ic) = junk(2);
  errgas_3rate(ic) = err.se(2);
end
quickgas_3rate = quickgas_3rate';
errgas_3rate = errgas_3rate';
subplot(121); plot(nwp_trends.era5_100_layertrends.gas_3(1 : p72avg.nlevs-1,ind),1 : p72avg.nlevs-1,'o-',quickgas_3rate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('frac O3(z) rates : old vs new');
subplot(122); plot(nwp_trends.era5_100_layertrends.gas_3(1 : p72avg.nlevs-1,ind) ./ quickgas_3rate(1 : p72avg.nlevs-1),1 : p72avg.nlevs-1); set(gca,'ydir','reverse'); title('frac O3(z) rates : old vs new');
  xlim([0.75 1.25]);
%disp('ret to continue'); pause;
pause(0.1)
warning on

obs = load(['../Output/Quantile16/test' num2str(ind) '.mat']);
obsx   = obs.rateset.rates;
calcsx = obs.oem.fit;

figure(1); clf
plot(fKc,sartatrend(:,iLonBin),'b.-',fKc,raaReconstruct(:,iLonBin),'r',fKc,quickBTtrend,'k',fKc,obsx,fKc,calcsx); 
  hl = legend('SARTA NWP 12monthx19yr trend','Reconstructed jac x dX/dt 12monthsx19yr','quicktrend','ACTUAL AIRS OBS','OEM FITS','location','best','fontsize',10);
  xlim([640 1640])
%disp('ret to continue'); pause;
pause(0.1)

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

if ~exist('yy2002')
  yy2002 = yy + (mm)/12;
end
figure(2); clf;  plot(yy2002,p72x.stemp,yy2002(usethesemonths),p72x.stemp(usethesemonths),'r.'); title('stemp timeseries'); xlim([min(yy2002) max(yy2002)]); plotaxis2;

icc = 1;
if icc > 0;
  boo = usethesemonths;
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

make_the_SSFW_trends_ln

%{
ax = load('reconstruct_SSFW_season1.mat'); quickBTx(:,1) = ax.quickBTtrend; raReconx(:,1) = ax.raReconstruct;
ax = load('reconstruct_SSFW_season2.mat'); quickBTx(:,2) = ax.quickBTtrend; raReconx(:,2) = ax.raReconstruct;
ax = load('reconstruct_SSFW_season3.mat'); quickBTx(:,3) = ax.quickBTtrend; raReconx(:,3) = ax.raReconstruct;
ax = load('reconstruct_SSFW_season4.mat'); quickBTx(:,4) = ax.quickBTtrend; raReconx(:,4) = ax.raReconstruct;
ax = load('reconstruct_allseasons.mat');
figure(1); clf; plot(ax.fKc,quickBTx,'b.-',ax.fKc,raReconx,'r-');
figure(2); clf; plot(ax.fKc,mean(quickBTx'),'b.-',ax.fKc,mean(raReconx'),'r-'); plotaxis2; xlim([640 1640]);

figure(3); clf; plot(ax.fKc,mean(quickBTx'),'b.-',ax.fKc,ax.sartatrendX,'r-'); plotaxis2; xlim([640 1640]); title('Verifying the sarta trends \newline (linear season vs annual 4term')
figure(4); clf; plot(ax.fKc,mean(raReconx'),'b.-',ax.fKc,ax.raReconstruct,'r-'); plotaxis2; xlim([640 1640]); title('Ooopsing the reconstruction \newline (linear season vs annual 4term)')
figure(5); clf; plot(ax.fKc,mean(quickBTx'),'b.-',ax.fKc,mean(raReconx'),'r.-',ax.fKc,ax.raReconstruct,'k-'); plotaxis2; xlim([640 1640]); 
%}
