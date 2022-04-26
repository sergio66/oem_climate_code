if ~exist('ccrtrend')
  load ccrtrends_2002_09_2021_08.mat
end

l1ctrends = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR_zonalavg/strow_rates_with_unc.mat');

figure(1); plot(wnum2378,squeeze(const(:,9,:)),'b',wnum2378,squeeze(ccrconst(:,9,:)),'r'); title('N month mean BTobs');
  plotaxis2; axis([640 1640 200 300])
figure(2); plot(wnum2378,squeeze(trend(:,9,:)),'b',wnum2378,squeeze(ccrtrend(:,9,:)),'r',l1ctrends.f,l1ctrends.airsobs,'k'); title('19 year, 12 month trend mean BTobs');
  plotaxis2; axis([640 1640 -0.1 +0.05])
  plotaxis2; axis([640 2780 -0.1 +0.05])
  hl = legend('CCR straight line fit','CCR Annual Cycle fit','LIC for STM','location','best','fontsize',10); ylabel('dBT/dt K/year'); xlabel('Wavenumber cm^{-1}')

xtrend  = reshape(trend,2378,18*18);
xbtrend = reshape(ccrtrend,2378,18*18);
figure(3); plot(wnum2378,nanmean(xtrend,2),'b',wnum2378,nanmean(xbtrend,2),'r',l1ctrends.f,l1ctrends.airsobs,'k'); title('19 year, 12 month trend mean BTobs');
  plotaxis2; axis([640 1640 -0.1 +0.05])
  plotaxis2; axis([640 2780 -0.1 +0.05])
  hl = legend('CCR straight line fit','CCR Annual Cycle fit','LIC for STM','location','best','fontsize',10); ylabel('dBT/dt K/year'); xlabel('Wavenumber cm^{-1}')

rlat19 = -90:10:+90;
rlon19 = (-90:10:+90)*2;

load llsmap5
figure(4); aslmap(4,rlat19,rlon19,smoothn((reshape(xbtrend(1291,:),18,18)'),1),[-90 +90],[-180 +180]); colormap(jet);  title('Evan Manning BT1231 CCR'); caxis([-1 +1]*0.15); colormap(llsmap5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare to actual 4608 tiles
%% read in data
moonoise = load('iType_4_convert_sergio_clearskygrid_obsonly_Q16.mat','b_err_desc');
b_err_desc = moonoise.b_err_desc; clear moonoise;
b_err_desc = permute(b_err_desc,[3 1 2]);
b_err_desc = reshape(b_err_desc,2645,72*64);

moonoise = load('iType_4_convert_sergio_clearskygrid_obsonly_Q16.mat','b_desc');
b_desc = moonoise.b_desc; clear moonoise;
b_desc = permute(b_desc,[3 1 2]);
b_desc = reshape(b_desc,2645,72*64);

moonoise = load('iType_4_convert_sergio_clearskygrid_obsonly_Q16.mat','h');
f = moonoise.h.vchan; clear moonoise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% map 2645 chans to 2378 chans
load /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;

figure(5); aslmap(5,rlat65,rlon73,smoothn((reshape(b_desc(1520,:),72,64)'),1),[-90 +90],[-180 +180]); colormap(jet);  title('UMBC BT1231 HottestQuant'); caxis([-1 +1]*0.15); colormap(llsmap5)

%% see CRODGERS_FAST_CLOUD/driver_stage1_getparam_setoutput.m
miaow = load('h2645structure.mat');
h = miaow.h; clear miaow
[Y,iX,iY] = intersect(h.ichan,1:2378);
iNotY = setdiff(1:2378,iY);
rates2378 = zeros(2378,4608);  rates2378(iY,:) = b_desc(iX,:);
unc2378   = zeros(2378,4608);  unc2378(iY,:)   = b_err_desc(iX,:);
airs2378 = wnum;
for ijunk = 1 : 4608
  junk = interp1(h.vchan(iX),b_desc(iX,ijunk),airs2378(iNotY));
    rates2378(iNotY,ijunk) = junk;
  junk = interp1(h.vchan(iX),b_err_desc(iX,ijunk),airs2378(iNotY));
    unc2378(iNotY,ijunk) = junk;
end

settings.iIgnoreChans_CH4 = -1;
settings.iIgnoreChans_N2O = -1;
settings.iIgnoreChans_SO2 = -1;
%chanset = jacobian.chanset;
chanset = 1 : 2378;
plotopt.iUpperWavenumLimit = 1620;
plotopt.rlon = meanvaluebin(rlon19);
plotopt.rlat = meanvaluebin(rlat19);
plotopt.scatterORaslmap = -1;

%% this compares 72x64 tiles FAKE FAKE FAKE
airsobs2378 = rates2378 + rand(size(rates2378))*0.01;
plot_spectral_region_chisqr(rates2378,0*rates2378,0*rates2378,airsobs2378,wnum,unc2378,-1,settings,plotopt);
figure(11); ylim([-1 +1]*0.1/2)
figure(12); ylim([-1 +1]*5)
for ii = 15:20; figure(ii); colormap jet; caxis([0 1]*0.5); end; figure(20); caxis([0 1]*0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% map 72x64 tiles to 18x18 CCR tiles
xmap = zeros(18,4);
xwgt = ones(18,4);
xnum = zeros(18,1);
for ii = 1 : 18
  boo = (1:4)+(ii-1)*4;
  xmap(ii,:) = boo;
  xnum(ii) = 4;
end

rlat65s = rlat65(1:64);
rlat65e = rlat65(2:65);
ymap = zeros(18,6);
ywgt = zeros(18,6);
ynum = zeros(18,1);
for ii = 1 : 18
  boo = [rlat19(ii) rlat19(ii+1)];
  miaow1 = find(rlat65e-eps > boo(1)); miaow1 = miaow1(1);
  miaow2 = find(rlat65s < boo(2)-eps); miaow2 = miaow2(end);
  booM = miaow1:miaow2;
  ynum(ii)            = length(booM);
  ymap(ii,1:ynum(ii)) = booM;
  ywgt(ii,1:ynum(ii)) = 1;

  booS = booM(1);
  booE = booM(end);
  if rlat65s(booS) < boo(1)
    wgt = abs(boo(1)-rlat65s(booS))/abs(rlat65s(booS)-rlat65e(booS));
    ywgt(ii,1) = 1-wgt;
  end
  if rlat65e(booE) > boo(end)
    wgt = abs(boo(end)-rlat65e(booE))/abs(rlat65s(booE)-rlat65e(booE));
    ywgt(ii,length(booM)) = 1-wgt;
  end
end

clear rates2378_18x18 unc2378_18x18
rates2378x = reshape(rates2378,2378,72,64);
unc2378x  = reshape(unc2378,2378,72,64);
for ii = 1 : 18
  for jj = 1 : 18
    cumwgt = 0;
    ny = ynum(jj);
    xxind = xmap(ii,:);
    xxwgt = xwgt(ii,:);
    yyind = ymap(jj,1:ny);
    yywgt = ywgt(jj,1:ny);

    dawgt = ones(length(yyind),length(xxind));
    for kk = 1 : length(yyind)
      dawgt(kk,:) = yywgt(kk);
    end

    boo = zeros(2378,1);
    woo = zeros(2378,1);    
    for kkx = 1 : length(xxind)
      for kky = 1 : length(yyind)
        boo = boo + rates2378x(:,xxind(kkx),yyind(kky)) * yywgt(kky);
        woo = woo + unc2378x(:,xxind(kkx),yyind(kky)) * yywgt(kky);      
        cumwgt = cumwgt + yywgt(kky);
      end
    end
    rates2378_18x18(:,ii,jj) = boo/cumwgt;
    unc2378_18x18(:,ii,jj)   = woo/cumwgt;
  end
end

xrates2378_18x18 = rates2378_18x18;
rates2378_18x18 = reshape(rates2378_18x18,2378,18*18);
unc2378_18x18 = reshape(unc2378_18x18,2378,18*18);
ccrtrend_18x18 = reshape(ccrtrend,2378,18*18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check that things look ok BT1231 trend
figure(1);  aslmap(1,rlat65,rlon73,smoothn((reshape(b_desc(1520,:),72,64)'),1),[-90 +90],[-180 +180]); colormap(jet);  title('L1c trends 72x64 BT1231'); cx = caxis;
figure(2);  aslmap(2,rlat65,rlon73,smoothn((reshape(rates2378x(1291,:),72,64)'),1),[-90 +90],[-180 +180]); colormap(jet);  title('L1c --> Lib trends 72x64 BT1231'); caxis(cx);
figure(3);  aslmap(3,rlat19,rlon19,smoothn((reshape(rates2378_18x18(1291,:),18,18)'),0.01),[-90 +90],[-180 +180]); colormap(jet);  title('L1c --> Lib trends 18x18 BT1231'); caxis(cx);
figure(4);  aslmap(4,rlat19,rlon19,smoothn((reshape(ccrtrend_18x18(1291,:),18,18)'),0.01),[-90 +90],[-180 +180]); colormap(jet);  title('CCR trends 18x18 BT1231'); caxis(cx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check that things look ok BT1419 trend
i1419 = find(l1ctrends.f >= 1419,1);
figure(1);  aslmap(1,rlat65,rlon73,smoothn((reshape(b_desc(i1419,:),72,64)'),1),[-90 +90],[-180 +180]); colormap(jet);  title('L1c trends 72x64 BT1419'); cx = caxis;
i1419 = find(wnum2378 >= 1419,1);
figure(2);  aslmap(2,rlat65,rlon73,smoothn((reshape(rates2378x(i1419,:),72,64)'),1),[-90 +90],[-180 +180]); colormap(jet);  title('L1c --> Lib trends 72x64 BT1419'); caxis(cx);
figure(3);  aslmap(3,rlat19,rlon19,smoothn((reshape(rates2378_18x18(i1419,:),18,18)'),0.01),[-90 +90],[-180 +180]); colormap(jet);  title('L1c --> Lib trends 18x18 BT1419'); caxis(cx);
figure(4);  aslmap(4,rlat19,rlon19,smoothn((reshape(ccrtrend_18x18(i1419,:),18,18)'),0.01),[-90 +90],[-180 +180]); colormap(jet);  title('CCR trends 18x18 BT1419'); caxis(cx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(21); plot(f,nanmean(b_desc,2),wnum,nanmean(rates2378_18x18,2))
figure(21); plot(wnum,rates2378_18x18(:,9*18+9),wnum,ccrtrend_18x18(:,9*18+9)); ylim([-1 +1]/10)
damean = mean(ccrtrend_18x18');
good = find(abs(damean) < 0.1);

% this compares 72x64 tiles --> 18x18 tiles versus CCR
plotopt.scatterORaslmap = -2;
plot_spectral_region_chisqr(rates2378_18x18(good,:),0*rates2378_18x18(good,:),0*rates2378_18x18(good,:),ccrtrend_18x18(good,:),...
                            wnum(good),unc2378_18x18(good,:),-1,settings,plotopt);
figure(11); ylim([-1 +1]*0.1/2)
figure(12); ylim([-1 +1]*5)
for ii = 15:20; figure(ii); colormap jet; caxis([0 1]*5/2); end; figure(20); caxis([0 1]/2);

figure(3); plot(wnum2378,nanmean(xtrend,2),'b',wnum2378,nanmean(xbtrend,2),'r',l1ctrends.f,l1ctrends.airsobs,'k'); title('19 year, 12 month trend mean BTobs');
  plotaxis2; axis([640 1640 -0.1 +0.05])
  plotaxis2; axis([640 2780 -0.1 +0.05])
  hl = legend('CCR straight line fit','CCR Annual Cycle fit','LIC for STM','location','best','fontsize',10); ylabel('dBT/dt K/year'); xlabel('Wavenumber cm^{-1}')

figure(3); plot(l1ctrends.f,l1ctrends.airsobs,'g.-',wnum2378,nanmean(rates2378_18x18,2),'b'); title('19 year, 12 month trend mean BTobs');
  plotaxis2; axis([640 1640 -0.1 +0.05])
  plotaxis2; axis([640 2780 -0.1 +0.05])
  hl = legend('2645 LIC, 72x64 for STM','2378 LIC, 18x18 for STM','location','best','fontsize',10); 
ylabel('dBT/dt K/year'); xlabel('Wavenumber cm^{-1}')
 xlim([640 1620])

figure(3); plot(wnum2378(good,:),nanmean(xbtrend(good,:),2),'r',l1ctrends.f,l1ctrends.airsobs,'k',wnum2378,nanmean(rates2378_18x18,2),'b'); title('19 year, 12 month trend mean BTobs');
  plotaxis2; axis([640 1640 -0.1 +0.05])
  plotaxis2; axis([640 2780 -0.1 +0.05])
  hl = legend('2378 CCR, 18x18 Annual Cycle fit','2645 LIC, 72x64 for STM','2378 LIC, 18x18 for STM','location','best','fontsize',10); 
ylabel('dBT/dt K/year'); xlabel('Wavenumber cm^{-1}')
 xlim([640 1620])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
