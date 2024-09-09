%[junk,wahind,~] = intersect(h.vchan,f_ind);
for ii = 1 : length(f_ind)
  boo = abs(h.vchan-f_ind(ii));
  wahind(ii) = find(boo == min(boo),1);
end
plot(h.vchan(wahind),nanmean(rad2bt(h.vchan(wahind),p.rcalc(wahind,:)),2))

%% /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/Readme_Anomaly_TWP_lat35_lon66
sartajac = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';
[hjunk,pjunk] = subset_rtp(h,p,[],[],(35-1)*72+66);
rtpwrite('junk.op.rtp',hjunk,ha,pjunk,pa);
addpath /home/sergio/KCARTA/MATLAB/
sartaer = ['!' sartajac ' fin=junk.op.rtp fout=junk.rp.rtp listj=-1'];  %% all

sartaer = ['!' sartajac ' fin=junk.op.rtp fout=junk.rp.rtp listj=200']; %% wgt
eval(sartaer); [wsarta,wgtfcn] = readsarta_jacV2('junk.rp.rtp_WGTFCN',200); wgtfcn = squeeze(wgtfcn)';

sartaer = ['!' sartajac ' fin=junk.op.rtp fout=junk.rp.rtp listj=1']; %% G1
eval(sartaer); [wsarta,g1jac] = readsarta_jacV2('junk.rp.rtp_jacG1',1); g1jac = squeeze(g1jac)';

sartaer = ['!' sartajac ' fin=junk.op.rtp fout=junk.rp.rtp listj=100']; %% T
eval(sartaer); [wsarta,tzjac] = readsarta_jacV2('junk.rp.rtp_jacTZ',100); tzjac = squeeze(tzjac)';

rmer = ['!rm junk.rp.rtp_WGTFCN']; eval(rmer)
rmer = ['!rm junk.rp.rtp_jacG1']; eval(rmer)
rmer = ['!rm junk.rp.rtp_jacTZ']; eval(rmer)

jett = jet(128); jett(1,:) = 1;
plays = meanvaluebin(p.plevs(:,(35-1)*72+66));
palts = meanvaluebin(p.palts(:,(35-1)*72+66))/1000;
plays = plays(1:97);
pcolor(h.vchan(wahind),plays,wgtfcn(wahind,1:97)'); shading interp; colorbar; xlim([645 1620]); colormap(jett); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([0 0.1]); ylim([10 1000])
pcolor(h.vchan(wahind),plays,g1jac(wahind,1:97)');  shading interp; colorbar; xlim([645 1620]); colormap(usa2); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 1]*0.1); ylim([10 1000])
pcolor(h.vchan(wahind),plays,tzjac(wahind,1:97)');  shading interp; colorbar; xlim([645 1620]); colormap(jett); 
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([0 1]*0.1); ylim([10 1000])

for ii = 1 : 2645
  boo = wgtfcn(ii,1:97); wgtpeak(ii) = find(boo == max(boo));
  boo = g1jac(ii,1:97);  g1peak(ii) = find(boo == min(boo));
  boo = tzjac(ii,1:97);  tzpeak(ii) = find(boo == max(boo));
end
wgtpeak = palts(wgtpeak);
g1peak  = palts(g1peak);
tzpeak  = palts(tzpeak);
plot(h.vchan,wgtpeak,'k',h.vchan,g1peak,'b',h.vchan,tzpeak,'r'); legend('Wgt','WV','Tz','location','best')
xlabel('Wavenumber'); ylabel('Peak Height [km]')
xlim([640 1640])
xlim([640 740])
xlim([840 1240])
xlim([1240 1640])
disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% look for good UT/LS and US WV channels 

i010mb = find(plays >= 10,1);
i050mb = find(plays >= 50,1);
i100mb = find(plays >= 100,1);
plot(wsarta,g1jac(:,[i010mb i050mb i100mb]));

clear utls_wv us_wv utls_wv_f_ind utls_wv_sarta_ind us_wv_f_ind us_wv_sarta_ind
utls_wv = find(g1jac(:,i050mb) < -0.04);
us_wv = find(g1jac(:,i010mb) > 0.03 & h.vchan >= 1400);
plot(h.vchan(utls_wv),g1jac(utls_wv,i050mb),'b',h.vchan(us_wv),g1jac(us_wv,i010mb),'r')

for ii = 1 : length(utls_wv)
  boo = abs(f_ind - wsarta(utls_wv(ii)));
  utls_wv_f_ind(ii) = find(boo == min(boo),1);
  boo = abs(f_ind(utls_wv_f_ind(ii)) - wsarta);
  utls_wv_sarta_ind(ii) = find(boo == min(boo),1);
end
% [f_ind(utls_wv_f_ind) wsarta(utls_wv_sarta_ind)]

for ii = 1 : length(us_wv)
  boo = abs(f_ind - wsarta(us_wv(ii)));
  us_wv_f_ind(ii) = find(boo == min(boo),1);
  boo = abs(f_ind(us_wv_f_ind(ii)) - wsarta);
  us_wv_sarta_ind(ii) = find(boo == min(boo),1);
end
% [f_ind(us_wv_f_ind) wsarta(us_wv_sarta_ind)]

plot(wsarta(us_wv),wsarta(us_wv)-f_ind(us_wv_f_ind),'bo-',wsarta(utls_wv),wsarta(utls_wv)-f_ind(utls_wv_f_ind),'rx-'); plotaxis2;
  xlabel('SARTA Wavenumber cm-1'); ylabel('SARTA - closest f\_ind'); legend('Upper Strat','Upper Trop/Lower Strat','location','best');
  set(gca,'fontsize',10);

plot(h.vchan(us_wv),g1jac(us_wv,i010mb),'bs',h.vchan(us_wv_sarta_ind),g1jac(us_wv_sarta_ind,i010mb),'cd',...,
     h.vchan(utls_wv),g1jac(utls_wv,i050mb),'ro',h.vchan(utls_wv_sarta_ind),g1jac(utls_wv_sarta_ind,i050mb),'mx','linewidth',2)
  plotaxis2;
  xlabel('SARTA Wavenumber cm-1'); ylabel('WV Jac'); 
  legend('US all SARTA chans','US Chosen Chans','UT/LS all SARTA chans','UT/LS Chosen Chans','location','best');
  set(gca,'fontsize',10);
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% look for good UT/LS and US Tz channels 

i010mb = find(plays >= 10,1);
i050mb = find(plays >= 50,1);
plot(wsarta,tzjac(:,[i010mb i050mb]));

clear utls_tz us_tz utls_tz_f_ind utls_tz_sarta_ind us_tz_f_ind us_tz_sarta_ind 
utls_tz = find(tzjac(:,i050mb) > 0.04 & h.vchan <= 700);
us_tz = find(tzjac(:,i010mb) > 0.04 & h.vchan <= 700);
plot(h.vchan(utls_tz),tzjac(utls_tz,i050mb),'b',h.vchan(us_tz),tzjac(us_tz,i010mb),'r')

clear utls_tz us_tz utls_tz_f_ind utls_tz_sarta_ind us_tz_f_ind us_tz_sarta_ind 
utls_tz = find(tzjac(wahind,i050mb) > 0.04 & h.vchan(wahind) <= 700);
us_tz = find(tzjac(wahind,i010mb) > 0.04 & h.vchan(wahind) <= 700);
utls_tz = wahind(utls_tz);
us_tz = wahind(us_tz);
plot(h.vchan(utls_tz),tzjac(utls_tz,i050mb),'b',h.vchan(us_tz),tzjac(us_tz,i010mb),'r')

for ii = 1 : length(utls_tz)
  boo = abs(f_ind - wsarta(utls_tz(ii)));
  utls_tz_f_ind(ii) = find(boo == min(boo),1);
  boo = abs(f_ind(utls_tz_f_ind(ii)) - wsarta);
  utls_tz_sarta_ind(ii) = find(boo == min(boo),1);
end
% [f_ind(utls_tz_f_ind) wsarta(utls_tz_sarta_ind)]

for ii = 1 : length(us_tz)
  boo = abs(f_ind - wsarta(us_tz(ii)));
  us_tz_f_ind(ii) = find(boo == min(boo),1);
  boo = abs(f_ind(us_tz_f_ind(ii)) - wsarta);
  us_tz_sarta_ind(ii) = find(boo == min(boo),1);
end
% [f_ind(us_tz_f_ind) wsarta(us_tz_sarta_ind)]

plot(wsarta(us_tz),wsarta(us_tz)-f_ind(us_tz_f_ind),'bo-',wsarta(utls_tz),wsarta(utls_tz)-f_ind(utls_tz_f_ind),'rx-'); plotaxis2;
  xlabel('SARTA Wavenumber cm-1'); ylabel('SARTA - closest f\_ind'); legend('Upper Strat','Upper Trop/Lower Strat','location','best');
  set(gca,'fontsize',10);

plot(h.vchan(us_tz),tzjac(us_tz,i010mb),'bs',h.vchan(us_tz_sarta_ind),tzjac(us_tz_sarta_ind,i010mb),'cd',...,
     h.vchan(utls_tz),tzjac(utls_tz,i050mb),'ro',h.vchan(utls_tz_sarta_ind),tzjac(utls_tz_sarta_ind,i050mb),'mx','linewidth',2)
  plotaxis2;
  xlabel('SARTA Wavenumber cm-1'); ylabel('Tz Jac'); 
  legend('US all SARTA chans','US Chosen Chans','UT/LS all SARTA chans','UT/LS Chosen Chans','location','best');
  set(gca,'fontsize',10);
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% i0667 = find(f_ind >= 0667.7,1);
% i1542 = find(f_ind >= 1542,1);

%% these are the best channels based on us_wv_f_ind
i1419 = find(f_ind >= 1419,1);
i1519 = find(f_ind >= 1519,1);
i1521 = find(f_ind >= 1521,1);
i1558 = find(f_ind >= 1558,1);
i1598 = find(f_ind >= 1598,1);

%% these are the best channels based on utls_wv_f_ind
i1419 = find(f_ind >= 1419,1);
i1465 = find(f_ind >= 1465,1);
i1541 = find(f_ind >= 1541,1);
i1554 = find(f_ind >= 1554,1);
i1568 = find(f_ind >= 1568,1);

%% these are the best channels based on utls_tz_f_ind
i0654 = find(f_ind >= 654.4,1);
i0677 = find(f_ind >= 677.0,1);
i0678 = find(f_ind >= 678.3,1);
i0680 = find(f_ind >= 680.1,1);
i0681 = find(f_ind >= 681.7,1);

%% these are the best channels based on us_tz_f_ind
i0658 = find(f_ind >= 658.1,1);
i0663 = find(f_ind >= 663.0,1);
i0669 = find(f_ind >= 669.3,1);
i0674 = find(f_ind >= 674.4,1);

%%%%%%%%%%

%utlswv0 = squeeze(btanomX(:,i1419,:) - btanomX(:,i1519,:));
%utlswv0 = squeeze(btanomX(:,i1519,:));
%utlswv0 = squeeze(btanomX(:,i1521,:));
%utlswv0 = squeeze(btanomX(:,i1542,:));
%utlswv0 = squeeze(btanomX(:,i1558,:));
%utlswv0 = squeeze(btanomX(:,i0667,:));
%%utlswv0 = squeeze(btanomX(:,i0667,:) - btanomX(:,i1558,:));
%utlswv0 = squeeze(btanomX(:,i1598,:) - btanomX(:,i1519,:));

%%%%%%%%%%

utlswv0 = btanomX(:,[i1419 i1519 i1521 i1558 i1598],:); utlswv0 = squeeze(nanmean(utlswv0,2));   strX = 'US WV';
utlswv0 = btanomX(:,[i1419 i1465 i1541 i1554 i1568],:); utlswv0 = squeeze(nanmean(utlswv0,2));   strX = 'UTLS WV';  %% shows something !!!!

utlswv0 = btanomX(:,[i0658 i0663 i0669 i0674],:);       utlswv0 = squeeze(nanmean(utlswv0,2));   strX = 'US Tz';    %% shows something
utlswv0 = btanomX(:,[i0654 i0677 i0678 i0680 i0681],:); utlswv0 = squeeze(nanmean(utlswv0,2));   strX = 'UTLS Tz';  %% shows something

for ii = 1 : 64
  P(ii,:) = polyfit(1:iNumAnomTimeSteps,utlswv0(ii,:),1);
  Q = polyval(P(ii,:),1:iNumAnomTimeSteps);
  utlswv(ii,:) = utlswv0(ii,:) - Q;
end
 
cosNlat = cos(rlat*pi/180) * ones(1,iNumAnomTimeSteps); 
utlswv = utlswv .* coslat;
%utlswv = squeeze(nanmean(reshape(utlswv,72,64,iNumAnomTimeSteps),1));

figure(31); clf; pcolor(yymm,rlat,utlswv);
figure(31); clf; pcolor(yymm,rlat,smoothn(utlswv,1));
shading interp; colormap(usa2); colorbar; caxis([-1 +1]*2); title(strX)
line([2002 2025],[-20.57 -20.57],'color','k','linewidth',2);  line([2022 2022],[-90 +90],'color','k','linewidth',2); %% hunga tonga

figure(32); plot(yymm,smooth(nanmean(utlswv,1),7)); plotaxis2(2022,0); title(strX)
disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(32); clf
BTmeanjunk = nanmean(rad2bt(h.vchan(wahind),p.rcalc(wahind,:)),2);
plot(h.vchan(wahind),BTmeanjunk,'.--')
xlim([650 680])
grid

iaTz1    = find(h.vchan(wahind) >= 667.7,1);

iaTz2(1) = find(h.vchan(wahind) >= 667.25,1);
iaTz2(2) = find(h.vchan(wahind) >= 669.20,1);

iaTz3(1) = find(h.vchan(wahind) >= 652.00,1);
iaTz3(2) = find(h.vchan(wahind) >= 658.00,1);
iaTz3(3) = find(h.vchan(wahind) >= 663.00,1);
iaTz3(4) = find(h.vchan(wahind) >= 674.40,1);

iaTz4(1) = find(h.vchan(wahind) >= 650.80,1);
iaTz4(2) = find(h.vchan(wahind) >= 653.20,1);
iaTz4(3) = find(h.vchan(wahind) >= 656.80,1);
iaTz4(4) = find(h.vchan(wahind) >= 659.30,1);
iaTz4(5) = find(h.vchan(wahind) >= 665.10,1);
iaTz4(6) = find(h.vchan(wahind) >= 673.10,1);
iaTz4(7) = find(h.vchan(wahind) >= 675.10,1);

BTmeanjunk = nanmean(rad2bt(h.vchan(wahind),p.rcalc(wahind,:)),2);
plot(h.vchan(wahind),BTmeanjunk,'.--')
xlim([650 680])
grid
hold on;
plot(h.vchan(wahind(iaTz1)),BTmeanjunk(iaTz1),'ro','linewidth',2)
plot(h.vchan(wahind(iaTz2)),BTmeanjunk(iaTz2),'bo','linewidth',2)
plot(h.vchan(wahind(iaTz3)),BTmeanjunk(iaTz3),'go','linewidth',2)
plot(h.vchan(wahind(iaTz4)),BTmeanjunk(iaTz4),'ko','linewidth',2)
hold off
xlabel('Wavenumber cm-1'); ylabel('BT mean over 64 latbins')
disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

%% https://doi.org/10.1038/s43247-024-01620-3     Nature 2024
%% Strong persistent cooling of the stratosphere after the Hunga eruption
%% Matthias Stocker 1 , Andrea K. Steiner 1, Florian Ladstädter 1, Ulrich Foelsche 1,2 & William Randall

bt_Tz1 = squeeze(btanomX(:,iaTz1,:));                                      bt_Tz1 = bt_Tz1.*cosNlat; %bt_Tz1 = squeeze(nanmean(reshape(bt_Tz1.*coslat,72,64,iNumAnomTimeSteps),1)); 
bt_Tz2 = squeeze(btanomX(:,iaTz2,:)); bt_Tz2 = squeeze(nanmean(bt_Tz2,2)); bt_Tz2 = bt_Tz2.*cosNlat; %bt_Tz2 = squeeze(nanmean(reshape(bt_Tz2.*coslat,72,64,iNumAnomTimeSteps),1)); 
bt_Tz3 = squeeze(btanomX(:,iaTz3,:)); bt_Tz3 = squeeze(nanmean(bt_Tz3,2)); bt_Tz3 = bt_Tz3.*cosNlat; %bt_Tz3 = squeeze(nanmean(reshape(bt_Tz3.*coslat,72,64,iNumAnomTimeSteps),1)); 
bt_Tz4 = squeeze(btanomX(:,iaTz4,:)); bt_Tz4 = squeeze(nanmean(bt_Tz4,2)); bt_Tz4 = bt_Tz4.*cosNlat; %bt_Tz4 = squeeze(nanmean(reshape(bt_Tz4.*coslat,72,64,iNumAnomTimeSteps),1)); 

ilat = 32; bt_show4 = [bt_Tz1(ilat,:); bt_Tz2(ilat,:); bt_Tz3(ilat,:); bt_Tz4(ilat,:)]; imagesc(yymm,1:4,bt_show4); colorbar; colormap(usa2); caxis([-1 +1]*4)
ilat = 32; bt_show4 = [bt_Tz1(ilat,:); bt_Tz2(ilat,:); bt_Tz3(ilat,:); bt_Tz4(ilat,:)]; pcolor(yymm,1:4,bt_show4); colorbar; colormap(usa2); caxis([-1 +1]*4)
  shading interp

ilat = find(rlat >= -30 & rlat < 10);
bt_show4 = [nanmean(bt_Tz2(ilat,:),1); nanmean(bt_Tz1(ilat,:),1); nanmean(bt_Tz3(ilat,:),1); nanmean(bt_Tz4(ilat,:),1)];

hgts = [tzpeak([iaTz2(1) iaTz1(1) iaTz3(1) iaTz4(1)]); 18];
%% note how I have ordered it as 2,1,3,4 to have [45 43 27 23] km

%% Figure 2a
for ii = 1 : 4
  Pjunk = polyfit(1:iNumAnomTimeSteps,bt_show4(ii,:),1);
  Q = polyval(Pjunk,1:iNumAnomTimeSteps);
  bt_show4_detrended(ii,:) = bt_show4(ii,:) - Q;
end
wagaboo = [0*ones(1,iNumAnomTimeSteps); bt_show4_detrended];
pcolor(yymm,hgts,smoothn(wagaboo,1)); colorbar; colormap(usa2); caxis([-1 +1]*2);  shading interp
%set(gca,'ydir','reverse');
line([2022+15/30/12 2022+15/30/12],[15 45],'color','k','linewidth',2)
xlim([2021.9 2024]); ylabel('Hgt [km]'); xlabel('Time')
disp('ret to continue'); pause

%% Figure 2b, T anomaly at 20 km
pcolor(yymm,rlat,bt_Tz4); 
pcolor(yymm,rlat,smoothn(bt_Tz4,1)); colorbar; colormap(usa2); caxis([-1 +1]*2);  shading interp
%set(gca,'ydir','reverse');
line([2022+15/30/12 2022+15/30/12],[15 45],'color','k','linewidth',2)
xlim([2021.9 2024]); ylabel('Latitude'); xlabel('Time')
line([2022+15/30/12 2022+15/30/12],[-1 +1]*50,'color','k','linewidth',2); ylim([-50 +50])
disp('ret to continue'); pause

%% Figure 2c, T anomaly at 28 km
pcolor(yymm,rlat,bt_Tz2); 
pcolor(yymm,rlat,smoothn(bt_Tz2,1)); colorbar; colormap(usa2); caxis([-1 +1]*2);  shading interp
%set(gca,'ydir','reverse');
line([2022+15/30/12 2022+15/30/12],[15 45],'color','k','linewidth',2)
xlim([2021.9 2024]); ylabel('Latitude'); xlabel('Time')
line([2022+15/30/12 2022+15/30/12],[-1 +1]*50,'color','k','linewidth',2); ylim([-50 +50])
disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iaWV2(1) = find(h.vchan(wahind) >= 1436,1);
iaWV2(2) = find(h.vchan(wahind) >= 1456,1);
iaWV2(3) = find(h.vchan(wahind) >= 1507,1);
iaWV2(4) = find(h.vchan(wahind) >= 1539,1);
iaWV2(5) = find(h.vchan(wahind) >= 1560,1);

iaWV1(1) = i1419;
iaWV1(2) = i1519;
iaWV1(3) = i1521;
iaWV1(4) = i1558;
iaWV1(5) = i1598;

BTmeanjunk = nanmean(rad2bt(h.vchan(wahind),p.rcalc(wahind,:)),2);
plot(h.vchan(wahind),BTmeanjunk,'.--')
xlim([1400 1600])
grid
hold on;
plot(h.vchan(wahind(iaWV1)),BTmeanjunk(iaWV1),'ro','linewidth',2)
plot(h.vchan(wahind(iaWV2)),BTmeanjunk(iaWV2),'bo','linewidth',2)
hold off
xlabel('Wavenumber cm-1'); ylabel('BT mean over 64 latbins')
disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

bt_WV2 = squeeze(btanomX(:,iaWV2,:)); bt_WV2 = squeeze(nanmean(bt_WV2,2)); bt_WV2 = bt_WV2.*cosNlat; %bt_WV2 = squeeze(nanmean(reshape(bt_WV2.*coslat,72,64,iNumAnomTimeSteps),1)); 
bt_WV1 = squeeze(btanomX(:,iaWV1,:)); bt_WV1 = squeeze(nanmean(bt_WV1,2)); bt_WV1 = bt_WV1.*cosNlat; %bt_WV1 = squeeze(nanmean(reshape(bt_WV1.*coslat,72,64,iNumAnomTimeSteps),1)); 

pcolor(yymm,rlat,bt_WV1); 
pcolor(yymm,rlat,smoothn(bt_WV1,1)); colorbar; colormap(usa2); caxis([-1 +1]*1);  shading interp
%set(gca,'ydir','reverse');
line([2022+15/30/12 2022+15/30/12],[15 45],'color','k','linewidth',2)
xlim([2021.9 2024]); ylabel('Latitude'); xlabel('Time')
line([2022+15/30/12 2022+15/30/12],[-1 +1]*50,'color','k','linewidth',2); ylim([-50 +50])
disp('ret to continue'); pause

pcolor(yymm,rlat,bt_WV2); 
pcolor(yymm,rlat,smoothn(bt_WV2,1)); colorbar; colormap(usa2); caxis([-1 +1]*1);  shading interp
%set(gca,'ydir','reverse');
line([2022+15/30/12 2022+15/30/12],[15 45],'color','k','linewidth',2)
xlim([2021.9 2024]); ylabel('Latitude'); xlabel('Time')
line([2022+15/30/12 2022+15/30/12],[-1 +1]*50,'color','k','linewidth',2); ylim([-50 +50])
disp('ret to continue'); pause

%%   dBT     dBT   dT      dBT  dq        dq    [ dBT       dBT  dT ] / dBT     anomalyt(6.7 um,t) - jac(T,6.7 um) Tonlyanomaly(15 um,T)
%%   ---  =  ---   --   +  ---  --   ==>  -- =  [ ---   -   ---  -- ] / --- =   -------------------------------------------------------
%%   dt       dT   dt       dq  dt        dt    [  dt       dT   dt ] / dq                       jac(q,6.7 um)
%
% also remember of jacg1 = q dBT/dq ==> dBT/dq = jacg1/q
bt_Q = (bt_WV2 - nanmean(nanmean(tzjac(iaWV2,i050mb),1),1) * bt_Tz2)/nanmean(nanmean(g1jac(iaWV2,i050mb),1),1)/100;
bt_Q = (bt_WV2 - nanmean(nansum(nansum(tzjac(iaWV2,1:i050mb),1)),1) * bt_Tz2)/nanmean(nansum(nansum(g1jac(iaWV2,1:i050mb),1)),1);
bt_Q = bt_Q - mean(bt_Q,2)*ones(1,iNumAnomTimeSteps);
pcolor(yymm,rlat,bt_Q); 
pcolor(yymm,rlat,smoothn(bt_Q,1)); colorbar; colormap(usa2); caxis([-1 +1]*1);  shading interp
%set(gca,'ydir','reverse');
line([2022+15/30/12 2022+15/30/12],[15 45],'color','k','linewidth',2)
xlim([2021.9 2024]); ylabel('Latitude'); xlabel('Time')
line([2022+15/30/12 2022+15/30/12],[-1 +1]*50,'color','k','linewidth',2); ylim([-50 +50])
disp('ret to continue'); pause

bt_Tz2_avg = sum(bt_Tz2.*cosNlat,1)./sum(cosNlat,1);
plot(yymm,bt_Tz2_avg); 
plot(yymm,smooth(bt_Tz2_avg,23)); xlim([2020 2025]); plotaxis2;
disp('ret to continue'); pause

bt_Q_avg = sum(bt_Q.*cosNlat,1)./sum(cosNlat,1);
plot(yymm,bt_Q_avg); 
plot(yymm,smooth(bt_Q_avg,23)); xlim([2020 2025]); plotaxis2;
disp('ret to continue'); pause
