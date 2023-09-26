addpath /asl/matlib/science/
addpath /home/sergio/MATLABCODE/TIME

if ~exist('btanom')
  load /asl/s1/sergio/JUNK/anomaly_ALL_QO3.mat
  daysSince2002 = change2days(yy,mm,dd,2002);
  load h2645structure.mat
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rlat65')
  rlat = load('/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/latB64.mat'); 
  rlat65 = rlat.latB2; rlat = 0.5*(rlat.latB2(1:end-1)+rlat.latB2(2:end));
  rlon73 = (1:73); rlon73 = -180 + (rlon73-1)*5;  rlon = (1:72); rlon = -177.5 + (rlon-1)*5;
  [Y,X] = meshgrid(rlat,rlon);
  [salti, landfrac] = usgs_deg10_dem(Y,X);
  
  lf = landfrac; lf = lf(:);
  %% XX = X'; XX = XX(:); %% MAN THIS IS CONFUSING??  ie do XX = X'; see pcolor below
  %% YY = Y'; YY = YY(:); %% MAN THIS IS CONFUSING??  ie do XX = X'; see pcolor below
  XX = X; XX = XX(:); %% MAN THIS IS CONFUSING BUT IT IS RIGHT  ie do not do XX = X'; see pcolor below
  YY = Y; YY = YY(:); %% MAN THIS IS CONFUSING BUT IT IS RIGHT  ie do not do XX = X'; see pcolor below
end

coslat = cos(YY*pi/180)*ones(1,454);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oni = load('ONI_sep2023.txt');
  [aaa,bbb] = size(oni);
  oniS = oni(1,1); oniE = oni(aaa,1);
  oni = oni(1:aaa,2:bbb); oni = oni'; oni = oni(:);
  onidd = 1:length(oni); onidd = (onidd-1)/12 + oniS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alone = load('anomaly_chID_1520_Q03.mat');
clear btavg* usethese

iCnt = 1;
usethese{iCnt} = 1 : 4608;
for ii = 1 : 2645
  moo = squeeze(btanom(:,ii,:));
  btavg1(ii,:)  = nansum(moo.*coslat)./nansum(coslat);
  xbtavg1(ii,:) = nanmean(moo.*coslat);
end
pcolor(2002+daysSince2002/365,h.vchan,btavg1); colorbar; shading flat;
pcolor(2002+daysSince2002/365,h.vchan,btavg1-xbtavg1); colorbar; shading flat;
plot(2002+daysSince2002/365,xbtavg1(1520,:),'bx-',2002+daysSince2002/365,nanmean(alone.btanom.*coslat,1),'linewidth',2); plotaxis2; xlim([2002 2023])
plot(2002+daysSince2002/365,btavg1(1520,:),'b',2002+daysSince2002/365,nanmean(alone.btanom.*coslat,1),'linewidth',2); plotaxis2; xlim([2002 2023])

xow = xbtavg1(1520,:); PX = polyfit(daysSince2002,xow,1); YPX = polyval(PX,daysSince2002);
wow = btavg1(1520,:);  PW = polyfit(daysSince2002,wow,1); YPW = polyval(PW,daysSince2002);
plot(onidd,oni,'k',2002+daysSince2002/365,10*smooth(xbtavg1(1520,:)-YPX,16*2),'b',2002+daysSince2002/365,10*smooth(btavg1(1520,:)-YPW,16*2),'r','linewidth',2); xlim([2002 2023])
plotaxis2; hl = legend('ONI','BT1231 anomaly vers1','BT1231 anomaly vers2','location','best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newLatGrid = [-90 -70 -50 -30 -10 +10 +30 +50 +70 +90];
newLatGrid = [-90:+10:+90];
newLatGrid = [-90 -75 -60 [-55:5:+55] +60 +75 +90];

for iCnt = 1 : length(newLatGrid)-1
  clear junk
  fprintf(1,'making avg %2i of %2i \n',iCnt, length(newLatGrid)-1);
  usethese{iCnt+1} = find(YY >= newLatGrid(iCnt) & YY < newLatGrid(iCnt+1));
  xusethese = usethese{iCnt+1};
  junkcoslat = coslat(xusethese,:);
  for ii = 1 : 2645
    moo = squeeze(btanom(:,ii,:));
    moo = moo(xusethese,:);
    junk(ii,:) = nansum(moo.*junkcoslat)./nansum(junkcoslat);
  end
  str = ['btavg' num2str(iCnt+1) ' = junk;'];
  eval(str);
end

%btavgAnomFinal = [btavg0 btavg1 btavg2 btavg3 btavg4 btavg5 btavg6 btavg7 btavg8 btavg9];
btavgAnomFinal = [];
for ii = 1 : length(newLatGrid)
  str = ['btavgAnomFinal = [btavgAnomFinal btavg' num2str(ii) '];'];
  eval(str);
end

disp('remember -0.06 K/yr trend for 792 Q branch  ==> 20 years we would span 0 to -0.06*20 = -1.2K, so on average the 792 Q branch anomaly should be -0.6 K')
disp('remember -0.06 K/yr trend for 792 Q branch  ==> 20 years we would span 0 to -0.06*20 = -1.2K, so on average the 792 Q branch anomaly should be -0.6 K')
disp('remember -0.06 K/yr trend for 792 Q branch  ==> 20 years we would span 0 to -0.06*20 = -1.2K, so on average the 792 Q branch anomaly should be -0.6 K')
for iCnt = 1 : length(newLatGrid)
  str = ['junk = btavg' num2str(iCnt) ';'];
  eval(str);
  figure(1); clf; pcolor(2002+daysSince2002/365,h.vchan,junk); colorbar; caxis([-1 +1]*2); shading flat; ylim([640 1640]); title(num2str(iCnt)); colormap(usa2);
  figure(2); clf; plot(h.vchan,nanmean(junk,2)); xlim([640 1640]); ylim([-1 +1]); title(num2str(iCnt)); plotaxis2;
  for ii = 1 : 2645
    P = polyfit(daysSince2002,junk(ii,:),1);
    quicktrend(ii) = P(1)*365;
  end
  figure(3); clf; plot(h.vchan,quicktrend); xlim([640 1640]); ylim([-1 +1]/10); ylabel('BT trend'); title(num2str(iCnt)); plotaxis2;
  trendavg(iCnt,:) = quicktrend;
  anomavg(iCnt,:) = nanmean(junk,2);
  pause(1);
end
figure(4); clf; pcolor(h.vchan,meanvaluebin(newLatGrid),anomavg(1:length(newLatGrid)-1,:)); colormap(usa2); colorbar; caxis([-1 +1]); xlim([640 1640]); shading flat; 
  title('mean spectral anomaly')
figure(5); clf; pcolor(h.vchan,meanvaluebin(newLatGrid),anomavg(1:length(newLatGrid)-1,:)/(iNumYears/2)); colormap(usa2); colorbar; caxis([-1 +1]/(iNumYears/2)); xlim([640 1640]); shading flat; 
  title(['mean spectral anomaly over ' num2str(iNumYears) ' years'])

iCnt = 1;
  str = ['junk = btavg' num2str(iCnt) ';'];
  eval(str);
  figure(1); clf; pcolor(2002+daysSince2002/365,h.vchan,junk); colorbar; caxis([-1 +1]*2); shading flat; ylim([640 1640]); title(num2str(iCnt)); colormap(usa2);
  figure(2); clf; plot(h.vchan,nanmean(junk,2)); xlim([640 1640]); ylim([-1 +1]); title(num2str(iCnt)); plotaxis2;
  for ii = 1 : 2645
    P = polyfit(daysSince2002,junk(ii,:),1);
    quicktrend(ii) = P(1)*365;
  end
  figure(3); clf; plot(h.vchan,quicktrend); xlim([640 1640]); ylim([-1 +1]/10); ylabel('BT trend'); title(num2str(iCnt)); plotaxis2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iSave = input('save -1/+1  : ');
if iSave > 0
  comment = 'see make_global_avg_and_9_averages.m';
  fout = ['anomaly_globalavg_and_' num2str(length(newLatGrid)-1) '_averages_timeseries_Q' num2str(iQuant,'%02d') '.mat'];
  saver = ['save -v7.3 ' fout ' btavgAnomFinal yy mm dd hh rtime comment usethese newLatGrid anomavg'];
  eval(saver);
end
