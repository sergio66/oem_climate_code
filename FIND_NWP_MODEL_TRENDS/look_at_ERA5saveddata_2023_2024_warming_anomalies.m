if ~exist('desc')
  [hERAI,ha,pERAI,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/summary_17years_all_lat_all_lon_2002_2019.rtp');
  desc = load('ERA5_atm_N_cld_data_2002_09_to_2024_08_desc.mat');
  desc.yymm = desc.all.yy + (desc.all.mm-1)/12;
  doy = change2days(desc.all.yy,desc.all.mm,ones(size(desc.all.yy))*15,2002);
  plevs = flipud(load('/home/sergio/MATLABCODE/airslevels.dat'));
  plays = meanvaluebin(plevs);
  do_XX_YY_from_X_Y
  coslat = cos(YY*pi/180);
  coslatall = ones(262,1) * cos(YY*pi/180);
end

if ~exist('asc')
  asc = load('ERA5_atm_N_cld_data_2002_09_to_2024_08_asc.mat');
  asc.yymm = asc.all.yy + (asc.all.mm-1)/12;
end

i300 = find(plays >= 300,1);
i500 = find(plays >= 500,1);
i800 = find(plays >= 800,1);

ind1 = 2000; 
ind1 = 2692; 
ind2 = 2514;

smN = 5;
smN = 5;
smN = 7;

ix = input('Enter 300,500,800 mb [default 500] ) : ');
if length(ix) == 0
  ix = 500;
end
if ix == 300
  ix = i300;
  str = '300 mb';
elseif ix == 500
  ix = i500;
  str = '500 mb';
elseif ix == 800
  ix = i800;
  str = '800 mb';
end

addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
%[x_anom b stats] = generic_compute_anomaly(dtime,r,k,iWhich,iNterms);
%[B,stats,btanomaly,radanomaly] = compute_anomaly_wrapper(k,x0,y0,N,f,iRad_or_OD,iDebug)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('testing the code : ')
figure(1); clf
plot(desc.yymm,35+smooth(squeeze(desc.all.ptemp(:,ix,ind1)),smN),'m',desc.yymm,smooth(squeeze(desc.all.stemp(:,ind1)),smN),'r.-',...
     desc.yymm,smooth(squeeze(desc.all.gas_1(:,ix,ind1)),smN)/4e20 + 300,'b',desc.yymm,smooth(squeeze(desc.all.RH(:,ix,ind1)),smN)/4+295,'c.-','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Tile ' num2str(ind1) ' at ' str]);

figure(2); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.gas_1(:,ix,ind1))/nanmean(squeeze(desc.all.gas_1(:,ix,ind1))),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.RH(:,ix,ind1)),4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.ptemp(:,ix,ind1)),4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.ptemp(:,ind1)),4,[],-1);
plot(desc.yymm,smooth(anomTZ,smN),'m',desc.yymm,smooth(anomST,smN),'r.-',desc.yymm,smooth(anomWV,smN)*5,'b',desc.yymm,smooth(anomRH,smN)/10,'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Tile ' num2str(ind1) ' at ' str]);

figure(3); clf
[xanomWV,~,~] = generic_compute_anomaly(doy,squeeze(desc.all.gas_1(:,ix,ind1))/nanmean(squeeze(desc.all.gas_1(:,ix,ind1))),1:length(doy),1,4);
[xanomRH,~,~] = generic_compute_anomaly(doy,squeeze(desc.all.RH(:,ix,ind1)),1:length(doy),1,4);
[xanomTZ,~,~] = generic_compute_anomaly(doy,squeeze(desc.all.ptemp(:,ix,ind1)),1:length(doy),1,4);
[xanomST,~,~] = generic_compute_anomaly(doy,squeeze(desc.all.ptemp(:,ind1)),1:length(doy),1,4);
 plot(desc.yymm,smooth(xanomTZ,smN),'m',desc.yymm,smooth(xanomST,smN),'r.-',desc.yymm,smooth(xanomWV,smN)*5,'b',desc.yymm,smooth(xanomRH,smN)/10,'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Tile ' num2str(ind1) ' at ' str]);

figure(4); clf
plot(desc.yymm,smooth(anomTZ-xanomTZ,smN),'m.-',desc.yymm,smooth(anomST-xanomST,smN),'m.-',desc.yymm,smooth(anomWV-xanomWV,smN),'b.-',desc.yymm,smooth(anomRH-xanomRH,smN),'c.-','linewidth',2);
xlim([2020 2025]); 
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Tile ' num2str(ind1) ' at ' str]);

figure(5); clf
plot(desc.yymm,smooth(squeeze(desc.all.RH(:,ix,ind1)),smN)-mean(squeeze(desc.all.RH(:,ix,ind1)))); title('RH')
xlim([2020 2025]); plotaxis2;

addpath /home/sergio/MATLABCODE/PLOTTER

figure(6); clf
plot_72x64_tiles([2000 2514 2692],ones(1,4608));

disp('RET to continue'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.gas_1(:,ix,ind1))/nanmean(squeeze(desc.all.gas_1(:,ix,ind1))),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.RH(:,ix,ind1)),4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.ptemp(:,ix,ind1)),4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.ptemp(:,ind1)),4,[],-1);
plot(desc.yymm,smooth(anomTZ,smN),'m',desc.yymm,smooth(anomST,smN),'r.-',desc.yymm,smooth(anomWV,smN)*5,'b',desc.yymm,smooth(anomRH,smN)/10,'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Tile ' num2str(ind1) ' at ' str]);

figure(2); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.gas_1(:,ix,ind2))/nanmean(squeeze(desc.all.gas_1(:,ix,ind2))),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.RH(:,ix,ind2)),4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.ptemp(:,ix,ind2)),4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(desc.all.ptemp(:,ind2)),4,[],-1);
plot(desc.yymm,smooth(anomTZ,smN),'m',desc.yymm,smooth(anomST,smN),'r.-',desc.yymm,smooth(anomWV,smN)*5,'b',desc.yymm,smooth(anomRH,smN)/10,'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Tile ' num2str(ind2) ' at ' str]);

figure(3); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(asc.all.gas_1(:,ix,ind1))/nanmean(squeeze(asc.all.gas_1(:,ix,ind1))),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(asc.all.RH(:,ix,ind1)),4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(asc.all.ptemp(:,ix,ind1)),4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(asc.all.ptemp(:,ind1)),4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN)/10,'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc Tile ' num2str(ind1) ' at ' str]);

figure(4); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(asc.all.gas_1(:,ix,ind2))/nanmean(squeeze(asc.all.gas_1(:,ix,ind2))),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(asc.all.RH(:,ix,ind2)),4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(asc.all.ptemp(:,ix,ind2)),4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,squeeze(asc.all.ptemp(:,ind2)),4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN)/10,'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc Tile ' num2str(ind2) ' at ' str]);

%%%%%%%%%%%%%%%%%%%%%%%%%
clear *blah* *anomST_time_*
xblahWV = squeeze(desc.all.gas_1(:,ix,:)); xblahRH = squeeze(desc.all.RH(:,ix,:)); xblahTZ = squeeze(desc.all.ptemp(:,ix,:)); xblahST = squeeze(desc.all.stemp(:,:)); 
%blahWV = sum(blahWV.*coslatall,2)./sum(coslatall,2); blahRH = sum(blahRH.*coslatall,2)./sum(coslatall,2); blahTZ = sum(blahTZ.*coslatall,2)./sum(coslatall,2); blahST = sum(blahST.*coslatall,2)./sum(coslatall,2);
for ii = 1 : 262
  moo = xblahWV(ii,:); blahWV(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahRH(ii,:); blahRH(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahTZ(ii,:); blahTZ(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahST(ii,:); blahST(ii) = nansum(moo.*coslat)/nansum(coslat);
end
figure(5); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,blahWV/nanmean(blahWV),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,blahRH,4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,blahTZ,4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,blahST,4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN),'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Tile Global Average at ' str]);

clear *blah* *anomST_time_*
xblahWV = squeeze(asc.all.gas_1(:,ix,:)); xblahRH = squeeze(asc.all.RH(:,ix,:)); xblahTZ = squeeze(asc.all.ptemp(:,ix,:)); xblahST = squeeze(asc.all.stemp(:,:)); 
%blahWV = sum(blahWV.*coslatall,2)./sum(coslatall,2); blahRH = sum(blahRH.*coslatall,2)./sum(coslatall,2); blahTZ = sum(blahTZ.*coslatall,2)./sum(coslatall,2); blahST = sum(blahST.*coslatall,2)./sum(coslatall,2);
for ii = 1 : 262
  moo = xblahWV(ii,:); blahWV(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahRH(ii,:); blahRH(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahTZ(ii,:); blahTZ(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahST(ii,:); blahST(ii) = nansum(moo.*coslat)/nansum(coslat);
end
figure(6); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,blahWV/nanmean(blahWV),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,blahRH,4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,blahTZ,4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,blahST,4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN),'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc Tile Global Average at ' str]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atlantic1 = find(pERAI.landfrac == 0 & pERAI.rlat > 10   & pERAI.rlat < 40 & pERAI.rlon > -90 & pERAI.rlon < 0);
atlantic2 = find(pERAI.landfrac == 0 & pERAI.rlat > -40 & pERAI.rlat <= 10 & pERAI.rlon > -75 & pERAI.rlon < 10);

atlantic1 = find(pERAI.landfrac == 0 & pERAI.rlat > 10   & pERAI.rlat < 15 & pERAI.rlon > -90 & pERAI.rlon < 0);
atlantic2 = find(pERAI.landfrac == 0 & pERAI.rlat > -15 & pERAI.rlat <= 10 & pERAI.rlon > -75 & pERAI.rlon < 10);

atlantic = union(atlantic1,atlantic2);
figure(7); clf; plot(pERAI.rlon,pERAI.rlat,'b.',pERAI.rlon(atlantic),pERAI.rlat(atlantic),'r.')
figure(7); clf; simplemap(pERAI.rlat,pERAI.rlon,pERAI.stemp); hold on; plot(pERAI.rlon(atlantic),pERAI.rlat(atlantic),'r.')'; hold off

clear *blah* *anomST_time_*
xblahWV = squeeze(desc.all.gas_1(:,ix,atlantic)); xblahRH = squeeze(desc.all.RH(:,ix,atlantic)); xblahTZ = squeeze(desc.all.ptemp(:,ix,atlantic)); xblahST = squeeze(desc.all.stemp(:,atlantic)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahTZ = nanmean(xblahTZ,2);  blahST = nanmean(xblahST,2); 
figure(7); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,blahWV/nanmean(blahWV),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,blahRH,4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,blahTZ,4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,blahST,4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN),'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Atlantic at ' str]);
for ii = 1 : size(xblahST,2)
  [~,~,anomST_time_D(ii,:)] = compute_anomaly_wrapper(1:length(doy),doy,xblahST(:,ii),4,[],-1);
end
% moo = nan(262,4608); moo(:,atlantic) = anomST_time_D';
% for ii = 1 : 262
%    scatter_coast(pERAI.rlon,pERAI.rlat,50,moo(ii,:)); title(num2str(ii)); caxis([-1 +1]); axis([-180 +180 -90 +90]); 
%   title([num2str(desc.all.yy(ii)) '/' num2str(desc.all.mm(ii))]);
%   pause(0.1); 
% end

clear *blah* *anomST_time_*
xblahWV = squeeze(asc.all.gas_1(:,ix,atlantic)); xblahRH = squeeze(asc.all.RH(:,ix,atlantic)); xblahTZ = squeeze(asc.all.ptemp(:,ix,atlantic)); xblahST = squeeze(asc.all.stemp(:,atlantic)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahTZ = nanmean(xblahTZ,2);  blahST = nanmean(xblahST,2); 
figure(8); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,blahWV/nanmean(blahWV),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,blahRH,4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,blahTZ,4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,blahST,4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN),'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc Atlantic at ' str]);
for ii = 1 : size(xblahST,2)
  [~,~,anomST_time_A(ii,:)] = compute_anomaly_wrapper(1:length(doy),doy,xblahST(:,ii),4,[],-1);
end
% moo = nan(262,4608); moo(:,atlantic) = anomST_time_A';
% for ii = 1 : 262
%   scatter_coast(pERAI.rlon,pERAI.rlat,50,moo(ii,:)); title(num2str(ii)); caxis([-1 +1]); axis([-180 +180 -90 +90]); 
%   title([num2str(desc.all.yy(ii)) '/' num2str(desc.all.mm(ii))]);
%   pause(0.1); 
% end

plot(asc.yymm,nanmean(anomST_time_A,1),asc.yymm,nanmean(anomST_time_D,1)); title('Atlantic anomaly (b) asc (r) desc');
plotaxis2;
pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%

twp1 = find(pERAI.landfrac == 0 & pERAI.rlat > -30 & pERAI.rlat < +30 & pERAI.rlon > +90 & pERAI.rlon <= 180);
twp2 = find(pERAI.landfrac == 0 & pERAI.rlat > -30 & pERAI.rlat < +30 & pERAI.rlon >= -180 & pERAI.rlon <= -75);
twp = union(twp1,twp2);

clear *blah* *anomST_time_*
xblahWV = squeeze(desc.all.gas_1(:,ix,twp)); xblahRH = squeeze(desc.all.RH(:,ix,twp)); xblahTZ = squeeze(desc.all.ptemp(:,ix,twp)); xblahST = squeeze(desc.all.stemp(:,twp)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahTZ = nanmean(xblahTZ,2);  blahST = nanmean(xblahST,2); 
figure(9); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,blahWV/nanmean(blahWV),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,blahRH,4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,blahTZ,4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,blahST,4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN),'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc TWP at ' str]);
for ii = 1 : size(xblahST,2)
  [~,~,anomST_time_D(ii,:)] = compute_anomaly_wrapper(1:length(doy),doy,xblahST(:,ii),4,[],-1);
end

clear *blah* *anomST_time_*
xblahWV = squeeze(asc.all.gas_1(:,ix,twp)); xblahRH = squeeze(asc.all.RH(:,ix,twp)); xblahTZ = squeeze(asc.all.ptemp(:,ix,twp)); xblahST = squeeze(asc.all.stemp(:,twp)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahTZ = nanmean(xblahTZ,2);  blahST = nanmean(xblahST,2); 
figure(10); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,blahWV/nanmean(blahWV),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,blahRH,4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,blahTZ,4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,blahST,4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN),'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc TWP at ' str]);
for ii = 1 : size(xblahST,2)
  [~,~,anomST_time_A(ii,:)] = compute_anomaly_wrapper(1:length(doy),doy,xblahST(:,ii),4,[],-1);
end

plot(asc.yymm,nanmean(anomST_time_A,1),asc.yymm,nanmean(anomST_time_D,1)); title('TWP anomaly (b) asc (r) desc');
plotaxis2;
pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%

tropical_land = find(pERAI.landfrac == 1 & pERAI.rlat > -30 & pERAI.rlat < +30);

xblahWV = squeeze(desc.all.gas_1(:,ix,tropical_land)); xblahRH = squeeze(desc.all.RH(:,ix,tropical_land)); xblahTZ = squeeze(desc.all.ptemp(:,ix,tropical_land)); xblahST = squeeze(desc.all.stemp(:,tropical_land)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahTZ = nanmean(xblahTZ,2);  blahST = nanmean(xblahST,2); 
figure(11); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,blahWV/nanmean(blahWV),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,blahRH,4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,blahTZ,4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,blahST,4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN),'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc TROPICAL\_LAND at ' str]);
for ii = 1 : size(xblahST,2)
  [~,~,anomST_time_D(ii,:)] = compute_anomaly_wrapper(1:length(doy),doy,xblahST(:,ii),4,[],-1);
end

xblahWV = squeeze(asc.all.gas_1(:,ix,tropical_land)); xblahRH = squeeze(asc.all.RH(:,ix,tropical_land)); xblahTZ = squeeze(asc.all.ptemp(:,ix,tropical_land)); xblahST = squeeze(asc.all.stemp(:,tropical_land)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahTZ = nanmean(xblahTZ,2);  blahST = nanmean(xblahST,2); 
figure(12); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,blahWV/nanmean(blahWV),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,blahRH,4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,blahTZ,4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,blahST,4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN),'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc TROPICAL\_LAND at ' str]);
for ii = 1 : size(xblahST,2)
  [~,~,anomST_time_A(ii,:)] = compute_anomaly_wrapper(1:length(doy),doy,xblahST(:,ii),4,[],-1);
end

plot(asc.yymm,nanmean(anomST_time_A,1),asc.yymm,nanmean(anomST_time_D,1)); title('Tropical Land anomaly (b) asc (r) desc');
plotaxis2;
pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coastlines = load('coast');

planet = find(pERAI.landfrac == 1 & pERAI.rlat > -30 & pERAI.rlat < +30);
planet = 1 : 4608;

clear *blah* *anomST_time_*
xblahWV = squeeze(desc.all.gas_1(:,ix,planet)); xblahRH = squeeze(desc.all.RH(:,ix,planet)); xblahTZ = squeeze(desc.all.ptemp(:,ix,planet)); xblahST = squeeze(desc.all.stemp(:,planet)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahTZ = nanmean(xblahTZ,2);  blahST = nanmean(xblahST,2); 
figure(13); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,blahWV/nanmean(blahWV),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,blahRH,4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,blahTZ,4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,blahST,4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN),'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc PLANET at ' str]);
for ii = 1 : size(xblahST,2)
  [~,~,anomST_time_D(ii,:)] = compute_anomaly_wrapper(1:length(doy),doy,xblahST(:,ii),4,[],-1);
end
pcolor(X,Y,reshape(anomST_time_D(:,262),72,64)); shading interp; colorbar; caxis([-1 +1]*20); hold on; plot(coastlines.long,coastlines.lat,'k'); hold off; title('ST DESC PLANET')

xblahWV = squeeze(asc.all.gas_1(:,ix,planet)); xblahRH = squeeze(asc.all.RH(:,ix,planet)); xblahTZ = squeeze(asc.all.ptemp(:,ix,planet)); xblahST = squeeze(asc.all.stemp(:,planet)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahTZ = nanmean(xblahTZ,2);  blahST = nanmean(xblahST,2); 
figure(14); clf
[~,~,anomWV] = compute_anomaly_wrapper(1:length(doy),doy,blahWV/nanmean(blahWV),4,[],-1);
[~,~,anomRH] = compute_anomaly_wrapper(1:length(doy),doy,blahRH,4,[],-1);
[~,~,anomTZ] = compute_anomaly_wrapper(1:length(doy),doy,blahTZ,4,[],-1);
[~,~,anomST] = compute_anomaly_wrapper(1:length(doy),doy,blahST,4,[],-1);
plot(asc.yymm,smooth(anomTZ,smN),'m',asc.yymm,smooth(anomST,smN),'r.-',asc.yymm,smooth(anomWV,smN)*5,'b',asc.yymm,smooth(anomRH,smN),'c.-','linewidth',2);
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc PLANET at ' str]);
for ii = 1 : size(xblahST,2)
  [~,~,anomST_time_A(ii,:)] = compute_anomaly_wrapper(1:length(doy),doy,xblahST(:,ii),4,[],-1);
end
pcolor(X,Y,reshape(anomST_time_A(:,262),72,64)); shading interp; colorbar; caxis([-1 +1]*20); hold on; plot(coastlines.long,coastlines.lat,'k'); hold off; title('SKT ASC PLANET')

%for ii = 1 : size(xblahST,2)
for ii = length(desc.all.mm) - 2*12 : length(desc.all.mm)
  figure(13); clf; pcolor(X,Y,reshape(anomST_time_D(:,ii),72,64)); shading interp; colorbar; caxis([-1 +1]*20); hold on; plot(coastlines.long,coastlines.lat,'k'); hold off; 
    title(['ST DESC PLANET ' num2str(desc.all.yy(ii)) '/' num2str(desc.all.mm(ii))])
  figure(14); clf; pcolor(X,Y,reshape(anomST_time_A(:,ii),72,64)); shading interp; colorbar; caxis([-1 +1]*20); hold on; plot(coastlines.long,coastlines.lat,'k'); hold off; 
    title(['SKT ASC PLANET ' num2str(desc.all.yy(ii)) '/' num2str(desc.all.mm(ii))])
  disp('RET to continue'); pause;
  pause(0.1)
end

plot(asc.yymm,nanmean(anomST_time_A,1),asc.yymm,nanmean(anomST_time_D,1)); title('Tropical Land anomaly (b) asc (r) desc');
plotaxis2;
pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
