if ~exist('desc')
  [hERAI,ha,pERAI,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/summary_17years_all_lat_all_lon_2002_2019.rtp');
  desc = load('ERA5_atm_N_cld_data_2002_09_to_2024_08_desc.mat');
  desc.yymm = desc.all.yy + (desc.all.mm-1)/12;
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
ind2 = 2514;
smN = 5;

ix = input('Enter 300,500,800 mb) : ');
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

figure(1); clf
plot(desc.yymm,35+smooth(squeeze(desc.all.ptemp(:,ix,ind1)),smN),'m',desc.yymm,smooth(squeeze(desc.all.stemp(:,ind1)),smN),'r',...
     desc.yymm,smooth(squeeze(desc.all.gas_1(:,ix,ind1)),smN)/4e20 + 300,'b',desc.yymm,smooth(squeeze(desc.all.RH(:,ix,ind1)),smN)/4+295,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Tile ' num2str(ind1) ' at ' str]);

figure(2); clf
plot(desc.yymm,35+smooth(squeeze(desc.all.ptemp(:,ix,ind2)),smN),'m',desc.yymm,smooth(squeeze(desc.all.stemp(:,ind2)),smN),'r',...
     desc.yymm,smooth(squeeze(desc.all.gas_1(:,ix,ind2)),smN)/4e20 + 300,'b',desc.yymm,smooth(squeeze(desc.all.RH(:,ix,ind2)),smN)/4+290,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Tile ' num2str(ind2) ' at ' str]);

figure(3); clf
plot(asc.yymm,35+smooth(squeeze(asc.all.ptemp(:,ix,ind1)),smN),'m',asc.yymm,smooth(squeeze(asc.all.stemp(:,ind1)),smN),'r',...
     asc.yymm,smooth(squeeze(asc.all.gas_1(:,ix,ind1)),smN)/4e20 + 300,'b',asc.yymm,smooth(squeeze(asc.all.RH(:,ix,ind1)),smN)/4+295,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc Tile ' num2str(ind1) ' at ' str]);

figure(4); clf
plot(asc.yymm,35+smooth(squeeze(asc.all.ptemp(:,ix,ind2)),smN),'m',asc.yymm,smooth(squeeze(asc.all.stemp(:,ind2)),smN),'r',...
    asc.yymm,smooth(squeeze(asc.all.gas_1(:,ix,ind2)),smN)/4e20 + 300,'b',asc.yymm,smooth(squeeze(asc.all.RH(:,ix,ind2)),smN)/4+290,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc Tile ' num2str(ind2) ' at ' str]);

%%%%%%%%%%%%%%%%%%%%%%%%%

xblahWV = squeeze(desc.all.gas_1(:,ix,:)); xblahRH = squeeze(desc.all.RH(:,ix,:)); xblahT = squeeze(desc.all.ptemp(:,ix,:)); xblahST = squeeze(desc.all.stemp(:,:)); 
%blahWV = sum(blahWV.*coslatall,2)./sum(coslatall,2); blahRH = sum(blahRH.*coslatall,2)./sum(coslatall,2); blahT = sum(blahT.*coslatall,2)./sum(coslatall,2); blahST = sum(blahST.*coslatall,2)./sum(coslatall,2);
for ii = 1 : 262
  moo = xblahWV(ii,:); blahWV(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahRH(ii,:); blahRH(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahT(ii,:); blahT(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahST(ii,:); blahST(ii) = nansum(moo.*coslat)/nansum(coslat);
end
figure(5); clf
plot(asc.yymm,smooth(blahT,smN)+20,'m',asc.yymm,smooth(blahST,smN)+10,'r',asc.yymm,smooth(blahWV,smN)/4e20 + 290,'b',asc.yymm,smooth(blahRH,smN)*5+15,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Tile Global Average at ' str]);

xblahWV = squeeze(asc.all.gas_1(:,ix,:)); xblahRH = squeeze(asc.all.RH(:,ix,:)); xblahT = squeeze(asc.all.ptemp(:,ix,:)); xblahST = squeeze(asc.all.stemp(:,:)); 
%blahWV = sum(blahWV.*coslatall,2)./sum(coslatall,2); blahRH = sum(blahRH.*coslatall,2)./sum(coslatall,2); blahT = sum(blahT.*coslatall,2)./sum(coslatall,2); blahST = sum(blahST.*coslatall,2)./sum(coslatall,2);
for ii = 1 : 262
  moo = xblahWV(ii,:); blahWV(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahRH(ii,:); blahRH(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahT(ii,:); blahT(ii) = nansum(moo.*coslat)/nansum(coslat);
  moo = xblahST(ii,:); blahST(ii) = nansum(moo.*coslat)/nansum(coslat);
end
figure(6); clf
plot(asc.yymm,smooth(blahT,smN)+20,'m',asc.yymm,smooth(blahST,smN)+10,'r',asc.yymm,smooth(blahWV,smN)/4e20 + 290,'b',asc.yymm,smooth(blahRH,smN)*5+15,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc Tile Global Average at ' str]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atlantic1 = find(pERAI.landfrac == 0 & pERAI.rlat > 10   & pERAI.rlat < 40 & pERAI.rlon > -90 & pERAI.rlon < 0);
atlantic2 = find(pERAI.landfrac == 0 & pERAI.rlat > -40 & pERAI.rlat <= 10 & pERAI.rlon > -75 & pERAI.rlon < 10);
atlantic = union(atlantic1,atlantic2);
figure(7); clf; plot(pERAI.rlon,pERAI.rlat,'b.',pERAI.rlon(atlantic),pERAI.rlat(atlantic),'r.')
figure(7); clf; simplemap(pERAI.rlat,pERAI.rlon,pERAI.stemp); hold on; plot(pERAI.rlon(atlantic),pERAI.rlat(atlantic),'r.')'; hold off

xblahWV = squeeze(desc.all.gas_1(:,ix,atlantic)); xblahRH = squeeze(desc.all.RH(:,ix,atlantic)); xblahT = squeeze(desc.all.ptemp(:,ix,atlantic)); xblahST = squeeze(desc.all.stemp(:,atlantic)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahT = nanmean(xblahT,2);  blahST = nanmean(xblahST,2); 
figure(7); clf
plot(desc.yymm,smooth(blahT,smN)+20,'m',desc.yymm,smooth(blahST,smN)+10,'r',desc.yymm,smooth(blahWV,smN)/4e20 + 290,'b',desc.yymm,smooth(blahRH,smN)*4+105,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc Atlantic at ' str]);

xblahWV = squeeze(asc.all.gas_1(:,ix,atlantic)); xblahRH = squeeze(asc.all.RH(:,ix,atlantic)); xblahT = squeeze(asc.all.ptemp(:,ix,atlantic)); xblahST = squeeze(asc.all.stemp(:,atlantic)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahT = nanmean(xblahT,2);  blahST = nanmean(xblahST,2); 
figure(8); clf
plot(asc.yymm,smooth(blahT,smN)+20,'m',asc.yymm,smooth(blahST,smN)+10,'r',asc.yymm,smooth(blahWV,smN)/4e20 + 290,'b',asc.yymm,smooth(blahRH,smN)*4+105,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc Atlantic at ' str]);

%%%%%%%%%%%%%%%%%%%%%%%%%

twp1 = find(pERAI.landfrac == 0 & pERAI.rlat > -30 & pERAI.rlat < +30 & pERAI.rlon > +90 & pERAI.rlon <= 180);
twp2 = find(pERAI.landfrac == 0 & pERAI.rlat > -30 & pERAI.rlat < +30 & pERAI.rlon >= -180 & pERAI.rlon <= -75);
twp = union(twp1,twp2);

xblahWV = squeeze(desc.all.gas_1(:,ix,twp)); xblahRH = squeeze(desc.all.RH(:,ix,twp)); xblahT = squeeze(desc.all.ptemp(:,ix,twp)); xblahST = squeeze(desc.all.stemp(:,twp)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahT = nanmean(xblahT,2);  blahST = nanmean(xblahST,2); 
figure(9); clf
plot(desc.yymm,smooth(blahT,smN)+20,'m',desc.yymm,smooth(blahST,smN)+10,'r',desc.yymm,smooth(blahWV,smN)/4e20 + 290,'b',desc.yymm,smooth(blahRH,smN)*4+105,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc TWP at ' str]);

xblahWV = squeeze(asc.all.gas_1(:,ix,twp)); xblahRH = squeeze(asc.all.RH(:,ix,twp)); xblahT = squeeze(asc.all.ptemp(:,ix,twp)); xblahST = squeeze(asc.all.stemp(:,twp)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahT = nanmean(xblahT,2);  blahST = nanmean(xblahST,2); 
figure(10); clf
plot(asc.yymm,smooth(blahT,smN)+20,'m',asc.yymm,smooth(blahST,smN)+10,'r',asc.yymm,smooth(blahWV,smN)/4e20 + 290,'b',asc.yymm,smooth(blahRH,smN)*4+105,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc TWP at ' str]);

%%%%%%%%%%%%%%%%%%%%%%%%%

tropical_land = find(pERAI.landfrac == 1 & pERAI.rlat > -30 & pERAI.rlat < +30);

xblahWV = squeeze(desc.all.gas_1(:,ix,tropical_land)); xblahRH = squeeze(desc.all.RH(:,ix,tropical_land)); xblahT = squeeze(desc.all.ptemp(:,ix,tropical_land)); xblahST = squeeze(desc.all.stemp(:,tropical_land)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahT = nanmean(xblahT,2);  blahST = nanmean(xblahST,2); 
figure(11); clf
plot(desc.yymm,smooth(blahT,smN)+20,'m',desc.yymm,smooth(blahST,smN)+10,'r',desc.yymm,smooth(blahWV,smN)/4e20 + 290,'b',desc.yymm,smooth(blahRH,smN)*4+105,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Desc TROPICAL\_LAND at ' str]);

xblahWV = squeeze(asc.all.gas_1(:,ix,tropical_land)); xblahRH = squeeze(asc.all.RH(:,ix,tropical_land)); xblahT = squeeze(asc.all.ptemp(:,ix,tropical_land)); xblahST = squeeze(asc.all.stemp(:,tropical_land)); 
blahWV = nanmean(xblahWV,2);  blahRH = nanmean(xblahRH,2);  blahT = nanmean(xblahT,2);  blahST = nanmean(xblahST,2); 
figure(12); clf
plot(asc.yymm,smooth(blahT,smN)+20,'m',asc.yymm,smooth(blahST,smN)+10,'r',asc.yymm,smooth(blahWV,smN)/4e20 + 290,'b',asc.yymm,smooth(blahRH,smN)*4+105,'c','linewidth',2)
xlim([2020 2025]); plotaxis2;
hl = legend('Tz(X mb)','STemp','WV(X mb)','RH(X mb)','location','best','fontsize',10); title(['Asc TROPICAL\_LAND at ' str]);

