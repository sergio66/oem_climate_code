addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/

clear all

iYorX = -1; %% keep trend in data
i2or3 = 3;  %% use climate index as is

iYorX = +1; %% remove trend from data
i2or3 = 3;  %% use climate index anomaly

anom_X = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/anomaly_chID_1520_Q03.mat');  %% BT 1231
%anom_X = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/anomaly_chID_2025_Q03.mat');  %% BT 1519
%anom_X = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/anomaly_chID_1861_Q03.mat');  %% BT 1419

miaow = anom_X.btanom;
for ii = 1 : 4608
  lala = miaow(ii,:);
  P = polyfit(1:454,lala,1);
  miaow(ii,:) = lala - polyval(P,1:454);
end

[xL, xEOFs, xEC, xerror] = detrend_time_series_make_EOF_v1(anom_X.btanom');
[yL, yEOFs, yEC, yerror] = detrend_time_series_make_EOF_v1(miaow');

iYorX
if iYorX < 0
  EC = xEC;
else
  EC = yEC;
end
figure(1); semilogy(1:454,xL,'b.-',1:454,yL,'r')

wah = ones(454,1) * max(abs(EC));
EC = EC ./ wah;

%% 30 days apart  
oni = load('oni_index.txt');
soi = load('soi_index.txt');
nao = load('nao_index.txt');
pdo = load('pdo_index.txt');

dtime = (1:length(oni))*30; [oni(:,3) b stats] = generic_compute_anomaly(dtime,oni(:,2));
dtime = (1:length(soi))*30; [soi(:,3) b stats] = generic_compute_anomaly(dtime,soi(:,2));
dtime = (1:length(nao))*30; [nao(:,3) b stats] = generic_compute_anomaly(dtime,nao(:,2));
dtime = (1:length(pdo))*30; [pdo(:,3) b stats] = generic_compute_anomaly(dtime,pdo(:,2));

thetime = (1:454);
%thetime = (thetime-1)/12 + 2002.75;
thetime = (thetime-1)*16/365 + 2002.75;

for ii = 1 : 12
  figure(1); plot(thetime,smooth(EC(:,ii),23)-polyval(polyfit(thetime,smooth(EC(:,ii),23),1),thetime)','r','linewidth',4);
  hold on; plot(oni(:,1),10*smooth(oni(:,i2or3),4),'b',soi(:,1),10*smooth(soi(:,i2or3),4),'c','linewidth',2); hold off
  xlim([2002 2023]); 
  title(num2str(ii))
  plotaxis2;
  hl = legend('EOF timeseies','ONI','SOI','location','best','fontsize',10);

  figure(2); plot(thetime,smooth(EC(:,ii),23)-polyval(polyfit(thetime,smooth(EC(:,ii),23),1),thetime)','r','linewidth',4);
  hold on; plot(nao(:,1),10*smooth(nao(:,i2or3),4),'b',pdo(:,1),10*smooth(pdo(:,i2or3),4),'c','linewidth',2); hold off
  xlim([2002 2023]); 
  title(num2str(ii))
  plotaxis2;
  hl = legend('EOF timeseies','NAO','PDO','location','best','fontsize',10);

  figure(3); plot(oni(:,1),10*smooth(oni(:,i2or3),4),soi(:,1),10*smooth(soi(:,i2or3),4),nao(:,1),10*smooth(nao(:,i2or3),4),pdo(:,1),10*smooth(pdo(:,i2or3),4),'linewidth',2); 
  hold on; plot(thetime,smooth(EC(:,ii),23)-polyval(polyfit(thetime,smooth(EC(:,ii),23),1),thetime)','k','linewidth',4); hold off
  xlim([2002 2023]); 
  title(num2str(ii))
  plotaxis2;
  hl = legend('ONI','SOI','NAO','PDO','EOF timeseies','location','best','fontsize',10);

  boo0 = smooth(EC(:,ii),23)-polyval(polyfit(thetime,smooth(EC(:,ii),23),1),thetime)';
  y1 = smooth(interp1(oni(:,1),oni(:,i2or3),thetime,[],'extrap'),23);
  y2 = smooth(interp1(soi(:,1),soi(:,i2or3),thetime,[],'extrap'),23);
  y3 = smooth(interp1(nao(:,1),nao(:,i2or3),thetime,[],'extrap'),23);
  y4 = smooth(interp1(pdo(:,1),pdo(:,i2or3),thetime,[],'extrap'),23);
  [xcf1,lags] = crosscorr(boo0,y1,NumLags=100);
  [xcf2,lags] = crosscorr(boo0,y2,NumLags=100);
  [xcf3,lags] = crosscorr(boo0,y3,NumLags=100);
  [xcf4,lags] = crosscorr(boo0,y4,NumLags=100);
  figure(4); 
  plot(lags*23/365,[xcf1 xcf2 xcf3 xcf4],'linewidth',2); plotaxis2;
  hl = legend('ONI','SOI','NAO','PDO','location','best','fontsize',10);
  title(num2str(ii))

  boo = find(abs(xcf1) == max(abs(xcf1)),1); themax(ii,1) = xcf1(boo); thelag(ii,1) = boo;
  boo = find(abs(xcf2) == max(abs(xcf2)),1); themax(ii,2) = xcf2(boo); thelag(ii,2) = boo;
  boo = find(abs(xcf3) == max(abs(xcf3)),1); themax(ii,3) = xcf3(boo); thelag(ii,3) = boo;
  boo = find(abs(xcf4) == max(abs(xcf4)),1); themax(ii,4) = xcf4(boo); thelag(ii,4) = boo;

  pause(0.1);
end

thelagg = thelag - (length(lags)+1)/2;
xstrstr = {'ONI','SOI','NAO','PDO'};
figure(1); clf; imagesc(themax); colorbar; colormap(usa2); caxis([-1 +1]); title('Corr Coef'); ylabel('EOF number'); 
  set(gca,'xtick',[1:4],'xticklabel',xstrstr,'fontsize',10);
figure(2); clf; imagesc(thelagg/23); colorbar; colormap(usa2); title('Time of Max Corr (years)');  ylabel('EOF number')
  set(gca,'xtick',[1:4],'xticklabel',xstrstr,'fontsize',10); caxis([-4 +4])
figure(3); clf;  plot(1:12,abs(themax),'linewidth',2); xlabel('EOF'); ylabel('max(abs(corr coef))'); grid
  hl = legend('ONI','SOI','NAO','PDO','location','best','fontsize',10);
