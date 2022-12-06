addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
[g2378,g2645] = compare_goodchans_2378_2645();   %% shows g is index into f-ABCD and NOT chanID

figure(1)
junk = squeeze(b_desc(:,:,1520));
scatter_coast(X(:),Y(:),50,junk(:)); pause(0.1);
caxis([-0.2 +0.2]);; colormap(usa2)
title(['BT1231 rate Quantile ' num2str(iQuantile,'%02d')])

figure(2); clf
junk = rad2bt(1231,squeeze(mean_rad(:,:,1520)));
scatter_coast(X(:),Y(:),50,junk(:));
colormap jet
title(['rad2bt(1231,rad1231) Quantile ' num2str(iQuantile,'%02d')])

figure(2); clf
junk = squeeze(mean_BT(:,:,1520));
scatter_coast(X(:),Y(:),50,junk(:));
colormap jet
title(['BT1231 Quantile ' num2str(iQuantile,'%02d')])

figure(3); clf
junk = squeeze(airs_noiseTtrue(:,:,1520));
junk = reshape(airs_noiseTtrue,72*64,2645)'/sqrt(120);    %% need sqrt(N) from about 12000 obs/tile/16 days .. so 1% of this is 120 ... noise goes down by sqrt(N)
plot(h.vchan,nanmean(reshape(b_err_desc,72*64,2645)',2),h.vchan(g2645),nanmean(junk(g2645,:),2));  
  ylim([0 0.3]); 
  hl = legend('from b_err','from 1/sqrt(N)','location','best','fontsize',10);
title(['Quantile ' num2str(iQuantile,'%02d')])

figure(4); clf
junkA = reshape(b_asc,4608,2645);
junkD = reshape(b_desc,4608,2645);
plot(h.vchan,nanmean(junkA,1),h.vchan,nanmean(junkD,1),h.vchan,nanstd(junkA,1),'--',h.vchan,nanstd(junkD,1),'--')
  ylim([-0.1 +0.1]*0.75)
  xlim([640 1640])
  hl = legend('mean ratesA','mean ratesD','std ratesA','std ratesD','location','best','fontsize',10);
title(['Quantile ' num2str(iQuantile,'%02d')])

load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2; 
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;
junk = load('h2645structure.mat');
f           = junk.h.vchan;

i1419 = find(f >= 1419,1);
i1231 = find(f >= 1231,1);
i0900 = find(f >= 0900,1);

load llsmap5
figure(5); clf;
moo = b_desc(:,:,i1231);
aslmap(5,rlat65,rlon73,smoothn(moo',1), [-90 +90],[-180 +180]); title(['dBT1231/dt Quantile ' num2str(iQuantile,'%02d')]); caxis([-1 +1]*0.15); colormap(llsmap5)

if exist('b_cal_desc')
  figure(6); clf;
  moo = b_cal_desc(:,:,i1231);
  aslmap(6,rlat65,rlon73,smoothn(moo',1), [-90 +90],[-180 +180]); title(['dBT1231/dt ERA5 calc']); caxis([-1 +1]*0.15); colormap(llsmap5)

  figure(7); clf;
  moo = b_err_desc(:,:,i1231);
  aslmap(7,rlat65,rlon73,smoothn(moo',1), [-90 +90],[-180 +180]); title(['unc dBT1231/dt Quantile ' num2str(iQuantile,'%02d')]); caxis([0 +1]*0.25); colormap(llsmap5)

  figure(8); clf;
  moo = b_cal_err_desc(:,:,i1231);
  aslmap(8,rlat65,rlon73,smoothn(moo',1), [-90 +90],[-180 +180]); title(['unc dBT1231/dt ERA5 calc']); caxis([0 +1]*0.25); colormap(llsmap5)

end
  