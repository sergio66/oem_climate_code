ratesQ00 = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q00.mat','rates');
stemp00 = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q00.mat','results');

ratesQ16 = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16.mat','rates');
stemp16 = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16.mat','results');

%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

addpath /home/sergio/MATLABCODE/matlib/science/            %% for usgs_deg10_dem.m that has correct paths
[salti, landfrac] = usgs_deg10_dem(Y(:),X(:));
Ylat = Y(:);
Xlon = X(:);

figure(1); clf; scatter_coast(Xlon,Ylat,50,stemp00.results(:,6)); caxis([-0.15 +0.15]); colormap(usa2); title('Mean Q00 dST/dt')
figure(2); clf; scatter_coast(Xlon,Ylat,50,stemp16.results(:,6)); caxis([-0.15 +0.15]); colormap(usa2); title('Hottest Q16 dST/dt')
tropics = find(abs(Ylat) < 30);
figure(3); plot(1:2645,nanmean(ratesQ00.rates,2),1:2645,nanmean(ratesQ16.rates,2)); hl = legend('Mean','Hottest'); title('Whole Planet');
figure(4); plot(1:2645,nanmean(ratesQ00.rates(:,tropics),2),1:2645,nanmean(ratesQ16.rates(:,tropics),2)); hl = legend('Mean','Hottest'); title('Tropics')
