meanwv = squeeze(nanmean(wv(iaTropics,:,:),1));
meantz = squeeze(nanmean(tz(iaTropics,:,:),1));
meano3 = squeeze(nanmean(o3(iaTropics,:,:),1));

%{
%https://www.mathworks.com/matlabcentral/fileexchange/49985-pcolor3
addpath /home/sergio/MATLABCODE/PLOTTER/PCOLOR3/pcolor3
pcolor3(okdates,latbins,playsN,tz);  caxis([-0.5 +0.5]); colormap(usa2); colorbar; title('T anom'); 
  set(gca,'zdir','reverse'); xlabel('time'); ylabel('latitude'); zlabel('pressure')
  axis([2002.75 2018.75 -90 +90 10 1000])
%}

%{
xdiff = wv;
xdiff = tz;
xdiff(abs(xdiff) > 1)=nan;                              % added line
h = slice(xdiff, [], [], 1:size(xdiff,3));
set(h, 'EdgeColor','none', 'FaceColor','interp')
alpha(.1)
colormap(usa2); colorbar; caxis([-0.5 +0.5]); set(gca,'zdir','reverse'); 
ylabel('latbin'); xlabel('TimeStep'); zlabel('press lay');
%}

clear smooth*

for ii = 1 : nlays
  smoothwv(:,ii) = smooth(meanwv(:,ii),2*5);
  smoothtz(:,ii) = smooth(meantz(:,ii),2*5);
  smootho3(:,ii) = smooth(meano3(:,ii),2*5);
end

addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/COLORMAP

for ii = 1 : 40
  for jj = 1 : nlays
    junk = polyfit(1:365,squeeze(wv(ii,:,jj)),1); trendwv(ii,jj) = junk(1)*365/16;  %% junk is per 16 days, so in 1 year (365 days) we convert!
    junk = polyfit(1:365,squeeze(tz(ii,:,jj)),1); trendtz(ii,jj) = junk(1)*365/16;
    junk = polyfit(1:365,squeeze(o3(ii,:,jj)),1); trendo3(ii,jj) = junk(1)*365/16;

    %[B, err, stats]=Math_tsfit_lin_robust(x,y,n);
    [B,err,stats] = Math_tsfit_lin_robust((1:365)*16,squeeze(wv(ii,:,jj)),4); trendwv2(ii,jj) = B(2); trendwv2_unc(ii,jj) = err(2);
    [B,err,stats] = Math_tsfit_lin_robust((1:365)*16,squeeze(tz(ii,:,jj)),4); trendtz2(ii,jj) = B(2); trendtz2_unc(ii,jj) = err(2);
    [B,err,stats] = Math_tsfit_lin_robust((1:365)*16,squeeze(o3(ii,:,jj)),4); trendo32(ii,jj) = B(2); trendo32_unc(ii,jj) = err(2);

    %[B, err, stats]=Math_tsfit_lin_robust(x,y,n)
    %[B,err,stats] = Math_tsfit_lin_robust((1:365)*16,smooth(squeeze(wv(ii,:,jj)),2*5*2),4); trendwv2(ii,jj) = B(2); trendwv2_unc(ii,jj) = err(2);
    %[B,err,stats] = Math_tsfit_lin_robust((1:365)*16,smooth(squeeze(tz(ii,:,jj)),2*5*2),4); trendtz2(ii,jj) = B(2); trendtz2_unc(ii,jj) = err(2);
    %[B,err,stats] = Math_tsfit_lin_robust((1:365)*16,smooth(squeeze(o3(ii,:,jj)),2*5*2),4); trendo32(ii,jj) = B(2); trendo32_unc(ii,jj) = err(2);

  end
end

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN1 = plevs(1:100)-plevs(2:101);
playsD1 = log(plevs(1:100)./plevs(2:101));
plays = playsN1./playsD1;
playsx = ones(1,103); playsx(4:103) = plays; playsx(3) = plays(4)*1.01; playsx(2) = playsx(3)*1.01; playsx(1) = playsx(2)*1.01;
for ii = 1 : nlays
  playsN(ii) = mean(playsx(a.jacobian.wvjaclays_used{ii}));
end
playsN = fliplr(playsN);

i400 = find(playsN <= 400); i400 = i400(end);
meanwv_400mb = squeeze(nanmean(wv(:,:,i400),3));
meantz_400mb = squeeze(nanmean(tz(:,:,i400),3));
meano3_400mb = squeeze(nanmean(o3(:,:,i400),3));
for ii = 1 : 40
  smoothwv_400mb(ii,:) = smooth(meanwv_400mb(ii,:),2*5);
  smoothtz_400mb(ii,:) = smooth(meantz_400mb(ii,:),2*5);
  smootho3_400mb(ii,:) = smooth(meano3_400mb(ii,:),2*5);
end
disp('no smoothing of 400 mb anoms')
for ii = 1 : 40
  smoothwv_400mb(ii,:) = meanwv_400mb(ii,:);
  smoothtz_400mb(ii,:) = meantz_400mb(ii,:);
  smootho3_400mb(ii,:) = meano3_400mb(ii,:);
end

quick_anom_T_WV_plots
