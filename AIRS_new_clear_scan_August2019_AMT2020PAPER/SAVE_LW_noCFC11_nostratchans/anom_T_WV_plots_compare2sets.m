meanwv = squeeze(nanmean(wv(iaTropics,:,:),1));
meantz = squeeze(nanmean(tz(iaTropics,:,:),1));
meano3 = squeeze(nanmean(o3(iaTropics,:,:),1));

xmeanwv = squeeze(nanmean(xwv(iaTropics,:,:),1));
xmeantz = squeeze(nanmean(xtz(iaTropics,:,:),1));
xmeano3 = squeeze(nanmean(xo3(iaTropics,:,:),1));

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

  xsmoothwv(:,ii) = smooth(xmeanwv(:,ii),2*5);
  xsmoothtz(:,ii) = smooth(xmeantz(:,ii),2*5);
  xsmootho3(:,ii) = smooth(xmeano3(:,ii),2*5);

end

addpath /home/sergio/MATLABCODE/FIND_TRENDS
addpath /home/sergio/MATLABCODE/COLORMAP

for ii = 1 : 40
  for jj = 1 : nlays
    %% junk is per 16 days, so in 1 year (365 days) we convert!

    junk = polyfit(1:365,squeeze(wv(ii,:,jj)),1); trendwv(ii,jj) = junk(1)*365/16;  
    junk = polyfit(1:365,squeeze(tz(ii,:,jj)),1); trendtz(ii,jj) = junk(1)*365/16;
    junk = polyfit(1:365,squeeze(o3(ii,:,jj)),1); trendo3(ii,jj) = junk(1)*365/16;

    %[B, err, stats]=Math_tsfit_lin_robust(x,y,n)
    B = Math_tsfit_lin_robust((1:365)*16,squeeze(wv(ii,:,jj)),4); trendwv2(ii,jj) = B(2);
    B = Math_tsfit_lin_robust((1:365)*16,squeeze(tz(ii,:,jj)),4); trendtz2(ii,jj) = B(2);
    B = Math_tsfit_lin_robust((1:365)*16,squeeze(o3(ii,:,jj)),4); trendo32(ii,jj) = B(2);

    %[B, err, stats]=Math_tsfit_lin_robust(x,y,n)
    %B = Math_tsfit_lin_robust((1:365)*16,smooth(squeeze(wv(ii,:,jj)),2*5*2),4); trendwv2(ii,jj) = B(2);
    %B = Math_tsfit_lin_robust((1:365)*16,smooth(squeeze(tz(ii,:,jj)),2*5*2),4); trendtz2(ii,jj) = B(2);
    %B = Math_tsfit_lin_robust((1:365)*16,smooth(squeeze(o3(ii,:,jj)),2*5*2),4); trendo32(ii,jj) = B(2);

    %%%%%%%%%%%%%%%%%%%%%%%%%

    junk = polyfit(1:365,squeeze(xwv(ii,:,jj)),1); xtrendwv(ii,jj) = junk(1)*365/16;  
    junk = polyfit(1:365,squeeze(xtz(ii,:,jj)),1); xtrendtz(ii,jj) = junk(1)*365/16;
    junk = polyfit(1:365,squeeze(xo3(ii,:,jj)),1); xtrendo3(ii,jj) = junk(1)*365/16;

    %[B, err, stats]=Math_tsfit_lin_robust(x,y,n)
    B = Math_tsfit_lin_robust((1:365)*16,squeeze(xwv(ii,:,jj)),4); xtrendwv2(ii,jj) = B(2);
    B = Math_tsfit_lin_robust((1:365)*16,squeeze(xtz(ii,:,jj)),4); xtrendtz2(ii,jj) = B(2);
    B = Math_tsfit_lin_robust((1:365)*16,squeeze(xo3(ii,:,jj)),4); xtrendo32(ii,jj) = B(2);

    %[B, err, stats]=Math_tsfit_lin_robust(x,y,n)
    %B = Math_tsfit_lin_robust((1:365)*16,smooth(squeeze(xwv(ii,:,jj)),2*5*2),4); xtrendwv2(ii,jj) = B(2);
    %B = Math_tsfit_lin_robust((1:365)*16,smooth(squeeze(xtz(ii,:,jj)),2*5*2),4); xtrendtz2(ii,jj) = B(2);
    %B = Math_tsfit_lin_robust((1:365)*16,smooth(squeeze(xo3(ii,:,jj)),2*5*2),4); xtrendo32(ii,jj) = B(2);

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

quick_anom_T_WV_plots_compare2sets
