junk = iNumAnomTimeSteps * iNumAnomTiles;
fprintf(1,'found %5i of %5i \n',sum(iaFound),length(iaFound))
figure(1); clf; waha = reshape(iaFound,iNumAnomTimeSteps,iNumAnomTiles);                                                   
  pcolor(daysSince2002,rlat,waha'); shading interp;   colorbar; set(gca,'ydir','normal');  
  title([num2str(sum(iaFound(:))) ' / ' num2str(junk)  ' = ' num2str(100*sum(iaFound(:))/junk) ' % made so far']);  
  xlabel('Longitude'); ylabel('Latitude'); colormap(jett); 

if ~exist('oni')
  oni = load('ONI_sep2023.txt');
  [aaa,bbb] = size(oni);
  oniS = oni(1,1); oniE = oni(aaa,1);
  oni = oni(1:aaa,2:bbb); oni = oni'; oni = oni(:);
  onidd = 1:length(oni); onidd = (onidd-1)/12 + oniS;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iNumYearSmooth = 5;
iNumYearSmooth = 3;
iNumYearSmooth = 1;
iNumYearSmooth = 0.25;

figure(2); clf
junk = results(1:iNumAnomTimeSteps,6); 
PX = polyfit(daysSince2002/365+2002,junk,1); YXval = polyval(PX,daysSince2002/365+2002); 
plot(onidd,oni/10,'k',daysSince2002/365+2002,smooth(junk,23*iNumYearSmooth),'b',daysSince2002/365+2002,smooth(junk-YXval',23*iNumYearSmooth),'r','linewidth',2); 
xlim([2002 2023]); ; plotaxis2; hl = legend('Ocean Nino Index','UMBC Global anomaly w/ trend','UMBC Global anomaly w/o trend','location','best','fontsize',10);
  title('Global Avg Stemp anomaly'); 

figure(3); clf; wahaT = resultsT(1:iNumAnomTimeSteps,:)'; 
  pcolor(daysSince2002/365+2002,pavg,wahaT);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Global AVG T(z,t)');      
  colormap(llsmap5); caxis([-1 +1]*2); set(gca,'yscale','log'); ylim([10 1000])
figure(4); clf; wahaWV = resultsWV(1:iNumAnomTimeSteps,:)';
  pcolor(daysSince2002/365+2002,pavg,wahaWV);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Global Avg WVfrac(z,t)'); 
  colormap(llsmap5); caxis([-1 +1]*0.1); ylim([100 1000])
figure(5); clf
  for ii = 1 : 49
    P = polyfit(daysSince2002/365,wahaT(ii,:),1); trendTavg(ii) = P(1);
    P = polyfit(daysSince2002/365,wahaWV(ii,:),1); trendWVavg(ii) = P(1);
  end
  semilogy(trendTavg,pavg,'r',trendWVavg,pavg,'b','linewidth',2); plotaxis2; 
  hl = legend('Global dT/dt','Global dWVfrac/dt','location','best'); set(gca,'ydir','reverse'); 
  ylim([10 1000])

  subplot(121); semilogy(trendTavg, pavg,'r','linewidth',2); title('Global dT/dt'); set(gca,'ydir','reverse'); plotaxis2; ylim([10 1000])  
  subplot(122); semilogy(trendWVavg,pavg,'b','linewidth',2); title('Global dWVfrac/dt'); set(gca,'ydir','reverse'); plotaxis2; ylim([10 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
pcolor(daysSince2002/365+2002,rlat(2:end),(reshape(results(iNumAnomTimeSteps+1:end,6),iNumAnomTimeSteps,iNumAnomTiles-1))'); colorbar; 
caxis([-1 +1]*2); colormap(llsmap5); shading interp; title('Anomaly Stemp(t,lat)')

figure(7); clf; waha = squeeze(nanmean(reshape(resultsT,iNumAnomTimeSteps,iNumAnomTiles,iNumLay),1)); waha = waha'; waha = waha(:,2:length(rlat));
  pcolor(rlat(2:end),pavg,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('<UMBC> T(z,lat)');      colormap(llsmap5); caxis([-1 +1]*5)
figure(8); clf; waha = squeeze(nanmean(reshape(resultsWV,iNumAnomTimeSteps,iNumAnomTiles,iNumLay),1)); waha = waha';waha = waha(:,2:length(rlat));       
  pcolor(rlat(2:end),pavg,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('<UMBC> WVfrac(z,lat)'); colormap(llsmap5); caxis([-1 +1]*0.5)

clear trendT trendWV
waha = reshape(resultsT,iNumAnomTimeSteps,iNumAnomTiles,iNumLay);
for jj = 1 : length(pavg)
  for ii = 1 : length(rlat)-1
    junk = squeeze(waha(:,ii+1,jj)); 
    P = nanpolyfit(daysSince2002/365,junk,1); trendT(ii,jj) = P(1);
  end
end
waha = reshape(resultsWV,iNumAnomTimeSteps,iNumAnomTiles,iNumLay);
for jj = 1 : length(pavg)
  for ii = 1 : length(rlat)-1
    junk = squeeze(waha(:,ii+1,jj)); 
    P = nanpolyfit(daysSince2002/365,junk,1); trendWV(ii,jj) = P(1);
  end
end
trendT  = trendT';
trendWV = trendWV';
figure(09); pcolor(rlat(2:end),pavg,trendT); title('dT/dt'); colormap(llsmap5); caxis([-1 +1]*0.15); 
  shading interp; colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
figure(10); pcolor(rlat(2:end),pavg,trendWV);title('dWVfrac/dt'); colormap(llsmap5); caxis([-1 +1]*0.015)
  shading interp; colorbar; set(gca,'ydir','reverse'); set(gca,'yscale','linear'); ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iaYears = 2002 : 2023;
iaYears = 2003 : 2022;
clear annualT annualWV
waha = reshape(resultsT,iNumAnomTimeSteps,iNumAnomTiles,iNumLay);
for jj = 1 : length(pavg)
  for ii = 1 : length(rlat)-1
    junk = squeeze(waha(:,ii+1,jj)); 
    annualT(ii,jj,:) = interp1(2002+daysSince2002/365,junk,iaYears,[],'extrap');
  end
end

waha = reshape(resultsWV,iNumAnomTimeSteps,iNumAnomTiles,iNumLay);
for jj = 1 : length(pavg)
  for ii = 1 : length(rlat)-1
    junk = squeeze(waha(:,ii+1,jj)); 
    annualWV(ii,jj,:) = interp1(2002+daysSince2002/365,junk,iaYears,[],'extrap');
  end
end

annualT0  = annualT;
annualWV0 = annualWV;

%figure(11); clf
%for ii = 1 : length(iaYears) 
%  pcolor(rlat(2:end),pavg,squeeze(annualT(:,:,ii))'); colorbar; caxis([-1 +1]*2); colormap(usa2); shading interp; set(gca,'ydir','reverse'); 
%  title(num2str(iaYears(ii)))
%  pause(0.1); 
%end

%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z] = ndgrid(double(rlat(2:end)),double(pavg),double(iaYears));

yslice = [-75:25:+75];
xslice = [100 200 500 800];
zslice = [1:5:length(iaYears)]; zslice = iaYears(zslice);

figure(11); clf
slice(Y,X,Z,double(annualT), xslice, yslice, zslice);
slice(Y,X,Z,double(annualT), [], [], zslice);
shading interp; colormap jet; colorbar; caxis([-2 +2]); colormap(usa2);
zlabel('Pav(mb)'); ylabel('Latitude'); zlabel('Time'); title('Tanomaly (K)')

figure(12); clf
slice(Y,X,Z,double(annualWV), xslice, yslice, zslice);
slice(Y,X,Z,double(annualWV), [], [], zslice);
shading interp; colormap jet; colorbar; caxis([-2 +2]/10); colormap(usa2);
zlabel('Pav(mb)'); ylabel('Latitude'); zlabel('Time'); title('WVanomaly (frac)')

%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z] = ndgrid(double(iaYears),double(rlat(2:end)),double(pavg));
annualWV = permute(annualWV,[3 1 2]);
annualT = permute(annualT,[3 1 2]);

yslice = [-75:25:+75];
zslice = [100 200 500 800];
xslice = [1:5:length(iaYears)]; xslice = iaYears(xslice);

figure(11); clf
h = slice(Y,X,Z,double(annualT), [], xslice, []);
ylabel('Time'); xlabel('Latitude'); zlabel('Pav (mb)'); title('Tanomaly (K)')
shading interp; colormap jet; colorbar; caxis([-1 +1]*1); colormap(usa2); set(gca,'zdir','reverse'); set(gca,'ydir','reverse');
ylim([2002 2023]); xlim([-90 +90]); zlim([100 1000])
set(gca,'zscale','log'); zlim([10 1000])
set(h,'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
% set transparency to correlate to the data values.
alpha('color');
%alpha('scaled');
colormap(usa2);

figure(12); clf;
h = slice(Y,X,Z,double(annualWV), [], xslice, []);
ylabel('Time'); xlabel('Latitude'); zlabel('Pav (mb)'); title('WVanomaly (frac)')
shading interp; colormap jet; colorbar; caxis([-1 +1]/5); colormap(usa2); set(gca,'zdir','reverse'); set(gca,'ydir','reverse');
ylim([2002 2023]); xlim([-90 +90]); zlim([100 1000])
set(h,'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
% set transparency to correlate to the data values.
alpha('color');
%alpha('scaled');
colormap(usa2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(14); clf
isosurface(Y,X,Z,annualT,-0.5)
isosurface(Y,X,Z,annualT,0)
isosurface(Y,X,Z,annualT,0.5)
colormap(llsmap5); caxis([-2 +2])
ylabel('Time'); xlabel('Latitude'); zlabel('Pav (mb)'); title('Tanomaly (K)')
set(gca,'zdir','reverse'); set(gca,'ydir','reverse');

figure(15); clf
lvls = [-2:0.25:+2];
h = contourslice(Y,X,Z,double(annualT), [], xslice, [], lvls);
ylabel('Time'); xlabel('Latitude'); zlabel('Pav (mb)'); title('Tanomaly (K)')
colormap(usa2); colorbar; caxis([-1 +1]*1); set(gca,'zdir','reverse'); set(gca,'ydir','reverse');
ylim([2002 2023]); xlim([-90 +90]); zlim([100 1000])
set(gca,'zscale','log'); zlim([10 1000])
view(3)
grid on
