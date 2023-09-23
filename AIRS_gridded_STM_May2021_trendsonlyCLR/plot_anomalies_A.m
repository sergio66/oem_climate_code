if ~exist('oni')
  oni = load('ONI_sep2023.txt');
  [aaa,bbb] = size(oni);
  oniS = oni(1,1); oniE = oni(aaa,1);
  oni = oni(1:aaa,2:bbb); oni = oni'; oni = oni(:);
  onidd = 1:length(oni); onidd = (onidd-1)/12 + oniS;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7)
junk = results(1:iNumAnomTimeSteps,6); 
PX = polyfit(daysSince2002/365+2002,junk,1); YXval = polyval(PX,daysSince2002/365+2002); 
plot(onidd,oni/10,'k',daysSince2002/365+2002,smooth(junk,16*2),'b',daysSince2002/365+2002,smooth(junk-YXval',16*2),'r','linewidth',2); 
xlim([2002 2023]); ; plotaxis2; hl = legend('Ocean Nino Index','UMBC Global anomaly w/ trend','UMBC Global anomaly w/o trend','location','best','fontsize',10);
  title('Global Avg Stemp anomaly'); 

figure(27); clf; waha = resultsT(1:iNumAnomTimeSteps,:)'; 
  pcolor(daysSince2002/365+2002,pavg,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Global AVG T(z,t)');      colormap(llsmap5); caxis([-1 +1]*2); set(gca,'yscale','log'); ylim([10 1000])
figure(28); clf; waha = resultsWV(1:iNumAnomTimeSteps,:)';
  pcolor(daysSince2002/365+2002,pavg,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Global Avg WVfrac(z,t)'); colormap(llsmap5); caxis([-1 +1]*0.1); ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
pcolor(daysSince2002/365+2002,rlat(2:end),(reshape(results(iNumAnomTimeSteps+1:end,6),iNumAnomTimeSteps,iNumAnomTiles-1))'); colorbar; caxis([-1 +1]*2); colormap(llsmap5); shading flat; title('Anomaly Stemp(t,lat)')

figure(29); clf; waha = squeeze(nanmean(reshape(resultsT,iNumAnomTimeSteps,iNumAnomTiles,iNumLay),1)); waha = waha'; waha = waha(:,2:length(rlat));
  pcolor(rlat(2:end),pavg,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('<UMBC> T(z,lat)');      colormap(llsmap5); caxis([-1 +1]*5)
figure(30); clf; waha = squeeze(nanmean(reshape(resultsWV,iNumAnomTimeSteps,iNumAnomTiles,iNumLay),1)); waha = waha';waha = waha(:,2:length(rlat));       
  pcolor(rlat(2:end),pavg,waha);  shading interp; colorbar; set(gca,'ydir','reverse'); title('<UMBC> WVfrac(z,lat)'); colormap(llsmap5); caxis([-1 +1]*0.5)
figure(31); clf; waha = reshape(iaFound,iNumAnomTimeSteps,iNumAnomTiles);                                                   
  pcolor(daysSince2002,rlat,waha'); shading flat;   colorbar; set(gca,'ydir','normal');  
  title([num2str(sum(iaFound(:))) ' / 4608 = ' num2str(100*sum(iaFound(:))/4608) ' % made so far']);  
  xlabel('Longitude'); ylabel('Latitude'); colormap(jett); 
