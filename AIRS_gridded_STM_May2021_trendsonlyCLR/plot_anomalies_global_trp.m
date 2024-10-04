junk = iNumAnomTimeSteps * iNumAnomTiles;
fprintf(1,'found %5i of %5i \n',sum(iaFound),length(iaFound))
figure(1); clf; waha = reshape(iaFound,iNumAnomTimeSteps,iNumAnomTiles);                                                   
  pcolor(yymm,1:2,waha'); shading interp;   colorbar; set(gca,'ydir','normal');  
  title([num2str(sum(iaFound(:))) ' / ' num2str(junk)  ' = ' num2str(100*sum(iaFound(:))/junk) ' % made so far']);  
  xlabel('Time'); ylabel('Latitude'); colormap(jett); 
  title('Anomaly Files made')

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = 10;

iFig = iFig + 1; figure(iFig); clf
junk = results(1:iNumAnomTimeSteps,6); 
plot(yymm,smooth(junk,15)); title('Global ST(t)'); plotaxis2;

PX = nanpolyfit(yymm,junk,1); YXval = polyval(PX,daysSince2002/365+2002); 
plot(onidd,oni/10,'k',yymm,smooth(junk,23*iNumYearSmooth),'b',...
     yymm,smooth(junk-YXval',23*iNumYearSmooth),'r','linewidth',2); 
  title('Global Avg Stemp anomaly'); ylim([-1 +1]*1)
xlim([2002 2025]); ; plotaxis2; hl = legend('Ocean Nino Index','UMBC Global anomaly w/ trend','UMBC Global anomaly w/o trend','location','best','fontsize',10);

iFig = iFig + 1; figure(iFig); clf; wahaT = resultsT(1:iNumAnomTimeSteps,:)'; 
  pcolor(yymm,pavg,wahaT);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Global AVG T(z,t)');      
  pcolor(yymm,pavg,smoothn(wahaT,1));  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Global AVG T(z,t)');      
  colormap(llsmap5); caxis([-1 +1]*2); set(gca,'yscale','log'); ylim([10 1000])

iFig = iFig + 1; figure(iFig); clf; wahaWV = resultsWV(1:iNumAnomTimeSteps,:)';
  pcolor(yymm,pavg,wahaWV);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Global Avg WVfrac(z,t)'); 
  pcolor(yymm,pavg,smoothn(wahaWV,1));  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Global Avg WVfrac(z,t)'); 
  colormap(llsmap5); caxis([-1 +1]*0.15); ylim([100 1000])
  ylim([50 1000]); set(gca,'yscale','log');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = iFig + 1; figure(iFig); clf
junk = results((1:iNumAnomTimeSteps)+iNumAnomTimeSteps,6); 
plot(yymm,smooth(junk,15)); title('Tropical ST(t)'); plotaxis2;

PX = nanpolyfit(yymm,junk,1); YXval = polyval(PX,daysSince2002/365+2002); 
plot(onidd,oni/10,'k',yymm,smooth(junk,23*iNumYearSmooth),'b',...
     yymm,smooth(junk-YXval',23*iNumYearSmooth),'r','linewidth',2); 
  title('Tropical Stemp anomaly'); ylim([-1 +1]*1)
xlim([2002 2025]); ; plotaxis2; hl = legend('Ocean Nino Index','UMBC Tropical anomaly w/ trend','UMBC Tropical anomaly w/o trend','location','best','fontsize',10);

iFig = iFig + 1; figure(iFig); clf; wahaT = resultsT((1:iNumAnomTimeSteps)+iNumAnomTimeSteps,:)'; 
  pcolor(yymm,pavg,wahaT);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Tropical AVG T(z,t)');      
  pcolor(yymm,pavg,smoothn(wahaT,1));  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Tropical AVG T(z,t)');      
  colormap(llsmap5); caxis([-1 +1]*2); set(gca,'yscale','log'); ylim([10 1000])

iFig = iFig + 1; figure(iFig); clf; wahaWV = resultsWV((1:iNumAnomTimeSteps)+iNumAnomTimeSteps,:)';
  pcolor(yymm,pavg,wahaWV);  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Tropical Avg WVfrac(z,t)'); 
  pcolor(yymm,pavg,smoothn(wahaWV,1));  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Tropical Avg WVfrac(z,t)'); 
  colormap(llsmap5); caxis([-1 +1]*0.15); ylim([100 1000])
  ylim([50 1000]); set(gca,'yscale','log');

