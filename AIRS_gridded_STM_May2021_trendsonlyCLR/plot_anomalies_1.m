junk = iNumAnomTimeSteps * iNumAnomTiles;
fprintf(1,'found %5i of %5i \n',sum(iaFound),length(iaFound))
figure(1); clf; waha = iaFound;
  plot(daysSince2002,waha); 
  title([num2str(sum(iaFound(:))) ' / ' num2str(junk)  ' = ' num2str(100*sum(iaFound(:))/junk) ' % made so far']);  

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
xlim([2002 2023]); ; plotaxis2; hl = legend('Ocean Nino Index','UMBC Tile anomaly w/ trend','UMBC Tile anomaly w/o trend','location','best','fontsize',10);
xlabel('Time'); ylabel('SurfTemp Anomaly (K)')
  title('Tile Avg Stemp anomaly'); 

figure(3); clf; wahaT = resultsT(1:iNumAnomTimeSteps,:)'; 
  pcolor(daysSince2002/365+2002,pavg,wahaT);             shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Tile AVG T(z,t)');      
  pcolor(daysSince2002/365+2002,pavg,smoothn(wahaT,10));  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Tile AVG T(z,t)');      
  colormap(llsmap5); caxis([-1 +1]*2); set(gca,'yscale','log'); ylim([10 1000])
  xlabel('Time'); ylabel('Pressure(mb)')
figure(4); clf; wahaWV = resultsWV(1:iNumAnomTimeSteps,:)';
  pcolor(daysSince2002/365+2002,pavg,wahaWV);             shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Tile Avg WVfrac(z,t)'); 
  pcolor(daysSince2002/365+2002,pavg,smoothn(wahaWV,10));  shading interp; colorbar; set(gca,'ydir','reverse'); title('UMBC Tile Avg WVfrac(z,t)'); 
  colormap(llsmap5); caxis([-1 +1]*0.2); ylim([100 1000])
  xlabel('Time'); ylabel('Pressure(mb)')

figure(5); clf
  for ii = 1 : 49
    P = polyfit(daysSince2002/365,wahaT(ii,:),1); trendTavg(ii) = P(1);
    P = polyfit(daysSince2002/365,wahaWV(ii,:),1); trendWVavg(ii) = P(1);
  end
  semilogy(trendTavg,pavg,'r',trendWVavg,pavg,'b','linewidth',2); plotaxis2; 
  hl = legend('Tile dT/dt','Tile dWVfrac/dt','location','best'); set(gca,'ydir','reverse'); 
  ylim([10 1000])

  subplot(121); semilogy(trendTavg, pavg,'r','linewidth',2); title('Tile dT/dt'); set(gca,'ydir','reverse'); plotaxis2; ylim([10 1000])  
  subplot(122); semilogy(trendWVavg,pavg,'b','linewidth',2); title('Tile dWVfrac/dt'); set(gca,'ydir','reverse'); plotaxis2; ylim([10 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%
era5   = load('../FIND_NWP_MODEL_TRENDS/ChrisHTrends/era5_newtrends_and_anom_2002_2022.mat');
airsL3 = load('../FIND_NWP_MODEL_TRENDS/ChrisHTrends/airsL3_newtrends_and_anom_2002_2022.mat');
giss   = load('../FIND_NWP_MODEL_TRENDS/ChrisHTrends/giss_trends_plus_anom_2002_2022.mat');
airsL1C = load('anomaly_tile_2515_timeseries_Q03.mat');
time_Giss_airsL3 = 2002.75 + (1:240)/12;

figure(6); clf
i67_35 = (35-1)*72 + 67;
junk = results(1:iNumAnomTimeSteps,6); 
plot(daysSince2002/365+2002,smooth(junk,23*iNumYearSmooth*2),'b',daysSince2002/365+2002,smooth(airsL1C.btavgAnomFinal(1520,:),23*iNumYearSmooth*2),'c',...
     time_Giss_airsL3,smooth(squeeze(giss.giss_anom4608(67,35,:)),23*iNumYearSmooth),'g',time_Giss_airsL3,smooth(squeeze(airsL3.airsL3_anom4608(67,35,:)),23*iNumYearSmooth),'r',...
      time_Giss_airsL3,smooth(era5.era5_anom4608(i67_35,:),23*iNumYearSmooth),'k','linewidth',2)
xlim([2002 2023]); ; plotaxis2; hl = legend('UMBC Tile anomaly','BT1231 AIRS L1C','GISS','AIRS L3','ERA5','location','best','fontsize',10);
  xlabel('Time'); ylabel('SurfTemp Anomaly(K)')
  title('Tile Avg Stemp anomaly'); 

plot(daysSince2002/365+2002,smooth(junk,23*iNumYearSmooth*2),'b',...
     time_Giss_airsL3,smooth(squeeze(giss.giss_anom4608(67,35,:)),23*iNumYearSmooth),'g',time_Giss_airsL3,smooth(squeeze(airsL3.airsL3_anom4608(67,35,:)),23*iNumYearSmooth),'r',...
      time_Giss_airsL3,smooth(era5.era5_anom4608(i67_35,:),23*iNumYearSmooth),'k','linewidth',2)
xlim([2002 2023]); ; plotaxis2; hl = legend('UMBC Tile anomaly','GISS','AIRS L3','ERA5','location','best','fontsize',10);
  xlabel('Time'); ylabel('SurfTemp Anomaly(K)')
  title('Tile Avg Stemp anomaly'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7);
i050 = find(pavg >= 050,1);
i200 = find(pavg >= 200,1);
i500 = find(pavg >= 500,1);
i800 = find(pavg >= 800,1);
junk050 = resultsT(1:iNumAnomTimeSteps,i050); 
junk200 = resultsT(1:iNumAnomTimeSteps,i200); 
junk500 = resultsT(1:iNumAnomTimeSteps,i500); 
junk800 = resultsT(1:iNumAnomTimeSteps,i800); 
plot(daysSince2002/365+2002,smooth(junk050,23*iNumYearSmooth*10),'b',...
     daysSince2002/365+2002,smooth(junk200,23*iNumYearSmooth*10),'g',...
     daysSince2002/365+2002,smooth(junk500,23*iNumYearSmooth*10),'r',...
     daysSince2002/365+2002,smooth(junk800,23*iNumYearSmooth*10),'k',...
'linewidth',2)
xlim([2002 2023]); ; plotaxis2; hl = legend('050 mb','200 mb','500 mb','800 mb','location','best','fontsize',10);
  xlabel('Time'); title('Tile Avg T(z) anomaly'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
addpath /asl/matlib/plotutils

dir0 = '/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';
figure(6); aslprint([dir0 'anomalytimeseries_Q05.pdf'])
%}
