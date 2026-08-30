addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/KCARTA/MATLAB

latbins = equal_area_spherical_bands(20);
xlat = (latbins(1:end-1)+latbins(2:end))/2;
tropics = find(abs(xlat) <= 30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yyStart = 2012; yyEnd = 2017;
yyStart = 2002; yyEnd = 2017;

yyStart = 2002; yyEnd = 2019;

iNewOrOld = +1;  %% steves new stats
iNewOrOld = -1;  %% steve/strow old stats

robs = [];
rclr = [];
rtime = [];

xrobs = [];
xrclr = [];
xrtime = [];

crobs = [];
crclr = [];
crtime = [];

iNewOrOld = +1;  %% steves new stats
if iNewOrOld > 0
  for yy0 = yyStart : yyEnd
    fprintf(1,'Steves new AIRS stats yy = %4i \n',yy0)
    fin = ['/asl/data/stats/airs/clear/rtp_airicrad_era_rad_kl_' num2str(yy0) '_clear_desc_ocean.mat'];
    a = load(fin);
    xrobs = cat(1,xrobs,a.robs);
    xrclr = cat(1,xrclr,a.rclr);
    xrtime = cat(1,xrtime,a.rtime_mean);
  end
end

iNewOrOld = -1;  %% steve/strow old stats 
if iNewOrOld < 0
  for ii = 1 : 40
    fprintf(1,'Strow/Steves old AIRS stats latbin = %4i \n',ii)

    %% day
    finX = ['/home/strow/Work/Airs/Stability/Data/L1c_v1/Asc/statlat' num2str(ii) '.mat'];
    finY = ['/home/strow/Work/Airs/Stability/Data/L1c_v1/Asc/nucal_statlat' num2str(ii) '.mat'];

    %% night
    %finX = ['/home/strow/Work/Airs/Stability/Data/Desc/statlat' num2str(ii) '.mat'];
    %finY = ['/home/strow/Work/Airs/Stability/Data/Desc/nucal_statlat' num2str(ii) '.mat'];
    finX = ['/home/strow/Work/Airs/Stability/Data/L1c_v1/Desc/statlat' num2str(ii) '.mat'];
    finY = ['/home/strow/Work/Airs/Stability/Data/L1c_v1/Desc/nucal_statlat' num2str(ii) '.mat'];

    aX = load(finX);
    aY = load(finY);
    aX.robs = aY.robs;
    robs(ii,:,:) = aX.robs;
    rclr(ii,:,:) = aX.rclr;
    rtime(ii,:) = aX.rtime_mean;
  end
  rtime = rtime';
  robs = permute(robs,[2 1 3]);
  rclr = permute(rclr,[2 1 3]);
end

for yy0 = 2012 : 2019
  fprintf(1,'Steves CRIS NSR stats yy = %4i \n',yy0)
  fin = ['/asl/data/stats/cris/clear/rtp_cris_lowres_rad_' num2str(yy0) '_clear_desc_ocean.mat'];
  a = load(fin);
  crobs = cat(1,crobs,a.robs);
  crclr = cat(1,crclr,a.rclr);
  crtime = cat(1,crtime,a.rtime_mean);
end

xcrobs = squeeze(nanmean(crobs,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whos robs xrobs rtime xrtome
addpath /home/sergio/MATLABCODE/TIME

[yy,mm,dd,hhjunk]  = tai2utcSergio(rtime(:,20));  daysSince2002 = change2days(yy,mm,dd,2002);
[yyx,mmx,ddx,hhx] = tai2utcSergio(xrtime(:,20)); xdaysSince2002 = change2days(yyx,mmx,ddx,2002);
[yyc,mmc,ddc,hhc] = tai2utcSergio(crtime(:,20)); cdaysSince2002 = change2days(yyc,mmc,ddc,2002);

load /home/sergio/MATLABCODE/QUICKTASKS_TELECON/AIRS2CRIS_and_IASI2CRIS/sergio_airs2cris_decon/hstructure_lores.mat
%% the breaks are at 717 and 1154
iC_chan_use = [[3:715] [720:1152] [1157:1315]];

%% this is for 1317 chan
iC_1231 = 737;
iC_900  = 403;
iC_790  = 227;
iC_792  = 230;

%% this is for 1305 chans used in the anomaly files
iC_1231 = find(hh.vchan(iC_chan_use) >= 1231,1);   iC_1231 = 731;
iC_900  = find(hh.vchan(iC_chan_use) >= 900,1);    iC_900 = 401;  
iC_790  = find(hh.vchan(iC_chan_use) >= 790,1);    iC_790 = 225;  
iC_792  = find(hh.vchan(iC_chan_use) >= 791.5,1);  iC_792 = 228;  
iC_745  = find(hh.vchan(iC_chan_use) >= 745,1);    iC_745 = 153;
iC_735  = find(hh.vchan(iC_chan_use) >= 735,1);    iC_735 = 137;

%% this is for 1317 chans used in the raw spectra stats files
iC_1231 = find(hh.vchan >= 1231,1);
iC_900 = find(hh.vchan >= 900,1);
iC_790  = find(hh.vchan >= 790,1); 
iC_792  = find(hh.vchan >= 791.5,1);
iC_745  = find(hh.vchan >= 745,1); 
iC_735  = find(hh.vchan >= 735,1); 

plot(squeeze(nanmean(xcrobs(:,tropics,:),1))','.-')
plot(hh.vchan,squeeze(nanmean(xcrobs(:,tropics,:),1))','.-')

%load_fairs
load f2645.mat; fairs = f2645;
i1231 = find(fairs >= 1231,1);   % 1520
i900 = find(fairs >= 900,1);     %  795
i745  = find(fairs >= 745,1);    %  349
i735  = find(fairs >= 735,1);    %  317
i790 = find(fairs >= 790.5,1);   %  486
i792 = find(fairs >= 791.5,1);   %  489

plot(hh.vchan,nanmean(squeeze(nanmean(xcrobs(:,tropics,:),1))),'b',fairs,nanmean(squeeze(nanmean(robs(:,tropics,:),1))),'r')

[~,iA,iC]   = intersect(daysSince2002,cdaysSince2002);
[~,iAx,iCx] = intersect(xdaysSince2002,cdaysSince2002);
plot(hh.vchan,rad2bt(hh.vchan,nanmean(squeeze(nanmean(xcrobs(iC,tropics,:),1)))),'b',...
     fairs,rad2bt(fairs,nanmean(squeeze(nanmean(robs(iA,tropics,:),1)))),'r',...
     fairs,rad2bt(fairs,nanmean(squeeze(nanmean(xrobs(iAx,tropics,:),1)))),'k')

meanBT1231_A  = rad2bt(1231,squeeze(nanmean(robs(iA,tropics,i1231),2)));
meanBT1231_C = rad2bt(1231,squeeze(nanmean(xcrobs(iC,tropics,iC_1231),2)));
meanBT1231_Ax = rad2bt(1231,squeeze(nanmean(xrobs(iAx,tropics,i1231),2)));
meanBT1231_Cx = rad2bt(1231,squeeze(nanmean(xcrobs(iCx,tropics,iC_1231),2)));

meanBT900_A = rad2bt(900,squeeze(nanmean(robs(iA,tropics,i900),2)));
meanBT900_C = rad2bt(900,squeeze(nanmean(xcrobs(iC,tropics,iC_900),2)));
meanBT900_Ax = rad2bt(900,squeeze(nanmean(xrobs(iAx,tropics,i900),2)));
meanBT900_Cx = rad2bt(900,squeeze(nanmean(xcrobs(iCx,tropics,iC_900),2)));

plot(daysSince2002(iA)/365+2002,meanBT1231_A,cdaysSince2002(iC)/365+2002,meanBT1231_C)
plot(daysSince2002(iA)/365+2002,meanBT900_A,cdaysSince2002(iC)/365+2002,meanBT900_C)

plot(daysSince2002(iA)/365+2002,smooth(meanBT1231_A-meanBT1231_C,365),'b',...
     daysSince2002(iA)/365+2002,smooth(meanBT900_A-meanBT900_C,365),'r')
hl = legend('1231 cm-1','900 cm-1','location','best'); 
ylabel('Tropical (AIRS L1C-CRIS NSR)'); xlabel('time');

plot(xdaysSince2002(iAx)/365+2002,smooth(meanBT1231_Ax-meanBT1231_Cx,365),'b',...
     xdaysSince2002(iAx)/365+2002,smooth(meanBT900_Ax-meanBT900_Cx,365),'r')
hl = legend('1231 cm-1','900 cm-1','location','best'); 
ylabel('Tropical (AIRS L1C-CRIS NSR)'); xlabel('time');
disp('ret to cont'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); 
plot(rtime(1:5549,20)-xrtime(1:5549,20)); title('difference in rtime for lat bins');
  nansum(rtime(1:5549,20)-xrtime(1:5549,20))

plot(xdaysSince2002-daysSince2002(1:5549))
  nansum(xdaysSince2002-daysSince2002(1:5549))

plot(squeeze(robs(1:5549,20,i1231))-squeeze(xrobs(1:5549,20,i1231)))
  nansum(squeeze(robs(1:5549,20,i1231))-squeeze(xrobs(1:5549,20,i1231)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(daysSince2002/365+2002,smooth(squeeze(robs(:,20,i900)),365),'b',...
     xdaysSince2002/365+2002,smooth(squeeze(xrobs(:,20,i900)),365),'r',...
     cdaysSince2002/365+2002,smooth(squeeze(xcrobs(:,20,iC_900)),365),'k',...
     cdaysSince2002/365+2002,smooth(squeeze(crobs(:,20,5,iC_900)),365),'g','linewidth',2);
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','CRIS center','location','best');
set(hl,'fontsize',10); title('rad 900');
grid

plot(daysSince2002/365+2002,rad2bt(900,smooth(squeeze(robs(:,20,i900)),365)),'b',...
     xdaysSince2002/365+2002,rad2bt(900,smooth(squeeze(xrobs(:,20,i900)),365)),'r',...
     cdaysSince2002/365+2002,rad2bt(900,smooth(squeeze(xcrobs(:,20,iC_900)),365)),'k',...
     cdaysSince2002/365+2002,rad2bt(900,smooth(squeeze(crobs(:,20,5,iC_900)),365)),'g','linewidth',2);
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','CRIS center','location','best');
set(hl,'fontsize',10); title('BT 900');
grid

plot(daysSince2002/365+2002,smooth(squeeze(robs(:,20,i1231)),365),'b',...
     xdaysSince2002/365+2002,smooth(squeeze(xrobs(:,20,i1231)),365),'r',...
     cdaysSince2002/365+2002,2+smooth(squeeze(xcrobs(:,20,iC_1231)),365),'k',...
     cdaysSince2002/365+2002,2+smooth(squeeze(crobs(:,20,5,iC_1231)),365),'g','linewidth',2);
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','CRIS center','location','best');
set(hl,'fontsize',10); title('rad 1231');
grid

plot(daysSince2002/365+2002,smooth(squeeze(robs(:,20,i745)),365),'b',...
     xdaysSince2002/365+2002,smooth(squeeze(xrobs(:,20,i745)),365),'r',...
     cdaysSince2002/365+2002,-5+smooth(squeeze(xcrobs(:,20,iC_745)),365),'k',...
     cdaysSince2002/365+2002,-5+smooth(squeeze(crobs(:,20,5,iC_745)),365),'g','linewidth',2);
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','CRIS center','location','best');
set(hl,'fontsize',10); title('rad 745');
grid

plot(daysSince2002/365+2002,smooth(squeeze(robs(:,20,i735)),365),'b',...
     xdaysSince2002/365+2002,smooth(squeeze(xrobs(:,20,i735)),365),'r',...
     cdaysSince2002/365+2002,+2+smooth(squeeze(xcrobs(:,20,iC_735)),365),'k',...
     cdaysSince2002/365+2002,+2+smooth(squeeze(crobs(:,20,5,iC_735)),365),'g','linewidth',2);
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','CRIS center','location','best');
set(hl,'fontsize',10); title('rad 735');
grid

plot(daysSince2002/365+2002,smooth(squeeze(robs(:,20,i790))-squeeze(robs(:,20,i792)),365),'b',...
     xdaysSince2002/365+2002,smooth(squeeze(xrobs(:,20,i790))-squeeze(xrobs(:,20,i792)),365),'r',...
     cdaysSince2002/365+2002,8+smooth(squeeze(xcrobs(:,20,iC_790)-xcrobs(:,20,iC_792)),365),'k',...
     cdaysSince2002/365+2002,8+smooth(squeeze(crobs(:,20,5,iC_790)-crobs(:,20,5,iC_792)),365),'g','linewidth',2);
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','CRIS center','location','best');
set(hl,'fontsize',10); title('rad 790-792');
grid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%note 2012,06,01 instead of 2012,05,01
radanom  = quickanom(squeeze(robs(:,20,i735)),yy,mm,dd,2012,06,01);
xradanom = quickanom(squeeze(xrobs(:,20,i735)),yyx,mmx,ddx,2012,06,01);
xcradanom = quickanom(squeeze(xcrobs(:,20,iC_735)),yyc,mmc,ddc,2012,06,01);

plot(daysSince2002/365+2002,smooth(radanom,180),'b',xdaysSince2002/365+2002,smooth(xradanom,180),'r',...
     cdaysSince2002/365+2002,smooth(xcradanom,180),'k','linewidth',2)
grid; axis([2002 2020 -1 +1])
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','location','best');
set(hl,'fontsize',10); title('rad 735 anom');

radanom   = quickanom(mean(squeeze(robs(:,tropics,i735)),2),yy,mm,dd,2012,06,01);
xradanom  = quickanom(mean(squeeze(xrobs(:,tropics,i735)),2),yyx,mmx,ddx,2012,06,01);
xcradanom = quickanom(mean(squeeze(xcrobs(:,tropics,iC_735)),2),yyc,mmc,ddc,2012,06,01);

plot(daysSince2002/365+2002,0.3+smooth(radanom,180),'b',xdaysSince2002/365+2002,0.3+smooth(xradanom,180),'r',...
     cdaysSince2002/365+2002,0.3+smooth(xcradanom,180),'k','linewidth',2)
grid; axis([2002 2020 -1 +1])
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','location','best');
set(hl,'fontsize',10); title('rad 735 anom');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radanom  = quickanom(squeeze(robs(:,20,i1231)),yy,mm,dd,2012,06,01);
xradanom = quickanom(squeeze(xrobs(:,20,i1231)),yyx,mmx,ddx,2012,06,01);
xcradanom = quickanom(squeeze(xcrobs(:,20,iC_1231)),yyc,mmc,ddc,2012,06,01);

plot(daysSince2002/365+2002,smooth(radanom,180),'b',xdaysSince2002/365+2002,smooth(xradanom,180),'r',...
     cdaysSince2002/365+2002,smooth(xcradanom,180),'k','linewidth',2)
grid; axis([2002 2020 -1 +1])
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','location','best');
set(hl,'fontsize',10); title('rad 1231 anom');

radanom   = quickanom(mean(squeeze(robs(:,tropics,i1231)),2),yy,mm,dd,2012,06,01);
xradanom  = quickanom(mean(squeeze(xrobs(:,tropics,i1231)),2),yyx,mmx,ddx,2012,06,01);
xcradanom = quickanom(mean(squeeze(xcrobs(:,tropics,iC_1231)),2),yyc,mmc,ddc,2012,06,01);

plot(daysSince2002/365+2002,-0.2+smooth(radanom,180),'b',xdaysSince2002/365+2002,-0.2+smooth(xradanom,180),'r',...
     cdaysSince2002/365+2002,0.0+smooth(xcradanom,180),'k','linewidth',2)
grid; axis([2002 2020 -0.5 +0.5])
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','location','best');
set(hl,'fontsize',10); title('rad 1231 anom');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radanom  = quickanom(squeeze(robs(:,20,i900)),yy,mm,dd,2012,06,01);
xradanom = quickanom(squeeze(xrobs(:,20,i900)),yyx,mmx,ddx,2012,06,01);
xcradanom = quickanom(squeeze(xcrobs(:,20,iC_900)),yyc,mmc,ddc,2012,06,01);

plot(daysSince2002/365+2002,smooth(radanom,180),'b',xdaysSince2002/365+2002,smooth(xradanom,180),'r',...
     cdaysSince2002/365+2002,smooth(xcradanom,180),'k','linewidth',2)
grid; axis([2002 2020 -1 +1])
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','location','best');
set(hl,'fontsize',10); title('rad 900 anom');

radanom   = quickanom(mean(squeeze(robs(:,tropics,i900)),2),yy,mm,dd,2012,06,01);
xradanom  = quickanom(mean(squeeze(xrobs(:,tropics,i900)),2),yyx,mmx,ddx,2012,06,01);
xcradanom = quickanom(mean(squeeze(xcrobs(:,tropics,iC_900)),2),yyc,mmc,ddc,2012,06,01);

plot(daysSince2002/365+2002,-0.2+smooth(radanom,180),'b',xdaysSince2002/365+2002,-0.2+smooth(xradanom,180),'r',...
     cdaysSince2002/365+2002,0.0+smooth(xcradanom,180),'k','linewidth',2)
grid; axis([2002 2020 -0.5 +0.5])
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','location','best');
set(hl,'fontsize',10); title('rad 900 anom');

radanom   = mean(squeeze(robs(:,tropics,i900)),2);
xradanom  = mean(squeeze(xrobs(:,tropics,i900)),2);
xcradanom = mean(squeeze(xcrobs(:,tropics,iC_900)),2);

plot(daysSince2002/365+2002,rad2bt(900,smooth(radanom,180)),'b',...
     xdaysSince2002/365+2002,rad2bt(900,smooth(xradanom,180)),'r',...
     cdaysSince2002/365+2002,rad2bt(900,smooth(xcradanom,180)),'k','linewidth',2)
grid; axis([2002 2020 294 296])
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','location','best');
set(hl,'fontsize',10); title('BT 900');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radanom  = quickanom(squeeze(robs(:,20,i790)-robs(:,20,i792)),yy,mm,dd,2012,06,01);
xradanom = quickanom(squeeze(xrobs(:,20,i790)-xrobs(:,20,i792)),yyx,mmx,ddx,2012,06,01);
xcradanom = quickanom(squeeze(xcrobs(:,20,iC_790)-crobs(:,20,i792)),yyc,mmc,ddc,2012,06,01);

%radanom  = quickanom(squeeze(robs(:,tropics,i790)-robs(:,tropics,i792)),yy,mm,dd,2012,06,01);
%xradanom = quickanom(squeeze(xrobs(:,tropics,i790)-xrobs(:,tropics,i792)),yyx,mmx,ddx,2012,06,01);
%xcradanom = quickanom(squeeze(xcrobs(:,tropics,iC_790)-crobs(:,tropics,iC_792)),yyc,mmc,ddc,2012,06,01);
%radanom = mean(radanom,2); xradanom = mean(xradanom,2); xcradanom = mean(xcradanom,2);

radanom   = quickanom(mean(squeeze(robs(:,tropics,i790)-robs(:,tropics,i792)),2),yy,mm,dd,2012,06,01);
xradanom  = quickanom(mean(squeeze(xrobs(:,tropics,i790)-xrobs(:,tropics,i792)),2),yyx,mmx,ddx,2012,06,01);
xcradanom = quickanom(mean(squeeze(xcrobs(:,tropics,iC_790)-xcrobs(:,tropics,iC_792)),2),yyc,mmc,ddc,2012,06,01);

plot(daysSince2002/365+2002,-0.25+smooth(radanom,180),'b',xdaysSince2002/365+2002,-0.25+smooth(xradanom,180),'r',...
     cdaysSince2002/365+2002,-0.1+smooth(xcradanom,180),'k','linewidth',2)
grid; axis([2002 2020 -1 +1])
hl = legend('Strow AIRS','Sergio/Steve redone AIRS','CRIS mean3x3','location','best');
set(hl,'fontsize',10); title('rad 790-792 anom');

