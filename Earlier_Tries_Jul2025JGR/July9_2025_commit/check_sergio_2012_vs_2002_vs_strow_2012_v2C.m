clear all
addpath /home/sergio/MATLABCODE/PLOTTER
latbins = equal_area_spherical_bands(20);
  tropics = find(abs(latbins) <= 30);
latbinsx = 0.5*(latbins(1:end-1) + latbins(2:end));
  tropics = find(abs(latbinsx) <= 30);

for ii = tropics(1) : tropics(end)
  fin = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/fit_robs_lat' num2str(ii) '.mat'];
  bonk = load(fin);
  fprintf(1,' 2002 Strow data, anomaly done by Strow ii=%2i fin=%s \n',ii,fin);
  %plot(smooth(1000*(bonk.all_bt_anom(:,i790)-bonk.all_bt_anom(:,i792)),365))
  strow_moo_obs(ii-10+1,:,:)  = bonk.all_bt_anom(1:5785,:);
  strow_moo_time(ii-10+1,:,:) = bonk.all_times(1:5785,:);
end

yyStart = 2012;
for ii = tropics(1) : tropics(end)
  fin = ['ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_' num2str(yyStart) '/sergio_latbin' num2str(ii) '.mat'];
  bonk = load(fin);
  fprintf(1,' 2002 Strow data, anomaly done by Sergio ii=%2i fin=%s \n',ii,fin);
  %plot(smooth(1000*(bonk.all_bt_anom(:,i790)-bonk.all_bt_anom(:,i792)),365))
  sergio_moo_obs(ii-10+1,:,:)  = bonk.all_bt_anom(1:1813,:);
  sergio_moo_time(ii-10+1,:,:) = bonk.all_times(1:1813,:);
end

yyStart = 2012;
for ii = tropics(1) : tropics(end)
  fin = ['/home/sergio/MATLABCODE/QUICKTASKS_TELECON/AIRS2CRIS_and_IASI2CRIS/sergio_airs2cris_decon/'];
  fin = [fin 'SergioRates/Desc_ocean_fits/fit_robs_lat' num2str(ii) '.mat'];
  bonk = load(fin);
  fprintf(1,' 2012 CRIS anomaly done by Sergio ii=%2i fin=%s \n',ii,fin);
  %plot(smooth(1000*(bonk.all_bt_anom(:,i790)-bonk.all_bt_anom(:,i792)),365))
  cris_moo_obs(ii-10+1,:,:)  = bonk.all_bt_anom(1:2502,:);
  cris_moo_time(ii-10+1,:,:) = bonk.all_times(1:2502,:);
end

%load_fairs
load f2645.mat; fairs = f2645;

i1231 = find(fairs >= 1231,1);
i790 = find(fairs >= 790.5,1);
i792 = find(fairs >= 791.5,1);

load /home/sergio/MATLABCODE/QUICKTASKS_TELECON/AIRS2CRIS_and_IASI2CRIS/sergio_airs2cris_decon/hstructure_lores.mat
%% the breaks are at 717 and 1154
iC_chan_use = [[3:715] [720:1152] [1157:1315]];

iC_1231 = 737;
iC_790  = 227;
iC_792  = 230;
iC_1231 = find(hh.vchan >= 1231,1);   iC_1231 = 731;
iC_790  = find(hh.vchan >= 790,1);    iC_790 = 225;
iC_792  = find(hh.vchan >= 791.5,1);  iC_792 = 228;

canom1231 = nanmean(cris_moo_obs(:,:,iC_1231),1);
canom790 = nanmean(cris_moo_obs(:,:,iC_790),1);
canom792 = nanmean(cris_moo_obs(:,:,iC_792),1);

xanom1231 = nanmean(sergio_moo_obs(:,:,i1231),1);
xanom790 = nanmean(sergio_moo_obs(:,:,i790),1);
xanom792 = nanmean(sergio_moo_obs(:,:,i792),1);

anom1231 = nanmean(strow_moo_obs(:,:,i1231),1);
anom790 = nanmean(strow_moo_obs(:,:,i790),1);
anom792 = nanmean(strow_moo_obs(:,:,i792),1);

figure(1)
plot(2012.31+(1:1813)/365,smooth(xanom1231,365),'b',2002.75+(1:5785)/365,0.02+smooth(anom1231,365),'r',...
     2012.31+(1:2502)/365,0.05+smooth(canom1231,365),'k','linewidth',2);
  title('BT 1231 anomaly'); hl = legend('Sergio','Strow','CRIS','location','best'); set(hl,'fontsize',10);
axis([2012 2019 -0.2 +0.4]); grid

figure(2)
plot(2012.31+(1:1813)/365,0.47+smooth(xanom790-xanom792,365),'b',...
     2002.75+(1:5785)/365,smooth(anom790-anom792,365),'r',...
     2012.31+(1:2502)/365,0.52+smooth(canom790-canom792,365),'k','linewidth',2);
  title('BT 790-792 anomaly'); hl = legend('Sergio','Strow','CRIS','location','best'); set(hl,'fontsize',10);
axis([2012 2019 0.4 +1.2]); grid

figure(3); 
addpath /home/sergio/MATLABCODE/PLOTMISC
plot(fairs,nanmean(squeeze(nanmean(sergio_moo_obs,2)),1),'b',...
     hh.vchan(iC_chan_use),nanmean(squeeze(nanmean(cris_moo_obs,2)),1),'k','linewidth',3)
hold on
shadedErrorBar(fairs,nanmean(squeeze(nanmean(sergio_moo_obs,2)),1),nanstd(squeeze(nanstd(sergio_moo_obs,1,2)),1),'b',0.25);
shadedErrorBar(hh.vchan(iC_chan_use),nanmean(squeeze(nanmean(cris_moo_obs,2)),1),nanstd(squeeze(nanstd(cris_moo_obs,1,2)),1),'k',0.25);
hold off
 axis([640 1640 -0.5 +0.5])
hl = legend('Sergio 2012-2017','CRIS 2012-2019','location','best');

figure(4); 
addpath /home/strow/MATLABCODE/PLOTMISC
plot(fairs,nanmean(squeeze(nanmean(strow_moo_obs,2)),1),'r',...
     hh.vchan(iC_chan_use),nanmean(squeeze(nanmean(cris_moo_obs,2)),1),'k','linewidth',3)
hold on
shadedErrorBar(fairs,nanmean(squeeze(nanmean(strow_moo_obs,2)),1),nanstd(squeeze(nanstd(strow_moo_obs,1,2)),1),'r',0.25);
shadedErrorBar(hh.vchan(iC_chan_use),nanmean(squeeze(nanmean(cris_moo_obs,2)),1),nanstd(squeeze(nanstd(cris_moo_obs,1,2)),1),'k',0.25);
hold off
 axis([640 1640 -0.5 +0.5])
hl = legend('Strow 2002-2017','CRIS 2012-2019','location','best');

figure(3); axis([640 840 -0.75 +0.75]); grid on
figure(4); axis([640 840 -0.75 +0.75]); grid on
