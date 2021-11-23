clear all

load_fairs

iNewOrOld = +1;  %% steves new stats
iNewOrOld = -1;  %% steve/strow old stats

i1231 = find(fairs >= 1231,1)
i790 = find(fairs >= 790.5,1)
i792 = find(fairs >= 791.5,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear bonk fin
yyStart = 2002;
for ii = 10 : 30
  fin = ['ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_' num2str(yyStart) '/sergio_latbin' num2str(ii) '.mat'];
  bonk = load(fin);
  fprintf(1,' 2002 Steve stats files, anomaly done by Sergio ii=%2i fin=%s \n',ii,fin);
  %plot(smooth(1000*(bonk.all_bt_anom(:,i790)-bonk.all_bt_anom(:,i792)),180))
  moo_obs(ii-10+1,:,:)  = bonk.all_bt_anom;
  moo_time(ii-10+1,:,:) = bonk.all_times;
end
avg_moo_time = nanmean(squeeze(moo_time(:,:,i1231)),1);  

anom1231 = squeeze(nanmean(moo_obs(:,:,i1231),1));
plot(2002.75 + (avg_moo_time-avg_moo_time(1))/365,(smooth(anom1231,180)))

anom790 = squeeze(nanmean(moo_obs(:,:,i790),1));
anom792 = squeeze(nanmean(moo_obs(:,:,i792),1));
plot(2002.75 + (avg_moo_time-avg_moo_time(1))/365,(smooth(anom790-anom792,180)))

%%%%%%%%%%%%%%%%%%%%%%%%%

clear bonk fin
yyStart = 2012;
for ii = 10 : 30
  fin = ['ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_' num2str(yyStart) '/sergio_latbin' num2str(ii) '.mat'];  
  bonk = load(fin);
  fprintf(1,' 2012 Steve stats files, anomaly done by Sergio ii=%2i fin=%s \n',ii,fin);
  %plot(smooth(1000*(bonk.all_bt_anom(:,i790)-bonk.all_bt_anom(:,i792)),180))
  xmoo_obs(ii-10+1,:,:)  = bonk.all_bt_anom;
  xmoo_time(ii-10+1,:,:) = bonk.all_times;
end
xavg_moo_time = nanmean(squeeze(xmoo_time(:,:,i1231)),1);  

xanom1231 = squeeze(nanmean(xmoo_obs(:,:,i1231),1));
plot(2012.42 + (xavg_moo_time-xavg_moo_time(1))/365,(smooth(xanom1231,180)))

xanom790 = squeeze(nanmean(xmoo_obs(:,:,i790),1));
xanom792 = squeeze(nanmean(xmoo_obs(:,:,i792),1));
plot(2012.42 + (xavg_moo_time-xavg_moo_time(1))/365,(smooth(xanom790-xanom792,180)))

plot(2002.75+(avg_moo_time-avg_moo_time(1))/365,(smooth(anom790-anom792,180)),'b',2012.42+(xavg_moo_time-xavg_moo_time(1))/365,0.5+(smooth(xanom790-xanom792,180)),'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear bonk fin strow_moo*
yyStart = 2002;
if iNewOrOld > 0
  for ii = 10 : 30
    fin = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/fit_robs_lat' num2str(ii) '.mat'];
    bonk = load(fin);
    fprintf(1,' 2002 Strow data, anomaly done by Strow ii=%2i fin=%s \n',ii,fin);
    %plot(smooth(1000*(bonk.all_bt_anom(:,i790)-bonk.all_bt_anom(:,i792)),180))
    strow_moo_obs(ii-10+1,:,:)  = bonk.all_bt_anom(1:5785,:);
    strow_moo_time(ii-10+1,:,:) = bonk.all_times(1:5785,:);
  end
else
  for ii = 10 : 30
    fin = ['ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_' num2str(yyStart) '/strow_latbin' num2str(ii) '.mat'];
    bonk = load(fin);
    fprintf(1,' 2002 Strow data, anomaly done by Sergio ii=%2i fin=%s \n',ii,fin);
    %plot(smooth(1000*(bonk.all_bt_anom(:,i790)-bonk.all_bt_anom(:,i792)),180))
    strow_moo_obs(ii-10+1,:,:)  = bonk.all_bt_anom(1:5424,:);
    strow_moo_time(ii-10+1,:,:) = bonk.all_times(1:5424,:);
  end
end
strow_avg_moo_time = nanmean(squeeze(strow_moo_time(:,:,i1231)),1);  

strow_anom1231 = squeeze(nanmean(strow_moo_obs(:,:,i1231),1));
plot(2002.75 + (strow_avg_moo_time-strow_avg_moo_time(1))/365,(smooth(strow_anom1231,180)))

figure(1)
strow_anom790 = squeeze(nanmean(strow_moo_obs(:,:,i790),1));
strow_anom792 = squeeze(nanmean(strow_moo_obs(:,:,i792),1));
plot(2002.75 + (avg_moo_time-avg_moo_time(1))/365,(smooth(anom1231,180)),'b',2002.75 + (strow_avg_moo_time-strow_avg_moo_time(1))/365,(smooth(strow_anom1231,180)),'r')
hold on; plot(2012.42+(xavg_moo_time-xavg_moo_time(1))/365,0.0+(smooth(xanom1231,180)),'kx-'); hold off
hl = legend('sergio 2002','strow 2002','sergio 2012','location','best'); set(hl,'fontsize',10)
title('BT1231 stemp anomaly')

figure(2)
plot(2002.75+(avg_moo_time-avg_moo_time(1))/365,(smooth(anom790-anom792,180)),'b',2002.75+(strow_avg_moo_time-strow_avg_moo_time(1))/365,(smooth(strow_anom790-strow_anom792,180)),'r')
hold on; plot(2012.42+(xavg_moo_time-xavg_moo_time(1))/365,0.5+(smooth(xanom790-xanom792,180)),'kx-'); hold off
hl = legend('sergio 2002','strow 2002','sergio 2012','location','best'); set(hl,'fontsize',10)
title('BT790-BT792 CO2 anomaly')
