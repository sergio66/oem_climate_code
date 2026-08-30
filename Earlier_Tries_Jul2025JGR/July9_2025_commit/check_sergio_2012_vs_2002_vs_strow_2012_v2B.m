clear all

for ii = 10 : 30
  fin = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/fit_robs_lat' num2str(ii) '.mat'];
  bonk = load(fin);
  fprintf(1,' 2002 Strow data, anomaly done by Strow ii=%2i fin=%s \n',ii,fin);
  %plot(smooth(1000*(bonk.all_bt_anom(:,i790)-bonk.all_bt_anom(:,i792)),180))
  strow_moo_obs(ii-10+1,:,:)  = bonk.all_bt_anom(1:5785,:);
  strow_moo_time(ii-10+1,:,:) = bonk.all_times(1:5785,:);
end

yyStart = 2012;
for ii = 10 : 30
  fin = ['ANOM_16dayavg_nonucal_2012_2018/Daily_Anomalies_' num2str(yyStart) '/sergio_latbin' num2str(ii) '.mat'];
  bonk = load(fin);
  fprintf(1,' 2002 Strow data, anomaly done by Sergio ii=%2i fin=%s \n',ii,fin);
  %plot(smooth(1000*(bonk.all_bt_anom(:,i790)-bonk.all_bt_anom(:,i792)),180))
  sergio_moo_obs(ii-10+1,:,:)  = bonk.all_bt_anom(1:1813,:);
  sergio_moo_time(ii-10+1,:,:) = bonk.all_times(1:1813,:);
end

%load_fairs
load f2645.mat; fairs = f2645;

i1231 = find(fairs >= 1231,1);
i790 = find(fairs >= 790.5,1);
i792 = find(fairs >= 791.5,1);

xanom1231 = nanmean(sergio_moo_obs(:,:,i1231),1);
xanom790 = nanmean(sergio_moo_obs(:,:,i790),1);
xanom792 = nanmean(sergio_moo_obs(:,:,i792),1);

anom1231 = nanmean(strow_moo_obs(:,:,i1231),1);
anom790 = nanmean(strow_moo_obs(:,:,i790),1);
anom792 = nanmean(strow_moo_obs(:,:,i792),1);

figure(1)
plot(2012.31+(1:1813)/365,smooth(xanom1231,180),'b',2002.75+(1:5785)/365,smooth(anom1231,180),'r');
  title('BT 1231 anomaly'); hl = legend('Sergio','Strow','location','best'); set(hl,'fontsize',10);

figure(2)
plot(2012.31+(1:1813)/365,0.49+smooth(xanom790-xanom792,180),'b',2002.75+(1:5785)/365,smooth(anom790-anom792,180),'r');
  title('BT 790-792 anomaly'); hl = legend('Sergio','Strow','location','best'); set(hl,'fontsize',10);
