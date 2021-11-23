a16 = load('ANOM/LatBin30/LonBin31/v1_Quantile16/desc_anaomaly_LatBin30_LonBin31_Quantile16_timesetps_001_412_V1.mat');
a08 = load('ANOM/LatBin30/LonBin31/v1_Quantile08/desc_anaomaly_LatBin30_LonBin31_Quantile08_timesetps_001_412_V1.mat');

figure(1); plot(1:2645,a08.Ball(:,2),1:2645,a16.Ball(:,2))
  hl = legend('Q08','Q16','location','best','fontsize',10);
figure(2); plot(1:2645,nanmean(a08.bt_anom_all,2),1:2645,nanmean(a16.bt_anom_all,2))
  hl = legend('Q08','Q16','location','best','fontsize',10);
  
