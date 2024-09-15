%-rw-rw-r-- 1 sergio pi_strow     554662 Jan  8 08:39 ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc_surf.mat
%-rw-rw-r-- 1 sergio pi_strow     555137 Jan  8 08:38 ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc_surf.mat
%-rw-rw-r-- 1 sergio pi_strow   18156117 Jan  7 13:21 ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc.mat
%-rw-rw-r-- 1 sergio pi_strow   18156874 Jan  7 13:10 ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc.mat
%-rw-r--r-- 1 sergio pi_strow   17806050 Oct 26  2022 ERA5_atm_data_2002_09_to_2022_08_trends_desc.mat

%-rw-rw-r-- 1 sergio pi_strow   18012919 Sep 18  2023 ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc_DJF.mat
%-rw-rw-r-- 1 sergio pi_strow   18065465 Sep 18  2023 ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc_MAM.mat
%-rw-rw-r-- 1 sergio pi_strow   17947419 Sep 18  2023 ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc_JJA.mat
%-rw-rw-r-- 1 sergio pi_strow   18069231 Sep 18  2023 ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc_SON.mat
%-rw-rw-r-- 1 sergio pi_strow 2106234673 Jan  8 08:33 ERA5_atm_N_cld_data_2002_09_to_2022_08_asc.mat
%-rw-rw-r-- 1 sergio pi_strow 2106331125 Jan  8 08:33 ERA5_atm_N_cld_data_2002_09_to_2022_08_desc.mat
%-rw-r--r-- 1 sergio pi_strow 2034640085 Oct 26  2022 ERA5_atm_data_2002_09_to_2022_08_desc.mat

disp('paper.tex : this is what Larrabee wants')
disp('paper.tex : this is what Larrabee wants')
disp('paper.tex : this is what Larrabee wants')
disp('paper.tex : \subsection{Test 2 : all ERA5 profiles within a tile : mean of the profiles versus hottest 10th quantile}')
disp('building up all 1/4 deg points over 3x5 grid ~ 240 to 273 points, at grid center; can do the mean over all points or the hottest 90th quantile')
%% /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/driver_load_in_summary_10percent_from_ERA5clearcalc.m --> /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/find_trends_summary_10percent_from_ERA5clearcalc.m
desc_era5_273files = load('/home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/trends_summary_10percent_from_ERA5clearcalc.mat');
oo = find(p.landfrac == 0 & abs(p.landfrac) < 60);
x1 = desc_era5_273files.skt_trend_mean(oo) - desc_era5_273files.skt_trend_Q90(oo);
x2 = desc_era5_273files.skt_trend_mean(oo) - desc_era5_273files.skt_trend_cntr(oo);
printarray([mean(x1) std(x1) mean(x2) std(x2)],'mean/std for loop over ERA 273 pts : meanSKT-Q90SKT and  meanSKT-cntrSKT')
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_era5_273files.skt_trend_mean,'273pts, mean stemp');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_era5_273files.skt_trend_Q90, '273pts, Q90  stemp');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_era5_273files.skt_trend_cntr,'273pts, center stemp');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_era5_273files.skt_trend_cntr - desc_era5_273files.skt_trend_Q90,'273pts, (center - Q90) stemp');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_era5_273files.skt_trend_cntr - desc_era5_273files.skt_trend_mean,'273pts, (center - mean) stemp');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_era5_273files.skt_trend_Q90 - desc_era5_273files.skt_trend_mean,'273pts, (Q90 - mean) stemp');

figure(1); clf; scatter_coast(p.rlon,p.rlat,50,desc_era5_273files.skt_trend_mean);             title('M : DESC : Mean over 252 points per tile'); caxis([-1 +1]*0.150); colormap(cmap)
figure(2); clf; scatter_coast(p.rlon,p.rlat,50,desc_era5_273files.skt_trend_Q90);              title('Q : DESC : Q90 hottest 10 percent');        caxis([-1 +1]*0.150); colormap(cmap)
figure(3); clf; scatter_coast(p.rlon,p.rlat,50,desc_era5_273files.skt_trend_cntr);             title('U : DESC : Center point (our usual)');      caxis([-1 +1]*0.150); colormap(cmap)
figure(4); clf; scatter_coast(p.rlon,p.rlat,50,desc_era5_273files.skt_trend_cntr - desc_era5_273files.skt_trend_Q90);   title('U-Q');  caxis([-1 +1]*0.025); colormap(cmap)
figure(5); clf; scatter_coast(p.rlon,p.rlat,50,desc_era5_273files.skt_trend_cntr - desc_era5_273files.skt_trend_mean);  title('U-M');  caxis([-1 +1]*0.025); colormap(cmap)
figure(6); clf; scatter_coast(p.rlon,p.rlat,50,desc_era5_273files.skt_trend_Q90 - desc_era5_273files.skt_trend_mean);   title('Q-M');  caxis([-1 +1]*0.025); colormap(cmap)

figure(1); clf; aslmap(1,rlat65,rlon73,smoothn(reshape(desc_era5_273files.skt_trend_mean,72,64)',1), [-90 +90],[-180 +180]);             title('M : DESC : Mean over 252 points per tile'); caxis([-1 +1]*0.150); colormap(cmap)
figure(2); clf; aslmap(2,rlat65,rlon73,smoothn(reshape(desc_era5_273files.skt_trend_Q90,72,64)',1), [-90 +90],[-180 +180]);              title('Q : DESC : Q90 hottest 10 percent');        caxis([-1 +1]*0.150); colormap(cmap)
figure(3); clf; aslmap(3,rlat65,rlon73,smoothn(reshape(desc_era5_273files.skt_trend_cntr,72,64)',1), [-90 +90],[-180 +180]);             title('U : DESC : Center point (our usual)');      caxis([-1 +1]*0.150); colormap(cmap)
figure(4); clf; aslmap(4,rlat65,rlon73,smoothn(reshape(desc_era5_273files.skt_trend_cntr - desc_era5_273files.skt_trend_Q90,72,64)',1), [-90 +90],[-180 +180]);   title('U-Q');  caxis([-1 +1]*0.025); colormap(cmap)
figure(5); clf; aslmap(5,rlat65,rlon73,smoothn(reshape(desc_era5_273files.skt_trend_cntr - desc_era5_273files.skt_trend_mean,72,64)',1), [-90 +90],[-180 +180]);  title('U-M');  caxis([-1 +1]*0.025); colormap(cmap)
figure(6); clf; aslmap(6,rlat65,rlon73,smoothn(reshape(desc_era5_273files.skt_trend_Q90 - desc_era5_273files.skt_trend_mean,72,64)',1), [-90 +90],[-180 +180]);   title('Q-M');  caxis([-1 +1]*0.025); colormap(cmap)

dbt = (-5 : 0.1:+5)*0.01/2; d1 = mean(diff(dbt));
ooA = find(p.landfrac == 0); ooB = find(p.landfrac == 0 & abs(p.rlat) < 60); 
figure(7); clf; 
  histA_UQ = histc(desc_era5_273files.skt_trend_cntr(ooA)-desc_era5_273files.skt_trend_Q90(ooA),dbt);
  histA_UM = histc(desc_era5_273files.skt_trend_cntr(ooA)-desc_era5_273files.skt_trend_mean(ooA),dbt);
  histA_QM = histc(desc_era5_273files.skt_trend_Q90(ooA)-desc_era5_273files.skt_trend_mean(ooA),dbt);
  histB_UQ = histc(desc_era5_273files.skt_trend_cntr(ooB)-desc_era5_273files.skt_trend_Q90(ooB),dbt);
  histB_UM = histc(desc_era5_273files.skt_trend_cntr(ooB)-desc_era5_273files.skt_trend_mean(ooB),dbt);
  histB_QM = histc(desc_era5_273files.skt_trend_Q90(ooB)-desc_era5_273files.skt_trend_mean(ooB),dbt);
  plot(dbt,histA_UQ,'b',dbt,histA_UM,'r',dbt,histA_QM,'g',dbt,histB_UQ,'c',dbt,histB_UM,'m',dbt,histB_QM,'k')
  hl = legend('all UQ','all UM','all QM','tropical UQ','tropical UM','tropical QM','location','best','fontsize',10); title('Ocean'); set(gca,'fontsize',10); grid
hold on
dbt = (-5 : 0.1:+5)*0.1; d2 = mean(diff(dbt));
  histA_U = histc(desc_era5_273files.skt_trend_cntr(ooA),dbt);
  histA_M = histc(desc_era5_273files.skt_trend_mean(ooA),dbt);
  histA_Q = histc(desc_era5_273files.skt_trend_Q90(ooA),dbt);
  plot(dbt,histA_U*d1/d2*6,'b',dbt,histA_M*d1/d2*6,'r',dbt,histA_Q*d1/d2*6,'g')
hold off
xlim([-1 +1]*0.2)

dbt = (-5 : 0.1:+5)*0.01/2; d1 = mean(diff(dbt));
ooA = find(p.landfrac == 1); ooB = find(p.landfrac == 0 & abs(p.rlat) < 60); 
figure(8); clf; 
  histA_UQ = histc(desc_era5_273files.skt_trend_cntr(ooA)-desc_era5_273files.skt_trend_Q90(ooA),dbt);
  histA_UM = histc(desc_era5_273files.skt_trend_cntr(ooA)-desc_era5_273files.skt_trend_mean(ooA),dbt);
  histA_QM = histc(desc_era5_273files.skt_trend_Q90(ooA)-desc_era5_273files.skt_trend_mean(ooA),dbt);
  histB_UQ = histc(desc_era5_273files.skt_trend_cntr(ooB)-desc_era5_273files.skt_trend_Q90(ooB),dbt);
  histB_UM = histc(desc_era5_273files.skt_trend_cntr(ooB)-desc_era5_273files.skt_trend_mean(ooB),dbt);
  histB_QM = histc(desc_era5_273files.skt_trend_Q90(ooB)-desc_era5_273files.skt_trend_mean(ooB),dbt);
  plot(dbt,histA_UQ,'b',dbt,histA_UM,'r',dbt,histA_QM,'g',dbt,histB_UQ,'c',dbt,histB_UM,'m',dbt,histB_QM,'k')
  hl = legend('all UQ','all UM','all QM','tropical UQ','tropical UM','tropical QM','location','best','fontsize',10); title('Land'); set(gca,'fontsize',10); grid
hold on
dbt = (-5 : 0.1:+5)*0.1; d2 = mean(diff(dbt));
  histA_U = histc(desc_era5_273files.skt_trend_cntr(ooA),dbt);
  histA_M = histc(desc_era5_273files.skt_trend_mean(ooA),dbt);
  histA_Q = histc(desc_era5_273files.skt_trend_Q90(ooA),dbt);
  plot(dbt,histA_U*d1/d2*15,'b',dbt,histA_M*d1/d2*15,'r',dbt,histA_Q*d1/d2*15,'g')
hold off
xlim([-1 +1]*0.2)

%%%%%

figure(9); clf
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_era5_273files.skt_trend_cntr - desc_era5_273files.skt_trend_Q90,'273pts, (center - Q90) stemp');
dbt_delta = (-5 : 0.1:+5)*0.01/2; d1 = mean(diff(dbt_delta));
ooA = find(p.landfrac == 0); ooB = find(p.landfrac == 0 & abs(p.rlat) < 60); 
histA_UQ_delta = histc(desc_era5_273files.skt_trend_cntr(ooA)-desc_era5_273files.skt_trend_Q90(ooA),dbt_delta);
  fit_UQ_delta = fit(dbt_delta',histA_UQ_delta','gauss2'); figure(9); clf; plot(fit_UQ_delta,dbt_delta,histA_UQ_delta)
  fit_UQ_delta = fit(dbt_delta',histA_UQ_delta','gauss1'); figure(9); clf; plot(fit_UQ_delta,dbt_delta,histA_UQ_delta)

stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_era5_273files.skt_trend_cntr,'273pts, (center) stemp');
dbt = (-5 : 0.1:+5)*0.1; d2 = mean(diff(dbt));
ooA = find(p.landfrac == 0); ooB = find(p.landfrac == 0 & abs(p.rlat) < 60); 
histA_UQ = histc(desc_era5_273files.skt_trend_cntr(ooA),dbt);
  fit_UQ = fit(dbt',histA_UQ','gauss2'); figure(9); clf; plot(fit_UQ,dbt,histA_UQ)
  fit_UQ = fit(dbt',histA_UQ','gauss1'); figure(9); clf; plot(fit_UQ,dbt,histA_UQ)

figure(9); clf; plot(dbt_delta,histA_UQ_delta/max(histA_UQ_delta),'b',dbt,histA_UQ/max(histA_UQ),'r')
xlim([-1 +1]*0.2); grid on; title('Normalized Ocean Histograms'); hl = legend('U-Q','U');

%%%%%

figure(10); clf
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_era5_273files.skt_trend_cntr - desc_era5_273files.skt_trend_Q90,'273pts, (center - Q90) stemp');
dbt_delta = (-5 : 0.1:+5)*0.01/2; d1 = mean(diff(dbt_delta));
ooA = find(p.landfrac == 1); ooB = find(p.landfrac == 1 & abs(p.rlat) < 60); 
histA_UQ_delta = histc(desc_era5_273files.skt_trend_cntr(ooA)-desc_era5_273files.skt_trend_Q90(ooA),dbt_delta);
  fit_UQ_delta = fit(dbt_delta',histA_UQ_delta','gauss2'); figure(10); clf; plot(fit_UQ_delta,dbt_delta,histA_UQ_delta)
  fit_UQ_delta = fit(dbt_delta',histA_UQ_delta','gauss1'); figure(10); clf; plot(fit_UQ_delta,dbt_delta,histA_UQ_delta)

stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_era5_273files.skt_trend_cntr,'273pts, (center) stemp');
dbt = (-5 : 0.1:+5)*0.1; d2 = mean(diff(dbt));
ooA = find(p.landfrac == 1); ooB = find(p.landfrac == 1 & abs(p.rlat) < 60); 
histA_UQ = histc(desc_era5_273files.skt_trend_cntr(ooA),dbt);
  fit_UQ = fit(dbt',histA_UQ','gauss2'); figure(10); clf; plot(fit_UQ,dbt,histA_UQ)
  fit_UQ = fit(dbt',histA_UQ','gauss1'); figure(10); clf; plot(fit_UQ,dbt,histA_UQ)

figure(10); clf; plot(dbt_delta,histA_UQ_delta/max(histA_UQ_delta),'b',dbt,histA_UQ/max(histA_UQ),'r')
xlim([-1 +1]*0.2); grid on; title('Normalized Land Histograms'); hl = legend('U-Q','U');

disp('reminder U = DESC usual tile center, UD = U');
disp('         Q = DESC quantile0.90 over whole tile');
disp('         M = DESC mean over whole tile');

disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                              >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
