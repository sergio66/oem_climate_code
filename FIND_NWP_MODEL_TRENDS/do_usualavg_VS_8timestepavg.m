%[sergio@taki-usr2 FIND_NWP_MODEL_TRENDS]$ ls -lt ERA5_atm_cld_data_2002_09_to_2022_08_*
%-rw-rw-r-- 1 sergio pi_strow   17764983 Jun 20 18:53 ERA5_atm_cld_data_2002_09_to_2022_08_trends_day_night.mat
%-rw-rw-r-- 1 sergio pi_strow     171390 Jun 20 12:13 ERA5_atm_cld_data_2002_09_to_2022_08_trends_day_night_surf.mat
%-rw-rw-r-- 1 sergio pi_strow   17789441 Jun 20 11:56 ERA5_atm_cld_data_2002_09_to_2022_08_trends_day_night_randompt.mat
%-rw-rw-r-- 1 sergio pi_strow     171478 Jun 20 06:07 ERA5_atm_cld_data_2002_09_to_2022_08_trends_day_night_randompt_surf.mat

%-rw-rw-r-- 1 sergio pi_strow 2070433308 Jun 20 12:10 ERA5_atm_cld_data_2002_09_to_2022_08_day_night.mat
%-rw-rw-r-- 1 sergio pi_strow 2075354799 Jun 20 05:54 ERA5_atm_cld_data_2002_09_to_2022_08_day_night_randompt.mat

disp('taking average over 8 timesteps : grid center vs random point : day+night')
disp('paper.tex : \subsection{Test 1 : Profile from random location within tile versus tile center}')
dn_grid      = load('ERA5_atm_cld_data_2002_09_to_2022_08_trends_day_night.mat');
dn_surf_grid = load('ERA5_atm_cld_data_2002_09_to_2022_08_trends_day_night_surf.mat');

dn_rand      = load('ERA5_atm_cld_data_2002_09_to_2022_08_trends_day_night_randompt.mat');
dn_surf_rand = load('ERA5_atm_cld_data_2002_09_to_2022_08_trends_day_night_randompt_surf.mat');

%%%%%%%%%%
%% this is day only
asc_grid      = load('ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc.mat');
asc_surf_grid = load('ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc_surf.mat');

%% this is night only
desc_grid      = load('ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc.mat');
desc_surf_grid = load('ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc_surf.mat');

%% this is the average
avg_asc_desc = 0.5*(asc_surf_grid.trend_stemp + desc_surf_grid.trend_stemp);
%%%%%%%%%%

stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,avg_asc_desc,'day U + night U, Usual regular grid');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,dn_surf_grid.trend_stemp,'day+night 8 timesetps, R8egular grid = Q');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,dn_surf_rand.trend_stemp,'day+night 8 timesteps, randoM8 grid = M');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,avg_asc_desc - dn_surf_grid.trend_stemp,'U - R8 grid');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,avg_asc_desc - dn_surf_rand.trend_stemp,'U - M8 grid');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,dn_surf_grid.trend_stemp - dn_surf_rand.trend_stemp,'R8 - M8 grid');

figure(1); clf; aslmap(1,rlat65,rlon73,smoothn(reshape(dn_surf_grid.trend_stemp,72,64)',1), [-90 +90],[-180 +180]);                             title('M : Random tile point, 8 timestep avg');         caxis([-1 +1]*0.150); colormap(cmap)
figure(2); clf; aslmap(2,rlat65,rlon73,smoothn(reshape(dn_surf_rand.trend_stemp,72,64)',1), [-90 +90],[-180 +180]);                             title('Q : Center point (our usual), 8 timestep avg');  caxis([-1 +1]*0.150); colormap(cmap)
figure(3); clf; aslmap(3,rlat65,rlon73,smoothn(reshape(avg_asc_desc,72,64)',1), [-90 +90],[-180 +180]);                                         title('U : DN : Center point (our usual) (UD+UN)/2');   caxis([-1 +1]*0.150); colormap(cmap)
figure(4); clf; aslmap(4,rlat65,rlon73,smoothn(reshape(avg_asc_desc - dn_surf_rand.trend_stemp,72,64)',1), [-90 +90],[-180 +180]);              title('U-Q');  caxis([-1 +1]*0.025); colormap(cmap)
figure(5); clf; aslmap(5,rlat65,rlon73,smoothn(reshape(avg_asc_desc - dn_surf_grid.trend_stemp,72,64)',1), [-90 +90],[-180 +180]);              title('U-M');  caxis([-1 +1]*0.025); colormap(cmap)
figure(6); clf; aslmap(6,rlat65,rlon73,smoothn(reshape(dn_surf_rand.trend_stemp - dn_surf_grid.trend_stemp,72,64)',1), [-90 +90],[-180 +180]);  title('Q-M');  caxis([-1 +1]*0.025); colormap(cmap)

dbt = (-5 : 0.1:+5)*0.01/2; d1 = mean(diff(dbt));
ooA = find(p.landfrac == 0); ooB = find(p.landfrac == 0 & abs(p.rlat) < 60); 
figure(7); clf; 
  histA_UQ = histc(avg_asc_desc(ooA)-dn_surf_rand.trend_stemp(ooA),dbt);
  histA_UM = histc(avg_asc_desc(ooA)-dn_surf_grid.trend_stemp(ooA),dbt);
  histA_QM = histc(dn_surf_rand.trend_stemp(ooA)-dn_surf_grid.trend_stemp(ooA),dbt);
  histB_UQ = histc(avg_asc_desc(ooB)-dn_surf_rand.trend_stemp(ooB),dbt);
  histB_UM = histc(avg_asc_desc(ooB)-dn_surf_grid.trend_stemp(ooB),dbt);
  histB_QM = histc(dn_surf_rand.trend_stemp(ooB)-dn_surf_grid.trend_stemp(ooB),dbt);
  plot(dbt,histA_UQ,'b',dbt,histA_UM,'r',dbt,histA_QM,'g',dbt,histB_UQ,'c',dbt,histB_UM,'m',dbt,histB_QM,'k')
  hl = legend('all UQ','all UM','all QM','tropical UQ','tropical UM','tropical QM','location','best','fontsize',10); title('Ocean'); set(gca,'fontsize',10); grid
hold on
dbt = (-5 : 0.1:+5)*0.1; d2 = mean(diff(dbt));
  histA_U = histc(avg_asc_desc(ooA),dbt);
  histA_M = histc(dn_surf_grid.trend_stemp(ooA),dbt);
  histA_Q = histc(dn_surf_rand.trend_stemp(ooA),dbt);
  plot(dbt,histA_U*d1/d2*6,'b',dbt,histA_M*d1/d2*6,'r',dbt,histA_Q*d1/d2*6,'g')
hold off
xlim([-1 +1]*0.2)

dbt = (-5 : 0.1:+5)*0.01/2; d1 = mean(diff(dbt));
ooA = find(p.landfrac == 1); ooB = find(p.landfrac == 0 & abs(p.rlat) < 60); 
figure(8); clf; 
  histA_UQ = histc(avg_asc_desc(ooA)-dn_surf_rand.trend_stemp(ooA),dbt);
  histA_UM = histc(avg_asc_desc(ooA)-dn_surf_grid.trend_stemp(ooA),dbt);
  histA_QM = histc(dn_surf_rand.trend_stemp(ooA)-dn_surf_grid.trend_stemp(ooA),dbt);
  histB_UQ = histc(avg_asc_desc(ooB)-dn_surf_rand.trend_stemp(ooB),dbt);
  histB_UM = histc(avg_asc_desc(ooB)-dn_surf_grid.trend_stemp(ooB),dbt);
  histB_QM = histc(dn_surf_rand.trend_stemp(ooB)-dn_surf_grid.trend_stemp(ooB),dbt);
  plot(dbt,histA_UQ,'b',dbt,histA_UM,'r',dbt,histA_QM,'g',dbt,histB_UQ,'c',dbt,histB_UM,'m',dbt,histB_QM,'k')
  hl = legend('all UQ','all UM','all QM','tropical UQ','tropical UM','tropical QM','location','best','fontsize',10); title('Land'); set(gca,'fontsize',10); grid
hold on
dbt = (-5 : 0.1:+5)*0.1; d2 = mean(diff(dbt));
  histA_U = histc(avg_asc_desc(ooA),dbt);
  histA_M = histc(dn_surf_grid.trend_stemp(ooA),dbt);
  histA_Q = histc(dn_surf_rand.trend_stemp(ooA),dbt);
  plot(dbt,histA_U*d1/d2*15,'b',dbt,histA_M*d1/d2*15,'r',dbt,histA_Q*d1/d2*15,'g')
hold off
xlim([-1 +1]*0.2)

%%%%%

figure(9); clf
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,avg_asc_desc - dn_surf_rand.trend_stemp,'273pts, (center - Q90) stemp');
dbt_delta = (-5 : 0.1:+5)*0.01/2; d1 = mean(diff(dbt_delta));
ooA = find(p.landfrac == 0); ooB = find(p.landfrac == 0 & abs(p.rlat) < 60); 
histA_UQ_delta = histc(avg_asc_desc(ooA)-dn_surf_rand.trend_stemp(ooA),dbt_delta);
  fit_UQ_delta = fit(dbt_delta',histA_UQ_delta','gauss2'); figure(9); clf; plot(fit_UQ_delta,dbt_delta,histA_UQ_delta)
  fit_UQ_delta = fit(dbt_delta',histA_UQ_delta','gauss1'); figure(9); clf; plot(fit_UQ_delta,dbt_delta,histA_UQ_delta)

stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,avg_asc_desc,'273pts, (center) stemp');
dbt = (-5 : 0.1:+5)*0.1; d2 = mean(diff(dbt));
ooA = find(p.landfrac == 0); ooB = find(p.landfrac == 0 & abs(p.rlat) < 60); 
histA_UQ = histc(avg_asc_desc(ooA),dbt);
  fit_UQ = fit(dbt',histA_UQ','gauss2'); figure(9); clf; plot(fit_UQ,dbt,histA_UQ)
  fit_UQ = fit(dbt',histA_UQ','gauss1'); figure(9); clf; plot(fit_UQ,dbt,histA_UQ)

figure(9); clf; plot(dbt_delta,histA_UQ_delta/max(histA_UQ_delta),'b',dbt,histA_UQ/max(histA_UQ),'r')
xlim([-1 +1]*0.2); grid on; title('Normalized Ocean Histograms'); hl = legend('U-Q','U');

%%%%%

figure(10); clf
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,avg_asc_desc - dn_surf_rand.trend_stemp,'273pts, (center - Q90) stemp');
dbt_delta = (-5 : 0.1:+5)*0.01/2; d1 = mean(diff(dbt_delta));
ooA = find(p.landfrac == 1); ooB = find(p.landfrac == 1 & abs(p.rlat) < 60); 
histA_UQ_delta = histc(avg_asc_desc(ooA)-dn_surf_rand.trend_stemp(ooA),dbt_delta);
  fit_UQ_delta = fit(dbt_delta',histA_UQ_delta','gauss2'); figure(10); clf; plot(fit_UQ_delta,dbt_delta,histA_UQ_delta)
  fit_UQ_delta = fit(dbt_delta',histA_UQ_delta','gauss1'); figure(10); clf; plot(fit_UQ_delta,dbt_delta,histA_UQ_delta)

stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,avg_asc_desc,'273pts, (center) stemp');
dbt = (-5 : 0.1:+5)*0.1; d2 = mean(diff(dbt));
ooA = find(p.landfrac == 1); ooB = find(p.landfrac == 1 & abs(p.rlat) < 60); 
histA_UQ = histc(avg_asc_desc(ooA),dbt);
  fit_UQ = fit(dbt',histA_UQ','gauss2'); figure(10); clf; plot(fit_UQ,dbt,histA_UQ)
  fit_UQ = fit(dbt',histA_UQ','gauss1'); figure(10); clf; plot(fit_UQ,dbt,histA_UQ)

figure(10); clf; plot(dbt_delta,histA_UQ_delta/max(histA_UQ_delta),'b',dbt,histA_UQ/max(histA_UQ),'r')
xlim([-1 +1]*0.2); grid on; title('Normalized Land Histograms'); hl = legend('U-Q','U');

disp('reminder U = our usual tile center, day/night mean (UD+UN)/2 = U');
disp('         Q = our usual tile center, average over 8 times');
disp('         M = random tile center,    average over 8 times');

disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                              >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
