addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
co2 = load('mauna_loa_co2_growth_rate.txt');
plot(co2(:,1),co2(:,2));
co2(:,3) = 350 + cumsum(co2(:,2));

yyS = 1956; yyE = 2092; ii = find(co2(:,1) >= yyS & co2(:,1) <= yyE);
[junk,err] = Math_tsfit_lin_robust(co2(ii,1)*365,co2(ii,3),0); 
fprintf(1,'mauna loa co2 rate between %4i - %4i = %8.6f +/- %8.6f ppm/yr \n',yyS,yyE,junk(2),err.se(2))

yyS = 2002; yyE = 2014; ii = find(co2(:,1) >= yyS & co2(:,1) <= yyE);
[junk,err] = Math_tsfit_lin_robust(co2(ii,1)*365,co2(ii,3),0); 
fprintf(1,'mauna loa co2 rate between %4i - %4i = %8.6f +/- %8.6f ppm/yr \n',yyS,yyE,junk(2),err.se(2))

yyS = 2002; yyE = 2021; ii = find(co2(:,1) >= yyS & co2(:,1) <= yyE);
[junk,err] = Math_tsfit_lin_robust(co2(ii,1)*365,co2(ii,3),0); 
fprintf(1,'mauna loa co2 rate between %4i - %4i = %8.6f +/- %8.6f ppm/yr \n',yyS,yyE,junk(2),err.se(2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
n2o = load('mauna_loa_n2o_growth_rate.txt');
plot(n2o(:,1),n2o(:,2));
n2o(:,3) = 350 + cumsum(n2o(:,2));

yyS = 1956; yyE = 2092; ii = find(n2o(:,1) >= yyS & n2o(:,1) <= yyE);
[junk,err] = Math_tsfit_lin_robust(n2o(ii,1)*365,n2o(ii,3),0); 
fprintf(1,'mauna loa n2o rate between %4i - %4i = %8.6f +/- %8.6f ppm/yr \n',yyS,yyE,junk(2),err.se(2))

yyS = 2002; yyE = 2014; ii = find(n2o(:,1) >= yyS & n2o(:,1) <= yyE);
[junk,err] = Math_tsfit_lin_robust(n2o(ii,1)*365,n2o(ii,3),0); 
fprintf(1,'mauna loa n2o rate between %4i - %4i = %8.6f +/- %8.6f ppm/yr \n',yyS,yyE,junk(2),err.se(2))

yyS = 2002; yyE = 2021; ii = find(n2o(:,1) >= yyS & n2o(:,1) <= yyE);
[junk,err] = Math_tsfit_lin_robust(n2o(ii,1)*365,n2o(ii,3),0); 
fprintf(1,'mauna loa n2o rate between %4i - %4i = %8.6f +/- %8.6f ppm/yr \n',yyS,yyE,junk(2),err.se(2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
ch4 = load('mauna_loa_ch4_growth_rate.txt');
plot(ch4(:,1),ch4(:,2));
ch4(:,3) = 350 + cumsum(ch4(:,2));

yyS = 1956; yyE = 2092; ii = find(ch4(:,1) >= yyS & ch4(:,1) <= yyE);
[junk,err] = Math_tsfit_lin_robust(ch4(ii,1)*365,ch4(ii,3),0); 
fprintf(1,'mauna loa ch4 rate between %4i - %4i = %8.6f +/- %8.6f ppm/yr \n',yyS,yyE,junk(2),err.se(2))

yyS = 2002; yyE = 2014; ii = find(ch4(:,1) >= yyS & ch4(:,1) <= yyE);
[junk,err] = Math_tsfit_lin_robust(ch4(ii,1)*365,ch4(ii,3),0); 
fprintf(1,'mauna loa ch4 rate between %4i - %4i = %8.6f +/- %8.6f ppm/yr \n',yyS,yyE,junk(2),err.se(2))

yyS = 2002; yyE = 2021; ii = find(ch4(:,1) >= yyS & ch4(:,1) <= yyE);
[junk,err] = Math_tsfit_lin_robust(ch4(ii,1)*365,ch4(ii,3),0); 
fprintf(1,'mauna loa ch4 rate between %4i - %4i = %8.6f +/- %8.6f ppm/yr \n',yyS,yyE,junk(2),err.se(2))
