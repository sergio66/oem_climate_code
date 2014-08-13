% Plot trace gas rates with in-situ and ECMWF comparison


%   CO2 Compared with global view


load /asl/s1/tangborn/trace_gas_insitu/co2_gv_rate_09_2002_08_2010;
load /asl/s1/tangborn/trace_gas_insitu/latbins; 


% Drop extreme latitudes 
for i=1:34
        co2_airs(i)=co2(i+1);
        co2_airs_errors(i)=co2_sigs(i+1);
        co2_airs_cal(i)=co2_cal(i+1)
        co2_airs_errors_cal(i)=co2_cal_sigs(i+1);
        lat_new(i)=save_lat(i+1);
        co2_airs_corrected(i)=co2_airs(i)-co2_airs_cal(i);
        co2_airs_errors_corrected(i)=co2_airs_errors(i)+co2_airs_errors_cal(i);
end

[dum n_gv]=size(co2_gv_rate);
ii=1;
for i=1:n_gv
   if co2_gv_hgt(i) > 6000   % Restrict GV comparison to above 6000 meters
      gv_error(ii)=co2_gv_error(i);
      gv_rate(ii)=co2_gv_rate(i);
      gv_lat(ii)=co2_gv_latitude(i);
      gv_hgt(ii)=co2_gv_hgt(i);
      ii=ii+1;
   end
end

figure
errorbar(gv_lat,gv_rate,gv_error,'r.')
hold on
errorbar(lat_new,co2_airs,co2_airs_errors,'.')
errorbar(lat_new,co2_airs_cal,co2_airs_errors_cal,'.k');
errorbar(lat_new,co2_airs_corrected,co2_airs_errors_corrected,'.g');
legend('Global View','AIRS','CALC','AIRS - CALC')
title('CO2 Rate')
xlabel('Latitude')
ylabel('CO2 rate (ppm/year)')


% CFC11 compared to AGAGE and NOAA HATS rates.



load /asl/s1/tangborn/trace_gas_insitu/agage_rate_09_2002_08_2012;
load /asl/s1/tangborn/trace_gas_insitu/noaahats_cfc11_rate_09_2002_08_2012;
load /asl/s1/tangborn/trace_gas_insitu/overocean_gsx_1day_clr_era_lays_spanday01_fatsummary_Nov02_2012_ALLcov;

% Remove first and lat latbin 

for i=1:34
        cfc11_airs(i)=cfc11(i+1);
        cfc11_airs_errors(i)=cfc11_sigs(i+1);
        cfc11_airs_cal(i)=cfc11_cal(i+1)
        cfc11_airs_errors_cal(i)=cfc11_cal_sigs(i+1);
        cfc11_airs_corrected(i)=cfc11_airs(i)-cfc11_airs_cal(i);
        cfc11_airs_errors_corrected(i)=cfc11_airs_errors(i)+cfc11_airs_errors_cal(i);
end


figure;
errorbar(hats_latitude,cfc11_hats_rate,cfc11_hats_error,'r.')
hold on
errorbar(agage_latitude,cfc11_agage_rate,cfc11_agage_error,'m.')
errorbar(lat_new,cfc11_airs,cfc11_airs_errors,'.')
%errorbar(lat_new,cfc11_airs_cal,cfc11_airs_errors_cal,'.k');
%errorbar(lat_new,cfc11_airs-cfc11_airs_cal,cfc11_airs_errors+cfc11_airs_errors_cal,'.g') 
legend('NOAA HATS','AGAGE','AIRS','CALC','AIRS - CALC'); 
title('CFC11 Rates')
xlabel('Latitude')
ylabel('CFC11 Rate (ppt/year)')

%  Compare N2O rates with AGAGE 


load /asl/s1/tangborn/trace_gas_insitu/agage_rate_09_2002_08_2012;

for i=1:34
        n2o_airs(i)=n2o(i+1);
        n2o_airs_errors(i)=n2o_sigs(i+1);
        n2o_airs_cal(i)=n2o_cal(i+1)
        n2o_airs_errors_cal(i)=n2o_cal_sigs(i+1);
        n2o_airs_corrected(i)=n2o_airs(i)-n2o_airs_cal(i);
        n2o_airs_errors_corrected(i)=n2o_airs_errors(i)+n2o_airs_errors_cal(i);
end

figure;
errorbar(agage_latitude,n2o_agage_rate,n2o_agage_error,'r.')
hold on
errorbar(lat_new,n2o_airs,n2o_airs_errors,'.')
errorbar(lat_new,n2o_airs_cal,n2o_airs_errors_cal,'.k');
errorbar(lat_new,n2o_airs_corrected,n2o_airs_errors_corrected,'.g');
legend('AGAGE','AIRS','Calc','AIRS - Calc');
title('N2O Rates')
xlabel('Latitude')
ylabel('N2O Rate (ppb/year)')



% Compare ozone with global view


load /asl/s1/tangborn/trace_gas_insitu/ozone_gv_rateNEW
load /asl/s1/tangborn/trace_gas_insitu/overocean_gsx_1day_clr_era_lays_spanday01_fatsummary_Nov02_2012_ALLcov;

% Remove first and last latbins 

for i=1:34
        ozone_airs(i)=o3(i+1);
        ozone_airs_errors(i)=o3_sigs(i+1);
end

[junk n_ozone_gv]=size(sitelat);
iozone=1; 
for i=1:n_ozone_gv
  if sitelat(i) < 61 && sitelat(i) > -61 
     lat_ozone(iozone)=sitelat(i); 
     gv_ozone_new(iozone)=ozone_gv_rate(i);  
     gv_error_new(iozone)=ozone_gv_error(i); 
     ozone_constant_new(iozone)=ozone_gv_constant(i); 
     iozone=iozone+1; 
  end 
end 

figure
errorbar(lat_ozone,gv_ozone_new./ozone_constant_new,gv_error_new./ozone_constant_new,'r.')
hold on
errorbar(lat_new,ozone_airs,ozone_airs_errors,'.')
legend('Global View','AIRS')
title('Ozone Rates')
xlabel('Latitude')
ylabel('Ozone Rate (frac/year)')


% Compare SST with ECMWF

load /asl/s1/rates/clear/Aug2013/overocean_gsx_1day_clr_era_lays_spanday01_profilerates_Nov02_2012_robust_span_09_2002_08_2012;
load /asl/s1/tangborn/trace_gas_insitu/overocean_gsx_1day_clr_era_lays_spanday01_fatsummary_Nov02_2012_ALLcov;

% remove first and last latbins 

for i=1:34
        sst_ecmwf(i)=stemprate(i+1); 
        sst_ecmwf_errors(i)=stempratestd(i+1); 
        sst_airs(i)=sst(i+1);
        sst_airs_errors(i)=sst_sigs(i+1);
        sst_airs_cal(i)=sst_cal(i+1)
        sst_airs_errors_cal(i)=sst_cal_sigs(i+1);
        sst_airs_corrected(i)=sst_airs(i)-sst_airs_cal(i);
        sst_airs_errors_corrected(i)=sst_airs_errors(i)+sst_airs_errors_cal(i);
end

% Plot Error bars 


figure;
errorbar(lat_new,sst_ecmwf,sst_ecmwf_errors,'r.')
hold on
errorbar(lat_new,sst_airs,sst_airs_errors,'.')
legend('ECMWF','AIRS')
title('SST Rates')
xlabel('Latitude')
ylabel('SST Rate (K/year)')

