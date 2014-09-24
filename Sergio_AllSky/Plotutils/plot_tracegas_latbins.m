% Plot trace gas rates with in-situ and ECMWF comparison
%---------------------------------------------------------------------------
% CO2 Compared with global view
load /asl/s1/tangborn/trace_gas_insitu/co2_gv_rate_09_2002_08_2010;
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
errorbar(lat,co2,co2_sigs,'.')
legend('Global View','AIRS');
title('CO2 Rate')
xlabel('Latitude')
ylabel('CO2 rate (ppm/year)')
%---------------------------------------------------------------------------
% Compare SST with ECMWF

% The next file has Sergio's older rates derived from ERA.  
% I think this *might* have some problems for us since he did 3-sigma
% filtering before the robust fit.  May not matter, or maybe he didn't 
% apply this filtering here.
% Variables: stemprate, stempratestd, waterrate, ozonerate, ptemprate, etc.
model_old = load('/asl/s1/rates/Clear/Profile_rates/overocean_gsx_1day_clr_era_lays_spanday01_profilerates_Nov02_2012_robust_span_09_2002_08_2012');

% This file has LLS fits to ERA using gdays set of days used for radiance rates
% Variables: gdays, sst (36x10), sst_rate (sst(:,2)), sst_rate_err
model_sst = load('/asl/s1/rates/Clear/lls_fit_robust_sst_out.mat');

figure;
% Old ERA SST rates: errorbar(lat,stemprate,stempratestd,'r.')
% hold on
errorbar(lat,sst,sst_sigs,'b-+')
hold on;
errorbar(lat,model_sst.sst_rate,model_sst.sst_rate_err,'r-+')
legend('AIRS','ECMWF')
title('SST Rates')
xlabel('Latitude')
ylabel('SST Rate (K/year)')
%---------------------------------------------------------------------------
% CFC11 compared to AGAGE and NOAA HATS rates.
load /asl/s1/tangborn/trace_gas_insitu/agage_rate_09_2002_08_2012;
load /asl/s1/tangborn/trace_gas_insitu/noaahats_cfc11_rate_09_2002_08_2012;

figure;
errorbar(hats_latitude,cfc11_hats_rate,cfc11_hats_error,'r.')
hold on
errorbar(agage_latitude,cfc11_agage_rate,cfc11_agage_error,'m.')
errorbar(lat,cfc11,cfc11_sigs,'.')
legend('NOAA HATS','AGAGE','AIRS')
title('CFC11 Rates')
xlabel('Latitude')
ylabel('CFC11 Rate (ppt/year)')
%---------------------------------------------------------------------------
%  Compare N2O rates with AGAGE 
load /asl/s1/tangborn/trace_gas_insitu/agage_rate_09_2002_08_2012;

figure;
errorbar(agage_latitude,n2o_agage_rate,n2o_agage_error,'r.')
hold on
errorbar(lat,n2o,n2o_sigs,'.')
legend('AGAGE','AIRS');
title('N2O Rates')
xlabel('Latitude')
ylabel('N2O Rate (ppb/year)')
%---------------------------------------------------------------------------
% Compare ozone with global view
load /asl/s1/tangborn/trace_gas_insitu/ozone_gv_rateNEW

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
errorbar(lat,o3,o3_sigs,'.')
legend('Global View','AIRS')
title('Ozone Rates')
xlabel('Latitude')
ylabel('Ozone Rate (frac/year)')
%---------------------------------------------------------------------------