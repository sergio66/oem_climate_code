%load Test/output_alllats_obs_cal_rate_offset_0p1K
load Test/alldriver_fit_obs.mat
for JOB=1:36
   obs.co2(JOB) = alldriver(JOB,1).oem.finalrates(1);
   obs.co2_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(1); 
   obs.o3(JOB) = alldriver(JOB,1).oem.finalrates(2); 
   obs.o3_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(2); 
   obs.n2o(JOB) = alldriver(JOB,1).oem.finalrates(3); 
   obs.n2o_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(3); 
   obs.ch4(JOB) = alldriver(JOB,1).oem.finalrates(4); 
   obs.ch4_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(4); 
   obs.cfc11(JOB) = alldriver(JOB,1).oem.finalrates(5); 
   obs.cfc11_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(5); 
   obs.sst(JOB) = alldriver(JOB,1).oem.finalrates(6); 
   obs.sst_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(6); 
end

load Test/alldriver_fit_cal.mat

for JOB=1:36
   cal.co2(JOB) = alldriver(JOB,1).oem.finalrates(1);
   cal.co2_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(1); 
   cal.o3(JOB) = alldriver(JOB,1).oem.finalrates(2); 
   cal.o3_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(2); 
   cal.n2o(JOB) = alldriver(JOB,1).oem.finalrates(3); 
   cal.n2o_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(3); 
   cal.ch4(JOB) = alldriver(JOB,1).oem.finalrates(4); 
   cal.ch4_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(4); 
   cal.cfc11(JOB) = alldriver(JOB,1).oem.finalrates(5); 
   cal.cfc11_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(5); 
   cal.sst(JOB) = alldriver(JOB,1).oem.finalrates(6); 
   cal.sst_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(6); 
end

load Test/alldriver_fit_bias.mat

for JOB=1:36
   bias.co2(JOB) = alldriver(JOB,1).oem.finalrates(1);
   bias.co2_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(1); 
   bias.o3(JOB) = alldriver(JOB,1).oem.finalrates(2); 
   bias.o3_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(2); 
   bias.n2o(JOB) = alldriver(JOB,1).oem.finalrates(3); 
   bias.n2o_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(3); 
   bias.ch4(JOB) = alldriver(JOB,1).oem.finalrates(4); 
   bias.ch4_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(4); 
   bias.cfc11(JOB) = alldriver(JOB,1).oem.finalrates(5); 
   bias.cfc11_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(5); 
   bias.sst(JOB) = alldriver(JOB,1).oem.finalrates(6); 
   bias.sst_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(6); 
end
% 
% Get in-situ CO2
load /asl/s1/tangborn/trace_gas_insitu/co2_gv_rate_09_2002_08_2010;
load /asl/s1/tangborn/trace_gas_insitu/latbins; 

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

% Get in-situ ERA SST

load /asl/s1/rates/Clear/Profile_rates/overocean_gsx_1day_clr_era_lays_spanday01_profilerates_Nov02_2012_robust_span_09_2002_08_2012;
load /asl/s1/tangborn/trace_gas_insitu/overocean_gsx_1day_clr_era_lays_spanday01_fatsummary_Nov02_2012_ALLcov;

% stemprate is grib sst rates
% stempratesstd are errors


