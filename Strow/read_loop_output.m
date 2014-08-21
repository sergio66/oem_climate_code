load Test/output_alllats_obs_cal_rate_offset_0p1K

for JOB=1:36
   off.co2(JOB) = alldriver(JOB,1).oem.finalrates(1);
   off.co2_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(1); 
   off.o3(JOB) = alldriver(JOB,1).oem.finalrates(2); 
   off.o3_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(2); 
   off.n2o(JOB) = alldriver(JOB,1).oem.finalrates(3); 
   off.n2o_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(3); 
   off.ch4(JOB) = alldriver(JOB,1).oem.finalrates(4); 
   off.ch4_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(4); 
   off.cfc11(JOB) = alldriver(JOB,1).oem.finalrates(5); 
   off.cfc11_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(5); 
   off.sst(JOB) = alldriver(JOB,1).oem.finalrates(6); 
   off.sst_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(6); 
end

load Test/output_alllats_obs_cal.mat

for JOB=1:36
   base.co2(JOB) = alldriver(JOB,1).oem.finalrates(1);
   base.co2_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(1); 
   base.o3(JOB) = alldriver(JOB,1).oem.finalrates(2); 
   base.o3_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(2); 
   base.n2o(JOB) = alldriver(JOB,1).oem.finalrates(3); 
   base.n2o_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(3); 
   base.ch4(JOB) = alldriver(JOB,1).oem.finalrates(4); 
   base.ch4_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(4); 
   base.cfc11(JOB) = alldriver(JOB,1).oem.finalrates(5); 
   base.cfc11_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(5); 
   base.sst(JOB) = alldriver(JOB,1).oem.finalrates(6); 
   base.sst_sigs(JOB) = alldriver(JOB,1).oem.finalsigs(6); 
end

for JOB=1:36
   calc.co2(JOB) = alldriver(JOB,2).oem.finalrates(1);
   calc.co2_sigs(JOB) = alldriver(JOB,2).oem.finalsigs(1); 
   calc.o3(JOB) = alldriver(JOB,2).oem.finalrates(2); 
   calc.o3_sigs(JOB) = alldriver(JOB,2).oem.finalsigs(2); 
   calc.n2o(JOB) = alldriver(JOB,2).oem.finalrates(3); 
   calc.n2o_sigs(JOB) = alldriver(JOB,2).oem.finalsigs(3); 
   calc.ch4(JOB) = alldriver(JOB,2).oem.finalrates(4); 
   calc.ch4_sigs(JOB) = alldriver(JOB,2).oem.finalsigs(4); 
   calc.cfc11(JOB) = alldriver(JOB,2).oem.finalrates(5); 
   calc.cfc11_sigs(JOB) = alldriver(JOB,2).oem.finalsigs(5); 
   calc.sst(JOB) = alldriver(JOB,2).oem.finalrates(6); 
   calc.sst_sigs(JOB) = alldriver(JOB,2).oem.finalsigs(6); 
end

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

load /asl/s1/rates/clear/Aug2013/overocean_gsx_1day_clr_era_lays_spanday01_profilerates_Nov02_2012_robust_span_09_2002_08_2012;
load /asl/s1/tangborn/trace_gas_insitu/overocean_gsx_1day_clr_era_lays_spanday01_fatsummary_Nov02_2012_ALLcov;

% stemprate is grib sst rates
% stempratesstd are errors


