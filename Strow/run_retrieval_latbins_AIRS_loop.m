%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS.m
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Select obs or cal. We do both to get "obs-cal" rates.

% Select the latitude bin
for JOB=1:36 
driver.iibin = JOB;
ix = driver.iibin;
%---------------------------------------------------------------------------
% Need oem_pkg
addpath ../../oem_pkg
%---------------------------------------------------------------------------
% Doing debug?
driver.debug = false;
driver.debug_dir = '../Debug';

% Open debug file if desired
if driver.debug
  writelog('open');
end;
%---------------------------------------------------------------------------
% Perform OEM fit?
driver.oem.dofit = true;
driver.lls.dofit = false;

% Oem loops?  Just one if linear.
driver.oem.nloop = 1;

% Debug plots inside rodgers?
driver.oem.doplots = false;
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Override many settings and add covariance matrix
[driver,aux] = strow_override_defaults_latbins_AIRS(driver);
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Do the retrieval
driver = retrieval(driver,aux);
%---------------------------------------------------------------------------
% Save retrieval output
driver.filename = ['Output/test' int2str(driver.iibin)];
%save(driver.filename,'-struct','driver');
%---------------------------------------------------------------------------
% Close debug file
if driver.debug
  writelog('close')
end

% Some simple output
fprintf('Scalar Retrievals from OEM\n')
fprintf('Lat Index %5.1f \n',JOB)
fprintf(1,'CO2   (ppm)   %5.3f  +- %5.3f \n',driver.oem.finalrates(1),driver.oem.finalsigs(1));
fprintf(1,'O3    (%%)     %5.3f  +- %5.3f \n',100*driver.oem.finalrates(2),100*driver.oem.finalsigs(2));
fprintf(1,'N2O   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(3),driver.oem.finalsigs(3));
fprintf(1,'CH4   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(4),driver.oem.finalsigs(4));
fprintf(1,'CFC11 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(5),driver.oem.finalsigs(5));
fprintf(1,'SST   (K)    %5.3f  +- %5.3f \n',driver.oem.finalrates(6),driver.oem.finalsigs(6));
%---------------------------------------------------------------------------
% Plot Results
addpath Plotutils
%plot_retrieval_latbins


% Save trace gas rates for all lat bins

if iob==1
   co2(JOB) = driver.oem.finalrates(1);
   co2_sigs(JOB) = driver.oem.finalsigs(1); 
   o3(JOB) = driver.oem.finalrates(2); 
   o3_sigs(JOB) = driver.oem.finalsigs(2); 
   n2o(JOB) = driver.oem.finalrates(3); 
   n2o_sigs(JOB) = driver.oem.finalsigs(3); 
   ch4(JOB) = driver.oem.finalrates(4); 
   ch4_sigs(JOB) = driver.oem.finalsigs(4); 
   cfc11(JOB) = driver.oem.finalrates(5); 
   cfc11_sigs(JOB) = driver.oem.finalsigs(5); 
   sst(JOB) = driver.oem.finalrates(6); 
   sst_sigs(JOB) = driver.oem.finalsigs(6); 
end 
if iob==2 
   co2_cal(JOB) = driver.oem.finalrates(1);
   co2_cal_sigs(JOB) = driver.oem.finalsigs(1);
   o3_cal(JOB) = driver.oem.finalrates(2);
   o3_cal_sigs(JOB) = driver.oem.finalsigs(2);
   n2o_cal(JOB) = driver.oem.finalrates(3);
   n2o_cal_sigs(JOB) = driver.oem.finalsigs(3);
   ch4_cal(JOB) = driver.oem.finalrates(4);
   ch4_cal_sigs(JOB) = driver.oem.finalsigs(4);
   cfc11_cal(JOB) = driver.oem.finalrates(5);
   cfc11_cal_sigs(JOB) = driver.oem.finalsigs(5);
   sst_cal(JOB) = driver.oem.finalrates(6);
   sst_cal_sigs(JOB) = driver.oem.finalsigs(6);
end 


alldriver(JOB,iob) = driver;

end % end of latbin loop  


% plot_tracegas_latbins; 
