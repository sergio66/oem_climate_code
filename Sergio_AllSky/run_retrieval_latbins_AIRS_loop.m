%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS_loop.m
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
addpath ../../oem_pkg
addpath Plotutils
load lat
load Data/airs_f
%---------------------------------------------------------------------------
% Doing debug?
   driver.debug = false;
   driver.debug_dir = '../Debug';

% Open debug file if desired
   if driver.debug
      writelog('open');
   end;
%---------------------------------------------------------------------------
% Loop over latitude bins
%---------------------------------------------------------------------------
for JOB=1:36 
   driver.iibin = JOB;
% Perform OEM fit?
   driver.oem.dofit = true;
   driver.lls.dofit = false;
% Oem loops?  Just one if linear.
   driver.oem.nloop = 1;
% Debug plots inside rodgers?
   driver.oem.doplots = false;
%---------------------------------------------------------------------------
% Override many settings and add covariance matrix
   [driver,aux] = strow_override_defaults_latbins_AIRS(driver);
%---------------------------------------------------------------------------
% Do the retrieval
   driver = retrieval(driver,aux);
%---------------------------------------------------------------------------
% Save retrieval output from this loop
   alld(JOB) = driver;
%---------------------------------------------------------------------------
% Some simple output
   fprintf('Scalar Retrievals from OEM\n')
   fprintf('Lat Index %5.1f  Lat: %5.1f \n',JOB,lat(JOB))
   fprintf(1,'CO2   (ppm)   %5.3f  +- %5.3f \n',driver.oem.finalrates(1),driver.oem.finalsigs(1));
   fprintf(1,'O3    (%%)     %5.3f  +- %5.3f \n',100*driver.oem.finalrates(2),100*driver.oem.finalsigs(2));
   fprintf(1,'N2O   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(3),driver.oem.finalsigs(3));
   fprintf(1,'CH4   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(4),driver.oem.finalsigs(4));
   fprintf(1,'CFC11 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(5),driver.oem.finalsigs(5));
   fprintf(1,'SST   (K)    %5.3f  +- %5.3f \n',driver.oem.finalrates(6),driver.oem.finalsigs(6));
   %---------------------------------------------------------------------------
% Pull interesting variable out for quick look
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

% Plot Results
   plot_retrieval_latbins
   disp('Hit return for next latitude')
   pause
   close all

end % end of latbin loop  
%---------------------------------------------------------------------------
% Close debug file
if driver.debug
   writelog('close')
end
%---------------------------------------------------------------------------
plot_tracegas_latbins; 

% To save everything
% save filename alld, or if include special vars
% save obs_no_ozone_chans  alld ...
%        co2 co2_sigs ...
%        o3 o3_sigs   ...
%        n2o n2o_sigs ...
%        ch4 ch4_sigs ...
%        cfc11 cfc11_sigs ...
%        sst sst_sigs