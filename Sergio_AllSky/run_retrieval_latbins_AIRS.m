%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS.m
% clustcmd -n 32 -q short run_retrieval_latbins_AIRS.m 1:36
%
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Select the latitude bin
driver.iibin = JOB;
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
[driver,aux] = strow_override_defaults_latbins_AIRS_6_6_97_97(driver);  %% 6 usual jacs, 6 cld jacs
%[driver,aux] = strow_override_defaults_latbins_AIRS_6_4_97_97(driver);  %% 6 usual jacs, 4 cld jacs
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Do the retrieval
driver = retrieval(driver,aux);
%---------------------------------------------------------------------------
% Save retrieval output
driver.filename = ['Output/test' int2str(driver.iibin)];
save(driver.filename,'-struct','driver');
%---------------------------------------------------------------------------
% Close debug file
if driver.debug
  writelog('close')
end
%---------------------------------------------------------------------------
% Some simple output
fprintf('Scalar Retrievals from OEM\n')
fprintf(1,'CO2   (ppm)   %5.3f  +- %5.3f \n',driver.oem.finalrates(1),driver.oem.finalsigs(1));
fprintf(1,'O3    (%%)     %5.3f  +- %5.3f \n',100*driver.oem.finalrates(2),100*driver.oem.finalsigs(2));
fprintf(1,'N2O   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(3),driver.oem.finalsigs(3));
fprintf(1,'CH4   (ppb)   %5.3f  +- %5.3f \n',driver.oem.finalrates(4),driver.oem.finalsigs(4));
fprintf(1,'CFC11 (ppt)  %5.3f  +- %5.3f \n',driver.oem.finalrates(5),driver.oem.finalsigs(5));
fprintf(1,'SST   (K)    %5.3f  +- %5.3f \n',driver.oem.finalrates(6),driver.oem.finalsigs(6));
%---------------------------------------------------------------------------
% Plot Results
addpath Plotutils
if exist('f','var')
  plot_retrieval_latbins
end
%---------------------------------------------------------------------------
co2(JOB) = driver.oem.finalrates(1);

