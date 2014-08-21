%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% run_retrieval_latbins_AIRS.m
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Select obs or cal. We do both to get "obs-cal" rates.

for iob=1:2 


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
% Raw rate data file
driver.rateset.datafile  = '/asl/s1/rates/clear/Aug2013/overocean_gsx_1day_clr_era_lays_spanday01_avgL1Brates_robust_Nov02_2012_span_09_2002_08_2012.mat';
%driver.rateset.datafile = '/home/sergio/MATLABCODE/RATES_TARO/ANOM/merra_linearratesfor_andy.mat';

% Fitting [obs][cal][bias], pick one  This is moved to outer loop for plotting
%driver.rateset.ocb_set  = 'obs';  

   if iob==1
     driver.rateset.ocb_set  = 'obs';
   end
   if iob==2
     driver.rateset.ocb_set  = 'cal';
   end


% Good channel set
load /asl/s1/rates/clear/good_chanset.mat 
driver.jacobian.chanset = chanset;

% Lag-1 correlation file; if using rate least-squares errors
driver.rateset.ncfile   = '../../oem_pkg/Test/all_lagcor.mat';

% Get rate data, do Q/A elsewhere
driver = get_rates(driver);

%---------------------------------------------------------------------------
% Jacobian file: f = 2378x1 and M_TS_jac_all = 36x2378x200
driver.jacobian.filename = '../../oem_pkg/Test/M_TS_jac_all.mat';
driver.jacobian.varname  = 'M_TS_jac_all';
driver.jacobian.scalar_i = 1:6;
driver.jacobian.water_i  = 7:103;
driver.jacobian.temp_i   = 104:200;
driver.jacobian.numlays  = 97;


% Get jacobians
jac             = load(driver.jacobian.filename);
aux.m_ts_jac    = squeeze(jac.M_TS_jac_all(ix,:,:));
driver.qrenorm  = jac.qrenorm;
[jac_max i_max(ix)] = max(abs(aux.m_ts_jac(1147,7:103)));
f = jac.f;
clear jac

% Test for removing CO2 from obs driver.rateset.rates=driver.rateset.rates - 1.0*aux.m_ts_jac(:,1); 

%---------------------------------------------------------------------------
% Apriori file
driver.oem.apriori_filename = 'apriori_lls';

% Load in apriori
xb = load(driver.oem.apriori_filename,'apriori');
xb = xb.apriori;
[mm,nn] = size(xb);
if nn > 1
  xb = xb(:,driver.ix);
end
% A Priori stored in aux.xb
aux.xb = xb./driver.qrenorm';
%---------------------------------------------------------------------------
% SARTA forward model and other "representation" errors
driver.oem.sarta_error = 0.0;
%---------------------------------------------------------------------------
% Override many settings and add covariance matrix
driver = strow_override_defaults_latbins_AIRS(driver);
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Do the retrieval
driver = retrieval(driver,aux);
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% Save retrieval output
driver.filename = ['Output/test' int2str(driver.iibin)];
%save(driver.filename,'-struct','driver');
% Close debug file
if driver.debug
  writelog('close')
end
%---------------------------------------------------------------------------
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

end  % end of iob loop 


plot_tracegas_latbins; 
