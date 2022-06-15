function [driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,iNlays_retrieve,topts);

settings.resetnorm2one = -1; %%% default, keep my scaling, set to +1 if you want to reset all to 1.00000000

settings.dataset = 1;       % (1) AIRS 16 year dataset (2) AIRS 11 year (IASI) (3) IASI2AIRS 11 year
settings.co2lays = 1;       % assume column jac for CO2, or 3 lays (gnd-500,500-trop,trop-TOA

settings.ocb_set = 0;       % 0 = obs, 1 = cal, -1 = bias
settings.numchan = 2645;    % L1b = 2378, L1c = 2645
settings.numchan = 1305;    % NSR
settings.chan_LW_SW = 0;    % 0 is 640 to 1640, 1 is 700 to 1640, -1 is 640 to 2740

settings.set_tracegas = -1;            %% do we leave apriori as 0 or set CO2/N2o/CH4/CFC to be 2.2, 1, 4.5, -1             
                                       %%   (-1 = no, +1 = yes)
                                       %%   note : if anomaly, adjust a priori for CO2,N2O,CH4,CFC
settings.offsetrates  = -1;            %% do we add a constant offset to the spectral rates (-1 = no, +1 = yes)
                                       %% if driver.rateset.ocb_set  == 'obs';
settings.addco2jacs   = -1;            %% do we add co2 jacs to the spectral rates (-1 = no, +1 = yes)
                                       %% if driver.rateset.ocb_set  == 'cal';
settings.obs_corr_matrix = -1;         %% just use nc_error (-1) or try to be fancy and use full cov matrix (+1)

settings.tie_sst_lowestlayer = +1;     %% tie together SST with lowest T(z)
settings.invtype         = 1;          %% pinv, see /home/sergio/MATLABCODE/oem_pkg/rodgers.m
settings.iNlays_retrieve = 97;         %% do all 97 layers
settings.descORasc = +1;               %% descending default
settings.iXJac = 0;                    %% const geo jacs, replace as needed CO2/CH4/N20;  
                                       %% +2 uses kCARTA varying geo/trace, +1 uses SARTA varying geo/trace jacs, 0 = constant kcarta jacs
settings.iDoStrowFiniteJac = -1;       %% -1 : do not change the time varying anomaly jacs                                done for all anomaly timesteps
                                       %% +1 stick to Sergio tracegas jacs = BT(1.001 X(t,latbin)) - BT(1.00 X(t,latbin)) interp in time
                                       %% +2 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin)))      interp in time .. 
                                       %% +3 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin)))      done for all anomaly timesteps DEFAULT
                                       %% +4 stick to Strow  tracegas jacs = BT(X(t,latbin)) - BT(2002 X(t,latbin)))      done for all anomaly timesteps with Age of Air for CO2
settings.iChSet = 1;                   %% +1 default, old chans (about 500)
                                       %% +2, new chans (about 400) with CFC11,CFC12      and weak WV, bad chans gone
                                       %% +3, new chans (about 400) w/o  CFC11 with CFC12 and weak WV, bad chans gone
settings.iFixTz_NoFit = -1;            %% -1 : do not fix Tz to ERA anomaly T(z,time) values, then fit for Tz
                                       %% +1 : do     fix Tz to ERA anomaly T(z,time) values, then keep Tz fixed (ie do not fit)
settings.iFixO3_NoFit = -1;            %% -1 : do not fix O3 to ERA anomaly O3(z,time) values, then fit for O3
                                       %% +1 : do     fix O3 to ERA anomaly O3(z,time) values, then keep O3 fixed (ie do not fit)
                                       %% +0 : do     fix O3 to zero        O3(z,time) values, then keep O3 fixed (ie do not fit)
settings.iFixTG_NoFit = -1;            %% -1 means retrieve all trace gases [CO2 N2O CH4 CFC11 CFC12]
                                       %% eg [4 5] means do not do CFC11,CFC12
                                       %% eg [4] means do not do CFC11

allowedparams = [{'ocb_set'},{'numchan'},{'chan_LW_SW'},{'iChSet'},{'set_tracegas'},{'offsetrates'},...
			    {'addco2jacs'},{'obs_corr_matrix'},{'invtype'},{'tie_sst_lowestlayer'},{'iNlays_retrieve'},...
                            {'descORasc'},{'dataset'},{'iXJac'},{'co2lays'},{'iDoStrowFiniteJac'},{'iFixTz_NoFit'},{'iFixO3_NoFit'},{'iFixTG_NoFit'},{'resetnorm2one'}];


%disp('settings before')
%settings

if nargin == 3
  optvar = fieldnames(topts);
  for i = 1 : length(optvar)
   if (length(intersect(allowedparams,optvar{i})) == 1)
     eval(sprintf('settings.%s = topts.%s;', optvar{i}, optvar{i}));
   else
     fprintf(1,'topts param not in allowed list ... %s \n',optvar{i});
     error('quitting ');
   end
 end
end

%disp('settings after')
%settings

aux.invtype = settings.invtype;

driver.topts = topts;

%---------------------------------------------------------------------------
% Which latitude bin
ix = driver.iibin;
%---------------------------------------------------------------------------
% Fitting [obs][cal][bias], pick one
if settings.ocb_set == -1
  driver.rateset.ocb_set  = 'bias';
  disp('bias')
elseif settings.ocb_set == 0
  driver.rateset.ocb_set  = 'obs';
  disp('obs')
elseif settings.ocb_set == +1
  driver.rateset.ocb_set  = 'cal';
  disp('cal')
elseif abs(settings.ocb_set) > 1
  settings.ocb_set
  error('incorrect settings.ocb_set')
end
%---------------------------------------------------------------------------
% Raw rate data file        
if settings.dataset == -1   
  disp('AIRS 16 year rates or anomalies, NO nu cal done')
  if settings.descORasc == +1 & driver.i16daytimestep < 0
    disp('doing descending latbin rates')
    error('oops not done')
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 0
    disp('doing descending OBS ANOMALY')
    driver.rateset.datafile = ['ANOM_16dayavg_nonucal/latbin_0dayavg_' num2str(driver.iibin) '.mat'];  
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 1
    disp('doing descending CAL ANOMALY')
    driver.rateset.datafile = ['ANOM_16dayavg_nonucal/latbin_0dayavg_' num2str(driver.iibin) '_cal.mat'];  
  end

elseif settings.dataset == 5   %% CRIS LORES DEFAULT
  disp('CRIS 6 years 2012/05 - 2018/04 6 year rates or anomalies')
  if settings.descORasc == +1 & driver.i16daytimestep < 0
    disp('doing descending latbin rates')
    driver.rateset.datafile  = 'convert_strowrates2oemrates_random_6_year_v32_clear_nucal.mat';
    if settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_random_6_year_v32_clear.mat';
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_random_6_year_v32_clear_nucal.mat';
    end
  elseif settings.descORasc == -1 & driver.i16daytimestep < 0
    disp('doing ascending latbin rates')
    driver.rateset.datafile  = 'convert_strowrates2oemrates_random_6_year_v32_clearasc_nucal.mat'; %% note bias/cal are from "usual" not nucal
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 0
    disp('doing descending OBS ANOMALY')
    driver.rateset.datafile = ['ANOM_16dayavg/latbin_0dayavg_' num2str(driver.iibin) '.mat'];  
    driver.rateset.datafile = ['ANOM_16dayavgDEBUG/latbin_0dayavg_' num2str(driver.iibin) '.mat'];  
    driver.rateset.datafile = ['ANOM_16dayavgDEBUG/sergio_latbin_0dayavg_' num2str(driver.iibin) '.mat'];  
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 1
    disp('doing descending CAL ANOMALY')
    driver.rateset.datafile = ['ANOM_16dayavg/latbin_0dayavg_' num2str(driver.iibin) '_cal.mat'];  
    driver.rateset.datafile = ['ANOM_16dayavgDEBUG/sergio_latbin_0dayavg_' num2str(driver.iibin) '_cal.mat'];  
  end

end

fprintf(1,' <<< driver.rateset.datafile >>> = %s \n',driver.rateset.datafile)

% Lag-1 correlation file; if using rate least-squares errors
driver.rateset.ncfile   = '../oem_pkg/Test/all_lagcor.mat';
driver.rateset.ncfile   = driver.rateset.datafile;

% Get rate data, do Q/A elsewhere
driver = get_rates(driver);  %% this gets spectral rates (driver.rateset.rates), and uncertainty (driver.rateset.unc_rates)
%keyboard_nowindow

%---------------------------------------------------------------------------
% Jacobian file: f = 2378x1 and M_TS_jac_all = 36x2378x200
% driver.jacobian.filename = '/home/sergio/MATLABCODE/oem_pkg/Test/M_TS_jac_all.mat';
%% clear sky

iXJac = settings.iXJac;
%if driver.i16daytimestep > 0
%  iXJac = 0; %% const geo kcarta jcs
%  iXJac = 1; %% varying geo sarta jacs
%  iXJac = 2; %% varying geo kcarta jacs
%end

if driver.i16daytimestep < 0
  if settings.descORasc == +1
    driver.jacobian.filename = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/JUNK_CRIS/kcarta_M_TS_jac_all_5_97_97_97_2235.mat';

    %% we can "fool" the code by using midpoint anomaly jac
    junk = num2str(180,'%03d');
    driver.jacobian.filename = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CLO_Anomaly137_16_12p8/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2235.mat']; 

    fprintf(1,'reading in constant kcarta jac file %s \n',driver.jacobian.filename)
  else
    %% for now assume same jacs
    driver.jacobian.filename = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/JUNK/kcarta_M_TS_jac_all_5_97_97_97_2235.mat';

    %% we can "fool" the code by using midpoint anomaly jac
    junk = num2str(180,'%03d');
    driver.jacobian.filename = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CLO_Anomaly137_16_12p8/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2235.mat']; 

    fprintf(1,'reading in constant kcarta jac file %s \n',driver.jacobian.filename)
  end
elseif driver.i16daytimestep > 0
  junk = num2str(driver.i16daytimestep,'%03d');
  %junk = num2str(150,'%03d');    %% testing
  if iXJac == 1
    %% sarta time vary jacs
    %asarta  = load('../MakeJacsSARTA/SARTA_AIRSL1c_CLO_Anomaly137_16/RESULTS/sarta_182_fixCFC_M_TS_jac_all_5_97_97_97_2235.mat');
    driver.jacobian.filename = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_CLO_Anomaly137_16_with_seasonal_OldSarta_largepert//RESULTS/sarta_' junk '_fixCFC_M_TS_jac_all_5_97_97_97_2235.mat']; %% old sarta
    driver.jacobian.filename = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_CLO_Anomaly137_16_no_seasonal_OldSarta_smallpert//RESULTS/sarta_' junk '_fixCFC_M_TS_jac_all_5_97_97_97_2235.mat'];   %% old sarta
    driver.jacobian.filename = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacsSARTA/SARTA_AIRSL1c_CLO_Anomaly137_16/RESULTS/sarta_' junk '_fixCFC_M_TS_jac_all_5_97_97_97_2235.mat'];                         %% new sarta
    fprintf(1,'iXJac == 1 reading in timestep sarta jac file %s \n',driver.jacobian.filename)

  elseif iXJac == 2  
    %% kcarta time vary jac
    %% after Dec 2019, the finite diff tracegas jacs used dQ = 0.001,dT = 0.01 GOOD tracegas profiles (basically the glatm.dat tracegas profiles for CO2/N2O/CH4 adjusted in time)
    %% the profiles do not have seasonal
    %% use with  together with iXJac = 2,iDoStrowFiniteJac = -1
    driver.jacobian.filename = ['/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/CLO_Anomaly137_16_12p8/RESULTS/kcarta_' junk '_M_TS_jac_all_5_97_97_97_2235.mat']; 

    fprintf(1,'iXJac == 2 reading in timestep kcarta jac file %s \n',driver.jacobian.filename)

  elseif iXJac == 0
    %% constant kcarta jacs
    driver.jacobian.filename = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeJacskCARTA_CLR/JUNK_CRIS/kcarta_M_TS_jac_all_5_97_97_97_2235.mat';
    fprintf(1,'iXJac == 0 reading in constant kcarta jac file %s \n',driver.jacobian.filename)
  end
end

%% driver.jacobian.filename

%% THIS IS DEFAULT -- 4 column trace gas (CO2/N2O/CH4/CFC11/CFC12), 1 stemp, (97x3) geo
driver.jacobian.varname  = 'M_TS_jac_all';
driver.jacobian.scalar_i = 1:6;
driver.jacobian.water_i  = 7:103;
driver.jacobian.temp_i   = 104:200;
driver.jacobian.ozone_i  = 201:297;
driver.jacobian.numlays  = 97;

% Get jacobians, and combine the 97 layer T(z)/WV(z)/O3(z) into N layers
jac               = load(driver.jacobian.filename);
jac.qrenorm(driver.jacobian.temp_i)          = 0.1;                                              %% bump up T jac
jac.M_TS_jac_all(:,:,driver.jacobian.temp_i) = 10*jac.M_TS_jac_all(:,:,driver.jacobian.temp_i);  %% bump up T jac

%% add in extra column for CFC12 >>>>>>
if iXJac == 0 | iXJac == 1
  m_ts_jac0_noCFC12 = squeeze(jac.M_TS_jac_all(ix,:,:));
  m_ts_jac0 = zeros(2235,297);
  m_ts_jac0(:,1:4)   = m_ts_jac0_noCFC12(:,1:4);
  m_ts_jac0(:,6:297) = m_ts_jac0_noCFC12(:,5:296);

  xm_ts_jac_coljac   = m_ts_jac0(:,1:6);
  if driver.i16daytimestep > 0  
    xm_ts_jac_coljac = replace_time_cfc12jac(xm_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3);
  else
    xm_ts_jac_coljac = replace_time_cfc12jac(xm_ts_jac_coljac,driver.iibin,floor(365/2),3);
  end
  m_ts_jac0(:,5) = xm_ts_jac_coljac(:,5);
  clear xm_ts_jac_coljac xm_ts_jac_coljac

  qrenormjunk = zeros(1,length(jac.qrenorm)+1);
  qrenormjunk(1:4)   = jac.qrenorm(1:4);
  qrenormjunk(5)     = jac.qrenorm(4);
  qrenormjunk(6:297) = jac.qrenorm(5:296);
  jac.qrenorm = qrenormjunk;
else
  m_ts_jac0 = squeeze(jac.M_TS_jac_all(ix,:,:));
  %% oops forgot to fix qrenorm, will do later
  %% qrenormjunk = zeros(1,length(jac.qrenorm)+1);
  %% qrenormjunk(1:4)   = jac.qrenorm(1:4);
  %% qrenormjunk(5)     = jac.qrenorm(4);
  %% qrenormjunk(6:297) = jac.qrenorm(5:296);
  %% jac.qrenorm = qrenormjunk;
end

m_ts_jac_coljac   = m_ts_jac0(:,1:6);
driver.qrenorm  = jac.qrenorm;       %% set this default

if iNlays_retrieve <= 60
  [m_ts_jac_wv,qWV,layWV]  = combinejaclays(m_ts_jac0,driver.jacobian.water_i,jac.qrenorm,iNlays_retrieve);
  [m_ts_jac_t,qT,layT]     = combinejaclays(m_ts_jac0,driver.jacobian.temp_i, jac.qrenorm,iNlays_retrieve);
  [m_ts_jac_o3,qO3,layO3]  = combinejaclays(m_ts_jac0,driver.jacobian.ozone_i,jac.qrenorm,iNlays_retrieve);
else
  fprintf(1,'setting iNlays_retrieve ( > 60) from %2i to 97 \n',iNlays_retrieve);
  iNlays_retrieve = 97;
end

%% replace CO2,N2O,CH4 jacs
if driver.i16daytimestep > 0  
  %% this is for 365 anomaly time steps
  %% put in time varying Jacobian, err no more need to do this??? well sarta has older CO2/CH4 but let's comment this for now
  iDoStrowFiniteJac = +2; %% testing Strows finite difference jacs CO2(t)-CO2(370) ...    6/24-27/2019 interp in time, used to be +1
  iDoStrowFiniteJac = +3; %% testing Strows finite difference jacs CO2(t)-CO2(370) ...    6/24-27/2019 at all anom timesteps
  iDoStrowFiniteJac = +4; %% testing Strows finite difference jacs CO2(t)-CO2(370) ...    6/24-27/2019 at all anom timesteps, age of air
  iDoStrowFiniteJac = +1; %% testing new Sergio finite diff jacs                          interp in time
  iDoStrowFiniteJac = -1; %% default, rely on time varying CO2/N20/CH4 jacs from kcarta,  done for all anom timsteps

  iDoStrowFiniteJac = settings.iDoStrowFiniteJac; %% from 6/29/2019

  if iXJac == 0 & iDoStrowFiniteJac == 1
    fprintf(1,'updating const kCARTA CO2/N2O/CH4 jacs with Sergio interpolated time varying jacs...\n');
    %% const kCARTA jacs, update the trace gases
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);

  elseif iXJac == 2 & iDoStrowFiniteJac == 2
    fprintf(1,'updating time varying kCARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs -1,+2 for testing 6/24-27/2019...\n');
    fprintf(1,'  note before July2, the tracegas profile (CO2/N2O/CH4) was US Std shoehorned and multiplied, so ppm was quite wonky except at 500 mb');
    fprintf(1,'  else turned off \n')
    %% const kCARTA jacs, update the trace gases
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% only used this on 6/26
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% added this on 6/27
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% added this on 6/27
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% added this on 8/21
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% added this on 8/29 ---> comment this out if iUgh == 4

  elseif iXJac == 2 & iDoStrowFiniteJac == 3
    fprintf(1,'updating time varying kCARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 6/24-27/2019...\n');
    fprintf(1,'  note before July2, the tracegas profile (CO2/N2O/CH4) was US Std shoehorned and multiplied, so ppm was quite wonky except at 500 mb');
    fprintf(1,'  note after July2, have improved the tracegas profile (CO2/N2O/CH4) so they are the same shape as glatm.dat');
    fprintf(1,'  else turned off \n')
    %% const kCARTA jacs, update the trace gases
    %% kcarta strow finite jacs
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% only used this on 6/26
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% added this on 6/27
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% added this on 6/27
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% added this on 8/21
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% added this on 8/29 ---> comment this out if iUgh == 4 <<<>>><<<>>><<<>>>

  elseif iXJac == 1 & iDoStrowFiniteJac == 3
    fprintf(1,'updating time varying SARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 6/24-27/2019...\n');
    fprintf(1,'  note before July2, the tracegas profile (CO2/N2O/CH4) was US Std shoehorned and multiplied, so ppm was quite wonky except at 500 mb');
    fprintf(1,'  note after July2, have improved the tracegas profile (CO2/N2O/CH4) so they are the same shape as glatm.dat');
    fprintf(1,'  else turned off \n')
    %% sarta strow finite jacs
    iVarType = -3; %% this uses SARTA  finitediff jacs, which I have shown are bad?? or good??
    iVarType = +3; %% this uses kCARTA finitediff jacs, which I have shown are good, just want to test
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/21
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/21 --> comment this out if iUgh == 4

  elseif iXJac == 2 & iDoStrowFiniteJac == 4
    fprintf(1,'updating time varying kCARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 12/08/2019 ... age of air\n');
    fprintf(1,'  note before July2, the tracegas profile (CO2/N2O/CH4) was US Std shoehorned and multiplied, so ppm was quite wonky except at 500 mb');
    fprintf(1,'  note after July2, have improved the tracegas profile (CO2/N2O/CH4) so they are the same shape as glatm.dat');
    fprintf(1,'  else turned off \n')
    %% const kCARTA jacs, update the trace gases
    %% kcarta strow finite jacs
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,4); %% only used this on 6/26
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,4); %% added this on 6/27
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,4); %% added this on 6/27
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,4); %% added this on 8/21
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,4); %% added this on 8/29 ---> comment this out if iUgh == 4 <<<>>><<<>>><<<>>>

  elseif iXJac == 1 & iDoStrowFiniteJac == 3
    fprintf(1,'updating time varying SARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 12/08/2019...\n');
    fprintf(1,'  note before July2, the tracegas profile (CO2/N2O/CH4) was US Std shoehorned and multiplied, so ppm was quite wonky except at 500 mb');
    fprintf(1,'  note after July2, have improved the tracegas profile (CO2/N2O/CH4) so they are the same shape as glatm.dat');
    fprintf(1,'  else turned off \n')
    %% sarta strow finite jacs
    iVarType = -3; %% this uses SARTA  finitediff jacs, which I have shown are bad?? or good??
    iVarType = +3; %% this uses kCARTA finitediff jacs, which I have shown are good, just want to test
    iVarType = +4; %% this uses kCARTA finitediff jacs, which I have shown are good, just want to test
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/21
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/21 --> comment this out if iUgh == 4
  end
else    
  %% this is for 1 average rate
  iTimeStepUse = 1;            %% put in constant Jacobian, at timestep 1 (2002/09)
  iTimeStepUse = 365;          %% put in constant Jacobian, at timestep 365 (2018/08)
  iTimeStepUse = floor(365/2); %% put in constant Jacobian, half way through (365/2 ==> 2009/09)
  if iXJac == 0
    fprintf(1,'updating CO2/N2O/CH4 jacs ...\n');
    %% const kCARTA jacs, update the trace gases
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,iTimeStepUse);
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,iTimeStepUse);
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,iTimeStepUse);
    m_ts_jac_coljac = replace_time_cfc11jac(m_ts_jac_coljac,driver.iibin,iTimeStepUse); %% added this on 8/21
    m_ts_jac_coljac = replace_time_cfc12jac(m_ts_jac_coljac,driver.iibin,iTimeStepUse); %% added this on 8/21
  end
end

if settings.co2lays == 3
  m_ts_jac_coljac = replace_time_co2_3layjac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep);
end

if iNlays_retrieve <= 60
  m_ts_jac = [m_ts_jac_coljac m_ts_jac_wv m_ts_jac_t m_ts_jac_o3];
else
  m_ts_jac = m_ts_jac0;
  qWV = (1:iNlays_retrieve) + 6;
end

bad = find(isnan(m_ts_jac) | isinf(m_ts_jac));
if length(bad) > 0
  fprintf('oopsy foound %5i NaN or Inf in jacobian, resetting to 0 \n',length(bad))
  m_ts_jac(bad) = 0;
end

iNlays_retrieve0 = iNlays_retrieve;
iNlays_retrieve = length(qWV);
ixlays = 1:iNlays_retrieve;
if settings.co2lays == 1
  driver.jacobian.scalar_i = 1:6;
  driver.jacobian.wvjaclays_offset = 6;
  if iNlays_retrieve <= 60
    driver.qrenorm  = [jac.qrenorm(1:6) qWV qT qO3];
  end
elseif settings.co2lays == 3
  driver.jacobian.scalar_i = 1:8;
  driver.jacobian.wvjaclays_offset = 8;
  if iNlays_retrieve <= 60
    driver.qrenorm  = [jac.qrenorm(1) jac.qrenorm(1) jac.qrenorm(1) jac.qrenorm(2:6) qWV qT qO3];
  end
end

iBruteForce_co2adjustjac = +1;
iBruteForce_co2adjustjac = -1;
if iBruteForce_co2adjustjac > 0
  disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CO2 jac ---> co2 jac * 0.85 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
  if settings.co2lays == 1
    m_ts_jac(:,1)  = m_ts_jac(:,1) * 0.85;
  elseif settings.co2lays == 3
    m_ts_jac(:,1:3)  = m_ts_jac(:,1:3) * 0.85;
  end
  disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CO2 jac ---> co2 jac * 0.85 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
end

if iNlays_retrieve <= 60
  driver.jacobian.water_i  = max(driver.jacobian.scalar_i) + ixlays;  
  driver.jacobian.temp_i   = max(driver.jacobian.water_i) + ixlays;  
  driver.jacobian.ozone_i  = max(driver.jacobian.temp_i) + ixlays;  
  driver.jacobian.numlays  = iNlays_retrieve;
  driver.jacobian.wvjaclays_used   = layWV;
end

aux.m_ts_jac    = m_ts_jac;
f = jac.f;

addpath /home/sergio/MATLABCODE/CRIS_Hi2Lo/
[fout,jacout] = translate_hi2lo(f,m_ts_jac);

fout = fout.vchan;

i1305 = load('/asl/matlib/cris/ch_std_from1317.mat');
f = fout(i1305.ch_std_i);
m_ts_jac = jacout(i1305.ch_std_i,:);
aux.m_ts_jac = m_ts_jac;

%whos f
%f1305 = f;
%save f1305.mat f1305
%error(';kjgs')

%keyboard_nowindow
%plot(f,m_ts_jac(:,1))

clear jac fout jacout
%---------------------------------------------------------------------------
% Good channel set
lenrates = length(driver.rateset.rates);

iChSet = 2; %% new chans
iChSet = 1; %% old chans (default)
iChSet = 3; %% new chans, but no CFC11
iChSet = topts.iChSet;

finitech = find(isfinite(driver.rateset.rates));
%ch = finitech(1:3:1000);

%keyboard_nowindow
%ch = find_the_oem_channels(f,lenrates,settings.numchan,settings.chan_LW_SW,iChSet);
load ../AIRS_new_clear_scan_August2019_AMT2020PAPER/f2645.mat
ch = find_the_oem_channels(f2645,2645,2645,settings.chan_LW_SW,iChSet);
ch = find_closest_airs2cris_chans(ch);
ch = intersect(finitech,ch);

guardchan = [1 3 716 717 718 719 1153 1154 1155 1156 1316 1317];   %% these are guard chans for 1 : 1317
guardchan = [1 713 714 1146 1147 1305];                            %% loading in eg ANOM_16dayavg/latbin_0dayavg_20.mat can see in latter times these are the NaN channels for 1:1305
ch = setdiff(ch,guardchan);

iUseMoreChans = -1;  %% use the AIRS list only, about 400 chans
iUseMoreChans = +2;  %% use all chans betwen 65- and 1800 cm=1 except for guard chans, abut 1100 chans
iUseMoreChans = +1;  %% gave really good results, June 4, 2022, about 657 chans

if iUseMoreChans > 0
  %% get great results
  f1305 = load('f1305.mat');
  ch1600 = find(f1305.f1305 > 1500 & f1305.f1305 < 1700);
  ch1600 = ch1600(1:1:length(ch1600));
  ch = union(ch,ch1600);
  
  ch700 = find(f1305.f1305 > 700 & f1305.f1305 < 800);
  ch700 = ch700(1:1:length(ch700));
  ch = union(ch,ch700);
  
  ch725 = find(f1305.f1305 > 725 & f1305.f1305 < 800);
  ch725 = find(f1305.f1305 > 655 & f1305.f1305 < 800);
  ch725 = ch725(1:1:length(ch725));
  ch = union(ch,ch725);
  
  if iUseMoreChans == 1
    chLWMW = find(f1305.f1305 > 725 & f1305.f1305 < 1800);
    chLWMW = find(f1305.f1305 > 650 & f1305.f1305 < 1800);
    ch = intersect(ch,chLWMW);
  elseif iUseMoreChans == 2
    chLWMW = find(f1305.f1305 > 725 & f1305.f1305 < 1800);
    chLWMW = find(f1305.f1305 > 650 & f1305.f1305 < 1800);
    ch = union(ch,chLWMW); 
    ch = setdiff(ch,guardchan);
  end
end

driver.topts.iChSet = iChSet;
driver.jacobian.chanset = ch;

%---------------------------------------------------------------------------
% Apriori file
%driver.oem.apriori_filename = 'apriori_lls';

% Load in apriori
%xb = load(driver.oem.apriori_filename,'apriori');
%xb = xb.apriori;

%xb = zeros(296,1);
xb = zeros(driver.jacobian.wvjaclays_offset + iNlays_retrieve*3,1);

if (abs(settings.set_tracegas) ~= 1) & (settings.set_tracegas ~= 2)
  settings.set_tracegas
  error('incorrect setting for overriding xb tracegas values');
end

if settings.iFixTG_NoFit(1) > 0
  fprintf(1,'oh oh you wanna get rid of a trace gas GID(1) = %2i \n',settings.iFixTG_NoFit(1))
  junk = settings.iFixTG_NoFit;
  badjunk = find(junk < 1 | junk > max(driver.jacobian.scalar_i) - 1);
  if length(badjunk) > 0 
    junk
    driver.jacobian.scalar_i
    error('can only thrown out a few of first 1 .. driver.jacobian.scalar_i - 1 gases!!!!')
  else
    numthrow = length(junk);
    fprintf(1,'throwing out %3i trace gases from the fit \n',numthrow);
%{
    aux.jacobian.scalar_i_allgases = 1:driver.jacobian.scalar_i;
    aux.jacobian.water_i_allgases  = driver.jacobian.water_i;
    aux.jacobian.temp_i_allgases   = driver.jacobian.temp_i;
    aux.jacobian.ozone_i_allgases  = driver.jacobian.ozone_i;

    driver.jacobian.scalar_i = 1:driver.jacobian.scalar_i - numthrow;
    driver.jacobian.water_i  = driver.jacobian.water_i - numthrow;
    driver.jacobian.temp_i  = driver.jacobian.temp_i - numthrow;
    driver.jacobian.ozone_i  = driver.jacobian.ozone_i - numthrow;
%}
  end
end
  

if settings.iFixTz_NoFit > 0 & strcmp(driver.rateset.ocb_set,'obs')
  iFixTz_NoFit = +1;    %%% LARABBEE LIKES THIS TURNED OFF ie keep spectra as is, just read in ERA anom and proceed
end
iZeroTVers = 0; %%% use my fit to sarta calcs, as a proxy to ERA T anomalies
iZeroTVers = 1; %%% use raw ERA T anomalies, and do the averaging here on the fly
iZeroTVers = 2; %%% use raw ERA T anomalies as saved in era_ptempanom.mat (see compare_era_anomaly_from_fit_and_model.m)
if exist('iFixTz_NoFit','var')
  %% from strow_override_defaults_latbins_AIRS_fewlays.m
  if iFixTz_NoFit > 0 & strcmp(driver.rateset.ocb_set,'obs')
    if iZeroTVers == 0
      izname = ['SAVE_BESTRUNv1/OutputAnomaly_CAL/' num2str(driver.iibin,'%02d') '/anomtest_timestep' num2str(driver.i16daytimestep) '.mat'];
      if exist(izname)
        fprintf(1,'WARNING setting dT(z,lat,t)/dt using CAL anom in %s \n',izname);       
        izt = load(izname);
        cal_T_rates = izt.oem.finalrates(driver.jacobian.temp_i);

        spectra_due_to_T_jac = zeros(size(driver.rateset.rates));
        renorm_cal_T_rates = cal_T_rates./driver.qrenorm(driver.jacobian.temp_i)';
        for iiT = 1 : length(driver.jacobian.temp_i)
          spectra_due_to_T_jac = spectra_due_to_T_jac + renorm_cal_T_rates(iiT)*m_ts_jac(:,driver.jacobian.temp_i(iiT));
        end
        load f2235.mat
        plot(f2235,spectra_due_to_T_jac,'b',f2235,driver.rateset.rates,'m.-',...
             f2235,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2235,driver.rateset.rates-spectra_due_to_T_jac,'r.-')
        plot(f2235,driver.rateset.rates,'m.-',f2235,driver.rateset.rates-spectra_due_to_T_jac,'r.-',f2235,m_ts_jac(:,1)*10+0.85,'k.-')
      else
        fprintf(1,'WARNING trying to setting dT(z,lat,t)/dt using CAL anom but %s DNE so set xb(T) = 0 \n',izname);
        cal_T_rates = zeros(size(driver.jacobian.temp_i));
        spectra_due_to_T_jac = zeros(size(driver.rateset.rates));
        xb(driver.jacobian.temp_i) = 0.0;
      end

    elseif iZeroTVers == 1
      era_model_file = ['/asl/s1/sergio/home/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Aug20_2019_Clr/'];
      era_model_file = [era_model_file '/Desc/16DayAvgNoS_withanom/latbin' num2str(driver.iibin,'%02d') '_16day_avg.rp.mat'];
      izname = era_model_file;
      x = ['SAVE_BESTRUNv1/OutputAnomaly_CAL/27/anomtest_timestep' num2str(123) '.mat'];
      xx = load(x);

      if exist(izname)
        fprintf(1,'WARNING setting dT(z,lat,t)/dt using raw ERA snom in %s \n',izname);       
        era_anom = load(era_model_file);
        figure(1); pcolor(era_anom.p16anomaly.ptempanomaly); shading flat; colorbar; colormap jet; caxis([-5 5])
        junk = era_anom.p16anomaly.ptempanomaly;
        for iiii = 1 : 20
          ixix = xx.jacobian.wvjaclays_used{iiii}-6;  %% 6 = xx.jacobian.scalar_i below
          if max(ixix) <= 98
            era_ptempanom(iiii,:) = mean(junk(ixix,:));
          end
        end        
        cal_T_rates = era_ptempanom(:,driver.i16daytimestep);

        spectra_due_to_T_jac = zeros(size(driver.rateset.rates));
        renorm_cal_T_rates = cal_T_rates./driver.qrenorm(driver.jacobian.temp_i)';
        for iiT = 1 : length(driver.jacobian.temp_i)
          spectra_due_to_T_jac = spectra_due_to_T_jac + renorm_cal_T_rates(iiT)*m_ts_jac(:,driver.jacobian.temp_i(iiT));
        end
        load f2235.mat
        plot(f2235,spectra_due_to_T_jac,'b',f2235,driver.rateset.rates,'m.-',...
             f2235,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2235,driver.rateset.rates-spectra_due_to_T_jac,'r.-')
        plot(f2235,driver.rateset.rates,'m.-',f2235,driver.rateset.rates-spectra_due_to_T_jac,'r.-',f2235,m_ts_jac(:,1)*10+0.85,'k.-')

      else
        fprintf(1,'WARNING trying to setting dT(z,lat,t)/dt using CAL anom but %s DNE so set xb(T) = 0 \n',izname);
        cal_T_rates = zeros(size(driver.jacobian.temp_i));
        spectra_due_to_T_jac = zeros(size(driver.rateset.rates));
        xb(driver.jacobian.temp_i) = 0.0;
      end

    elseif iZeroTVers == 2
      era_model_file = 'era_ptempanom.mat';
      izname = era_model_file;
      if exist(izname)
        fprintf(1,'WARNING setting dT(z,lat,t)/dt using raw ERA anom in %s \n',izname);       
        era_anom = load(era_model_file);
        junk = era_anom.era_ptempanom;
        era_ptempanom = squeeze(junk(driver.iibin,:,:));
        cal_T_rates = era_ptempanom(:,driver.i16daytimestep);

        spectra_due_to_T_jac = zeros(size(driver.rateset.rates));
        renorm_cal_T_rates = cal_T_rates./driver.qrenorm(driver.jacobian.temp_i)';
        for iiT = 1 : length(driver.jacobian.temp_i)
          spectra_due_to_T_jac = spectra_due_to_T_jac + renorm_cal_T_rates(iiT)*m_ts_jac(:,driver.jacobian.temp_i(iiT));
        end
        load f2235.mat
        plot(f2235,spectra_due_to_T_jac,'b',f2235,driver.rateset.rates,'m.-',...
             f2235,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2235,driver.rateset.rates-spectra_due_to_T_jac,'r.-')
        plot(f2235,driver.rateset.rates,'m.-',f2235,driver.rateset.rates-spectra_due_to_T_jac,'r.-',f2235,m_ts_jac(:,1)*10+0.85,'k.-')

      else
        fprintf(1,'WARNING trying to setting dT(z,lat,t)/dt using CAL anom but %s DNE so set xb(T) = 0 \n',izname);
        cal_T_rates = zeros(size(driver.jacobian.temp_i));
        spectra_due_to_T_jac = zeros(size(driver.rateset.rates));
        xb(driver.jacobian.temp_i) = 0.0;
      end

    end  %% if iZeroTVers == 0 ! if iZeroTVers == 1 | if iZeroTVers == 2
  end    %% if iFixTz_NoFit > 0 & strcmp(driver.rateset.ocb_set,'obs')
end      %% if exist('iFixTz_NoFit','var')

%%%%%%%%%%%%%%%%%%%%%%%%%

if settings.iFixO3_NoFit >= 0 & (strcmp(driver.rateset.ocb_set,'obs') | strcmp(driver.rateset.ocb_set,'cal'))
  iFixO3_NoFit = settings.iFixO3_NoFit;    %%% LARABBEE LIKES THIS TURNED OFF ie keep spectra as is, just read in ERA anom and proceed
end

iZeroO3Vers = 0; %%% use my fit to sarta calcs, as a proxy to ERA T anomalies
iZeroO3Vers = 1; %%% use raw ERA T anomalies, and do the averaging here on the fly
iZeroO3Vers = 2; %%% use raw ERA T anomalies as saved in era_ptempanom.mat (see compare_era_anomaly_from_fit_and_model.m)

if exist('iFixO3_NoFit','var')
  %% from strow_override_defaults_latbins_AIRS_fewlays.m
  if iFixO3_NoFit >= 0 & (strcmp(driver.rateset.ocb_set,'obs') | strcmp(driver.rateset.ocb_set,'cal'))
    if iZeroO3Vers == 0
      izname = ['SAVE_BESTRUNv1/OutputAnomaly_CAL/' num2str(driver.iibin,'%02d') '/anomtest_timestep' num2str(driver.i16daytimestep) '.mat'];
      if exist(izname)
        fprintf(1,'WARNING setting dO3(z,lat,t)/dt using CAL anom in %s \n',izname);       
        izt = load(izname);
        cal_O3_rates = izt.oem.finalrates(driver.jacobian.ozone_i);

        spectra_due_to_O3_jac = zeros(size(driver.rateset.rates));
        renorm_cal_O3_rates = cal_O3_rates./driver.qrenorm(driver.jacobian.ozone_i)';
        for iiO3 = 1 : length(driver.jacobian.ozone_i)
          spectra_due_to_O3_jac = spectra_due_to_O3_jac + renorm_cal_O3_rates(iiT)*m_ts_jac(:,driver.jacobian.ozone_i(iiO3));
        end
        load f2235.mat
        plot(f2235,spectra_due_to_O3_jac,'b',f2235,driver.rateset.rates,'m.-',...
             f2235,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2235,driver.rateset.rates-spectra_due_to_O3_jac,'r.-')
        plot(f2235,driver.rateset.rates,'m.-',f2235,driver.rateset.rates-spectra_due_to_O3_jac,'r.-',f2235,m_ts_jac(:,1)*10+0.85,'k.-')
      else
        fprintf(1,'WARNING trying to setting dO3(z,lat,t)/dt using CAL anom but %s DNE so set xb(O3) = 0 \n',izname);
        cal_O3_rates = zeros(size(driver.jacobian.ozone_i));
        spectra_due_to_O3_jac = zeros(size(driver.rateset.rates));
        xb(driver.jacobian.ozone_i) = 0.0;
      end

    elseif iZeroO3Vers == 1
      era_model_file = ['/asl/s1/sergio/home/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Aug20_2019_Clr/'];
      era_model_file = [era_model_file '/Desc/16DayAvgNoS_withanom/latbin' num2str(driver.iibin,'%02d') '_16day_avg.rp.mat'];
      izname = era_model_file;
      x = ['SAVE_BESTRUNv1/OutputAnomaly_CAL/27/anomtest_timestep' num2str(123) '.mat'];
      xx = load(x);

      if exist(izname)
        fprintf(1,'WARNING setting dO3(z,lat,t)/dt using raw ERA snom in %s \n',izname);       
        era_anom = load(era_model_file);
        figure(1); pcolor(era_anom.p16anomaly.gas_3anomaly); shading flat; colorbar; colormap jet; caxis([-5 5])
        junk = era_anom.p16anomaly.ozoneanomaly;
        for iiii = 1 : 20
          ixix = xx.jacobian.wvjaclays_used{iiii}-6;  %% 6 = xx.jacobian.scalar_i below
          if max(ixix) <= 98
            era_ozoneanom(iiii,:) = mean(junk(ixix,:));
          end
        end        
        cal_ozone_rates = era_ozoneanom(:,driver.i16daytimestep);

        spectra_due_to_ozone_jac = zeros(size(driver.rateset.rates));
        renorm_cal_O3_rates = cal_O3_rates./driver.qrenorm(driver.jacobian.ozone_i)';
        for iiO3 = 1 : length(driver.jacobian.ozone_i)
          spectra_due_to_O3_jac = spectra_due_to_O3_jac + renorm_cal_O3_rates(iiT)*m_ts_jac(:,driver.jacobian.ozone_i(iiO3));
        end
        load f2235.mat
        plot(f2235,spectra_due_to_O3_jac,'b',f2235,driver.rateset.rates,'m.-',...
             f2235,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2235,driver.rateset.rates-spectra_due_to_O3_jac,'r.-')
        plot(f2235,driver.rateset.rates,'m.-',f2235,driver.rateset.rates-spectra_due_to_O3_jac,'r.-',f2235,m_ts_jac(:,1)*10+0.85,'k.-')

      else
        fprintf(1,'WARNING trying to setting dO3(z,lat,t)/dt using CAL anom but %s DNE so set xb(O3) = 0 \n',izname);
        cal_O3_rates = zeros(size(driver.jacobian.ozone_i));
        spectra_due_to_O3_jac = zeros(size(driver.rateset.rates));
        xb(driver.jacobian.ozone_i) = 0.0;
      end

    elseif iZeroO3Vers == 2
      era_model_file = 'era_gas_3anom.mat';
      izname = era_model_file;
      if exist(izname)
        fprintf(1,'WARNING setting dO3(z,lat,t)/dt using raw ERA anom in %s \n',izname);       
        era_anom = load(era_model_file);
        junk = era_anom.era_gas_3anom;
        era_ozoneanom = squeeze(junk(driver.iibin,:,:));
        cal_O3_rates = era_ozoneanom(:,driver.i16daytimestep);

        spectra_due_to_O3_jac = zeros(size(driver.rateset.rates));
        renorm_cal_O3_rates = cal_O3_rates./driver.qrenorm(driver.jacobian.ozone_i)';
        for iiT = 1 : length(driver.jacobian.ozone_i)
          spectra_due_to_O3_jac = spectra_due_to_O3_jac + renorm_cal_O3_rates(iiT)*m_ts_jac(:,driver.jacobian.ozone_i(iiT));
        end
        load f2235.mat
        plot(f2235,spectra_due_to_O3_jac,'b',f2235,driver.rateset.rates,'m.-',...
             f2235,sum(m_ts_jac(:,driver.jacobian.ozone_i)')*100,'k.-',f2235,driver.rateset.rates-spectra_due_to_O3_jac,'r.-')
        plot(f2235,driver.rateset.rates,'m.-',f2235,driver.rateset.rates-spectra_due_to_O3_jac,'r.-',f2235,m_ts_jac(:,1)*10+0.85,'k.-')

      else
        fprintf(1,'WARNING trying to setting dT(z,lat,t)/dt using CAL anom but %s DNE so set xb(T) = 0 \n',izname);
        cal_O3_rates = zeros(size(driver.jacobian.ozone_i));
        spectra_due_to_O3_jac = zeros(size(driver.rateset.rates));
        xb(driver.jacobian.ozone_i) = 0.0;
      end

    end  %% if iZeroO3Vers == 0 ! if iZeroO3Vers == 1 | if iZeroO3Vers == 2
    if iFixO3_NoFit == 0
      disp('hmm in SW no need to worry about O3, zero all this O3 stuff ...')
      spectra_due_to_O3_jac = zeros(size(spectra_due_to_O3_jac));
      cal_O3_rates = zeros(size(cal_O3_rates));
    end 
  end    %% if iFixO3_NoFit > 0 & strcmp(driver.rateset.ocb_set,'obs')
end      %% if exist('iFixO3_NoFit','var')

%%%%%%%%%%%%%%%%%%%%%%%%%

%% driver.settings = settings;   %% no need to do this since topts == settings == topts

if settings.set_tracegas == +1 & driver.i16daytimestep < 0 & settings.ocb_set == 0
  disp('setting constant rates for tracegas apriori : CO2 = 2.2  CH4 = 4.5 N2O = 0.8 CFC = -1.25')
  if settings.co2lays == 1
    xb(1) = 2.2;  % Set CO2 apriori
    xb(2) = 1;
    xb(3) = 4.5;
    xb(4) = -1.0;
    xb(5) = -1.0;

    xb(1) = 2.2 * 1;    % Set CO2 apriori
    xb(2) = 0.8 * 1;    % set N2O 
    xb(3) = 4.5 * 1;    % set CH4
    xb(4) = -1.25 * 0;  % set CFC11, before Aug 23 the mult was 1
    xb(5) = -1.25 * 0;  % set CFC12, before Aug 23 the mult was 1

  elseif settings.co2lays == 3
    xb(1) = 2.2 * 1;        % Set CO2 apriori lower trop
    xb(2) = 2.2 * 1;        % Set CO2 apriori mid trop
    xb(3) = 2.2 * 1;        % Set CO2 apriori strat

    xb(4) = 0.8 * 1;        % set N2O 
    xb(5) = 4.5 * 1;        % set CH4
    xb(6) = -1.25 * 0;  % set CFC11, before Aug 23 the mult was 1
    xb(7) = -1.25 * 0;  % set CFC12, before Aug 23 the mult was 1
  end

elseif settings.set_tracegas == +1 & driver.i16daytimestep > 0 & settings.ocb_set == 0
  junk = 365/16; %% days per timestep 
  junk = (driver.i16daytimestep-1)/junk;
  str = ['setting time varying rates for tracegas apriori : CO2 = 2.2  CH4 = 4.5 N2O = 0.8 CFC = -1.25 for ' num2str(junk) ' years']; 
  disp(str);
  if settings.co2lays == 1
    deltaT = 365/16; %% days per timestep

    %% default all this while getting good results
    xb(1) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori
    xb(2) = 0.8 * (driver.i16daytimestep-1)/deltaT * 1.0;    % set N2O 
    xb(3) = 4.5 * (driver.i16daytimestep-1)/deltaT * 1.0;    % set CH4
    xb(4) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;    % set CFC11, before Aug 23 the mult was 1
    xb(5) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;    % set CFC12, before Aug 23 the mult was 1
    xb(4) = -1.25 * (driver.i16daytimestep-1)/deltaT * 1;    % set CFC11, before May 2022  the mult was 1
    xb(5) = -1.25 * (driver.i16daytimestep-1)/deltaT * 1;    % set CFC12, before May 2022 the mult was 1

  elseif settings.co2lays == 3
    deltaT = 365/16; %% days per timestep

    xb(1) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori lower trop
    xb(2) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori mid trop
    xb(3) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori strat

    xb(4) = 0.8 * (driver.i16daytimestep-1)/deltaT * 1.0;        % set N2O 
    xb(5) = 4.5 * (driver.i16daytimestep-1)/deltaT * 1.0;        % set CH4
    xb(6) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;  % set CFC11, before Aug 23 the mult was 1
    xb(7) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;  % set CFC12, before Aug 23 the mult was 1
    xb(6) = -1.25 * (driver.i16daytimestep-1)/deltaT * 1;  % set CFC11, before May 2022 the mult was 1
    xb(7) = -1.25 * (driver.i16daytimestep-1)/deltaT * 1;  % set CFC12, before May 2022 the mult was 1
  end

elseif settings.set_tracegas == +2 & driver.i16daytimestep > 1 & settings.ocb_set == 0
  junk = 365/16; %% days per timestep 
  junk = (driver.i16daytimestep-1);
  str = ['setting bootstrap time varying rates for tracegas apriori : CO2 = 2.2  CH4 = 4.5 N2O = 0.8 CFC = -1.25 based on previous timestep'];
  disp(str);
  fminus = ['OutputAnomaly_OBS/' num2str(driver.iibin,'%02d') '/anomtest_timestep' num2str(driver.i16daytimestep-1) '.mat'];
  fprintf(1,'looking for %s to fill in co2/n2o/ch4/cfc11/cfc12 rates ... \n',fminus)
  if exist(fminus)
    prev = load(fminus);
    if settings.co2lays == 1
      deltaT = 365/16; %% days per timestep

      %% default all this while getting good results
      xb0(1) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori
      xb0(2) = 0.8 * (driver.i16daytimestep-1)/deltaT * 1.0;    % set N2O 
      xb0(3) = 4.5 * (driver.i16daytimestep-1)/deltaT * 1.0;    % set CH4
      xb0(4) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;    % set CFC11, before Aug 23 the mult was 1
      xb0(5) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;    % set CFC12, before Aug 23 the mult was 1
      xb0(4) = -1.25 * (driver.i16daytimestep-1)/deltaT * 1;    % set CFC11, before May 2022  the mult was 1
      xb0(5) = -1.25 * (driver.i16daytimestep-1)/deltaT * 1;    % set CFC12, before May 2022 the mult was 1

      %% overwrite!!!
      xb(1) = prev.oem.finalrates(1);
      xb(2) = prev.oem.finalrates(2);
      xb(3) = prev.oem.finalrates(3);
      xb(4) = prev.oem.finalrates(4);
      xb(5) = prev.oem.finalrates(5);

      fprintf(1,'changed CO2   apriori from %8.6f to %8.6f ppv \n',xb0(1),xb(1))
      fprintf(1,'changed N2O   apriori from %8.6f to %8.6f ppb \n',xb0(2),xb(2))
      fprintf(1,'changed CH4   apriori from %8.6f to %8.6f ppb \n',xb0(3),xb(3))
      fprintf(1,'changed CRF11 apriori from %8.6f to %8.6f ppt \n',xb0(4),xb(4))
      fprintf(1,'changed CFC12 apriori from %8.6f to %8.6f ppt \n',xb0(5),xb(5))

    else
      error('not doing this')
    end
  else
    error('oops previous file DNE');
  end
end

[mm,nn] = size(xb);
if nn > 1
  xb = xb(:,driver.iibin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('iFixTz_NoFit','var') & strcmp(driver.rateset.ocb_set,'obs')

  disp(' >>>>> need to account for cal_T_rates by eg subbing in T(z,t) from fits to cal or raw ERA anoms >>>>> ')
  xb(driver.jacobian.temp_i) = cal_T_rates;

  aux.FixTz_NoFit = cal_T_rates;
  aux.orig_water_i = driver.jacobian.water_i;
  aux.orig_temp_i = driver.jacobian.temp_i;
  aux.orig_ozone_i = driver.jacobian.ozone_i;

  noT = setdiff(1:length(xb),driver.jacobian.temp_i);

  driver.jacobian.ozone_i =  driver.jacobian.temp_i;
  driver.jacobian.temp_i = [];

  m_ts_jac = m_ts_jac(:,noT);
  aux.m_ts_jac    = m_ts_jac;

  xb = xb(noT);
  driver.rateset.rates = driver.rateset.rates - spectra_due_to_T_jac;   %% LARRABEE DOES NOT WANT THIS
  driver.qrenorm = driver.qrenorm(noT);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('iFixO3_NoFit','var') & (strcmp(driver.rateset.ocb_set,'obs') | strcmp(driver.rateset.ocb_set,'cal'))

  disp(' >>>>> need to account for cal_O3_rates by eg subbing in O3(z,t) from fits to cal or raw ERA anoms >>>>> ')
  xb(driver.jacobian.ozone_i) = cal_O3_rates;

  aux.FixO3_NoFit = cal_O3_rates;
  aux.orig_water_i = driver.jacobian.water_i;
  aux.orig_temp_i = driver.jacobian.temp_i;
  aux.orig_ozone_i = driver.jacobian.ozone_i;

  noO3 = setdiff(1:length(xb),driver.jacobian.ozone_i);

  driver.jacobian.ozone_i =  [];

  m_ts_jac = m_ts_jac(:,noO3);
  aux.m_ts_jac    = m_ts_jac;

  xb = xb(noO3);
  driver.rateset.rates = driver.rateset.rates - spectra_due_to_O3_jac;   %% LARRABEE DOES NOT WANT THIS
  driver.qrenorm = driver.qrenorm(noO3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% A Priori stored in aux.xb
aux.xb = xb./driver.qrenorm';

%---------------------------------------------------------------------------
% SARTA forward model and other "representation" errors
driver.oem.sarta_error = 0.0;
driver.oem.xb = xb;  %% note this is un-normalized xb

%---------------------------------------------------------------------------
if abs(settings.offsetrates) ~= 1
  settings.offsetrates
  error('incorrect setting for overriding obs rates by constant');
elseif settings.offsetrates > 0 & driver.i16daytimestep < 0
   disp('Offset linear rates (obs or cal or bias) by constant to get sensitivity to AIRS BT drift')
   driver.rateset.rates = driver.rateset.rates + 0.01;
elseif settings.offsetrates > 0 & driver.i16daytimestep > 0 & settings.ocb_set == 0 & nansum(abs(driver.rateset.rates)) > 0
   disp('Offset anomaly (obs) by constant 0.01/year to get sensitivity to AIRS BT drift')
   oktimes = load('ok157times.mat');
   oktimes = oktimes.okdates(driver.i16daytimestep);   %% need to get the CO2 jac at this time!!!!
   fprintf(1,'  --> oktimes = %8.6f which is %8.6f years away from 2012.356 \n',oktimes,oktimes - 2012.3562)
   driver.rateset.rates = driver.rateset.rates + 0.01 * (oktimes - 2012.3562);
elseif settings.offsetrates > 0 & driver.i16daytimestep > 0 & settings.ocb_set == 1 & nansum(abs(driver.rateset.rates)) > 0
   disp('Offset anomaly (cal) by constant 0.01/year to get sensitivity to AIRS BT drift')
   oktimes = load('ok157times.mat');
   oktimes = oktimes.okdates(driver.i16daytimestep);   %% need to get the CO2 jac at this time!!!!
   fprintf(1,'  --> oktimes = %8.6f which is %8.6f years away from 2012.3562 \n',oktimes,oktimes - 2012.3562)
   disp('BUT LARRABEE SAID NOT TO DO THIS SO DONT ADD THIS OFFEST')
   %%%% driver.rateset.rates = driver.rateset.rates + 0.01 * (oktimes - 2012.3562);
end

%------------------------------------------------------------------------
if abs(settings.addco2jacs) ~= 1
  settings.addco2jacs
  error('incorrect setting for adding co2jacs to ERA calcrates');
end
if strcmp(driver.rateset.ocb_set,'cal') & settings.addco2jacs > 0
  disp('Offset rates to add in CO2 jacs to ERA rates, to pretend there is CO2')
  jac = load('co2_kcartajac.mat');
  %  haha = driver.rateset.rates;
  %  baba = jac.co2jac(ix,:)';
  %  whos haha baba
  driver.rateset.rates = driver.rateset.rates + jac.co2jac(ix,:)';
end

%---------------------------------------------------------------------------
% Modify rates with lag-1 correlation errors or add to above
if driver.i16daytimestep < 0
  nc_cor = nc_rates(driver);
  driver.rateset.unc_rates = driver.rateset.unc_rates.*nc_cor;   %% THIS IS AN ARRAY
end

finitech = find(isfinite(driver.rateset.rates));
finitechnoise = find(isfinite(driver.rateset.unc_rates));

%driver.rateset.unc_rates = 0.001 * ones(size(driver.rateset.unc_rates));   %% THIS IS AN ARRAY

%%%%%%%%%%%%%%%%%%%%%%%%% >>>>>>
%{
OLD!!!!
if settings.obs_corr_matrix > 0
  addpath /home/sergio/MATLABCODE
  thecov = load('/home/sergio/MATLABCODE/oem_pkg_run/Simulate_Calcs/thecov_clear');
  [f2235,i2235] = map_2834_to_2235;
  junk = thecov.thecov; junk = junk(i2235,i2235);

  %% see https://www.mathworks.com/help/finance/corr2cov.html
  %%junk = corr2cov(driver.rateset.unc_rates,junk);
  %% new, fast
  junk0 = junk;
  junk0 = diag(driver.rateset.unc_rates) * junk0 * diag(driver.rateset.unc_rates);

  junk = junk0;  
  fprintf(1,'numchans, rank, mean rank, condition number of obs cov matrix = %6i %6i %8.6e %8.6e \n',nn,rank(junk),rank(junk)/length(junk),cond(junk))

  driver.rateset.unc_rates = junk;  %% THIS IS A SQUARE MATRTIX
end
%}

%% NEW
iCovSe_OffDiag = +1;  %% send in off-diag matrix strow
iCovSe_OffDiag = 0;   %% send in diag matrix sergio which should be equivalent to iCovSe_OffDiag = -1 ???????
iCovSe_OffDiag = -1;  %% do usual in rodgers.m          <<< DEFAULT >>>

iCovSe_OffDiag = -1;  %% do usual in rodgers.m          <<< DEFAULT >>>
iCovSe_OffDiag = 0;   %% send in diag matrix sergio which should be equivalent to iCovSe_OffDiag = -1 ???????!
iCovSe_OffDiag = +1;  %% send in off-diag matrix strow
iCovSe_OffDiag = +2;  %% send in off-diag matrix serio
iCovSe_OffDiag = topts.obs_corr_matrix;

fprintf(1,'strow_override_defaults_latbins_CRIS_fewlays.m : iCovSe_OffDiag = %2i \n',iCovSe_OffDiag)
if iCovSe_OffDiag >= 0
  %% now turn driver.rateset.unc_rates into cov matrix
  %% var(i) = sig(i) sig(i)   where sig(i) = std dev (i)
  %% rij = cov(ij)/sqrt(var(i) var(j))
  %% so cov(ij) = rij sqrt(var(i) var(j))  = rij sqrt(sig(i) sig(i)  sig(j) sig(j)) = rij sig(i) sig(j) 
  %%    cov(ii) = 1.0 sqrt(var(i) var(i))  = var(i) = sig(i) sig(i)
  
  junk = driver.rateset.unc_rates;  %% sig(i)
  junk0 = junk' * junk;             %% sig(i) sig(j)
  junknew = diag(diag(junk0));      %% sig(i) sig(i) on diag, everything else zeros == "AIRS" normal cov matrix

  junkP1 = diag(junk0,+1);
  junkM1 = diag(junk0,-1);
  
  if iCovSe_OffDiag == 2
    % Characterization of the Observational Covariance Matrix of Hyper-Spectral Infrared Satellite Sensors Directly from Measured Earth Views
    % Carmine Serio 1 , Guido Masiello 1,∗ , Pietro Mastro 1 and David C. Tobin 2
    junknew = junknew + -0.3*(diag(junkP1,1) + diag(junkM1,-1)); %% see Fig 11
  elseif iCovSe_OffDiag == 1
    %% strow slide https://www.star.nesdis.noaa.gov/star/documents/meetings/2016JPSSAnnual/S4/S4_04_cris_strow.pdf
    junknew = junknew + 0.02*(diag(junkP1,1) + diag(junkM1,-1)); %% see Slide 7
  end

  driver.rateset.unc_rates = junknew;
end

%%%%%%%%%%%%%%%%%%%%%%%%% >>>>>>

% Modify with estimated error in freq + regress errors 
%driver.rateset.unc_rates = ones(2378,1)*0.001 +driver.rateset.unc_rates.*nc_cor;
%load corr_errors
%load 15yr_corr_errors_2314chans

%---------------------------------------------------------------------------
% Do rate Q/A (empty for now)
%---------------------------------------------------------------------------

if settings.resetnorm2one == +1
  %% we already cleared jac
  %aux

  disp('resetting all jabian norms to 1 ONE UNO UN MOJA')
  boo1 = driver.qrenorm;
  boo2 = aux.xb; 
  %whos boo1 boo2
  %[driver.qrenorm'  aux.xb]

  for ii = 1 : length(boo1)
    aux.m_ts_jac(:,ii) = aux.m_ts_jac(:,ii)/driver.qrenorm(ii);
  end;
  aux.xb = aux.xb./driver.qrenorm';
  driver.qrenorm = ones(size(driver.qrenorm));
end

build_cov_matrices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%