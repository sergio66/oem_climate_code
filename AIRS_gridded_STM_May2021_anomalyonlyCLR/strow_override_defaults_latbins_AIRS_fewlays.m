function [driver,aux] = strow_override_defaults_latbins_AIRS_fewlays(driver,iNlays_retrieve,topts);

settings.resetnorm2one = -1; %%% default, keep my scaling, set to +1 if you want to reset all to 1.00000000

settings.dataset = 1;       % (1) AIRS 16 year dataset (2) AIRS 11 year (IASI) (3) IASI2AIRS 11 year
settings.co2lays = 1;       % assume column jac for CO2, or 3 lays (gnd-500,500-trop,trop-TOA

settings.ocb_set = 0;       % 0 = obs, 1 = cal, -1 = bias
settings.numchan = 2645;    % L1b = 2378, L1c = 2645
settings.chan_LW_SW = 0;    % 0 is 640 to 1640, 1 is 700 to 1640, -1 is 640 to 2740

settings.set_tracegas = -1;            %% do we leave apriori as 0 or set CO2/N2o/CH4/CFC to be 2.2, 4.5, -1             
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
settings.UMBCvsERAjac = -1;            %% do not adjust jacobian, based on handful of clear sky retrieval days
                                       %% [+1] means to adjust

allowedparams = [{'ocb_set'},{'numchan'},{'chan_LW_SW'},{'iChSet'},{'set_tracegas'},{'offsetrates'},...
			    {'addco2jacs'},{'obs_corr_matrix'},{'invtype'},{'tie_sst_lowestlayer'},{'iNlays_retrieve'},...
                            {'descORasc'},{'dataset'},{'iXJac'},{'co2lays'},{'iDoStrowFiniteJac'},{'iFixTz_NoFit'},{'iFixO3_NoFit'},{'iFixTG_NoFit'},{'resetnorm2one'},{'UMBCvsERAjac'}];

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
% Which tile
ix = driver.iTile;
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
  error('oops not done')

elseif settings.dataset == 1   %% DEFAULT
  disp('AIRS 16 year rates or anomalies, nu cal done in there')
  if settings.descORasc == +1 & driver.i16daytimestep < 0
    disp('doing gridded rates')
    driver.rateset.datafile  = [];
    if settings.ocb_set == 0  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsonly_transpose.mat';  %% ah puts strides, as hoped/expected ie WRONG
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsonly.mat';           
      driver.rateset.datafile  = 'convert_sergio_clearskygrid_obsonly.mat';           
      driver.rateset.datafile  = 'oops look at ../AIRS_gridded_May2021_trendsonlyCLR/';
      error('look at ../AIRS_gridded_May2021_trendsonlyCLR/')
    elseif settings.ocb_set == +1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
      driver.rateset.datafile  = 'convert_strowrates2oemrates_allskygrid_obsNcalcs.mat';                 
      driver.rateset.datafile  = 'oops look at ../AIRS_gridded_May2021_trendsonlyCLR/';
      error('look at ../AIRS_gridded_May2021_trendsonlyCLR/')
    elseif settings.ocb_set == -1  & driver.i16daytimestep < 0
      driver.rateset.datafile  = 'AHAH';
    end
  elseif settings.descORasc == -1 & driver.i16daytimestep < 0
    disp('doing ascending latbin rates')
    driver.rateset.datafile  = [];
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 0
    disp('doing descending OBS ANOMALY')
    driver.rateset.datafile = [];
    driver.rateset.datafile = ['ANOM/LatBin' num2str(driver.iLat,'%02d') '/LonBin' num2str(driver.iLon,'%02d') '/v1_Quantile' num2str(driver.iQuantile,'%02d') '/'];
    driver.rateset.datafile = [driver.rateset.datafile 'desc_anaomaly_LatBin' num2str(driver.iLat,'%02d') '_LonBin' num2str(driver.iLon,'%02d') '_Quantile' num2str(driver.iQuantile,'%02d') '_timesetps_001_412_V1.mat'];
  elseif driver.i16daytimestep > 0 & settings.ocb_set == 1
    disp('doing descending CAL ANOMALY')
    driver.rateset.datafile = [];
  end
end

fprintf(1,' <<< driver.rateset.datafile >>> = %s \n',driver.rateset.datafile)

% Lag-1 correlation file; if using rate least-squares errors
driver.rateset.ncfile   = '../oem_pkg/Test/all_lagcor.mat';
driver.rateset.ncfile   = driver.rateset.datafile;

driver
% Get rate data, do Q/A elsewhere
driver = get_rates(driver);  %% this gets spectral rates (driver.rateset.rates), and uncertainty (driver.rateset.unc_rates)

%---------------------------------------------------------------------------
% Jacobian file: f = 2378x1 and M_TS_jac_all = 36x2378x200

iXJac = settings.iXJac;
%if driver.i16daytimestep > 0
%  iXJac = 0; %% const geo kcarta jcs, default for trends
%  iXJac = 1; %% varying geo sarta jacs
%  iXJac = 2; %% varying geo kcarta jacs, default for anomaly
%end

if driver.i16daytimestep < 0
  iMidPoint_TimeStepUse = 1;            %% put in constant Jacobian, at timestep 1 (2002/09)
  iMidPoint_TimeStepUse = 388;          %% put in constant Jacobian, at timestep 365 (2018/08)
  iMidPoint_TimeStepUse = floor(388/2); %% put in constant Jacobian, half way through (365/2 ==> 2009/09)

  iMidPoint_TimeStepUse = 270; 
  iMidPoint_TimeStepUse = 335; 
  iMidPoint_TimeStepUse = 194; %% default "fool" the code by using midpoint anomaly jac

  if settings.descORasc == -1
    driver.jacobian.filename = [AHAJAC];
    fprintf(1,'reading in constant kcarta jac file %s \n',driver.jacobian.filename)

  elseif settings.descORasc == +1
    %% we can "fool" the code by using midpoint anomaly jac
    %junk = num2str(iMidPoint_TimeStepUse,'%03d');
    %driver.jacobian.filename = ['AHA']; 

    %% for now assume same jacs
    AHA = '/asl/s1/sergio/rtp/MakeAvgProfs2002_2020/Retrieval/LatBin65/SubsetJacLatbin/';
    AHA = '/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/';
    %% figure out which latbin
    %% latbin 1 has jacs/points for indices 1:72 + (1-1)*72
    %% latbin 2 has jacs/points for indices 1:72 + (2-1)*72
    driver.jac_latbin         = floor((driver.iibin-1)/72)+1;
    driver.jac_indexINSIDEbin = driver.iibin - (driver.jac_latbin-1)*72;     %% so this should be lonbin

    iKCARTAorSARTA = +1;
    if iKCARTAorSARTA < 0
      %AHA = [AHA '/subjacLatBin' num2str(driver.jac_latbin,'%02i') '.mat'];
      AHA = [AHA '/clr_subjacLatBin' num2str(driver.jac_latbin,'%02i') '.mat'];
    else
      %AHA = [AHA '/kcarta_subjacLatBin' num2str(driver.jac_latbin,'%02i') '.mat'];
      AHA = [AHA '/kcarta_clr_subjacLatBin_newSARTA_' num2str(driver.jac_latbin,'%02i') '.mat'];
    end

    driver.jacobian.filename = AHA;
    clear AHA
    fprintf(1,'reading in constant kcarta jac file %s \n',driver.jacobian.filename)
  end

elseif driver.i16daytimestep > 0
  junk = num2str(driver.i16daytimestep,'%03d');
  if iXJac == 1
    %% sarta time vary jacs
    driver.jacobian.filename = [];
    fprintf(1,'iXJac == 1 reading in timestep sarta jac file %s \n',driver.jacobian.filename)

  elseif iXJac == 2  
    %% kcarta time vary jac
    driver.jacobian.filename = [];
    fprintf(1,'iXJac == 2 reading in timestep kcarta jac file %s \n',driver.jacobian.filename)

  elseif iXJac == 0
    %% constant kcarta jacs
    driver.jacobian.filename = [];
    fprintf(1,'iXJac == 0 reading in constant kcarta jac file %s \n',driver.jacobian.filename)
  end
end

% Get jacobians, and combine the 97 layer T(z)/WV(z)/O3(z) into N layers
driver
%[m_ts_jac0,nlays,qrenorm]  = get_jac(driver.jacobian.filename,driver.jac_indexINSIDEbin);
[m_ts_jac0,nlays,qrenorm,freq2645]  = get_jac_fast(driver.jacobian.filename,driver.iTile,driver.iLon,driver.iLat);
m_ts_jac0 = double(m_ts_jac0);

%% THIS IS DEFAULT -- 4 column trace gas (CO2/N2O/CH4/Cld1/CLd2), 1 stemp, (97x3) geo
driver.jacobian.varname  = 'M_TS_jac_all';
driver.jacobian.scalar_i = 1:6;  %% ST, CO2/N2O/CH4, Cld1, Cld2
driver.jacobian.water_i  = (1:nlays)+6+(nlays-1)*0;
driver.jacobian.temp_i   = (1:nlays)+6+(nlays-1)*1+1;
driver.jacobian.ozone_i  = (1:nlays)+6+(nlays-1)*2+2;
driver.jacobian.numlays  = nlays;

%driver.jacobian.fat = m_ts_jac0;
driver.jacobian.scalar_values = m_ts_jac0(:,driver.jacobian.scalar_i);
driver.jacobian.colWV_values   = nansum(m_ts_jac0(:,driver.jacobian.water_i),2);
driver.jacobian.colTz_values   = nansum(m_ts_jac0(:,driver.jacobian.temp_i), 2);
driver.jacobian.colO3_values   = nansum(m_ts_jac0(:,driver.jacobian.ozone_i),2);
figure(7); imagesc(m_ts_jac0(:,driver.jacobian.water_i)'); colorbar
figure(8); imagesc(m_ts_jac0(:,driver.jacobian.temp_i)'); colorbar
figure(9); imagesc(m_ts_jac0(:,driver.jacobian.ozone_i)'); colorbar
size(m_ts_jac0)
%[1 6 min(driver.jacobian.water_i) max(driver.jacobian.water_i) min(driver.jacobian.temp_i) max(driver.jacobian.temp_i) min(driver.jacobian.ozone_i) max(driver.jacobian.ozone_i)]
%disp('here 1'); pause

m_ts_jac_coljac = m_ts_jac0(:,driver.jacobian.scalar_i);

driver.qrenorm  = qrenorm;       %% set this default
jac.qrenorm     = qrenorm;
junk = load('h2645structure.mat');
jac.f           = junk.h.vchan;

if iNlays_retrieve <= 60
  [m_ts_jac_wv,qWV,layWV]  = combinejaclays(m_ts_jac0,driver.jacobian.water_i,qrenorm,iNlays_retrieve);
  [m_ts_jac_t,qT,layT]     = combinejaclays(m_ts_jac0,driver.jacobian.temp_i, qrenorm,iNlays_retrieve);
  [m_ts_jac_o3,qO3,layO3]  = combinejaclays(m_ts_jac0,driver.jacobian.ozone_i,qrenorm,iNlays_retrieve);
else
  fprintf(1,'setting iNlays_retrieve ( > 60) from %2i to 97 \n',iNlays_retrieve);
  m_ts_jac_wv = m_ts_jac0;
  iNlays_retrieve = 97;
end
%figure(10); imagesc(m_ts_jac_wv'); colorbar
%figure(11); imagesc(m_ts_jac_t'); colorbar
%figure(12); imagesc(m_ts_jac_o3'); colorbar
%disp('here 2'); pause

%% replace CO2,N2O,CH4 jacs
if driver.i16daytimestep > 0  
  %{
  %% this is for 388 anomaly time steps
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
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,1);
  elseif iXJac == 2 & iDoStrowFiniteJac == 2
    fprintf(1,'updating time varying kCARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs -1,+2 for testing 6/24-27/2019...\n');
    fprintf(1,'  note before July2, the tracegas profile (CO2/N2O/CH4) was US Std shoehorned and multiplied, so ppm was quite wonky except at 500 mb');
    fprintf(1,'else turned off \n')
    %% const kCARTA jacs, update the trace gases
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,2); %% only used this on 6/26
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,2); %% added this on 6/27
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,2); %% added this on 6/27
  elseif iXJac == 2 & iDoStrowFiniteJac == 3
    fprintf(1,'updating time varying kCARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 6/24-27/2019...\n');
    %% kcarta strow finite jacs
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,3); %% only used this on 6/26
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,3); %% added this on 6/27
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,3); %% added this on 6/27
  elseif iXJac == 1 & iDoStrowFiniteJac == 3
    fprintf(1,'updating time varying SARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 6/24-27/2019...\n');
    %% sarta strow finite jacs
    iVarType = -3; %% this uses SARTA  finitediff jacs, which I have shown are bad?? or good??
    iVarType = +3; %% this uses kCARTA finitediff jacs, which I have shown are good, just want to test
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,iVarType); %% added this on 8/3
  elseif iXJac == 2 & iDoStrowFiniteJac == 4
    fprintf(1,'updating time varying kCARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 12/08/2019 ... age of air\n');
    %% const kCARTA jacs, update the trace gases
    %% kcarta strow finite jacs
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,4); %% only used this on 6/26
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,4); %% added this on 6/27
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,4); %% added this on 6/27
  elseif iXJac == 1 & iDoStrowFiniteJac == 3
    fprintf(1,'updating time varying SARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 12/08/2019...\n');
    fprintf(1,'else turned off \n')
    %% sarta strow finite jacs
    iVarType = -3; %% this uses SARTA  finitediff jacs, which I have shown are bad?? or good??
    iVarType = +3; %% this uses kCARTA finitediff jacs, which I have shown are good, just want to test
    iVarType = +4; %% this uses kCARTA finitediff jacs, which I have shown are good, just want to test
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep,iVarType); %% added this on 8/3
  end
  %}
elseif driver.i16daytimestep < 0
  %% this is for 1 average rate
  if iXJac == 0
    fprintf(1,'not updating CO2/N2O/CH4 jacs for TRENDS ... keeping same jacs, give better results\n');
  elseif iXJac == 1 | iXJac == 2
    fprintf(1,'updating CO2/N2O/CH4 jacs for TRENDS ...\n');
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iTile,iMidPoint_TimeStepUse,3);
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iTile,iMidPoint_TimeStepUse,3);
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iTile,iMidPoint_TimeStepUse,3);
  end
end

if settings.co2lays == 3
  m_ts_jac_coljac = replace_time_co2_3layjac(m_ts_jac_coljac,driver.iTile,driver.i16daytimestep);
end

if iNlays_retrieve <= 60
  m_ts_jac = [m_ts_jac_coljac m_ts_jac_wv m_ts_jac_t m_ts_jac_o3];
else
  m_ts_jac = m_ts_jac0;
  qWV = (1:iNlays_retrieve) + 6;
end

%whos m_ts_jac*
%figure(13); colormap jet; imagesc(log10(abs(m_ts_jac'))); caxis([-12 0]); colorbar; 
%disp('here 3b'); pause

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

if settings.UMBCvsERAjac == 1 & settings.iDoStrowFiniteJac > 0
  disp('adjusting first 5 (tracegas ONLY) jacobian based on clear sky retrievals')
  addpath /home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/NewJacs_forAMT2020
  m_ts_jac = add_deltaJ_to_jacs(m_ts_jac,driver.iTile,-1);
end

aux.m_ts_jac = m_ts_jac;
aux.f        = jac.f;

f = jac.f;
clear jac
%---------------------------------------------------------------------------
% Good channel set
lenrates = length(driver.rateset.rates);

iChSet = 2; %% new chans
iChSet = 1; %% old chans (default)
iChSet = 3; %% new chans, but no CFC11
iChSet = topts.iChSet;

ch = find_the_oem_channels(f,lenrates,settings.numchan,settings.chan_LW_SW,iChSet);

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
  disp('oh oh you wanna get rid of a trace gas')
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
      izname = ['SAVE_BESTRUNv1/OutputAnomaly_CAL/' num2str(driver.iTile,'%02d') '/anomtest_timestep' num2str(driver.i16daytimestep) '.mat'];
      if exist(izname)
        fprintf(1,'WARNING setting dT(z,lat,t)/dt using CAL anom in %s \n',izname);       
        izt = load(izname);
        cal_T_rates = izt.oem.finalrates(driver.jacobian.temp_i);

        spectra_due_to_T_jac = zeros(size(driver.rateset.rates));
        renorm_cal_T_rates = cal_T_rates./driver.qrenorm(driver.jacobian.temp_i)';
        for iiT = 1 : length(driver.jacobian.temp_i)
          spectra_due_to_T_jac = spectra_due_to_T_jac + renorm_cal_T_rates(iiT)*m_ts_jac(:,driver.jacobian.temp_i(iiT));
        end
        load f2645.mat
        plot(f2645,spectra_due_to_T_jac,'b',f2645,driver.rateset.rates,'m.-',...
             f2645,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2645,driver.rateset.rates-spectra_due_to_T_jac,'r.-')
        plot(f2645,driver.rateset.rates,'m.-',f2645,driver.rateset.rates-spectra_due_to_T_jac,'r.-',f2645,m_ts_jac(:,1)*10+0.85,'k.-')
      else
        fprintf(1,'WARNING trying to setting dT(z,lat,t)/dt using CAL anom but %s DNE so set xb(T) = 0 \n',izname);
        cal_T_rates = zeros(size(driver.jacobian.temp_i));
        spectra_due_to_T_jac = zeros(size(driver.rateset.rates));
        xb(driver.jacobian.temp_i) = 0.0;
      end

    elseif iZeroTVers == 1
      era_model_file = ['/asl/s1/sergio/home/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Aug20_2019_Clr/'];
      era_model_file = [era_model_file '/Desc/16DayAvgNoS_withanom/latbin' num2str(driver.iTile,'%02d') '_16day_avg.rp.mat'];
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
        load f2645.mat
        plot(f2645,spectra_due_to_T_jac,'b',f2645,driver.rateset.rates,'m.-',...
             f2645,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2645,driver.rateset.rates-spectra_due_to_T_jac,'r.-')
        plot(f2645,driver.rateset.rates,'m.-',f2645,driver.rateset.rates-spectra_due_to_T_jac,'r.-',f2645,m_ts_jac(:,1)*10+0.85,'k.-')

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
        era_ptempanom = squeeze(junk(driver.iTile,:,:));
        cal_T_rates = era_ptempanom(:,driver.i16daytimestep);

        spectra_due_to_T_jac = zeros(size(driver.rateset.rates));
        renorm_cal_T_rates = cal_T_rates./driver.qrenorm(driver.jacobian.temp_i)';
        for iiT = 1 : length(driver.jacobian.temp_i)
          spectra_due_to_T_jac = spectra_due_to_T_jac + renorm_cal_T_rates(iiT)*m_ts_jac(:,driver.jacobian.temp_i(iiT));
        end
        load f2645.mat
        plot(f2645,spectra_due_to_T_jac,'b',f2645,driver.rateset.rates,'m.-',...
             f2645,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2645,driver.rateset.rates-spectra_due_to_T_jac,'r.-')
        plot(f2645,driver.rateset.rates,'m.-',f2645,driver.rateset.rates-spectra_due_to_T_jac,'r.-',f2645,m_ts_jac(:,1)*10+0.85,'k.-')

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
      izname = ['SAVE_BESTRUNv1/OutputAnomaly_CAL/' num2str(driver.iTile,'%02d') '/anomtest_timestep' num2str(driver.i16daytimestep) '.mat'];
      if exist(izname)
        fprintf(1,'WARNING setting dO3(z,lat,t)/dt using CAL anom in %s \n',izname);       
        izt = load(izname);
        cal_O3_rates = izt.oem.finalrates(driver.jacobian.ozone_i);

        spectra_due_to_O3_jac = zeros(size(driver.rateset.rates));
        renorm_cal_O3_rates = cal_O3_rates./driver.qrenorm(driver.jacobian.ozone_i)';
        for iiO3 = 1 : length(driver.jacobian.ozone_i)
          spectra_due_to_O3_jac = spectra_due_to_O3_jac + renorm_cal_O3_rates(iiT)*m_ts_jac(:,driver.jacobian.ozone_i(iiO3));
        end
        load f2645.mat
        plot(f2645,spectra_due_to_O3_jac,'b',f2645,driver.rateset.rates,'m.-',...
             f2645,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2645,driver.rateset.rates-spectra_due_to_O3_jac,'r.-')
        plot(f2645,driver.rateset.rates,'m.-',f2645,driver.rateset.rates-spectra_due_to_O3_jac,'r.-',f2645,m_ts_jac(:,1)*10+0.85,'k.-')
      else
        fprintf(1,'WARNING trying to setting dO3(z,lat,t)/dt using CAL anom but %s DNE so set xb(O3) = 0 \n',izname);
        cal_O3_rates = zeros(size(driver.jacobian.ozone_i));
        spectra_due_to_O3_jac = zeros(size(driver.rateset.rates));
        xb(driver.jacobian.ozone_i) = 0.0;
      end

    elseif iZeroO3Vers == 1
      era_model_file = ['/asl/s1/sergio/home/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeProfs/LATS40_avg_made_Aug20_2019_Clr/'];
      era_model_file = [era_model_file '/Desc/16DayAvgNoS_withanom/latbin' num2str(driver.iTile,'%02d') '_16day_avg.rp.mat'];
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
        load f2645.mat
        plot(f2645,spectra_due_to_O3_jac,'b',f2645,driver.rateset.rates,'m.-',...
             f2645,sum(m_ts_jac(:,driver.jacobian.temp_i)')*100,'k.-',f2645,driver.rateset.rates-spectra_due_to_O3_jac,'r.-')
        plot(f2645,driver.rateset.rates,'m.-',f2645,driver.rateset.rates-spectra_due_to_O3_jac,'r.-',f2645,m_ts_jac(:,1)*10+0.85,'k.-')

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
        era_ozoneanom = squeeze(junk(driver.iTile,:,:));
        cal_O3_rates = era_ozoneanom(:,driver.i16daytimestep);

        spectra_due_to_O3_jac = zeros(size(driver.rateset.rates));
        renorm_cal_O3_rates = cal_O3_rates./driver.qrenorm(driver.jacobian.ozone_i)';
        for iiT = 1 : length(driver.jacobian.ozone_i)
          spectra_due_to_O3_jac = spectra_due_to_O3_jac + renorm_cal_O3_rates(iiT)*m_ts_jac(:,driver.jacobian.ozone_i(iiT));
        end
        load f2645.mat
        plot(f2645,spectra_due_to_O3_jac,'b',f2645,driver.rateset.rates,'m.-',...
             f2645,sum(m_ts_jac(:,driver.jacobian.ozone_i)')*100,'k.-',f2645,driver.rateset.rates-spectra_due_to_O3_jac,'r.-')
        plot(f2645,driver.rateset.rates,'m.-',f2645,driver.rateset.rates-spectra_due_to_O3_jac,'r.-',f2645,m_ts_jac(:,1)*10+0.85,'k.-')

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

if settings.set_tracegas == +1 & driver.i16daytimestep < 0
  disp('setting constant rates for tracegas apriori : CO2 = 2.2  CH4 = 4.5 N2O = 0.8')
  if settings.co2lays == 1
    xb(1) = 2.2 * 1;  % Set CO2 apriori
    xb(2) = 0.8 * 1;
    xb(3) = 4.5 * 1;
    xb(4) = 0.0; %% clouds, so dunno value
    xb(5) = 0.0; %% clouds, so dunno value


  elseif settings.co2lays == 3
    xb(1) = 2.2 * 1;        % Set CO2 apriori lower trop
    xb(2) = 2.2 * 1;        % Set CO2 apriori mid trop
    xb(3) = 2.2 * 1;        % Set CO2 apriori strat

    xb(4) = 0.8 * 1;        % set N2O 
    xb(5) = 4.5 * 1;        % set CH4
    xb(6) = 0.0; %% CFC11
    xb(7) = 0.0; %% CFC12

  end

elseif settings.set_tracegas == +1 & driver.i16daytimestep > 0
  %% 16 day timesteps so 365/16 = 22.81 per year
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

  elseif settings.co2lays == 3
    deltaT = 365/16; %% days per timestep

    xb(1) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori lower trop
    xb(2) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori mid trop
    xb(3) = 2.2 * (driver.i16daytimestep-1)/deltaT * 1.0;    % Set CO2 apriori strat

    xb(4) = 0.8 * (driver.i16daytimestep-1)/deltaT * 1.0;        % set N2O 
    xb(5) = 4.5 * (driver.i16daytimestep-1)/deltaT * 1.0;        % set CH4
    xb(6) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;  % set CFC11, before Aug 23 the mult was 1
    xb(7) = -1.25 * (driver.i16daytimestep-1)/deltaT * 0;  % set CFC12, before Aug 23 the mult was 1
  end

elseif settings.set_tracegas == +2 & driver.i16daytimestep > 1
  junk = 365/16; %% days per timestep 
  junk = (driver.i16daytimestep-1);
  str = ['setting bootstrap time varying rates for tracegas apriori : CO2 = 2.2  CH4 = 4.5 N2O = 0.8 CFC = -1.25 based on previous timestep'];
  disp(str);
  fminus = ['OutputAnomaly_OBS/' num2str(driver.iTile,'%02d') '/anomtest_timestep' num2str(driver.i16daytimestep-1) '.mat'];
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

      %% overwrite!!!
      xb(1) = prev.oem.finalrates(1);
      xb(2) = prev.oem.finalrates(2);
      xb(3) = prev.oem.finalrates(3);
      xb(4) = prev.oem.finalrates(4);
      xb(5) = prev.oem.finalrates(5);

      fprintf(1,'changed CO2   apriori from %8.6f to %8.6f ppv \n',xb0(1),xb(1))
      fprintf(1,'changed N2O   apriori from %8.6f to %8.6f ppb \n',xb0(2),xb(2))
      fprintf(1,'changed CH4   apriori from %8.6f to %8.6f ppb \n',xb0(3),xb(3))
      fprintf(1,'changed CFC11 apriori from %8.6f to %8.6f ppt \n',xb0(4),xb(4))
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
  xb = xb(:,driver.iTile);
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
elseif settings.offsetrates > 0 & driver.i16daytimestep > 0 & settings.ocb_set == 0 & sum(abs(driver.rateset.rates)) > 0
   disp('Offset anomaly (obs) by constant 0.01/year to get sensitivity to AIRS BT drift')
   oktimes = load('ok365times.mat');
   oktimes = oktimes.okdates(driver.i16daytimestep);   %% need to get the CO2 jac at this time!!!!
   fprintf(1,'  --> oktimes = %8.6f which is %8.6f years away from 2002.75 \n',oktimes,oktimes - 2002.75)
   driver.rateset.rates = driver.rateset.rates + 0.01 * (oktimes - 2002.75);
elseif settings.offsetrates > 0 & driver.i16daytimestep > 0 & settings.ocb_set == 1 & sum(abs(driver.rateset.rates)) > 0
   disp('Offset anomaly (cal) by constant 0.01/year to get sensitivity to AIRS BT drift')
   oktimes = load('ok365times.mat');
   oktimes = oktimes.okdates(driver.i16daytimestep);   %% need to get the CO2 jac at this time!!!!
   fprintf(1,'  --> oktimes = %8.6f which is %8.6f years away from 2002.75 \n',oktimes,oktimes - 2002.75)
   driver.rateset.rates = driver.rateset.rates + 0.01 * (oktimes - 2002.75);
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
  driver.rateset.unc_rates = driver.rateset.unc_rates.*nc_cor;    %% THIS IS AN ARRAY
  %driver.rateset.unc_rates = driver.rateset.unc_rates *sqrt(1/8); %% this accounts for counting .... 
end

%driver.rateset.unc_rates = 0.001 * ones(size(driver.rateset.unc_rates));   %% THIS IS AN ARRAY

%%%%%%%%%%%%%%%%%%%%%%%%% >>>>>>
if settings.obs_corr_matrix > 0
  addpath /home/sergio/MATLABCODE
  thecov = load('/home/sergio/MATLABCODE/oem_pkg_run/Simulate_Calcs/thecov_clear');
  [f2645,i2645] = map_2834_to_2645;
  junk = thecov.thecov; junk = junk(i2645,i2645);

  %% see https://www.mathworks.com/help/finance/corr2cov.html
  %%junk = corr2cov(driver.rateset.unc_rates,junk);
  %% new, fast
  junk0 = junk;
  junk0 = diag(driver.rateset.unc_rates) * junk0 * diag(driver.rateset.unc_rates);

  %% new,slow
  %[mm,nn] = size(junk);
  %for ii = 1 : mm
  %  for jj = 1 : nn
  %    junk(ii,jj) = junk(ii,jj) * driver.rateset.unc_rates(ii) * driver.rateset.unc_rates(jj);
  %  end
  %end
  %sum(sum(junk-junk0))/nn/nn
  junk = junk0;  
  fprintf(1,'numchans, rank, mean rank, condition number of obs cov matrix = %6i %6i %8.6e %8.6e \n',nn,rank(junk),rank(junk)/length(junk),cond(junk))

  %% old
  %plot(f2645,driver.rateset.unc_rates.^2,'r',thecov.fairs(i2645),diag(junk))
  %%keyboard_nowindow
  %slope = (2645+1);
  %xind = 1:2645; diagind=slope*(xind-1)+1;
  %junk(diagind) = driver.rateset.unc_rates.^2;

  driver.rateset.unc_rates = junk;  %% THIS IS A SQUARE MATRTIX
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

%figure(13); colormap jet; imagesc(log10(abs(aux.m_ts_jac'))); caxis([-12 0]); colorbar 
%disp('here 3'); pause

%driver.jacobian.thin = aux.m_ts_jac;
%keyboard_nowindow