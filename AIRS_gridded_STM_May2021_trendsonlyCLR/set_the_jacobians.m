% Get jacobians, and combine the 97 layer T(z)/WV(z)/O3(z) into N layers
%[m_ts_jac0,nlays,qrenorm]  = get_jac(driver.jacobian.filename,driver.jac_indexINSIDEbin,iVersJac);
if iVersJac == 2012 | iVersJac == 2019
  %% [m_ts_jac0,nlays,qrenorm,freq2645]  = get_jac_fast(driver.jacobian.filename,driver.iibin,driver.iLon,driver.iLat,2021);   %% I have not made jacs for this time period
  [m_ts_jac0,nlays,qrenorm,freq2645]  = get_jac_fast(driver.jacobian.filename,driver.iibin,driver.iLon,driver.iLat,iVersJac,topts);
else
  [m_ts_jac0,nlays,qrenorm,freq2645]  = get_jac_fast(driver.jacobian.filename,driver.iibin,driver.iLon,driver.iLat,iVersJac,topts);
end
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
  playsx = load('/home/sergio/MATLABCODE/airslevels.dat');
  playsx = flipud(plevs2plays(playsx));
  for kkk = 1 : length(layWV)
    iavg = layWV{kkk}-max(driver.jacobian.scalar_i);
    plays(kkk) = mean(playsx(iavg));
  end
else
  fprintf(1,'setting iNlays_retrieve ( > 60) from %2i to 97 \n',iNlays_retrieve);
  m_ts_jac_wv = m_ts_jac0;
  iNlays_retrieve = 97;
  plays = load('/home/sergio/MATLABCODE/airslevels.dat');
  plays = flipud(plevs2plays(plays));
  plays = plays(4:100);
end
%figure(10); imagesc(m_ts_jac_wv'); colorbar
%figure(11); imagesc(m_ts_jac_t'); colorbar
%figure(12); imagesc(m_ts_jac_o3'); colorbar
%disp('here 2'); pause

%% replace CO2,N2O,CH4 jacs
if driver.i16daytimestep > 0  
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
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,1);
  elseif iXJac == 2 & iDoStrowFiniteJac == 2
    fprintf(1,'updating time varying kCARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs -1,+2 for testing 6/24-27/2019...\n');
    fprintf(1,'  note before July2, the tracegas profile (CO2/N2O/CH4) was US Std shoehorned and multiplied, so ppm was quite wonky except at 500 mb');
    fprintf(1,'else turned off \n')
    %% const kCARTA jacs, update the trace gases
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% only used this on 6/26
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% added this on 6/27
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,2); %% added this on 6/27
  elseif iXJac == 2 & iDoStrowFiniteJac == 3
    fprintf(1,'updating time varying kCARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 6/24-27/2019...\n');
    %% kcarta strow finite jacs
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% only used this on 6/26
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% added this on 6/27
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,3); %% added this on 6/27
  elseif iXJac == 1 & iDoStrowFiniteJac == 3
    fprintf(1,'updating time varying SARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 6/24-27/2019...\n');
    %% sarta strow finite jacs
    iVarType = -3; %% this uses SARTA  finitediff jacs, which I have shown are bad?? or good??
    iVarType = +3; %% this uses kCARTA finitediff jacs, which I have shown are good, just want to test
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
  elseif iXJac == 2 & iDoStrowFiniteJac == 4
    fprintf(1,'updating time varying kCARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 12/08/2019 ... age of air\n');
    %% const kCARTA jacs, update the trace gases
    %% kcarta strow finite jacs
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,4); %% only used this on 6/26
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,4); %% added this on 6/27
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,4); %% added this on 6/27
  elseif iXJac == 1 & iDoStrowFiniteJac == 3
    fprintf(1,'updating time varying SARTA CO2/N2O/CH4 jacs with interpolated BIGSTEP time varying jacs3 for testing 12/08/2019...\n');
    fprintf(1,'else turned off \n')
    %% sarta strow finite jacs
    iVarType = -3; %% this uses SARTA  finitediff jacs, which I have shown are bad?? or good??
    iVarType = +3; %% this uses kCARTA finitediff jacs, which I have shown are good, just want to test
    iVarType = +4; %% this uses kCARTA finitediff jacs, which I have shown are good, just want to test
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,driver.i16daytimestep,iVarType); %% added this on 8/3
  end
elseif driver.i16daytimestep < 0
  %% this is for 1 average rate
  if iXJac == 0
    fprintf(1,'not updating CO2/N2O/CH4 jacs for TRENDS ... keeping same jacs, give better results\n');
  elseif iXJac == 1 | iXJac == 2
    fprintf(1,'updating CO2/N2O/CH4 jacs for TRENDS ...\n');
    m_ts_jac_coljac = replace_time_co2jac(m_ts_jac_coljac,driver.iibin,iMidPoint_TimeStepUse,3);
    m_ts_jac_coljac = replace_time_n2ojac(m_ts_jac_coljac,driver.iibin,iMidPoint_TimeStepUse,3);
    m_ts_jac_coljac = replace_time_ch4jac(m_ts_jac_coljac,driver.iibin,iMidPoint_TimeStepUse,3);
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
  m_ts_jac = add_deltaJ_to_jacs(m_ts_jac,driver.iibin,-1);
end
