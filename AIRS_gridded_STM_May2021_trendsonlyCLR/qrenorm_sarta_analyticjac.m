% m_ts_jac0 = [jac2; jac4; jac6; 0*jac2; 0*jac2; jacST; jacWV; jacTZ; jacOZ]';
% nlays = iaNumLay;
% [mm,nn] = size(m_ts_jac0);
% qrenorm = ones(1,mm);
% freq2645 = w;
% colo3 = nansum(jac3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iVers = 0;   %% qrenorm = ones                     EASY BASIC
iVers = 1;   %% qrenorm is wiertd : jac/100/ppmv   WRONG
iVers = 2;   %% try to duplicate get_jac_fast.m    DEFAULT

nlays = iaNumLay;

m_ts_jac0 = [jac2; jac4; jac6; 0*jac2; 0*jac2; jacST; jacWV; jacTZ; jacOZ]';
[mm,nn] = size(m_ts_jac0);
qrenorm = ones(1,mm);
qrenorm = ones(1,mm);
      
if iVers == 1
  %% see [xm_ts_jac0,xnlays,xqrenorm,xfreq2645,~,xprofilejunk]  = get_jac_fast(driver.jacobian.filename,driver.iibin,driver.iLon,driver.iLat,iVersJac,iOldORNew,topts)
  qrenorm(1) = 2.2;
  qrenorm(2) = 1.0;
  qrenorm(2) = 0.3;
  qrenorm(3) = 5.0;
  qrenorm(4) = 1.0;
  qrenorm(5) = 1.0;
  qrenorm(6) = 0.1;
  
  qrenorm(7:end) = 0.01;
  scale = 0.01;
  
  m_ts_jac0 = [jac2/100/qrenorm(1); jac4/100/qrenorm(2); jac6/100/qrenorm(3); 0*jac2; 0*jac2; jacST/100/qrenorm(6);];
  m_ts_jac0 = [m_ts_jac0;  jacWV*scale; jacTZ*scale; jacOZ*scale;]';
  
elseif iVers == 2
  [~,zsarta.subjac.ppmv2] = layers2ppmv(havg,pavg,1,2);
  [~,zsarta.subjac.ppmv4] = layers2ppmv(havg,pavg,1,4);
  [~,zsarta.subjac.ppmv6] = layers2ppmv(havg,pavg,1,6);
  qrenorm(1) = 2.2/zsarta.subjac.ppmv2;            %% 2.2 ppmv/400 ppmv       
  qrenorm(2) = 1.0/(zsarta.subjac.ppmv4*1000);     %% 1 ppb/300 ppb
  qrenorm(3) = 5.0/(zsarta.subjac.ppmv6*1000);     %% 5 ppb/1840 ppb
  qrenorm(4) = 1.0;
  qrenorm(5) = 1.0;
  qrenorm(6) = 0.1;

  qrenorm(7:end) = 0.01;
  scale = 0.01;
  
  m_ts_jac0 = [jac2*qrenorm(1); jac4*qrenorm(2); jac6*qrenorm(3); 0*jac2; 0*jac2; jacST*qrenorm(6);];
  m_ts_jac0 = [m_ts_jac0;  jacWV*scale; jacTZ*scale; jacOZ*scale;]';

  %% but really
  qrenorm(1) = 2.2;         %% 2.2 ppmv/400 ppm
  qrenorm(2) = 1.0;         %% 1 ppb/300 ppb
  qrenorm(3) = 5.0;         %% 5 ppm/1840 ppb
  qrenorm(4) = 1.0;         %% 1 ppt/300 ppt
  qrenorm(5) = 1.0;         %% 1 ppt/600 ppt
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq2645 = w;
colo3 = nansum(xjac3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
