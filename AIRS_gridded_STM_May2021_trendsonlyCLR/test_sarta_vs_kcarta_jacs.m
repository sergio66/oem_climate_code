%% testing the code : see    set_the_jacobians.m    from lines 30-47
%% SARTA jacs on the fly
%%    [m_ts_jac0,nlays,qrenorm,freq2645,~,profilejunk]  = get_sarta_analyticjac_fast(settings.iNumYears,driver.iibin,driver.iLon,driver.iLat,topts,driver);

% before I did qrenorm in get_sarta_analyticjac_fast.m --> sarta_analytic_jac_trend.m
% figure(1); clf;
% plot(freq2645,m_ts_jac0(:,1)/100/2.2,freq2645,xm_ts_jac0(:,1)); xlim([645 2645])
% plot(freq2645,m_ts_jac0(:,1)/100/2.2,freq2645,xm_ts_jac0(:,1)); xlim([645 1645])
% plot(freq2645,m_ts_jac0(:,3)/100/5,freq2645,xm_ts_jac0(:,3)); xlim([645 1645])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% reading in pre-computed kCARTA jacs
[xm_ts_jac0,xnlays,xqrenorm,xfreq2645,~,xprofilejunk]  = get_jac_fast(driver.jacobian.filename,driver.iibin,driver.iLon,driver.iLat,iVersJac,iOldORNew,topts);
%% after reading analytic kCARTA dBT/d(logQ), this "get_jac_fast.m" routine does
%%  qrenorm(1) = 2.2/kcarta.subjac.ppmv2;            %% 2.2 ppmv/400 ppmv       
%%  qrenorm(2) = 1.0/(kcarta.subjac.ppmv4*1000);     %% 1 ppb/300 ppb
%%             = 1.0/(0.32*1000)
%%             = 0.0031
%%  qrenorm(3) = 5.0/(kcarta.subjac.ppmv6*1000);     %% 5 ppb/1840 ppb
%% then
%%   m_ts_jac = (ones(2645,1) * qrenorm) .* m_ts_jac_fast;  %%% SO THIS IS NORMALIZED JAC!!!!!!!!!!!!!!!!! JTRUE * RENORM
%% So effectively, 
%%  qrenorm(1) = 2.2;         %% 2.2 ppmv/400 ppm
%%  qrenorm(2) = 1.0;         %% 1 ppb/300 ppb
%%  qrenorm(3) = 5.0;         %% 5 ppm/1840 ppb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% after I did qrenorm in get_sarta_analyticjac_fast.m --> sarta_analytic_jac_trend.m

figure(1); clf; ixx = 1; plot(freq2645,m_ts_jac0(:,ixx),freq2645,xm_ts_jac0(:,ixx)); xlim([645 1645]); title('CO2'); legend('SARTA','KCARTA','location','best')
figure(2); clf; ixx = 2; plot(freq2645,m_ts_jac0(:,ixx),freq2645,xm_ts_jac0(:,ixx)); xlim([645 1645]); title('N2O'); legend('SARTA','KCARTA','location','best')
figure(3); clf; ixx = 3; plot(freq2645,m_ts_jac0(:,ixx),freq2645,xm_ts_jac0(:,ixx)); xlim([645 1645]); title('CH4'); legend('SARTA','KCARTA','location','best')
figure(4); clf; ixx = 6; plot(freq2645,m_ts_jac0(:,ixx),freq2645,xm_ts_jac0(:,ixx)); xlim([645 1645]); title('SKT'); legend('SARTA','KCARTA','location','best')
figure(5); clf; ixx = (1:nlays); ixx = 6 + ixx + 0*nlays; plot(freq2645,sum(m_ts_jac0(:,ixx),2),freq2645,sum(xm_ts_jac0(:,ixx),2));   xlim([645 1645]); title('WV'); legend('SARTA','KCARTA','location','best')
figure(6); clf; ixx = (1:nlays); ixx = 6 + ixx + 1*nlays; plot(freq2645,sum(m_ts_jac0(:,ixx),2),freq2645,sum(xm_ts_jac0(:,ixx),2));   xlim([645 1645]); title('Tz'); legend('SARTA','KCARTA','location','best')
figure(7); clf; ixx = (1:nlays); ixx = 6 + ixx + 2*nlays; plot(freq2645,sum(m_ts_jac0(:,ixx),2),freq2645,sum(xm_ts_jac0(:,ixx),2));   xlim([645 1645]); title('O3'); legend('SARTA','KCARTA','location','best')

printarray(qrenorm([1:6 7 110 200]), 'sarta qrenom 1:6, WV,T,O3')
printarray(xqrenorm([1:6 7 110 200]),'kcarta qrenom 1:6, WV,T,O3')
