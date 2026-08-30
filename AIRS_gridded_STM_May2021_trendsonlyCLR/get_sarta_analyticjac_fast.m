function [xm_ts_jac0,xnlays,xqrenorm,xfreq2645,xcolo3,xprofilejunk]  = get_sarta_analyticjac_fast(iNumYears,iibin,iLon,iLat,topts,driver);

%% see read_fileMean17years.m for h,p
read_fileMean17years

[xm_ts_jac0,xnlays,xqrenorm,xfreq2645,xcolo3,xprofilejunk] = sarta_analytic_jac_trend(driver,hMean17years,ha,pMean17years,pa,iibin);


