% Statistical significance of seasonal warming/cooling trends
% Josef Ludeschera Armin Bundea and Hans Joachim Schellnhuber
% E2998–E3003 | PNAS | Published online March 27, 2017
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5393220/pdf/pnas.201700838.pdf
%
% also see JOURNAL OF  GEOPHYSICAL RESEARCH, VOL. 105, NO.  D6,  PAGES  7337-7356,  MARCH 27, 2000 
% Statistical significance of trends and trend  differences 
% in  layer-average atmospheric temperature  time  series 
% B. D.  Santer, • T.  M.  L. Wigley, 2 J. S. Boyle, • D.  J. Gaffen, 3 J. J. Hnilo, • 
% D.  Nychka, 2 D.  E.  Parker, 4 and K. E.  Taylor
%
% see ../../oem_pkg_run_sergio_AuxJacs/MakeProfs/fit_robust_one_lat.m for lag1
% Get lag-1 correlation (ignoring that we don't have all days)
%for i = 1:nf
%   y = squeeze(all_bt_resid(:,i));
%   k = remove_nan(y);
%   if length(k) > 100
%      l = xcorr(y(k),1,'coeff');
%      lag(i) = l(1);
%   else
%      lat(i) = NaN;
%   end
%end
%
% also see strow_override_defaults_latbins_AIRS_fewlays.m
% if driver.i16daytimestep < 0
%   nc_cor = nc_rates(driver);
%   driver.rateset.unc_rates = driver.rateset.unc_rates.*nc_cor;    %% THIS IS AN ARRAY
%   %driver.rateset.unc_rates = driver.rateset.unc_rates *sqrt(1/8); %% this accounts for counting ....
% end
%
% where from nc_rates.m
%    nc = (1+lagcor_obs_anom(iibin,:))./(1-lagcor_obs_anom(iibin,:));
%    nc = sqrt(nc);
%    nc_cor = real(nc');

% https://www.statology.org/standard-error-of-regression-slope/
% t-statistic = 1.678 The p-value that corresponds to this test (score vs number of hours studied)
%  statistic is 0.107. Since this p-value is not less than .05, this
%  indicates that ‘hours studied’ does not have a statistically
%  significant relationship with final exam score.
% 25 students in class, fitting linear polynmial so d of f = 2; so pvalue = 2*(1-tcdf(1.678,25-2)) = 0.1069

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iDo = -1;
if iDo > 0
  x = 1:100;
  y = 3*x + 5 + 100*randn(size(x));
  y = x.*x + 3*x + 5 + 100*randn(size(x));
  y = 0.1*x + 5 + 100*randn(size(x));    %%%% we are trying to see if the slope could be 0!!!!!!! H0 = slope == 0
  P = polyfit(x,y,1)
  ynew = polyval(P,x);
  [B,err] = Math_tsfit_lin_robust(x,y,0);
  slopeB = B(2)/365;
  interceptB = B(1);
  PB = [slopeB interceptB]
  
  plot(err.resid,y-ynew,err.resid,err.resid)
  
  dof = length(x)-2;
  var_res = sum(err.resid.^2)/(length(x)-2);   %% Eqn 4 of Santer = variance of residuals
  denom = (x-mean(x));
    denom = denom/365; %% remember change to years, to compare Math_tsfit
  denom = sum(denom.*denom);
  std_error = sqrt(var_res/denom);
  [err.se(2) std_error]
  
  %% https://www.mathworks.com/help/stats/f-statistic-and-t-statistic.html
  % Let n be your sample size
  % Let v be your degrees of freedom = samplesize - number of params fitted, in the case of a straight line = 2
  % Then:
  % pvalues = 2*(1-tcdf(abs(t),n-v))
  % and of course n-v = n-(n-2) = 2
  [err.t(2) B(2)/err.se(2)]                           %% now figure out t statistic = trend/se_trend = coeff estimate / standard error of estimate
  [err.p(2) 2*(1-tcdf(abs(err.t(2)),length(x)-dof))]  %% and the p value
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../../oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /home/sergio/MATLABCODE/COLORMAP
load('/home/sergio/MATLABCODE/COLORMAP/LLS/llsmap5.mat')

if ~exist('era5')
  %load /asl/s1/sergio/JUNK//gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX100_50fatlayers_CLIMCAPS_MERRA2_AMIP6_feedback.mat
  load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX3_50fatlayers_AIRSL3_ERA5_CMIP6_feedback.mat
end

iX = 4;
tvalue = 0.025;
tvalue = 0.050;  %% ok for ERA5, bad for UMBC
tvalue = 0.7500;

moo = load('latB64.mat');
lat64 = meanvaluebin(moo.latB2);
moo = load('/home/sergio/MATLABCODE/airslevels.dat');
plays = plevs2plays(flipud(moo));

iERA5orUMBC = input('Enter (-1) UMBC or (+1) ERA5 : ');
if length(iERA5orUMBC) == 0
  iERA5orUMBC = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%

iX = iX + 1;
if iERA5orUMBC == 1
  xTtrend = reshape(era5.trend_ptemp,100,72,64); xTtrend = squeeze(nanmean(xTtrend,2));
else
  xTtrend = reshape(deltaT,101,72,64); xTtrend = squeeze(nanmean(xTtrend,2));
  xTtrend = xTtrend(1:100,:);
end
deltaTlat = xTtrend';
figure(1); pcolor(deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(2); pcolor(xTtrend); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(3); pcolor(xTtrend-deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); %sum(sum(xTtrend-deltaTlat'))

if iERA5orUMBC == 1
  xTunc = reshape(era5.trend_ptemp_err,100,72,64); 
else
  xTunc = reshape(deltaTunc,101,72,64);
  xTunc = xTunc(1:100,:,:);
end
[lag,ratio,xTunc] = find_72_lag(xTunc);
xTunc = xTunc; 
figure(2); pcolor(xTunc); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02)
xTunc = xTunc .* ratio;
xTunc = xTunc';

%normcdf([-1.96 0 +1.96])
[zscore,pvalueT] = find_z_score(deltaTlat,xTunc,ratio);
figure(3); pcolor(pvalueT'); shading interp; colorbar; set(gca,'ydir','reverse'); 

figure(iX); pcolor(lat64,plays(1:100),deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dT/dt K/yr');
figure(2); pcolor(xTunc'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dT_{\sigma}/dt K/yr');
figure(3); pcolor(pvalueT'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueT');

ppvalueT = pvalueT; ppvalueT(ppvalueT > tvalue) = NaN;
ppvalueT = zeros(size(pvalueT)); ppvalueT(pvalueT < tvalue) = 1;
figure(4); pcolor(ppvalueT'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueT');

figure(iX); hold on
for jj = 1 : 100
  for ii = 1 : 64
    if ppvalueT(ii,jj) == 1
      hold on; plot(lat64(ii),plays(jj),'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iX); hold off; set(gca,'yscale','log'); ylim([1 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iX = iX + 1;
if iERA5orUMBC == 1
  xRHtrend = reshape(era5.trend_RH,100,72,64); xRHtrend = squeeze(nanmean(xRHtrend,2));
else
  xRHtrend = reshape(deltaRH,100,72,64);xRHtrend = squeeze(nanmean(xRHtrend,2));
end
%deltaRHlat = xRHtrend';
%figure(1); pcolor(deltaRHlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(2); pcolor(xRHtrend); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
%figure(3); pcolor(xRHtrend-deltaRHlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); %nansum(nansum(xRHtrend-deltaRHlat'))

if iERA5orUMBC == 1
  xRHunc = reshape(era5.trend_RH_err,100,72,64); 
else
  xRHunc = reshape(deltaRHunc,100,72,64);
end
[lag,ratio,xRHunc] = find_72_lag(xRHunc);
figure(2); pcolor(xRHunc); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02)
xRHunc = xRHunc .* ratio;
xRHunc = xRHunc';

[zscore,pvalueRH] = find_z_score(deltaRHlat,xRHunc,ratio);
figure(3); pcolor(pvalueRH'); shading interp; colorbar; set(gca,'ydir','reverse'); 

figure(iX); pcolor(lat64,plays(1:100),deltaRHlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dRH/dt pc/yr');
figure(2); pcolor(xRHunc'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dRH_{\sigma}/dt pc/yr');
figure(3); pcolor(pvalueRH'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueRH');

ppvalueRH = pvalueRH; ppvalueRH(ppvalueRH > tvalue) = NaN;
ppvalueRH = zeros(size(pvalueRH)); ppvalueRH(pvalueRH < tvalue) = 1;
figure(4); pcolor(ppvalueRH'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueRH');

figure(iX); hold on
for jj = 1 : 100
  for ii = 1 : 64
    if ppvalueRH(ii,jj) == 1
      hold on; plot(lat64(ii),plays(jj),'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iX); hold off; set(gca,'yscale','log'); ylim([1 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iX = iX + 1;
if iERA5orUMBC == 1
  xWVtrend = reshape(era5.trend_gas_1,100,72,64); xWVtrend = squeeze(nanmean(xWVtrend,2));
else
  xWVtrend = reshape(fracWV,101,72,64); xWVtrend = squeeze(nanmean(xWVtrend,2));
  xWVtrend = xWVtrend(1:100,:);
end
%fracWVlat = xWVtrend';
%figure(1); pcolor(fracWVlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01)
figure(2); pcolor(xWVtrend); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01)
%figure(3); pcolor(xWVtrend-fracWVlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01); nansum(nansum(xWVtrend-fracWVlat'))
xWVtrend = xWVtrend';

if iERA5orUMBC == 1
  xWVunc = reshape(era5.trend_gas_1_err,100,72,64); 
else
  xWVunc = reshape(fracWVunc,101,72,64);
  xWVunc = xWVunc(1:100,:,:);
end
[lag,ratio,xWVunc] = find_72_lag(xWVunc);
figure(2); pcolor(xWVunc); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02)
xWVunc = xWVunc .* ratio;
xWVunc = xWVunc';

[zscore,pvalueWV] = find_z_score(xWVtrend,xWVunc,ratio);
figure(3); pcolor(pvalueWV'); shading interp; colorbar; set(gca,'ydir','reverse'); 

figure(iX); pcolor(lat64,plays(1:100),xWVtrend'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dWV/dt frac/yr');
figure(2); pcolor(xWVunc'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dWV_{\sigma}/dt frac/yr');
figure(3); pcolor(pvalueWV'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueWV');

ppvalueWV = pvalueWV; ppvalueWV(ppvalueWV > tvalue) = NaN;
ppvalueWV = zeros(size(pvalueWV)); ppvalueWV(pvalueWV < tvalue) = 1;
figure(4); pcolor(ppvalueWV'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueWV');

figure(iX); hold on
for jj = 1 : 100
  for ii = 1 : 64
    if ppvalueWV(ii,jj) == 1
      hold on; plot(lat64(ii),plays(jj),'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iX); hold off; set(gca,'yscale','log'); ylim([1 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
