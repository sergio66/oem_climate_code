function y = atmospheric_amplification(airsL3,iCosWgt,iPlot)

% Global Changes in Water Vapor 1979–2020
% Richard P. Allan1 , Kate M. Willett2 , Viju O. John3 , and Tim Trent4
% JGR Atmospheres RESEARCH ARTICLE 10.1029/2022JD036728
% The troposphere is resolved across seven pressure levels (300, 400, 500, 600, 700, 850, and 925 hPa). 

%% see eg atmospheric_amplification_plots2.m

%% this is same as atmospheric_amplification_plots2.m except we use the "nwp_spectral_trends_cmip6_era5_airsL3_umbc" variable instead of eg cmip6,era5,airsL3 vars
%% and vice versa

%% and here we do smoothing
epsx = 1e-2;
epsx = 1e-3; %% default
epsx = 1e-4;

%% wants airsL3.fields of this size
%%    stemprate: [1x4608 double]
%%       RHrate: [100x4608 double]
%%    ptemprate: [100x4608 double]
%%    waterrate: [100x4608 double]

if nargin == 1
  iPlot = -1;
  iCosWgt = +1;
elseif nargin == 2
  iPlot = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iCosWgt > 0
  load latB64.mat
  rlat65 = latB2; rlon73 = -180 : 5 : +180;
  rlon = -180 : 5 : +180;  rlat = latB2; 
  rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
  rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
  [Y,X] = meshgrid(rlat,rlon);
  X = X; Y = Y;
  YY = Y(:)'; 
  YY = cos(YY*pi/180);
elseif iCosWgt == -1
  YY = ones(1,4608);; %% everything equally weighted
end

stemprate = ones(100,1) * (airsL3.stemprate .* YY); 
kamask = ones(size(stemprate));
kamask(abs(stemprate) < epsx) = NaN;
kamask100 = ones(100,1) * YY;

delta = airsL3.ptemprate .* kamask100;
  boo = delta./stemprate;
  boo = boo.*kamask;
  boo = reshape(boo,100,72,64); 
  boo = squeeze(nanmean(boo,2));
  boo = boo(1:97,:);
if iPlot > 0
  figure(6); clf
  pcolor(boo); colorbar; shading interp; caxis([-1 +1]*2); set(gca,'ydir','reverse'); 
  plotaxis2;
  clf; plot(nanmean(boo,2),airsL3.pavg(1:97),'color','k','linewidth',3); plotaxis2; set(gca,'ydir','reverse'); 
  title('T amplification dT/dST')
  ylim([100 1000])
  xlim([-1 +1])
end
  y.T_amp     = nanmean(boo,2);
  y.T_amp_sig = nanstd(boo,[],2);

delta = airsL3.RHrate .* kamask100;
  boo = delta./stemprate;
  boo = boo.*kamask;
  boo = reshape(boo,100,72,64); 
  boo = squeeze(nanmean(boo,2));
  boo = boo(1:97,:);
if iPlot > 0
  figure(7); clf
  pcolor(boo); colorbar; shading interp; caxis([-1 +1]*2); set(gca,'ydir','reverse');
  plotaxis2;
  clf; plot(nanmean(boo,2),airsL3.pavg(1:97),'color','k','linewidth',3); plotaxis2; set(gca,'ydir','reverse'); 
  title('RH amplification dRH/dST')
  ylim([100 1000])
  xlim([-1 +1]*10)
end
  y.RH_amp = nanmean(boo,2);
  y.RH_amp_sig = nanstd(boo,[],2);

delta = airsL3.waterrate .* kamask100;
  boo = delta./stemprate;
  boo = boo.*kamask;
  boo = reshape(boo,100,72,64); 
  boo = squeeze(nanmean(boo,2))*100;
  boo = boo(1:97,:);
if iPlot > 0
  pcolor(boo); colorbar; shading interp; caxis([-1 +1]*2); set(gca,'ydir','reverse');
  plotaxis2;
  clf; plot(nanmean(boo,2),airsL3.pavg(1:97),'color','k','linewidth',3); plotaxis2; set(gca,'ydir','reverse');
  title('WV amplification 100*dWVfrac/dST')
  ylim([100 1000])
  xlim([-1 +1]*20)
end
  y.WVfrac_percent_amp = nanmean(boo,2);
  y.WVfrac_percent_amp_sig = nanstd(boo,[],2);
