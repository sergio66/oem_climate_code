addpath /home/sergio/MATLABCODE/COLORMAP

if ~exist('deltaT')
  load('/home/sergio/MATLABCODE/COLORMAP/LLS/llsmap5.mat')
  %load /asl/s1/sergio/JUNK//gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX100_50fatlayers_CLIMCAPS_MERRA2_AMIP6_feedback.mat
  load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX3_50fatlayers_AIRSL3_ERA5_CMIP6_feedback.mat
end

airslevels = flipud(load('/home/sergio/MATLABCODE/airslevels.dat'));
airslays = meanvaluebin(airslevels);

iX = 4;
tvalue = 0.050;
tvalue = 0.025;
alphavalue = 0.05;

if ~exist('iOffset')
  iOffset = 0;
end

fprintf(1,'tvalue     = %8.6f \n',tvalue')
fprintf(1,'alphavalue = %8.6f \n',alphavalue')

%%%%%%%%%%%%%%%%%%%%%%%%%

iX = iX + 1;
xTtrend = reshape(deltaT,101,72,64); xTtrend = squeeze(nanmean(xTtrend,2)); % whos xTtrend
if ~exist('deltaTlat')
   deltaTlat = xTtrend';
end
figure(iOffset+1); pcolor(deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(iOffset+2); pcolor(xTtrend); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(iOffset+3); pcolor(xTtrend-deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); sum(sum(xTtrend-deltaTlat'))

N72 = 72;
N72 = 1;
xTunc = reshape(deltaTunc,101,72,64); xTunc = squeeze(nanmean(xTunc,2))/sqrt(N72); % whos xTunc
figure(iOffset+2); pcolor(xTunc); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02)
xTunc = xTunc';

normcdf([-1.96 0 +1.96])
zscore = (abs(deltaTlat) - 0)./xTunc;
pvalueT = 1-normcdf(zscore);
figure(iOffset+3); pcolor(pvalueT'); shading interp; colorbar; set(gca,'ydir','reverse'); 

figure(iOffset+iX); pcolor(deltaTlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dT/dt K/yr');
figure(iOffset+2); pcolor(xTunc'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dT_{\sigma}/dt K/yr');
figure(iOffset+3); pcolor(pvalueT'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueT');

ppvalueT = pvalueT; ppvalueT(ppvalueT > tvalue) = NaN;
ppvalueT = zeros(size(pvalueT)); ppvalueT(pvalueT < tvalue) = 1;
figure(iOffset+4); pcolor(ppvalueT'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueT');

figure(iOffset+iX); hold on
for jj = 1 : 101
  for ii = 1 : 64
    if ppvalueT(ii,jj) == 1
      hold on; plot(ii,jj,'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iOffset+iX); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%

%% https://www.mathworks.com/help/stats/hypothesis-testing.html
bonk    = deltaTlat(:);
bonkunc = abs(xTunc(:));
for ii = 1 : length(bonk)
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii));
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',0.01);
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',0.05);
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',0.10);
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',alphavalue);
end

figure(iOffset+iX+1); pcolor(rlat,airslays,deltaTlat(:,1:100)'); shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); caxis([-1 +1]*0.101); xlabel('Latitude'); ylabel('Pressure [mb]'); 
colorbar('horizontal'); % title('dT/dt K/yr');
figure(iOffset+iX+1); hold on
ht = reshape(h,64,101);
for jj = 1 : 100
  for ii = 1 : 64
    if ht(ii,jj) == 1
      hold on; plot(rlat(ii),airslays(jj),'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iOffset+iX+1); hold off
set(gca,'yscale','log'); ylim([10 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iX = iX + 2;
xRHtrend = reshape(deltaRH,100,72,64); xRHtrend = squeeze(nanmean(xRHtrend,2)); % whos xRHtrend
if ~exist('deltaRHlat')
   deltaRHlat = xRHtrend';
end
figure(iOffset+1); pcolor(deltaRHlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(iOffset+2); pcolor(xRHtrend); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125)
figure(iOffset+3); pcolor(xRHtrend-deltaRHlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); nansum(nansum(xRHtrend-deltaRHlat'))

xRHunc = reshape(deltaRHunc,100,72,64); xRHunc = squeeze(nanmean(xRHunc,2))/sqrt(N72); % whos xRHunc
figure(iOffset+2); pcolor(xRHunc); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02)
xRHunc = xRHunc';

normcdf([-1.96 0 +1.96])
zscore = (abs(deltaRHlat) - 0)./xRHunc;
pvalueRH = 1-normcdf(zscore);
figure(iOffset+3); pcolor(pvalueRH'); shading interp; colorbar; set(gca,'ydir','reverse'); 

figure(iOffset+iX); pcolor(deltaRHlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.125); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dRH/dt pc/yr');
figure(iOffset+2); pcolor(xRHunc'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dRH_{\sigma}/dt pc/yr');
figure(iOffset+3); pcolor(pvalueRH'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueRH');

ppvalueRH = pvalueRH; ppvalueRH(ppvalueRH > tvalue) = NaN;
ppvalueRH = zeros(size(pvalueRH)); ppvalueRH(pvalueRH < tvalue) = 1;
figure(iOffset+4); pcolor(ppvalueRH'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueRH');

figure(iOffset+iX); hold on
for jj = 1 : 100
  for ii = 1 : 64
    if ppvalueRH(ii,jj) == 1
      hold on; plot(ii,jj,'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iOffset+iX); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%

%% https://www.mathworks.com/help/stats/hypothesis-testing.html
bonk    = deltaRHlat(:);
bonkunc = abs(xRHunc(:));
for ii = 1 : length(bonk)
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii));
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',0.01);
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',0.05);
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',0.10);
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',alphavalue);
end

figure(iOffset+iX+1); pcolor(rlat,airslays,deltaRHlat(:,1:100)'); shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); caxis([-1 +1]*0.126); xlabel('Latitude'); ylabel('Pressure [mb]'); 
colorbar('horizontal'); % title('dRH/dt percentyr');
figure(iOffset+iX+1); hold on
ht = reshape(h,64,101);
for jj = 1 : 100
  for ii = 1 : 64
    if ht(ii,jj) == 1
      hold on; plot(rlat(ii),airslays(jj),'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iOffset+iX+1); hold off
ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iX = iX + 2;
xWVtrend = reshape(fracWV,101,72,64); xWVtrend = squeeze(nanmean(xWVtrend,2)); % whos xWVtrend
if ~exist('fracWVlat')
   fracWVlat = xWVtrend';
end
figure(iOffset+1); pcolor(fracWVlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01)
figure(iOffset+2); pcolor(xWVtrend); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01)
figure(iOffset+3); pcolor(xWVtrend-fracWVlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01); nansum(nansum(xWVtrend-fracWVlat'))

xWVunc = reshape(fracWVunc,101,72,64); xWVunc = squeeze(nanmean(xWVunc,2))/sqrt(N72); % whos xWVunc
figure(iOffset+2); pcolor(xWVunc); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02)
xWVunc = xWVunc';

normcdf([-1.96 0 +1.96])
zscore = (abs(fracWVlat) - 0)./xWVunc;
pvalueWV = 1-normcdf(zscore);
figure(iOffset+3); pcolor(pvalueWV'); shading interp; colorbar; set(gca,'ydir','reverse'); 

figure(iOffset+iX); pcolor(fracWVlat'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([-1 +1]*0.01); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dWV/dt frac/yr');
figure(iOffset+2); pcolor(xWVunc'); shading interp; colormap(llsmap5); colorbar; set(gca,'ydir','reverse'); caxis([0 +1]*0.02); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('dWV_{\sigma}/dt frac/yr');
figure(iOffset+3); pcolor(pvalueWV'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueWV');

ppvalueWV = pvalueWV; ppvalueWV(ppvalueWV > tvalue) = NaN;
ppvalueWV = zeros(size(pvalueWV)); ppvalueWV(pvalueWV < tvalue) = 1;
figure(iOffset+4); pcolor(ppvalueWV'); shading interp; colorbar; set(gca,'ydir','reverse'); xlabel('Latitude'); ylabel('AIRS PressLayer'); title('PvalueWV');

figure(iOffset+iX); hold on
for jj = 1 : 101
  for ii = 1 : 64
    if ppvalueWV(ii,jj) == 1
      hold on; plot(ii,jj,'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iOffset+iX); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%

%% https://www.mathworks.com/help/stats/hypothesis-testing.html
bonk    = fracWVlat(:);
bonkunc = abs(xWVunc(:));
for ii = 1 : length(bonk)
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii));
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',0.01);
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',0.05);
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',0.10);
  [h(ii),pvalue(ii),ci(ii,:)] = ztest(bonk(ii),0,bonkunc(ii),'Alpha',alphavalue);
end

figure(iOffset+iX+1); pcolor(rlat,airslays,fracWVlat(:,1:100)'); shading interp; colormap(llsmap5); set(gca,'ydir','reverse'); caxis([-1 +1]*0.0101); xlabel('Latitude'); ylabel('Pressure [mb]'); 
colorbar('horizontal'); % title('dWV/dt frac/yr');
figure(iOffset+iX+1); hold on
ht = reshape(h,64,101);
for jj = 1 : 100
  for ii = 1 : 64
    if ht(ii,jj) == 1
      hold on; plot(rlat(ii),airslays(jj),'k.','Markersize',4,'linewidth',4);
    end
  end
end
figure(iOffset+iX+1); hold off
ylim([100 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
