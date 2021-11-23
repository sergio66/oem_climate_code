addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TROPOPAUSE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/matlib/science/
addpath /asl/matlib/aslutil
addpath /asl/matlib/maps
addpath /asl/matlib/h4tools

load llsmap5
%if length(llsmap5) == 64
%  llsmap5 = llsmap5(2:end,:);
%end

iJunk = input('clear 50 figs???? (-1 = no = default) ');
if length(iJunk) == 0
  iJunk = -1;
end
if iJunk > 0
  for ii = 1 : 50
    figure(ii); colormap jet
  end
else
  figure(6);  clf; %% this puts d/dt ST plot with correct title
  figure(28); clf; %% this puts d/dt ST plot with correct title
  figure(29); clf; %% this puts d/dt ST plot with correct title
  figure(30); clf; %% this puts d/dt ST plot with correct title
end

iNumLay = 33;
iNumLay = 20;

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
playsN = plevs(1:100)-plevs(2:101);
playsD = log(plevs(1:100)./plevs(2:101));
plays = playsN./playsD;
plays = flipud(plays);

clear pavg

iNorD = input('night (+1) or day (-1)  (+1 = night = default) ');
if length(iNorD) == 0
  iNorD = +1;
end
iDorA = iNorD;

%dataset = input('Enter (+1) Strow 2002/09-2020/08 Q1-16 (-1) Sergio 2002/09-2020/08 Q1-16  (2) Sergio 2002/09-2021/08  Q1-16 (3) Sergio 2002/09-2021/08 Extreme (-3) Sergio 2002/09-2021/08 Mean : [2 = Default] : ');
dataset = 2;
if length(dataset) == 0
  dataset = 2;
end

if dataset == -3
  iQuantile = 00;
elseif dataset ~= 3
  iQuantile = 16;  %% AIRS STM 2021, hottest
  iQuantile = 08;  %% 
  iQuantile = input('Which quantile -1 for extremes [(1--16) (99 for orig Q16, done for AIRS STM)]   [16 = Default] : ');
  if length(iQuantile) == 0
    iQuantile = 16;
  end
else
  iQuantile = [];
end

if dataset == 3
  data_trends = load(['iType_3_extreme_convert_sergio_clearskygrid_obsonly.mat']);
elseif dataset == -3
  data_trends = load(['iType_-3_mean_convert_sergio_clearskygrid_obsonly.mat']);
elseif dataset == 1
  if iQuantile >= 1 & iQuantile <= 16
    data_trends = load(['convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
  elseif iQuantile == 99
    data_trends = load(['convert_sergio_clearskygrid_obsonly_Q16.mat']);
  end
elseif dataset == -1
  if iQuantile >= 1 & iQuantile <= 16
    data_trends = load(['iType_-1_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
  elseif iQuantile == 99
    data_trends = load(['XYZconvert_sergio_clearskygrid_obsonly_Q16.mat']);
  end
elseif dataset == 2
  if iQuantile >= 1 & iQuantile <= 16
    data_trends = load(['iType_2_convert_sergio_clearskygrid_obsonly_Q' num2str(iQuantile,'%02d') '.mat']);
  elseif iQuantile == 99
    data_trends = load(['XYZconvert_sergio_clearskygrid_obsonly_Q16.mat']);
  end
else
  dataset
  error('huh unknown dataset')
end

if ~exist('iaFound')
  clear results*
  iaFound = zeros(1,64);
  existfname = zeros(1,64);

  cdofs = nan(64,66);
  lencdofs = nan(1,64);
  thedofs = nan(1,64);

  results = nan(64,6);
  resultsWV = nan(64,iNumLay);
  resultsT  = nan(64,iNumLay);
  resultsO3 = nan(64,iNumLay);

  resultsunc = nan(64,6);
  resultsWVunc = nan(64,iNumLay);
  resultsTunc  = nan(64,iNumLay);
  resultsO3unc = nan(64,iNumLay);

  rates = nan(2645,64);
  fits  = nan(2645,64);
end

iWarning = 0;
clear fname

for ii = 1 : 64
  if dataset <= 2
    if iNorD > 0
      fname = ['/asl/s1/sergio/Tiles4608/Output_WORKS_May18_2021_Great_AIRS_STM/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here after July 2021
      fname = ['Output/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021, and fornew test comparisons
    elseif iNorD < 0
      fname = ['Output_Day/Quantile' num2str(iQuantile,'%02d') '/test' num2str(ii) '.mat']; %% stored here before before July 2021
    end
  elseif dataset == 3
    if iNorD > 0
      fname = ['Output/Extreme/test' num2str(ii) '.mat']; %% stored here before before July 2021, and fornew test comparisons
    elseif iNorD < 0
      fname = ['Output_Day/Extreme/test' num2str(ii) '.mat']; %% stored here before before July 2021
    end
  elseif dataset == -3
    if iNorD > 0
      fname = ['Output/Quantile00/test' num2str(ii) '.mat']; %% stored here before before July 2021, and fornew test comparisons
    elseif iNorD < 0
      fname = ['Output_Day/Quantile00/test' num2str(ii) '.mat']; %% stored here before before July 2021
    end
  end
  existfname(ii) = exist(fname);
  %fprintf(1,'%5i    %s \n',existfname(ii),fname)
  if exist(fname) > 0 & iaFound(ii) == 0

    iExist = +1;
    loader = ['load ' fname];
    eval(loader);
    iaFound(ii) = +1;
    results(ii,1:6)    = oem.finalrates(1:6);
    resultsunc(ii,1:6) = oem.finalsigs(1:6);
    [mmn,nn] = size(oem.ak_water);

    if length(jacobian.water_i) >= 18
      for iii = 1 : length(jacobian.water_i)
        %iavg = jacobian.wvjaclays_used{iNumLay}-6;
        iavg = jacobian.wvjaclays_used{iii}-jacobian.wvjaclays_offset;
        pavg(iii) = mean(plays(iavg));
      end
    end

    nlays(ii) = nn;
    nn0 = min(nn,iNumLay);
    if nn0 == nn
      thedofs(ii) = oem.dofs;
      lencdofs(ii) = length(oem.cdofs);
      cdofs(ii,1:length(oem.cdofs)) = oem.cdofs;
      resultsWV(ii,1:nn) = oem.finalrates((1:nn)+6+nn*0);
      resultsT(ii,1:nn)  = oem.finalrates((1:nn)+6+nn*1);
      resultsO3(ii,1:nn) = oem.finalrates((1:nn)+6+nn*2);

      resultsWVunc(ii,1:nn) = oem.finalsigs((1:nn)+6+nn*0);
      resultsTunc(ii,1:nn)  = oem.finalsigs((1:nn)+6+nn*1);
      resultsO3unc(ii,1:nn) = oem.finalsigs((1:nn)+6+nn*2);

    else
      iWarning = iWarning + 1
      iaWarning(iWarning) = ii;
      wah = oem.finalrates((1:nn)+6+nn*0);   resultsWV(ii,1:nn0) = wah(1:nn0);
      wah = oem.finalrates((1:nn)+6+nn*1);   resultsT(ii,1:nn0)  = wah(1:nn0);
      wah = oem.finalrates((1:nn)+6+nn*2);   resultsO3(ii,1:nn0) = wah(1:nn0);

      wah = oem.finalsigs((1:nn)+6+nn*0);   resultsWVunc(ii,1:nn0) = wah(1:nn0);
      wah = oem.finalsigs((1:nn)+6+nn*1);   resultsTunc(ii,1:nn0)  = wah(1:nn0);
      wah = oem.finalsigs((1:nn)+6+nn*2);   resultsO3unc(ii,1:nn0) = wah(1:nn0);
    end

    rates(:,ii) = rateset.rates;
    fits(:,ii)  = oem.fit';
  else
    iExist = -1;
    iaFound(ii) = 0;
    results(ii,1:6) = NaN;
    resultsWV(ii,:) = NaN;
    resultsT(ii,:)  = NaN;
    resultsO3(ii,:) = NaN;

    resultsunc(ii,1:6) = NaN;
    resultsWVunc(ii,:) = NaN;
    resultsTunc(ii,:)  = NaN;
    resultsO3unc(ii,:) = NaN;

    rates(:,ii) = NaN;
    fits(:,ii) = NaN;
  end
  if mod(ii,10) == 0
    fprintf(1,'+');
  else
    fprintf(1,'.');
  end
end
fprintf(1,'\n');
resultsunc = real(resultsunc);
resultsWVunc = real(resultsWVunc);
resultsTunc  = real(resultsTunc);
resultsO3unc = real(resultsO3unc);

display('Trying to look at atributes of the last file that should have been read in .... ')
lser = ['!ls -lt ' fname]; eval(lser)

fprintf(1,'found %4i of %4i \n',sum(iaFound),64)

%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 12;  figure(ii); clf; end;

figure(1); clf; semilogy(nanmean(resultsT,1),pavg,'bx-',nanmean(resultsTunc,1),pavg,'c.-');   set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('T and \sigma T')
figure(2); clf; semilogy(nanmean(resultsWV,1),pavg,'bx-',nanmean(resultsWVunc,1),pavg,'c.-'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('WV and \sigma WV')
figure(3); clf; semilogy(nanmean(resultsO3,1),pavg,'bx-',nanmean(resultsO3unc,1),pavg,'c.-'); set(gca,'ydir','reverse'); plotaxis2; ylim([1 1000]); title('O3 and \sigma O3')

%%%%%%%%%%%%%%%%%%%%%%%%%

%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2; 
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

addpath /home/sergio/MATLABCODE/matlib/science/            %% for usgs_deg10_dem.m that has correct paths
[salti, landfrac] = usgs_deg10_dem(Y(:),X(:));
Ylat = Y(:);
Xlon = X(:);
% save landfrac_mask4608.mat landfrac Ylat Xlon rlat65 rlon73 rlon rlat

junk = load('h2645structure.mat');
f           = junk.h.vchan;

i1419 = find(f >= 1419,1);
i1231 = find(f >= 1231,1);
i0900 = find(f >= 0900,1);

figure(1); clf; pcolor(rlat,pavg,resultsWV'); title('dWV/dt'); ylim([100 1000]); shading interp; caxis([-0.15 +0.15]*1/10); colorbar; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); 
figure(2); clf; pcolor(rlat,pavg,resultsT');  title('dTz/dt'); ylim([010 1000]); shading interp; caxis([-0.15 +0.15]*1);    colorbar; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log')
%figure(3); clf; pcolor(rlat,pavg,resultsRH'); title('dRH/dt'); ylim([010 1000]); shading interp; caxis([-0.25 +0.25]*1);    colorbar; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log')
figure(4); clf; pcolor(rlat,pavg,resultsO3'); title('dO3/dt'); ylim([0.1 1000]); shading interp; caxis([-0.1 +0.1]*1/20);   colorbar; colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log')

find_T_WV_RH_03_trends

figure(5); clf; plot(rlat,results(:,1)); title('dCO2/dt');
figure(6); clf; plot(rlat,results(:,2)); title('dN2O/dt');
figure(7); clf; plot(rlat,results(:,3)); title('dCH4/dt');
figure(8); clf; plot(rlat,results(:,6)); title('dST/dt');

figure(9); plot(rlat,thedofs); title('DOFs')
figure(10); boo = find(lencdofs == 66); sumc = mean(cdofs(boo,:),1); 
  plot(cumsum(sumc)); fl = ceil(sum(sumc)); line([6 6],[0 fl],'color','k'); line([26 26],[0 fl],'color','k'); line([46 46],[0 fl],'color','k');
  text(2,10,'TG','fontsize',10);   text(16,10,'WV','fontsize',10);   text(36,10,'T','fontsize',10);   text(56,10,'O3','fontsize',10); xlim([0 66])
for ii = 1 : 20
  ind = (1:5) + (ii-1)*5;
  pjunk20(ii) = mean(plays(ind));
end
cssumc = cumsum(sumc);
pflip20 = fliplr(pjunk20);
wvsumc     = sumc(07:26);         semilogy(wvsumc,pjunk20); set(gca,'ydir','reverse'); ylim([0.01 1000])
wvsumcflip = fliplr(sumc(07:26)); semilogy(wvsumc,pjunk20,'o-',wvsumcflip,pflip20,cumsum(wvsumcflip),pflip20); set(gca,'ydir','reverse'); ylim([0.01 1000])
wvsumcflip = cumsum(fliplr(sumc(07:26)));
tsumcflip  = cumsum(fliplr(sumc(27:46)));
o3sumcflip = cumsum(fliplr(sumc(47:66)));
figure(11); semilogy(wvsumcflip,pflip20,tsumcflip,pflip20,o3sumcflip,pflip20,'linewidth',2)
  set(gca,'ydir','reverse'); ylim([0.050 1000]); hl = legend('WV','T','O3','location','best'); grid; xlabel('DOF'); ylabel('P(mb)')

figure(12); plot(f,nanmean(rates'),f,nanmean(fits')); plotaxis2; xlim([640 1640]); hl = legend('AIRS obs','UMBC fits','location','best'); grid; 

disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../FIND_NWP_MODEL_TRENDS/
plot_ERA_ERA5_AIRSL3_CMIP6_trends

%make_UMBC_ERA5_AIRSL3_CMIP6_trends_64
make_UMBC_ERA5_AIRSL3_CMIP6_trends_64_withunc

do_spectral_closure
figure(1); figure(2); figure(3); figure(21)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
