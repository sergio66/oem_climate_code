function [xall,x200,xtrop,xUTLS] = find_umbc_nwp_profile_trends_mmw_strat(h,p,pavg,nwp_spectral_trends,resultsT,resultsWV,resultsO3,resultsST,maskLFmatr,iERAorCMIP6);

addpath /home/sergio/MATLABCODE/TROPOPAUSE/

if nargin == 9
  iERAorCMIP6 = 1;  % use ERA titles
end

if iERAorCMIP6 == 1
  strERAorCMIP6 = 'ERA';
elseif iERAorCMIP6 == 2
  strERAorCMIP6 = 'CMIP6';
  nwp_spectral_trends.era_100_layertrends = nwp_spectral_trends.cmip6_100_layertrends;
end

load('llsmap5.mat');
%load /home/motteler/shome/obs_stats/airs_tiling/latB64.mat
load latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2; 
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

if ~isfield(p,'plays')
  playsN = p.plevs(1:100,:)-p.plevs(2:101,:);
  playsD = log(p.plevs(1:100,:)./p.plevs(2:101,:));
  p.plays = zeros(size(p.plevs));
  p.plays(1:100,:) = playsN./playsD;
end

pert = p;
for ii = 1 : length(p.stemp)
  %%% CHECK THESW ARE FINITE
  
  nlays = p.nlevs(ii)-1;
  playsjunk = p.plays(1:nlays,ii);

  roo = interp1(log(pavg),resultsT(ii,:),log(playsjunk),[],'extrap');
  ohoh = find(~isfinite(roo));
  if length(ohoh) > 0
    fprintf(1,' OHOH profile %4i dT(z) has %3i not finite members, set to 0 \n',ii,length(ohoh));
    roo(ohoh) = 0;
  end
  pert.ptemp(1:nlays,ii) =  pert.ptemp(1:nlays,ii) + roo;

  roo = interp1(log(pavg),resultsWV(ii,:),log(playsjunk),[],'extrap');
  ohoh = find(~isfinite(roo));
  if length(ohoh) > 0
    fprintf(1,' OHOH profile %4i dWV(z) has %3i not finite members, set to 0 \n',ii,length(ohoh));
    roo(ohoh) = 0;
  end
  pert.gas_1(1:nlays,ii) =  pert.gas_1(1:nlays,ii) .* (1+roo);

  roo = interp1(log(pavg),resultsO3(ii,:),log(playsjunk),[],'extrap');
  ohoh = find(~isfinite(roo));
  if length(ohoh) > 0
    fprintf(1,' OHOH profile %4i dO3(z) has %3i not finite members, set to 0 \n',ii,length(ohoh));
    roo(ohoh) = 0;
  end
  pert.gas_3(1:nlays,ii) =  pert.gas_3(1:nlays,ii) .* (1+roo);

  if ~isfinite(resultsST(ii))
    fprintf(1,' OHOH profile %4i dST is not finite, set to 0 \n',ii);
  else
    pert.stemp(ii) = pert.stemp(ii) + resultsST(ii);
  end

  pert.gas_2(1:nlays,ii) =  pert.gas_2(1:nlays,ii) .* (1+2.2/385);
  pert.gas_4(1:nlays,ii) =  pert.gas_4(1:nlays,ii) .* (1+0.8/300);
  pert.gas_6(1:nlays,ii) =  pert.gas_6(1:nlays,ii) .* (1+4.5/1700);

end

pera = p;
pera.stemp = pera.stemp + nwp_spectral_trends.era_100_layertrends.stemp;
pera.ptemp(1:100,:) = pera.ptemp(1:100,:) + nwp_spectral_trends.era_100_layertrends.ptemp;
pera.gas_1(1:100,:) = pera.gas_1(1:100,:) .* (1 + nwp_spectral_trends.era_100_layertrends.gas_1);
pera.gas_3(1:100,:) = pera.gas_3(1:100,:) .* (1 + nwp_spectral_trends.era_100_layertrends.gas_3);
  pera.gas_2 =  pera.gas_2 * (1+2.2/385);
  pera.gas_4 =  pera.gas_4 * (1+0.8/300);
  pera.gas_6 =  pera.gas_6 * (1+4.5/1700);

pera5 = p;
pera5.stemp = pera5.stemp + nwp_spectral_trends.era5_100_layertrends.stemp;
pera5.ptemp(1:100,:) = pera5.ptemp(1:100,:) + nwp_spectral_trends.era5_100_layertrends.ptemp;
pera5.gas_1(1:100,:) = pera5.gas_1(1:100,:) .* (1 + nwp_spectral_trends.era5_100_layertrends.gas_1);
pera5.gas_3(1:100,:) = pera5.gas_3(1:100,:) .* (1 + nwp_spectral_trends.era5_100_layertrends.gas_3);
  pera5.gas_2 =  pera5.gas_2 * (1+2.2/385);
  pera5.gas_4 =  pera5.gas_4 * (1+0.8/300);
  pera5.gas_6 =  pera5.gas_6 * (1+4.5/1700);

pairsL3 = p;
pairsL3.stemp = pairsL3.stemp + nwp_spectral_trends.airsL3_100_layertrends.stemp;
pairsL3.ptemp(1:100,:) = pairsL3.ptemp(1:100,:) + nwp_spectral_trends.airsL3_100_layertrends.ptemp;
pairsL3.gas_1(1:100,:) = pairsL3.gas_1(1:100,:) .* (1 + nwp_spectral_trends.airsL3_100_layertrends.gas_1);
pairsL3.gas_3(1:100,:) = pairsL3.gas_3(1:100,:) .* (1 + nwp_spectral_trends.airsL3_100_layertrends.gas_3);
  pairsL3.gas_2 =  pairsL3.gas_2 * (1+2.2/385);
  pairsL3.gas_4 =  pairsL3.gas_4 * (1+0.8/300);
  pairsL3.gas_6 =  pairsL3.gas_6 * (1+4.5/1700);

xall.mmw0      = mmwater_rtp(h,p);
xall.mmwUMBC   = mmwater_rtp(h,pert);
xall.mmwERA    = mmwater_rtp(h,pera);
xall.mmwERA5   = mmwater_rtp(h,pera5);
xall.mmwAIRSL3 = mmwater_rtp(h,pairsL3);

figure(1); scatter_coast(p.rlon,p.rlat,50,xall.mmw0); title('mmw0 to surface'); colormap jet
figure(2); scatter_coast(p.rlon,p.rlat,50,xall.mmwUMBC-xall.mmw0);   title('dmmwUMBC/dt');   colormap(usa2); caxis([-0.25 +0.25]);
figure(3); scatter_coast(p.rlon,p.rlat,50,xall.mmwERA-xall.mmw0);    title('dmmwERA/dt');    colormap(usa2); caxis([-0.25 +0.25])
figure(4); scatter_coast(p.rlon,p.rlat,50,xall.mmwERA5-xall.mmw0);   title('dmmwERA5/dt');   colormap(usa2); caxis([-0.25 +0.25])
figure(5); scatter_coast(p.rlon,p.rlat,50,xall.mmwAIRSL3-xall.mmw0); title('dmmwAIRSL3/dt'); colormap(usa2); caxis([-0.25 +0.25])
if iERAorCMIP6 == 2
  figure(3); title('d(mmwCMIP6)/dt');
end

figure(1); aslmap(1,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xall.mmw0,72,64)'),1),            [-90 +90],[-180 +180]); colormap(jet);      title('mmw to surface');  caxis([0 70])
figure(2); aslmap(2,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xall.mmwUMBC-xall.mmw0,72,64)'),1),  [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwUMBC)/dt');   caxis([-0.25 +0.25])
figure(3); aslmap(3,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xall.mmwERA-xall.mmw0,72,64)'),1),   [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwERA)/dt');    caxis([-0.25 +0.25])
figure(4); aslmap(4,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xall.mmwERA5-xall.mmw0,72,64)'),1),  [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwERA5)/dt');   caxis([-0.25 +0.25])
figure(5); aslmap(5,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xall.mmwAIRSL3-xall.mmw0,72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwAIRSL3)/dt'); caxis([-0.25 +0.25])
if iERAorCMIP6 == 2
  figure(3); title('d(mmwCMIP6)/dt');
end

figure(6); clf; 
  %junk = maskLFmatr.*smoothn((reshape(xall.mmw0,72,64)'),1);             plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xall.mmwUMBC-xall.mmw0,72,64)'),1);    plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xall.mmwAIRSL3-xall.mmw0,72,64)'),1);  plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xall.mmwERA-xall.mmw0,72,64)'),1);     plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xall.mmwERA5-xall.mmw0,72,64)'),1);    plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  hold off; plotaxis2;
  hl = legend('UMBC','AIRSL3',strERAorCMIP6,'ERA5','location','best','fontsize',10); xlabel('Latitude'); ylabel('Mean ColumnWater d(mmw)/dt'); title('To Surface')

%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to continue to colum water to 200 mb'); pause
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
x200.mmw0      = mmwater_rtp(h,p,200);
x200.mmwUMBC   = mmwater_rtp(h,pert,200);
x200.mmwERA    = mmwater_rtp(h,pera,200);
x200.mmwERA5   = mmwater_rtp(h,pera5,200);
x200.mmwAIRSL3 = mmwater_rtp(h,pairsL3,200);

figure(1); scatter_coast(p.rlon,p.rlat,50,x200.mmw0); title('mmw0 to 200 mb'); colormap jet
figure(2); scatter_coast(p.rlon,p.rlat,50,x200.mmwUMBC-x200.mmw0);   title('dmmwUMBC/dt');   colormap(usa2); caxis([-0.00025 +0.00025]);
figure(3); scatter_coast(p.rlon,p.rlat,50,x200.mmwERA-x200.mmw0);    title('dmmwERA/dt');    colormap(usa2); caxis([-0.00025 +0.00025])
figure(4); scatter_coast(p.rlon,p.rlat,50,x200.mmwERA5-x200.mmw0);   title('dmmwERA5/dt');   colormap(usa2); caxis([-0.00025 +0.00025])
figure(5); scatter_coast(p.rlon,p.rlat,50,x200.mmwAIRSL3-x200.mmw0); title('dmmwAIRSL3/dt'); colormap(usa2); caxis([-0.00025 +0.00025])
if iERAorCMIP6 == 2
  figure(3); title('d(mmwCMIP6)/dt');
end

figure(1); aslmap(1,rlat65,rlon73,maskLFmatr.*smoothn((reshape(x200.mmw0,72,64)'),1),          [-90 +90],[-180 +180]); colormap(jet);      title('mmw to 200 mb');  caxis([0 0.05])
figure(2); aslmap(2,rlat65,rlon73,maskLFmatr.*smoothn((reshape(x200.mmwUMBC-x200.mmw0,72,64)'),1),  [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwUMBC)/dt');   caxis([-0.00025 +0.00025])
figure(3); aslmap(3,rlat65,rlon73,maskLFmatr.*smoothn((reshape(x200.mmwERA-x200.mmw0,72,64)'),1),   [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwERA)/dt');    caxis([-0.00025 +0.00025])
figure(4); aslmap(4,rlat65,rlon73,maskLFmatr.*smoothn((reshape(x200.mmwERA5-x200.mmw0,72,64)'),1),  [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwERA5)/dt');   caxis([-0.00025 +0.00025])
figure(5); aslmap(5,rlat65,rlon73,maskLFmatr.*smoothn((reshape(x200.mmwAIRSL3-x200.mmw0,72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwAIRSL3)/dt'); caxis([-0.00025 +0.00025])
if iERAorCMIP6 == 2
  figure(3); title('d(mmwCMIP6)/dt');
end

figure(6); clf; 
  %junk = maskLFmatr.*smoothn((reshape(x200.mmw0,72,64)'),1);            plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(x200.mmwUMBC-x200.mmw0,72,64)'),1);    plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(x200.mmwAIRSL3-x200.mmw0,72,64)'),1);  plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(x200.mmwERA-x200.mmw0,72,64)'),1);     plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(x200.mmwERA5-x200.mmw0,72,64)'),1);    plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  hold off; plotaxis2;
  hl = legend('UMBC','AIRSL3',strERAorCMIP6,'ERA5','location','best','fontsize',10); xlabel('Latitude'); ylabel('Mean ColumnWater d(mmw)/dt'); title('To 200 mb')

%%%%%%%%%%%%%%%%%%%%%%%%%

disp('ret to continue to colum water to tropopause'); pause
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
xtrop.mmw0      = mmwater_rtp(h,p,-1);
xtrop.mmwUMBC   = mmwater_rtp(h,pert,-1);
xtrop.mmwERA    = mmwater_rtp(h,pera,-1);
xtrop.mmwERA5   = mmwater_rtp(h,pera5,-1);
xtrop.mmwAIRSL3 = mmwater_rtp(h,pairsL3,-1);

figure(1); scatter_coast(p.rlon,p.rlat,50,xtrop.mmw0); title('mmw0 to tropopause'); colormap jet
figure(2); scatter_coast(p.rlon,p.rlat,50,xtrop.mmwUMBC-xtrop.mmw0);   title('dmmwUMBC/dt');   colormap(usa2); caxis([-0.05 +0.05]);
figure(3); scatter_coast(p.rlon,p.rlat,50,xtrop.mmwERA-xtrop.mmw0);    title('dmmwERA/dt');    colormap(usa2); caxis([-0.05 +0.05])
figure(4); scatter_coast(p.rlon,p.rlat,50,xtrop.mmwERA5-xtrop.mmw0);   title('dmmwERA5/dt');   colormap(usa2); caxis([-0.05 +0.05])
figure(5); scatter_coast(p.rlon,p.rlat,50,xtrop.mmwAIRSL3-xtrop.mmw0); title('dmmwAIRSL3/dt'); colormap(usa2); caxis([-0.05 +0.05])
if iERAorCMIP6 == 2
  figure(3); title('d(mmwCMIP6)/dt');
end

figure(1); aslmap(1,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xtrop.mmw0,72,64)'),1),          [-90 +90],[-180 +180]); colormap(jet);      title('mmw to tropopause'); caxis([0 0.01])
figure(2); aslmap(2,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xtrop.mmwUMBC-xtrop.mmw0,72,64)'),1),  [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwUMBC)/dt');     caxis([-0.01 +0.01]/100)
figure(3); aslmap(3,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xtrop.mmwERA-xtrop.mmw0,72,64)'),1),   [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwERA)/dt');      caxis([-0.01 +0.01]/100)
figure(4); aslmap(4,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xtrop.mmwERA5-xtrop.mmw0,72,64)'),1),  [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwERA5)/dt');     caxis([-0.01 +0.01]/100)
figure(5); aslmap(5,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xtrop.mmwAIRSL3-xtrop.mmw0,72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwAIRSL3)/dt');   caxis([-0.01 +0.01]/100)
if iERAorCMIP6 == 2
  figure(3); title('d(mmwCMIP6)/dt');
end

figure(6); clf; 
  %junk = maskLFmatr.*smoothn((reshape(xtrop.mmw0,72,64)'),1);            plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xtrop.mmwUMBC-xtrop.mmw0,72,64)'),1);    plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xtrop.mmwAIRSL3-xtrop.mmw0,72,64)'),1);  plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xtrop.mmwERA-xtrop.mmw0,72,64)'),1);     plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xtrop.mmwERA5-xtrop.mmw0,72,64)'),1);    plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  hold off; plotaxis2;
  hl = legend('UMBC','AIRSL3',strERAorCMIP6,'ERA5','location','best','fontsize',10); xlabel('Latitude'); ylabel('Mean ColumnWater d(mmw)/dt'); title('To Tropopause')

%%%%%%%%%%%%%%%%%%%%%%%%%

% https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2011RG000355 THE
% EXTRATROPICAL UPPER TROPOSPHERE AND LOWER STRATOSPHERE A. Gettelman,
% P. Hoor, L. L. Pan, W. J. Randel, M. I. Hegglin, T. Birner First
% published: 09 August 2011
% https://doi.org/10.1029/2011RG000355Citations: 224 
%
% The extratropical upper troposphere and lower stratosphere (Ex-UTLS)
% is a transition region between the stratosphere and the
% troposphere. The Ex-UTLS includes the tropopause, a strong static
% stability gradient and dynamic barrier to transport.
%
% The upper troposphere and lower stratosphere (UTLS) is a coupling
% layer in the atmosphere. It can be broadly defined as the region ±5
% km around the tropopause, the traditional boundary between the
% troposphere (from the Greek τρequation imageπω or “to turn over”)
% and the stratosphere (from the Latin stratus, “to spread out”). The
% UTLS is a consequence of the transition between the troposphere and
% stratosphere, and processes in the region may alter both the
% troposphere and stratosphere. Stratosphere-troposphere exchange
% (STE) across the tropopause is an important bidirectional process
% for influencing the chemistry of the upper troposphere and lower
% stratosphere [Holton et al., 1995]. STE is important for
% understanding tropospheric ozone (O3) concentrations that affect air
% quality. But the UTLS is important for more than just STE.

disp('ret to continue to colum water to upper troposphere '); pause
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
xUTLS.mmw0      = mmwater_rtp(h,p,-2);
xUTLS.mmwUMBC   = mmwater_rtp(h,pert,-2);
xUTLS.mmwERA    = mmwater_rtp(h,pera,-2);
xUTLS.mmwERA5   = mmwater_rtp(h,pera5,-2);
xUTLS.mmwAIRSL3 = mmwater_rtp(h,pairsL3,-2);

figure(1); scatter_coast(p.rlon,p.rlat,50,xUTLS.mmw0); title('mmw0 to UT'); colormap jet
figure(2); scatter_coast(p.rlon,p.rlat,50,xUTLS.mmwUMBC-xUTLS.mmw0);   title('dmmwUMBC/dt');   colormap(usa2); caxis([-0.05 +0.05]);
figure(3); scatter_coast(p.rlon,p.rlat,50,xUTLS.mmwERA-xUTLS.mmw0);    title('dmmwERA/dt');    colormap(usa2); caxis([-0.05 +0.05])
figure(4); scatter_coast(p.rlon,p.rlat,50,xUTLS.mmwERA5-xUTLS.mmw0);   title('dmmwERA5/dt');   colormap(usa2); caxis([-0.05 +0.05])
figure(5); scatter_coast(p.rlon,p.rlat,50,xUTLS.mmwAIRSL3-xUTLS.mmw0); title('dmmwAIRSL3/dt'); colormap(usa2); caxis([-0.05 +0.05])
if iERAorCMIP6 == 2
  figure(3); title('d(mmwCMIP6)/dt');
end

figure(1); aslmap(1,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xUTLS.mmw0,72,64)'),1),          [-90 +90],[-180 +180]); colormap(jet);      title('mmw to UT');       caxis([0 0.1])
figure(2); aslmap(2,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xUTLS.mmwUMBC-xUTLS.mmw0,72,64)'),1),  [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwUMBC)/dt');   caxis([-0.05 +0.05]/10)
figure(3); aslmap(3,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xUTLS.mmwERA-xUTLS.mmw0,72,64)'),1),   [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwERA)/dt');    caxis([-0.05 +0.05]/10)
figure(4); aslmap(4,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xUTLS.mmwERA5-xUTLS.mmw0,72,64)'),1),  [-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwERA5)/dt');   caxis([-0.05 +0.05]/10)
figure(5); aslmap(5,rlat65,rlon73,maskLFmatr.*smoothn((reshape(xUTLS.mmwAIRSL3-xUTLS.mmw0,72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  title('d(mmwAIRSL3)/dt'); caxis([-0.05 +0.05]/10)
if iERAorCMIP6 == 2
  figure(3); title('d(mmwCMIP6)/dt');
end

figure(6); clf; 
  %junk = maskLFmatr.*smoothn((reshape(xUTLS.mmw0,72,64)'),1);            plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xUTLS.mmwUMBC-xUTLS.mmw0,72,64)'),1);    plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xUTLS.mmwAIRSL3-xUTLS.mmw0,72,64)'),1);  plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xUTLS.mmwERA-xUTLS.mmw0,72,64)'),1);     plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  junk = maskLFmatr.*smoothn((reshape(xUTLS.mmwERA5-xUTLS.mmw0,72,64)'),1);    plot(rlat,nanmean(junk,2),'linewidth',2); hold on;
  hold off; plotaxis2;
  hl = legend('UMBC','AIRSL3',strERAorCMIP6,'ERA5','location','best','fontsize',10); xlabel('Latitude'); ylabel('Mean ColumnWater d(mmw)/dt'); title('To UpperTrop')

disp('ret to exit find_umbc_nwp_profile_trends_mmw_strat.m'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
