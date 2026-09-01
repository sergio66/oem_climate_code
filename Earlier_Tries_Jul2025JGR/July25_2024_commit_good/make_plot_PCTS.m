%% clear all
addpath /home/sergio/MATLABCODE/
addpath /asl/matlib/h4tools

cd /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR

if ~exist('results')
  load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX100_50fatlayers_CLIMCAPS_MERRA2_AMIP6_feedback.mat
end

%% see compute_feedbacks_umbc_ecRad.m

indSST    = results(:,6)';
globalSST = nanmean(results(:,6));

px = p;
px.gas_2 = px.gas_2*(1+2.2/400);
px.gas_4 = px.gas_4*(1+1/300);
px.gas_6 = px.gas_6*(1+5/1800);
px.stemp = px.stemp + indSST;
deltaTJUNK = deltaT; bad = find(isnan(deltaTJUNK)); deltaTJUNK(bad) = 0;
px.ptemp = px.ptemp + deltaTJUNK;
fracJUNK = fracWV; bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
  px.gas_1 = px.gas_1 .* (1 + fracJUNK);
fracJUNK = fracO3; bad = find(isnan(fracJUNK)); fracJUNK(bad) = 0;
  px.gas_3 = px.gas_3 .* (1 + fracJUNK);

lats = load('latB64.mat');
tropical_ocean = find(abs(p.rlat) <= 30 & p.landfrac == 0);
tropical_land  = find(abs(p.rlat) <= 30 & p.landfrac == 1);
midlat_ocean   = find(abs(p.rlat) > 30  & abs(p.rlat) <= 60 & p.landfrac == 0);
midlat_land    = find(abs(p.rlat) > 30  & abs(p.rlat) <= 60 & p.landfrac == 1);
polar_ocean    = find(abs(p.rlat) > 60  & p.landfrac == 0);
polar_land     = find(abs(p.rlat) > 60  & p.landfrac == 1);

whos polar_ocean polar_land midlat_ocean midlat_land tropical_ocean tropical_land
100 * (length(polar_ocean) + length(midlat_ocean) + length(tropical_ocean))/4608
100 * (length(tropical_land) + length(tropical_ocean))/4608
100 * (length(midlat_land) + length(midlat_ocean))/4608
100 * (length(polar_land) + length(polar_ocean))/4608

[hnew,pnew]  = replicate_rtp_headprof(h,p,3000,6);
[hnew,pnewx] = replicate_rtp_headprof(h,p,3000,6);

for ind = 1 : 6
  if ind == 1
    oo = tropical_ocean; 
  elseif ind == 2
    oo = tropical_land; 
  elseif ind == 3
    oo = midlat_ocean; 
  elseif ind == 4
    oo = midlat_land; 
  elseif ind == 5
    oo = polar_ocean; 
  elseif ind == 6
    oo = polar_land; 
  end
  pnew.stemp(ind)     = mean(p.stemp(oo));
  pnew.ptemp(:,ind)   = mean(p.ptemp(:,oo),2);
  pnew.gas_1(:,ind)   =  mean(p.gas_1(:,oo),2);
  pnew.gas_2(:,ind)   = mean(p.gas_2(:,oo),2);
  pnew.gas_3(:,ind)   = mean(p.gas_3(:,oo),2);
  pnew.gas_4(:,ind)   = mean(p.gas_4(:,oo),2);
  pnew.gas_5(:,ind)   = mean(p.gas_5(:,oo),2);
  pnew.gas_6(:,ind)   = mean(p.gas_6(:,oo),2);
  pnew.gas_9(:,ind)   = mean(p.gas_9(:,oo),2);
  pnew.gas_12(:,ind)  = mean(p.gas_12(:,oo),2);
  pnew.bttrend(:,ind) = mean(rates(:,oo),2);
  pnew.nemis(ind)     = mean(p.nemis(oo));
  pnew.efreq(:,ind)   = mean(p.efreq(:,oo),2);
  pnew.emis(:,ind)    = mean(p.emis(:,oo),2);
  pnew.rlat(ind)      = mean(abs(p.rlat(oo)));
  pnew.rlon(ind)      = mean(p.rlon(oo));
  pnew.plat(ind)      = mean(abs(p.plat(oo)));
  pnew.plon(ind)      = mean(p.plon(oo));
  
  pnewx.stemp(ind)     = mean(px.stemp(oo));
  pnewx.ptemp(:,ind)   = mean(px.ptemp(:,oo),2);
  pnewx.gas_1(:,ind)   = mean(px.gas_1(:,oo),2);
  pnewx.gas_2(:,ind)   = mean(px.gas_2(:,oo),2);
  pnewx.gas_3(:,ind)   = mean(px.gas_3(:,oo),2);
  pnewx.gas_4(:,ind)   = mean(px.gas_4(:,oo),2);
  pnewx.gas_5(:,ind)   = mean(px.gas_5(:,oo),2);
  pnewx.gas_6(:,ind)   = mean(px.gas_6(:,oo),2);
  pnewx.gas_9(:,ind)   = mean(px.gas_9(:,oo),2);
  pnewx.gas_12(:,ind)  = mean(px.gas_12(:,oo),2);
  pnewx.bttrend(:,ind) = mean(rates(:,oo),2);
  pnewx.nemis(ind)     = mean(px.nemis(oo));
  pnewx.efreq(:,ind)   = mean(px.efreq(:,oo),2);
  pnewx.emis(:,ind)    = mean(px.emis(:,oo),2);
  pnewx.rlat(ind)      = mean(abs(px.rlat(oo)));
  pnewx.rlon(ind)      = mean(px.rlon(oo));
  pnewx.plat(ind)      = mean(abs(px.plat(oo)));
  pnewx.plon(ind)      = mean(px.plon(oo));
end

[hall,pall] = cat_rtp(hnew,pnew,hnew,pnewx);
pall = rmfield(pall,'rcalc');
pall = rmfield(pall,'sarta_rclearcalc');
pall = rmfield(pall,'mmw');
pall = rmfield(pall,'Tw');
pall = rmfield(pall,'Tdew');
pall = rmfield(pall,'RH');
pall = rmfield(pall,'TwSurf');
pall = rmfield(pall,'TdewSurf');
pall = rmfield(pall,'RHSurf');
pall.cfrac = 0 * pall.cfrac;
pall.cngwat = 0 * pall.cngwat;
pall.cprtop = -9999 * ones(size(pall.cfrac));
pall.cprbot = -9999 * ones(size(pall.cfrac));
pall.ctype = -9999 * ones(size(pall.cfrac));
pall.cfrac2 = 0 * pall.cfrac;
pall.cngwat2 = 0 * pall.cngwat;
pall.cprtop2 = -9999 * ones(size(pall.cfrac));
pall.cprbot2 = -9999 * ones(size(pall.cfrac));
pall.ctype2 = -9999 * ones(size(pall.cfrac));

if ~exist('make_plot_PCTS.op.rtp')
  rtpwrite('make_plot_PCTS.op.rtp',hall,[],pall,[]);
end
figure(1); plot(f,pall.bttrend,'c',f,nanmean(pall.bttrend,2),'b'); xlim([640 1640])
  title('AIRS obs rates')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now read in the kCARTA rads
for ii = 1 : 12
  fname = ['KCARTA_SpectralTrends_AllBand_PCTS/allkcbands_prof' num2str(ii) '.mat'];
  a = load(fname);
  kcf = a.wall;
  kcrad(:,ii) = a.dall;
end
boo = find(kcf >= 605);

mu12 = ones(size(pall.rlat)); mu12 = cos(pall.rlat*pi/180);
mu06 = ones(size(pnew.rlat)); mu06 = cos(pnew.rlat*pi/180);
meanStemp = nansum(pall.stemp .* mu12)/sum(mu12);
deltaStemp = nansum((pnewx.stemp - pnew.stemp) .* mu06)/sum(mu06);
mu06_V = ones(length(fcr),1) * mu06;
mu12_V = ones(length(fcr),1) * mu12;
mu_2645_06 = ones(2645,1) * cos(pnew.rlat*pi/180);
mu_2645_12 = ones(2645,1) * cos(pall.rlat*pi/180);

[fc,qc] = convolve_airs(kcf(boo),kcrad(boo,:),double(hall.ichan)); tc = rad2bt(fc,qc); tcrate = tc(:,7:12)-tc(:,1:6);
figure(2); plot(hall.vchan,tcrate,'m',hall.vchan,nansum(mu_2645_06.*tcrate,2)/sum(mu06),'r'); xlim([640 1640])
  title('kCARTA calc rates')

figure(3); plot(f,pnew.bttrend,'c',f,nansum(mu_2645_06.*pnew.bttrend,2)/sum(mu06),'b',hall.vchan,tcrate,'m',hall.vchan,nansum(mu_2645_06.*tcrate,2)/sum(mu06),'r'); xlim([640 1640])
figure(3); plot(f,nansum(mu_2645_06.*pnew.bttrend,2)/sum(mu06),'b',hall.vchan,nansum(mu_2645_06.*tcrate,2)/sum(mu06),'r'); xlim([640 1640])

rad_trend = kcrad(:,7:12)-kcrad(:,1:6);
[fcr,qcr] = quickconvolve(kcf,kcrad,1,1);  tcr = rad2bt(fcr,qcr);
rad_trend = qcr(:,7:12)-qcr(:,1:6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4); 
scr_siz = get(gcf); a0 = scr_siz.Position;
xlimits = [0 max(f)];
ylimits1 = [200 300];
ylimits2 = [-0.1 +0.05];

ta = tiledlayout(2,1,'TileSpacing','Tight', 'Padding','Tight');
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

plotoptions.xstr  = 'Wavenumber cm-1';
plotoptions.ystr1 = 'BT (K)';
plotoptions.ystr2 = 'd(BT)/dt (K/yr-1)';
xstr =  plotoptions.xstr;
ystr1 = plotoptions.ystr1;
ystr2 = plotoptions.ystr2;

tafov(1) = nexttile;
plot(fcr,tcr(:,1:6),'c',fcr,nansum(mu06_V .* tcr(:,1:6),2)/sum(mu06),'b',fcr,ones(size(fcr)) * meanStemp,'k'); 
  line([min(hall.vchan)  min(hall.vchan)+00],[200 300],'color',[1 1 1]*0.6,'linewidth',2,'linestyle','--')
  %line([max(hall.vchan)  max(hall.vchan)+1],[200 300],'color',[1 1 1]*0.6,'linewidth',2,'linestyle','--')
  plotaxis2;
  ylabel(ystr1,'Fontsize',16)
  xlim(xlimits); ylim(ylimits1);
X0 = 0; W = 650; Y0 = 200; H = 100; color = [1 1 1]*0.6; transperancy = 0.3;
hh = patch([X0 X0+W X0+W X0 X0],[Y0 Y0 Y0+H Y0+H Y0],[0 0 0 0 0],'FaceColor',color);
set(hh,'FaceAlpha',transperancy);  %%[0 0.5 1 = transperant inbetween opaque

addpath /home/sergio/MATLABCODE/NANROUTINES/
[f2 y21] = plotAIRS_NANchans(hall.vchan,nansum(mu_2645_06.*tcrate,2)/sum(mu06));
[f2 y22] = plotAIRS_NANchans(f,nansum(mu_2645_06.*pnew.bttrend(:,1:6),2)/sum(mu06));
tafov(2) = nexttile;
plot(fcr,tcr(:,7:12)-tcr(:,1:6),'c',hall.vchan,nansum(mu_2645_06.*tcrate,2)/sum(mu06),'b',f,nansum(mu_2645_06.*pnew.bttrend,2)/sum(mu06),'r'); 
plot(fcr,tcr(:,7:12)-tcr(:,1:6),'c',f2,y21,'b',f2,y22,'r',fcr,ones(size(fcr)) * deltaStemp,'k'); 
  line([min(hall.vchan)  min(hall.vchan)+0],[-0.1 +0.05],'color',[1 1 1]*0.6,'linewidth',2,'linestyle','--')
  %line([max(hall.vchan)  max(hall.vchan)+01],[-0.1 +0.05],'color',[1 1 1]*0.6,'linewidth',2,'linestyle','--')
  plotaxis2;
  ylabel(ystr2,'Fontsize',16)
  xlabel(xstr,'Fontsize',16)
  xlim(xlimits); ylim(ylimits2);
X0 = 0; W = 650; Y0 = -0.1; H = 0.25; color = [1 1 1]*0.6; transperancy = 0.3;
hh = patch([X0 X0+W X0+W X0 X0],[Y0 Y0 Y0+H Y0+H Y0],[0 0 0 0 0],'FaceColor',color);
set(hh,'FaceAlpha',transperancy);  %%[0 0.5 1 = transperant inbetween opaque

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;

% Get rid of all extra space I can
%ta.Padding = 'none';
%ta.TileSpacing = 'none';

% Remove all xtick labels except for 3rd row
for ii = [1]
   tafov(ii).XTickLabel = '';
   tafov(ii).XLabel.String = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

joao = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/areaweight_19years_trop_ocean.txt');
joao = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_for_TileTrends/areaweight_19years.txt');

figure(5); 
scr_siz = get(gcf); a0 = scr_siz.Position;
xlimits = [0 max(f)];
ylimits1 = [0 150];
ylimits2 = [-0.1 +0.05];

ta = tiledlayout(2,1,'TileSpacing','Tight', 'Padding','Tight');
ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

plotoptions.xstr  = 'Wavenumber cm-1';
plotoptions.ystr1 = 'Radiance';
plotoptions.ystr2 = 'd(Radiance)/dt (yr-1)';
xstr =  plotoptions.xstr;
ystr1 = plotoptions.ystr1;
ystr2 = plotoptions.ystr2;

tafov(1) = nexttile;
plot(fcr,qcr(:,1:6),'c',fcr,nansum(mu06_V.*qcr(:,1:6),2)/sum(mu06),'b',fcr,ttorad(fcr,meanStemp),'k')
  line([min(hall.vchan)  min(hall.vchan)+00],[0 150],'color',[1 1 1]*0.6,'linewidth',2,'linestyle','--')
  %line([max(hall.vchan)  max(hall.vchan)_01],[0 150],'color',[1 1 1]*0.6,'linewidth',2,'linestyle','--')
  plotaxis2;
  ylabel(ystr1,'Fontsize',16)
  xlim(xlimits); ylim(ylimits1);
X0 = 0; W = 650; Y0 = 0; H = 150; color = [1 1 1]*0.6; transperancy = 0.3;
hh = patch([X0 X0+W X0+W X0 X0],[Y0 Y0 Y0+H Y0+H Y0],[0 0 0 0 0],'FaceColor',color);
set(hh,'FaceAlpha',transperancy);  %%[0 0.5 1 = transperant inbetween opaque

tafov(2) = nexttile;
plot(fcr,rad_trend,'c',fcr,nansum(mu06_V.*rad_trend,2)/sum(mu06),'b',joao(:,1),joao(:,5),'r',...
     fcr,ttorad(fcr,nansum(mu06.*pnewx.stemp)/sum(mu06))-ttorad(fcr,nansum(mu06.*pnew.stemp)/sum(mu06)),'k'); 
  line([min(hall.vchan)  min(hall.vchan)+00],[-0.1 +0.05],'color',[1 1 1]*0.6,'linewidth',2,'linestyle','--')
  %line([max(hall.vchan)  max(hall.vchan)+01],[-0.1 +0.05],'color',[1 1 1]*0.6,'linewidth',2,'linestyle','--')
  plotaxis2;
  ylabel(ystr2,'Fontsize',16)
  xlabel(xstr,'Fontsize',16)
  xlim(xlimits); ylim(ylimits2);
X0 = 0; W = 650; Y0 = -0.1; H = 0.25; color = [1 1 1]*0.6; transperancy = 0.3;
hh = patch([X0 X0+W X0+W X0 X0],[Y0 Y0 Y0+H Y0+H Y0],[0 0 0 0 0],'FaceColor',color);
set(hh,'FaceAlpha',transperancy);  %%[0 0.5 1 = transperant inbetween opaque

tafov(1).FontSize = 10;
tafov(2).FontSize = 10;

% Get rid of all extra space I can
%ta.Padding = 'none';
%ta.TileSpacing = 'none';

% Remove all xtick labels except for 3rd row
for ii = [1]
   tafov(ii).XTickLabel = '';
   Tafov(ii).XLabel.String = [];
end

i650 = find(fcr >= 650);
boo = nansum(mu06_V.*qcr(:,1:6),2)/sum(mu06);
wooboo = nansum(mu06_V(i650,:).*qcr(i650,1:6),2)/sum(mu06);
fprintf(1,'ratio of OLR(650:3000) : OLR(0:3000) = %8.4f \n',sum(wooboo)/sum(boo))

boo1 = nansum(mu06_V.*qcr(:,1:6),2)/sum(mu06);
boo2 = nansum(mu06_V.*qcr(:,7:12),2)/sum(mu06);
fprintf(1,'OLR = %8.4f deltaOLR = %8.6f W/m2 \n',[sum(boo1)*pi/1000 sum(boo2-boo1)*pi/1000])
i550 = find(fcr <= 550);
i800a = find(fcr >= 550 & fcr <= 800);
i800b = find(fcr > 800);
junk = [sum(boo2-boo1) sum(boo2(i550)-boo1(i550)) sum(boo2(i800a)-boo1(i800a)) sum(boo2(i800b)-boo1(i800b))] *pi/1000;
fprintf(1,'deltaOLR all = %8.6f W/m2  0<f<550 = %8.6f 550<f<800 = %8.6f f>800 = %8.6f\n',junk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6);
if ~exist('airsL3') & ~exist('era5') & ~exist('cmip6')
  error('if needed, run "addpath ../FIND_NWP_MODEL_TRENDS" and "driver_get_the_model_trends"  --- this next line sfrom that code')
end
figure(6); plot(fchanx,nanmean(airsL3.airsL3_spectral_rates'),'b',fchanx,nanmean(era5.era5_spectral_rates'),'r',fchanx,nanmean(cmip6.cmip6_spectral_rates'),'k');
plot(fcr,tcr(:,7:12)-tcr(:,1:6),'c',hall.vchan,nansum(mu_2645_06.*tcrate,2)/sum(mu06),'b',f,nansum(mu_2645_06.*pnew.bttrend,2)/sum(mu06),'r'); 

figure(6); plot(f,nansum(mu_2645_06.*pnew.bttrend,2)/sum(mu06),'k',...
                fchanx,nanmean(airsL3.airsL3_spectral_rates'),'b',fchanx,nanmean(era5.era5_spectral_rates'),'r',fchanx,nanmean(cmip6.cmip6_spectral_rates'),'g');
xlim([640 1640]); ylim([-0.075 +0.05])
plotaxis2;
hl = legend('AIRS obs','AIRS L3','ERA-5','CMIP6','location','best','fontsize',10); xlabel('Wavenumber cm-1'); ylabel('dBT/dt K/yr')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
addpath /asl/matlib/plotutils
figure(4); aslprint('/home/sergio/PAPERS/SPECIALTALKS/Princeton2022/NEWFIGS/kcarta_0_3000__BTtrends.pdf');
figure(5); aslprint('/home/sergio/PAPERS/SPECIALTALKS/Princeton2022/NEWFIGS/kcarta_0_3000__radtrends.pdf');
figure(6); aslprint('/home/sergio/PAPERS/SPECIALTALKS/Princeton2022/NEWFIGS/model_vs_obs_BTtrends.pdf');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see eg driver_olr_fluxchanges_UMBC.m
load /asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX100_50fatlayers_AIRSL3_ERA5_CMIP6_feedback_olr_ceres_UMBC.mat
plot(olr.outflux.umbc_trends.rlat,olr.outflux.umbc_trends.umbc_all)
[mean(olr.outflux.umbc_trends.umbc_all) mean(olr.outflux.umbc_trends.umbc_all(abs(olr.outflux.umbc_trends.rlat)<30))]
