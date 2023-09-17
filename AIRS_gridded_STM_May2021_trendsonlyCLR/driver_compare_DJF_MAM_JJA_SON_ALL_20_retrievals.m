%% cp /home/sergio/PAPERS/SUBMITPAPERS/trends/Matlab/contrast_05_10_15_20.m compare_05_10_15_20_spectral_rates.m

addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/plotutils
addpath /asl/matlib/maps

load llsmap5
clear plotoptions

if ~exist('aALL')
  dir0 = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_GULP_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2';

  aDJF = load([dir0 '_DJF.mat']);
  aMAM = load([dir0 '_MAM.mat']);
  aJJA = load([dir0 '_JJA.mat']);
  aSON = load([dir0 '_SON.mat']);
  aALL = load([dir0 '.mat']);
  
  loader = ['load h2645structure.mat'];
  eval(loader)

  [h,ha,p,pa] = rtpread(['summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_PERTv1.rtp']);
end
  
rlat = p.rlat;
rlat = cos(rlat*pi/180);

load latB64.mat
rlon73 = -180 : 5 : +180;  rlat65 = latB2;
rlon = 0.5*(rlon73(1:end-1)+rlon73(2:end));
rlat = 0.5*(rlat65(1:end-1)+rlat65(2:end));
[Y,X] = meshgrid(rlat,rlon);
X = X; Y = Y;

plays = aDJF.pert.plays(:,1000);
pavgLAY = aDJF.pert.plays;
iFig = 6;
rlat100x64 = ones(100,1)*rlat';
plays100x64 = plays(1:100)*ones(1,64);

plotoptions.str11 = 'DJF';
plotoptions.str12 = 'MAM';
plotoptions.str21 = 'JJA';
plotoptions.str22 = 'SON';
plotoptions.xstr = ' ';        plotoptions.ystr = ' ';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = 0;
iFig = iFig + 1; mooDJF = squeeze(aDJF.results(:,6)); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooDJF,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('DJF BT1231 trend')
iFig = iFig + 1; mooMAM = squeeze(aMAM.results(:,6)); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooMAM,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('MAM BT1231 trend')
iFig = iFig + 1; mooJJA = squeeze(aJJA.results(:,6)); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooJJA,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('JJA BT1231 trend')
iFig = iFig + 1; mooSON = squeeze(aSON.results(:,6)); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooSON,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('SON BT1231 trend')
iFig = iFig + 1; mooALL = squeeze(aALL.results(:,6)); aslmap(iFig,rlat65,rlon73,smoothn((reshape(mooALL,72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('ALL BT1231 trend')
aXYZ.results = (aDJF.results + aMAM.results + aJJA.results + aSON.results)/4;
iFig = iFig + 1; mooALL = squeeze(aXYZ.results(:,6)); aslmap(iFig,rlat65,rlon73,smoothn((reshape(aXYZ.results(:,6),72,64)') ,1), [-90 +90],[-180 +180]); colormap(usa2); caxis([-1 +1]*0.15); title('<SEASON> BT1231 trend')

iFig = iFig + 1; 
plotoptions.yLinearOrLog = +1;
plotoptions.cx = [-1 +1]*0.151; plotoptions.maintitle = 'dST/dt'; plotoptions.plotcolors = llsmap5;
aslmap_2x2tiledlayout(mooDJF,mooMAM,mooJJA,mooSON,iFig,plotoptions);
disp('ret to continue'); pause

%iFig = iFig + 1; 
%  mooALL = squeeze(aXYZ.results(:,6))./squeeze(aALL.results(:,6)); 
%    figure(iFig); simplemap_thickcoast(Y,X,mooALL-1,5); colormap(usa2); caxis([-1 +1]*0.5); title('1- ALL<4>/ALL BT1231 trend')
%  mooALL = squeeze(aXYZ.results(:,6)) - squeeze(aALL.results(:,6)); 
%    figure(iFig); simplemap_thickcoast(Y,X,mooALL,5); colormap(usa2); caxis([-1 +1]*0.015); title('1- ALL<4>/ALL BT1231 trend')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = iFig + 1; figure(iFig); clf; 
  moo = aDJF.deltaT; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('DJF dT/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aMAM.deltaT; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('MAM dT/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aJJA.deltaT; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('JJA dT/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aSON.deltaT; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('SON dT/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aALL.deltaT; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('ALL dT/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = 1/4*(aDJF.deltaT + aMAM.deltaT + aJJA.deltaT + aSON.deltaT); moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('<SEASONAL> dT/dt');

iFig = iFig + 1; 
plotoptions.yLinearOrLog = -1;
plotoptions.yReverseDir  = +1;
plotoptions.yLimits = [10 1000];
plotoptions.xaxis_metric = 'sine';

plotoptions.maintitle = 'dT/dt trends [K/yr]';
plotoptions.cx = [-1 +1]*0.151;
miaow11 = squeeze(nanmean(reshape(aDJF.deltaT,101,72,64),2)); 
miaow12 = squeeze(nanmean(reshape(aMAM.deltaT,101,72,64),2)); 
miaow21 = squeeze(nanmean(reshape(aJJA.deltaT,101,72,64),2)); 
miaow22 = squeeze(nanmean(reshape(aSON.deltaT,101,72,64),2)); 
profile_plots_2x2tiledlayout(rlat,pavgLAY(1:97,1000),miaow11(1:97,:),miaow12(1:97,:),miaow21(1:97,:),miaow22(1:97,:),iFig,plotoptions);

disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = iFig + 1; figure(iFig); clf; 
  moo = aDJF.fracWV; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('DJF d(fracWV)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aMAM.fracWV; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('MAM d(fracWV)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aJJA.fracWV; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('JJA d(fracWV)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aSON.fracWV; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('SON d(fracWV)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aALL.fracWV; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('ALL d(fracWV)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = 1/4*(aDJF.fracWV + aMAM.fracWV + aJJA.fracWV + aSON.fracWV); moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','linear'); caxis([-1 +1]*0.15/10); shading interp; ylim([10 1000]); colormap(usa2); colorbar; title('<SEASONAL> d(fracWV)/dt');

iFig = iFig + 1; 
plotoptions.yLinearOrLog = +1;
plotoptions.yReverseDir  = +1;
plotoptions.yLimits = [100 1000];
plotoptions.xaxis_metric = 'sine';

plotoptions.maintitle = 'dfracWV/dt trends [1/yr]';
plotoptions.cx = [-1 +1]*0.151/10;
miaow11 = squeeze(nanmean(reshape(aDJF.fracWV,101,72,64),2)); 
miaow12 = squeeze(nanmean(reshape(aMAM.fracWV,101,72,64),2)); 
miaow21 = squeeze(nanmean(reshape(aJJA.fracWV,101,72,64),2)); 
miaow22 = squeeze(nanmean(reshape(aSON.fracWV,101,72,64),2)); 
profile_plots_2x2tiledlayout(rlat,pavgLAY(1:97,1000),miaow11(1:97,:),miaow12(1:97,:),miaow21(1:97,:),miaow22(1:97,:),iFig,plotoptions);

disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = iFig + 1; figure(iFig); clf; 
  moo = aDJF.fracO3; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('DJF d(fracO3)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aMAM.fracO3; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('MAM d(fracO3)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aJJA.fracO3; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('JJA d(fracO3)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aSON.fracO3; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('SON d(fracO3)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = aALL.fracO3; moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('ALL d(fracO3)/dt');
iFig = iFig + 1; figure(iFig); clf; 
  moo = 1/4*(aDJF.fracO3 + aMAM.fracO3 + aJJA.fracO3 + aSON.fracO3); moo = squeeze(nanmean(reshape(moo,101,72,64),2)); 
  pcolor(rlat,pavgLAY(1:97,1000),smoothn(moo(1:97,:),1)); set(gca,'ydir','reverse'); set(gca,'yscale','log'); caxis([-1 +1]*0.15/10); shading interp; ylim([0.01 100]); colormap(usa2); colorbar; title('<SEASONAL> d(fracO3)/dt');

iFig = iFig + 1; 
plotoptions.yLinearOrLog = -1;
plotoptions.yReverseDir  = +1;
plotoptions.yLimits = [0.010 100];
plotoptions.xaxis_metric = 'sine';

plotoptions.maintitle = 'dfracO3/dt trends [1/yr]';
plotoptions.cx = [-1 +1]*0.151/10;
miaow11 = squeeze(nanmean(reshape(aDJF.fracO3,101,72,64),2)); 
miaow12 = squeeze(nanmean(reshape(aMAM.fracO3,101,72,64),2)); 
miaow21 = squeeze(nanmean(reshape(aJJA.fracO3,101,72,64),2)); 
miaow22 = squeeze(nanmean(reshape(aSON.fracO3,101,72,64),2)); 
profile_plots_2x2tiledlayout(rlat,pavgLAY(1:97,1000),miaow11(1:97,:),miaow12(1:97,:),miaow21(1:97,:),miaow22(1:97,:),iFig,plotoptions);


disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iPrint = -1;
if iPrint > 0
  figure(07); aslprint('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/skt_trends_lat_lon_4panels20_years_seasonal.pdf');
  figure(14); aslprint('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/tz_trends_zonal_p_4panels20_years_seasonal.pdf')
  figure(21); aslprint('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/fracWV_trends_zonal_p_4panels20_years_seasonal.pdf');
  figure(28); aslprint('/home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/fracO3_trends_zonal_p_4panels20_years_seasonal.pdf');
end
