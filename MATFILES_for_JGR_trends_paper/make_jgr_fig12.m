load fig12.mat

iFig = 1;
figure(iFig); clf
   profile_plots_2x1x2tiledlayout_tall(fig12_WVtrends.rlat,fig12_WVtrends.plays100,fig12_WVtrends.umbc,fig12_WVtrends.airsL3,fig12_WVtrends.climcaps,fig12_WVtrends.merra2,fig12_WVtrends.era5,iFig,fig12_WVtrends.plotoptions2x1x2);
