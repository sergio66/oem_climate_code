load fig11.mat

iFig = 1;
figure(iFig); clf
   profile_plots_2x1x2tiledlayout_tall(fig11_Ttrends.rlat,fig11_Ttrends.plays100,fig11_Ttrends.umbc,fig11_Ttrends.airsL3,fig11_Ttrends.climcaps,fig11_Ttrends.merra2,fig11_Ttrends.era5,iFig,fig11_Ttrends.plotoptions2x1x2);
