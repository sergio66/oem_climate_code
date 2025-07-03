load fig10.mat

figure(1); clf

plot(fig10_colwvtrend.rlat,fig10_colwvtrend.umbc,'k',fig10_colwvtrend.rlat,fig10_colwvtrend.airsL3,'b',fig10_colwvtrend.rlat,fig10_colwvtrend.climcaps,'g',...
     fig10_colwvtrend.rlat,fig10_colwvtrend.era5,'r',fig10_colwvtrend.rlat,fig10_colwvtrend.merra2,'m','linewidth',2);

ylim([-1 +2]*0.05);
xlim([-1 +1]*90);
plotaxis2; hl = legend('AIRS\_RT','AIRS L3','CLIMCAPS L3','ERA5','MERRA2','location','best','fontsize',8);
xlabel('Latitude [deg]'); ylabel('d mmw/dt [mm/yr]');
