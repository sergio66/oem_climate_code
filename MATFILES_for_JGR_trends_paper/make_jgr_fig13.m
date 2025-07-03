load fig13.mat

load llsmap5

  aslmap(1,fig13_400mb.rlat65,fig13_400mb.rlon73,fig13_400mb.umbc_wv_400,[-90 +90],[-180 +180]); caxis([-1 +1]*0.015); colormap(llsmap5);
    set(gca,'fontsize',18)
    text(-0.25,-1.625,'yr^{-1}','fontsize',14)
  aslmap(2,fig13_400mb.rlat65,fig13_400mb.rlon73,fig13_400mb.era5_wv_400,[-90 +90],[-180 +180]); caxis([-1 +1]*0.015); colormap(llsmap5);
    set(gca,'fontsize',18)
    text(-0.25,-1.625,'yr^{-1}','fontsize',14)
