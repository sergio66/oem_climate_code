a20_0 = load('iType_9_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat');
a20_X = load('iType_20_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat');

a23 = load('iType_19_iQAX_3_convert_sergio_clearskygrid_obsonly_Q03.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_XX_YY_from_X_Y

figure(1); clf; pcolor(a23.h.vchan,nanmean(a23.Y,1),squeeze(nanmean(a23.b_desc,1)));   shading interp; colormap(usa2); xlim([645 1645]); colorbar; caxis([-01 +1]*0.1); title('23 years 2002/09-2025/08');
figure(2); clf; pcolor(a23.h.vchan,nanmean(a23.Y,1),squeeze(nanmean(a20_0.b_desc,1))); shading interp; colormap(usa2); xlim([645 1645]); colorbar; caxis([-01 +1]*0.1); title('20 years 2002/09-2022/08');
figure(3); clf; pcolor(a23.h.vchan,nanmean(a23.Y,1),squeeze(nanmean(a20_X.b_desc,1))); shading interp; colormap(usa2); xlim([645 1645]); colorbar; caxis([-01 +1]*0.1); title('20 years 2003/01-2022/12');

x23   = nanmean(squeeze(nanmean(a23.b_desc,1)),1);
x20_0 = nanmean(squeeze(nanmean(a20_0.b_desc,1)),1);
x20_X = nanmean(squeeze(nanmean(a20_X.b_desc,1)),1);

figure(4); clf;
  plot(a23.h.vchan,x20_0,'b',a23.h.vchan,x20_X,'g',a23.h.vchan,x23,'r'); xlim([645 1645])
  plot(a23.h.vchan,x20_X - x20_0,'g',a23.h.vchan,x23 - x20_0,'r',a23.h.vchan,x20_0,'k.-'); xlim([645 1645]) ; plotaxis2;
    title('X-20 usual'); legend('20ryan','23usual','location','best')

