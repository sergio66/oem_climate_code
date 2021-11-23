llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat');

figure(1); clf
figname = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/SAVE_BESTRUNv1/Figs/era_400mb_tz_anom_nosmooth.fig';
%figname = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/SAVE_BESTRUNv1/Figs/era_400mb_tz_anom_nosmooth_AK.fig';
[tera,yera,cera] = get_fig_image_data(figname); 
shading interp; caxis([-5 +5]); colorbar; colormap(llsmap4.llsmap4);
title('ERA')

%%%%%%%%%%

figure(2); clf
figname = 'SAVE_LW_noCFC11_noN2O/Figs/obs_tz_anom_400mb.fig';
[t2,y2,c2] = get_fig_image_data(figname);
shading interp; caxis([-5 +5]); colorbar; colormap(llsmap4.llsmap4);
title('cov and reg,obs')

figure(2); clf
figname = 'SAVE_LW_noCFC11_noN2O/Figs/cal_tz_anom_400mb.fig';
[cal_t2,cal_y2,cal_c2] = get_fig_image_data(figname);
shading interp; caxis([-5 +5]); colorbar; colormap(llsmap4.llsmap4);
title('cov and reg,obs')

figure(3); clf
figname = 'SAVE_LW_noCFC11_noN2O_covx10/Figs/obs_tz_anom_400mb.fig';
[t2_10,y2_10,c2_10] = get_fig_image_data(figname);
shading interp; caxis([-5 +5]); colorbar; colormap(llsmap4.llsmap4);
title('cov and regx10,obs')

figure(3); clf
figname = 'SAVE_LW_noCFC11_noN2O_covx10/Figs/cal_tz_anom_400mb.fig';
[cal_t2_10,cal_y2_10,cal_c2_10] = get_fig_image_data(figname);
shading interp; caxis([-5 +5]); colorbar; colormap(llsmap4.llsmap4);
title('cov and regx10,cal')

figure(4); clf
figname = 'SAVE_LW_noCFC11_noN2O_cov1_noreg/Figs/obs_tz_anom_400mb.fig';
[t3,y3,c3] = get_fig_image_data(figname);
shading interp; caxis([-5 +5]); colorbar; colormap(llsmap4.llsmap4);
title('no reg, cov x1')

figure(5); clf
figname = 'SAVE_LW_noCFC11_noN2O_cov10_noreg/Figs/obs_tz_anom_400mb.fig';
[t4,y4,c4] = get_fig_image_data(figname);
shading interp; caxis([-5 +5]); colorbar; colormap(llsmap4.llsmap4);
title('no reg, cov x10')

figure(6); clf
figname = 'SAVE_LW_noCFC11_noN2O_reg/Figs/obs_tz_anom_400mb.fig';
[t5,y5,c5] = get_fig_image_data(figname);
shading interp; caxis([-5 +5]); colorbar; colormap(llsmap4.llsmap4);
title('no cov, reg x1')

figure(6); clf
figname = 'SAVE_LW_noCFC11_noN2O_qrenorm=1/Figs/obs_tz_anom_400mb.fig';
[t6,y6,c6] = get_fig_image_data(figname);
shading interp; caxis([-5 +5]); colorbar; colormap(llsmap4.llsmap4);
title('qrenorm = 1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat50 = find(abs(yera) <= 50);
dTMaxR = 10.0; %% find correlation within these values of ERA T
dTMaxR = 1.50; %% find correlation within these values of ERA T
dTMaxR = 1.00; %% find correlation within these values of ERA T
oera = cera(lat50,:); oera = oera(:);
oo2 = c2(lat50,:);       oo2 = oo2(:);              ix = find(isfinite(oera) & abs(oera) < dTMaxR & isfinite(oo2));        R2 = corrcoef(oera(ix),oo2(ix));
coo2 = cal_c2(lat50,:);  coo2 = coo2(:);            ix = find(isfinite(oera) & abs(oera) < dTMaxR & isfinite(coo2));       cR2 = corrcoef(oera(ix),coo2(ix));
oo2_10 = c2_10(lat50,:);      oo2_10 = oo2_10(:);   ix = find(isfinite(oera) & abs(oera) < dTMaxR & isfinite(oo2_10));     R2_10 = corrcoef(oera(ix),oo2_10(ix));
coo2_10 = cal_c2_10(lat50,:); coo2_10 = coo2_10(:); ix = find(isfinite(oera) & abs(oera) < dTMaxR & isfinite(coo2_10));    cR2_10 = corrcoef(oera(ix),coo2_10(ix));
oo3 = c3(lat50,:);       oo3 = oo3(:);              ix = find(isfinite(oera) & abs(oera) < dTMaxR & isfinite(oo3));        R3 = corrcoef(oera(ix),oo3(ix));
oo4 = c4(lat50,:);       oo4 = oo4(:);              ix = find(isfinite(oera) & abs(oera) < dTMaxR & isfinite(oo4));        R4 = corrcoef(oera(ix),oo4(ix));
oo5 = c5(lat50,:);       oo5 = oo5(:);              ix = find(isfinite(oera) & abs(oera) < dTMaxR & isfinite(oo5));        R5 = corrcoef(oera(ix),oo5(ix));
oo6 = c6(lat50,:);       oo6 = oo6(:);              ix = find(isfinite(oera) & abs(oera) < dTMaxR & isfinite(oo6));        R6 = corrcoef(oera(ix),oo6(ix));

[Y,I] = sort(oera);
figure(7); plot(oera(I),oo2(I),'ro',oera(I),oo3(I),'bx',oera(I),oo4(I),'g.'); 
axis([-5 +4 -5 +5]); grid
xlabel('ERA'); ylabel('retr'); hl = legend('reg+cov','cov only, x1','cov only, x10','location','best');

set(0,'DefaultLegendAutoUpdate','off')
dT = -4:0.1:+4;
delta = 0.001*ones(size(dT));

figure(8);
addpath /home/sergio/MATLABCODE/SHOWSTATS
[nz,ny,nx,nmean2,nstd2]         = myhist2d(oera(I),oo2(I),dT,dT);       pcolor(nx,ny,nz); shading flat; title('2'); colorbar; hold on; plot(dT,nmean2,'k','linewidth',2); hold off; disp('ret'); pause;
[nz,ny,nx,cnmean2,cnstd2]       = myhist2d(oera(I),coo2(I),dT,dT);      pcolor(nx,ny,nz); shading flat; title('cal 2'); colorbar; hold on; plot(dT,nmean2,'k','linewidth',2); hold off; disp('ret'); pause;
[nz,ny,nx,nmean2_10,nstd2_10]   = myhist2d(oera(I),oo2_10(I),dT,dT);    pcolor(nx,ny,nz); shading flat; title('2-10'); colorbar; hold on; plot(dT,nmean2_10,'k','linewidth',2); hold off; disp('ret'); pause;
[nz,ny,nx,cnmean2_10,cnstd2_10] = myhist2d(oera(I),coo2_10(I),dT,dT);   pcolor(nx,ny,nz); shading flat; title('cal 2-10'); colorbar; hold on; plot(dT,nmean2_10,'k','linewidth',2); hold off; disp('ret'); pause;
[nz,ny,nx,nmean3,nstd3]         = myhist2d(oera(I),oo3(I),dT,dT);       pcolor(nx,ny,nz); shading flat; title('3'); colorbar; hold on; plot(dT,nmean3,'k','linewidth',2); hold off; disp('ret'); pause;
[nz,ny,nx,nmean4,nstd4]         = myhist2d(oera(I),oo4(I),dT,dT);       pcolor(nx,ny,nz); shading flat; title('4'); colorbar; hold on; plot(dT,nmean4,'k','linewidth',2); hold off; disp('ret'); pause;
[nz,ny,nx,nmean5,nstd5]         = myhist2d(oera(I),oo5(I),dT,dT);       pcolor(nx,ny,nz); shading flat; title('5'); colorbar; hold on; plot(dT,nmean5,'k','linewidth',2); hold off; disp('ret'); pause;
[nz,ny,nx,nmean6,nstd6]         = myhist2d(oera(I),oo6(I),dT,dT);       pcolor(nx,ny,nz); shading flat; title('6'); colorbar; hold on; plot(dT,nmean6,'k','linewidth',2); hold off; disp('ret'); pause;

clf
plot(dT-2*delta,nmean2,'ro',dT-1*delta,nmean2_10,'bs',dT+2*delta,cnmean2,'mx',dT+1*delta,cnmean2_10,'cs','linewidth',2)
hl = legend('reg+cov','reg+covx10','cal reg+cov','cal reg+covx10','location','best');
hold on; errorbar(dT-2*delta,nmean2,nstd2,'color','r'); 
hold on; errorbar(dT-1*delta,nmean2_10,nstd2_10,'color','b'); 
hold on; errorbar(dT+2*delta,cnmean2,cnstd2,'color','m'); 
hold on; errorbar(dT+1*delta,cnmean2_10,cnstd2_10,'color','c'); 
hold off; grid
xlabel('ERA'); ylabel('retr'); 
line([-4 +4],[-4 +4],'color','k');


clf
plot(dT-2*delta,nmean2,'ro',dT-1*delta,nmean2_10,'ms',dT+0*delta,nmean3,'bx',dT+1*delta,nmean4,'g.-',dT+2*delta,nmean5,'co',dT+2*delta,nmean6,'kd','linewidth',2)
hl = legend('reg+cov','reg+covx10','cov only, x1','cov only, x10','reg only, x1','qrenorm 1','location','best');
hold on; errorbar(dT-2*delta,nmean2,nstd2,'color','r'); 
hold on; errorbar(dT-1*delta,nmean2_10,nstd2_10,'color','m'); 
hold on; errorbar(dT-0*delta,nmean3,nstd3,'color','b'); 
hold on; errorbar(dT+1*delta,nmean4,nstd4,'color','g'); 
hold on; errorbar(dT+2*delta,nmean5,nstd5,'color','c'); 
hold on; errorbar(dT+2*delta,nmean6,nstd6,'color','k'); 
hold off; grid
xlabel('ERA'); ylabel('retr'); 
line([-4 +4],[-4 +4],'color','k');

junk = [sqrt(R2(2,1)) sqrt(R2_10(2,1)) sqrt(R3(2,1)) sqrt(R4(2,1)) sqrt(R5(2,1)) sqrt(R6(2,1))];
fprintf(1,'the corr coeffs are %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f\n',junk);

