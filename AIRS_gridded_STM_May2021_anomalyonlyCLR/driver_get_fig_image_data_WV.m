llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat');

figure(1); clf
figname = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/SAVE_BESTRUNv1/Figs/era_400mb_wv_anom_nosmooth.fig';
%figname = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/SAVE_BESTRUNv1/Figs/era_400mb_wv_anom_nosmooth_AK.fig';
[tera,yera,cera] = get_fig_image_data(figname); 
shading interp; caxis([-1 +1]); colorbar; colormap(llsmap4.llsmap4);
title('ERA')

%%%%%%%%%%

figure(2); clf
figname = 'SAVE_LW_noCFC11_noN2O_Bootstrap/Figs/obs_wv_anom_400mb.fig';
[t2,y2,c2] = get_fig_image_data(figname);
shading interp; caxis([-1 +1]); colorbar; colormap(llsmap4.llsmap4);
title('bootstrap,obs')

figure(3); clf
figname = 'SAVE_LW_noCFC11_noN2O_cov100/Figs/obs_wv_anom_400mb.fig';
[t3,y3,c3] = get_fig_image_data(figname);
shading interp; caxis([-1 +1]); colorbar; colormap(llsmap4.llsmap4);
title('cov100+reg, obs')

figure(4); clf
figname = 'SAVE_LW_noCFC11_noN2O_covx10/Figs/obs_wv_anom_400mb.fig';
figname = 'SAVE_LW_noCFC11_noN2O_covx10/Figs/cal_wv_anom_400mb.fig';
[t4,y4,c4] = get_fig_image_data(figname);
shading interp; caxis([-1 +1]); colorbar; colormap(llsmap4.llsmap4);
title('cov 10+reg, obs')
title('cov 10+reg, cal')

figure(4); clf
figname = 'SAVE_LW_noCFC11_noN2O/Figs/obs_wv_anom_400mb.fig';
[t5,y5,c5] = get_fig_image_data(figname);
shading interp; caxis([-1 +1]); colorbar; colormap(llsmap4.llsmap4);
title('cov+reg, obs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lat50 = find(abs(yera) <= 50);
dWVMaxR = 1.00; %% find correlation within these values of ERA WV
oera = cera(lat50,:); oera = oera(:);
oo2 = c2(lat50,:);       oo2 = oo2(:);              ix = find(isfinite(oera) & abs(oera) < dWVMaxR & isfinite(oo2));        R2 = corrcoef(oera(ix),oo2(ix));
oo3 = c3(lat50,:);       oo3 = oo3(:);              ix = find(isfinite(oera) & abs(oera) < dWVMaxR & isfinite(oo3));        R3 = corrcoef(oera(ix),oo3(ix));
oo4 = c4(lat50,:);       oo4 = oo4(:);              ix = find(isfinite(oera) & abs(oera) < dWVMaxR & isfinite(oo4));        R4 = corrcoef(oera(ix),oo4(ix));
oo5 = c5(lat50,:);       oo5 = oo5(:);              ix = find(isfinite(oera) & abs(oera) < dWVMaxR & isfinite(oo5));        R5 = corrcoef(oera(ix),oo5(ix));

[Y,I] = sort(oera);
figure(7); plot(oera(I),oo2(I),'ro',oera(I),oo3(I),'bx',oera(I),oo4(I),'g.',oera(I),oo5(I),'ks'); 
axis([-1 +1 -1 +1]); grid
xlabel('ERA'); ylabel('retr'); hl = legend('bootstrap','cov100+reg','cov10+reg CAL','cov+reg','location','best');

set(0,'DefaultLegendAutoUpdate','off')
dWV = -1:0.05:+1;
delta = 0.001*ones(size(dWV));

figure(8);
addpath /home/sergio/MATLABCODE/SHOWSTATS
[nz,ny,nx,nmean2,nstd2]         = myhist2d(oera(I),oo2(I),dWV,dWV);       pcolor(nx,ny,nz); shading flat; title('2'); colorbar; hold on; plot(dWV,nmean2,'k','linewidth',2); hold off; disp('ret'); pause;
[nz,ny,nx,nmean3,nstd3]         = myhist2d(oera(I),oo3(I),dWV,dWV);       pcolor(nx,ny,nz); shading flat; title('3'); colorbar; hold on; plot(dWV,nmean3,'k','linewidth',2); hold off; disp('ret'); pause;
[nz,ny,nx,nmean4,nstd4]         = myhist2d(oera(I),oo4(I),dWV,dWV);       pcolor(nx,ny,nz); shading flat; title('4'); colorbar; hold on; plot(dWV,nmean4,'k','linewidth',2); hold off; disp('ret'); pause;
[nz,ny,nx,nmean5,nstd5]         = myhist2d(oera(I),oo5(I),dWV,dWV);       pcolor(nx,ny,nz); shading flat; title('5'); colorbar; hold on; plot(dWV,nmean5,'k','linewidth',2); hold off; disp('ret'); pause;

clf
plot(dWV-2*delta,nmean2,'ro',dWV-1*delta,nmean3,'mx',dWV+1*delta,nmean4,'bo',dWV+2*delta,nmean5,'cs','linewidth',2)
hl = legend('bootstrap','cov100+reg','cov10+reg CAL','cov+reg','location','best');
hold on; errorbar(dWV-2*delta,nmean2,nstd2,'color','r'); 
hold on; errorbar(dWV-1*delta,nmean3,nstd3,'color','m'); 
hold on; errorbar(dWV+1*delta,nmean4,nstd4,'color','b'); 
hold on; errorbar(dWV+2*delta,nmean5,nstd5,'color','c'); 
hold off; grid
xlabel('ERA'); ylabel('retr'); 
line([-1 +1],[-1 +1],'color','k');

junk = [sqrt(R2(2,1)) sqrt(R3(2,1)) sqrt(R4(2,1)) sqrt(R5(2,1)) ];
fprintf(1,'the corr coeffs are %8.6f %8.6f %8.6f %8.6f\n',junk);

figure(9); plot(dWV,histc(oera,dWV),'k',dWV,histc(oo2,dWV),'r',dWV,histc(oo3,dWV),'m',dWV,histc(oo4,dWV),'b',dWV,histc(oo5,dWV),'c','linewidth',2);
hl = legend('ERA','bootstrap','cov100+reg','cov10+reg CAL','cov+reg','location','best');
