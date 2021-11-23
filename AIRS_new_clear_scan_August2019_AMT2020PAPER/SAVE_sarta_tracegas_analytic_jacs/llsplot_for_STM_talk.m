nchan = 2314;
%nchan = 2378;

plotsigs = true;
dotitle = false;

addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE//PLOTTER

load /home/sergio/MATLABCODE/COLORMAP/LLS/llsmap5

%dirname = 'Simulation_lessTdamp';
%dirname = 'Obs_simul_lessTdamp_params';
%dirname = 'Obs_lessTdamp_realerrors';
%dirname = 'Obs_simul_lessTdamp_params_fixedCO2';
dirname = 'No_mec_no_sw';

% for ii = 1 : 6
%   fig_screen(ii);
% end

plevs = load('/home/sergio/MATLABCODE//airslevels.dat');
pN = plevs(1:end-1)-plevs(2:end);
pD = log(plevs(1:end-1)./plevs(2:end));
plays = pN./pD;
plays = flipud(plays(4:100));

for ii = 1 : 40
  filex = ['/asl/s1/strow/Airs_random_fits/Output/' dirname '/test' num2str(ii) '.mat'];
  filex = ['/asl/s1/strow/Airs_random_fits/No_mec_no_sw//test' num2str(ii) '.mat'];
  a = load(filex);

  chanset(:,ii) = zeros(1,nchan);
  g = a.jacobian.chanset;
  chanset(g,ii) = 1;

%{
  traceNstemp(:,ii) = a.oem.finalrates(1:11);
  wv_ret(:,ii) = a.oem.finalrates(12:108);
  temp_ret(:,ii) = a.oem.finalrates(109:205);
  ozone_ret(:,ii) = a.oem.finalrates(206:302);

  traceNstemp_sigs(:,ii) = a.oem.finalsigs(1:11);
  wv_ret_sigs(:,ii) = a.oem.finalsigs(12:108);
  temp_ret_sigs(:,ii) = a.oem.finalsigs(109:205);
  ozone_ret_sigs(:,ii) = a.oem.finalsigs(206:302);
%}

  traceNstemp(:,ii) = a.oem.finalrates(1:9);
  wv_ret(:,ii) = a.oem.finalrates(10:106);
  temp_ret(:,ii) = a.oem.finalrates(107:203);
  ozone_ret(:,ii) = a.oem.finalrates(204:300);

  traceNstemp_sigs(:,ii) = a.oem.finalsigs(1:9);
  wv_ret_sigs(:,ii) = a.oem.finalsigs(10:106);
  temp_ret_sigs(:,ii) = a.oem.finalsigs(107:203);
  ozone_ret_sigs(:,ii) = a.oem.finalsigs(204:300);

  input_rates(:,ii) = a.rateset.rates;
  fit_to_rates(:,ii) = a.oem.fit;
 
end

addpath /home/sergio/MATLABCODE//PLOTTER
addpath /home/sergio/MATLABCODE//COLORMAP
lats = equal_area_spherical_bands(20);
latx = 0.5*(lats(1:end-1) + lats(2:end));

%llsmap4 = load('/home/sergio/MATLABCODE//COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);

addpath MATLABCODE/oem_pkg_run/cris_stability/Fits/Plotutils
load labels
load xtickl

if plotsigs

figure(1); clf
%  sm_wv_ret = smoothn(wv_ret,4);
  pcolor(1:40,plays,100*wv_ret_sigs); shading flat
  caxis([0 0.04]);
  colorbar; 
  if dotitle title('WV(lat,z) Unc') ; end;
  colormap('default');
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000])
  set(gca,'YTick',ypt);
  set(gca,'YTickLabel',ytickl);
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);
  xlabel('Latitude');ylabel('Pressure (mbar)');

figure(2); clf
%  sm_temp_ret = smoothn(temp_ret,4);
  pcolor(1:40,plays,10*temp_ret_sigs); shading flat; 
  colormap('default');
  caxis([0 0.02])
  colorbar; 
  if dotitle title('T(lat,z) Unc'); end;
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  xlabel('Latitude');ylabel('Pressure (mbar)');

  figure(3); clf
%  sm_ozone_ret = smoothn(ozone_ret,4);
  pcolor(1:40,plays,100*ozone_ret_sigs); shading flat; 
  colormap('default')
  caxis([0.15 0.2]);
  %caxis([-0.02  +0.02]);  
  colorbar; 
  if dotitle title('O3(lat,z)'); end;
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  xlabel('Latitude');ylabel('Pressure (mbar)');
  
end
   figure(4); clf
  pcolor(1:40,plays,100*wv_ret); shading flat
  caxis([-1  +1]); colorbar; 
  if dotitle title('WV(lat,z)'); end;
  colormap(llsmap5);
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000])
  set(gca,'YTick',ypt);
  set(gca,'YTickLabel',ytickl);
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);
  shading interp;
  xlabel('Latitude');ylabel('Pressure (mbar)');

figure(5); clf
  pcolor(1:40,plays,10*temp_ret); shading flat; colormap(llsmap5);
  caxis([-01.5 +01.5]); colorbar; 
  if dotitle title('T(lat,z)');end
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp;
  xlabel('Latitude');ylabel('Pressure (mbar)');
  
  figure(6); clf
  pcolor(1:40,plays,100*ozone_ret); shading flat; colormap(llsmap5)
  caxis([-2  2]);  colorbar; 
  if dotitle title('O3(lat,z)'); end;
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp;  
  xlabel('Latitude');ylabel('Pressure (mbar)');

  figure(7); clf
%   plot(latx,traceNstemp(1:4,:),latx,traceNstemp(5,:)*1,'g','linewidth',2); grid
%   hl = legend('CO2','N2O','CH4','CFC11','stemp*1'); set(hl,'fontsize',10);
subplot(211);  plot(latx,traceNstemp(5,:)*1,'linewidth',2); grid
  hl = legend('T Surf'); set(hl,'fontsize',12);
  xlim([-80 80]);
  set(gca,'fontsize',16);
  %xlabel('Latitude');
  ylabel('\Delta T in K/Year')
set(gca,'xtick',[])
set(gca,'xticklabel',[])

%figure(8); clf
subplot(212)
  plot(latx,traceNstemp(6:7,:),'--',latx,traceNstemp(8:9,:),'-','linewidth',2); grid
%  hl = legend('cng1','cng2','csz1','csz2'); 
  hl = legend('Ice Cloud (g/m^2/year)','Water Cloud (g/m^2/year)','Ice Size micron/year','Water Size micron/year'); 
  set(hl,'fontsize',12);
  xlim([-80 80]);
  xlabel('Latitude');
  ylabel('Cloud params/yr')
  
g = a.jacobian.chanset;
fairs = instr_chans;
figure(9); clf
  plot(fairs(g),nanmean(input_rates(g,:)'-fit_to_rates(g,:)'),...
       fairs(g),nanstd(input_rates(g,:)'-fit_to_rates(g,:)'),...
       fairs(g),nanmean(input_rates(g,:)')); grid
hl = legend('Fit Bias','Fit Std','Obs Rate')
addpath MATLABCODE/oem_pkg_run/AIRS_new_random_scan_Aug2018/Plotutils/WORKS_Sept16_2018/
xwn;
ylabel('\Delta B(T) in K/Year')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iYes = input('show ERA rates (with AK?) : ');
if iYes > 0
  compute_geo_rates_ak
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
m = mean(chanset');
s = std(chanset');
woo = find(m == 1 & s == 0); 

comment = 'see Plotutils/plot_all_latbins.m';
save good_chans.mat comment woo
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure(4); set(gca,'fontsize',16); aslprint([savedir '/umbc_wv_rates_from_spectralrates.png'],1);
figure(5); set(gca,'fontsize',16); aslprint([savedir '/umbc_T_rates_from_spectralrates.png'],1);
figure(6); set(gca,'fontsize',16); aslprint([savedir '/umbc_o3_rates_from_spectralrates.png'],1);
figure(7); set(gca,'fontsize',16); aslprint([savedir '/umbc_STnCLD_rates_from_spectralrates.pdf']);
%}
