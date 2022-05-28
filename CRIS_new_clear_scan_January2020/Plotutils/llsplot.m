plotsigs = true;

addpath ~/MATLABCODE
addpath ~/MATLABCODE/PLOTTER

%dirname = 'Simulation_lessTdamp';
%dirname = 'Obs_simul_lessTdamp_params';
%dirname = 'Obs_lessTdamp_realerrors';
%dirname = 'Obs_simul_lessTdamp_params_fixedCO2';
dirname = 'Test2';

% for ii = 1 : 6
%   fig_screen(ii);
% end

plevs = load('~/MATLABCODE/airslevels.dat');
pN = plevs(1:end-1)-plevs(2:end);
pD = log(plevs(1:end-1)./plevs(2:end));
plays = pN./pD;
plays = flipud(plays(4:100));

for ii = 1 : 40
  filex = ['Output/' dirname '/test' num2str(ii) '.mat'];
  filex = ['Output/test' num2str(ii) '.mat'];  
  a = load(filex);

  chanset(:,ii) = zeros(1,2378);
  g = a.jacobian.chanset;
  chanset(g,ii) = 1;
 
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

addpath ~/MATLABCODE/PLOTTER
addpath ~/MATLABCODE/COLORMAP
lats = equal_area_spherical_bands(20);
latx = 0.5*(lats(1:end-1) + lats(2:end));

llsmap4 = load('~/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);

load labels
load xtickl

if plotsigs

figure(1); clf
  pcolor(1:40,plays,wv_ret); shading flat
  caxis([-0.01  +0.01]); colorbar; title('WV(lat,z)') ; colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000])
  set(gca,'YTick',ypt);
  set(gca,'YTickLabel',ytickl);
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);

figure(2); clf
  pcolor(1:40,plays,temp_ret); shading flat; colormap(llsmap4.llsmap4);
  caxis([-0.15 +0.15]); colorbar; title('T(lat,z)')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  

  figure(3); clf
  pcolor(1:40,plays,ozone_ret); shading flat; colormap(llsmap4.llsmap4)
  caxis([-0.02  +0.02]);  colorbar; title('O3(lat,z)')  
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  
else
   figure(1); clf
  pcolor(1:40,plays,wv_ret); shading flat
  caxis([-0.01  +0.01]); colorbar; title('WV(lat,z)') ; colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([100 1000])
  set(gca,'YTick',ypt);
  set(gca,'YTickLabel',ytickl);
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);

figure(2); clf
  pcolor(1:40,plays,temp_ret); shading flat; colormap(llsmap4.llsmap4);
  caxis([-0.15 +0.15]); colorbar; title('T(lat,z)')
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  

  figure(3); clf
  pcolor(1:40,plays,ozone_ret); shading flat; colormap(llsmap4.llsmap4)
  caxis([-0.02  +0.02]);  colorbar; title('O3(lat,z)')  
  set(gca,'ydir','reverse'); set(gca,'yscale','log'); ylim([10 1000])
  set(gca,'XTick',xpt);
  set(gca,'XTickLabel',xtickl);
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  
end
  
figure(4); clf
  plot(latx,traceNstemp(1:5,:),latx,traceNstemp(6,:)*1,'g','linewidth',2); grid
  hl = legend('CO2','N2O','CH4','CFC11','CFC12','stemp*1'); set(hl,'fontsize',10);

figure(5); clf
  plot(latx,traceNstemp(6:7,:),'--',latx,traceNstemp(8:9,:),'-','linewidth',2); grid
  hl = legend('cng1','cng2','csz1','csz2'); set(hl,'fontsize',10);

g = a.jacobian.chanset;
fairs = instr_chans;
figure(6); clf
  plot(fairs(g),nanmean(input_rates(g,:)'-fit_to_rates(g,:)'),'b',...
       fairs(g),nanstd(input_rates(g,:)'-fit_to_rates(g,:)'),'r',...
       fairs(g),nanmean(input_rates(g,:)'),'k','linewidth',2); grid
  title('(b) bias (r) std (k) signal over all lats')

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
