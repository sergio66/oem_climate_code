addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER

for ii = 1 : 6
  fig_screen(ii);
end

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
pN = plevs(1:end-1)-plevs(2:end);
pD = log(plevs(1:end-1)./plevs(2:end));
plays = pN./pD;
plays = flipud(plays(4:100));

for ii = 1 : 40
  filex = ['Output/test' num2str(ii) '.mat'];
  a = load(filex);

  chanset(:,ii) = zeros(1,2645);
  g = a.jacobian.chanset;
  chanset(g,ii) = 1;
 
  traceNstemp(:,ii) = a.oem.finalrates(1:6);
  wv_ret(:,ii)      = a.oem.finalrates(7:103);
  temp_ret(:,ii)    = a.oem.finalrates(104:200);
  ozone_ret(:,ii)   = a.oem.finalrates(201:297);

  traceNstemp_sigs(:,ii) = a.oem.finalsigs(1:6);
  wv_ret_sigs(:,ii)      = a.oem.finalsigs(7:103);
  temp_ret_sigs(:,ii)    = a.oem.finalsigs(104:200);
  ozone_ret_sigs(:,ii)   = a.oem.finalsigs(201:297);

  input_rates(:,ii) = a.rateset.rates;
  fit_to_rates(:,ii) = a.oem.fit;
end

addpath /home/sergio/MATLABCODE/PLOTTER
lats = -90 : 5 : +90;
lats = equal_area_spherical_bands(20);
latx = 0.5*(lats(1:end-1) + lats(2:end));
%latx = latx(1:36);

addpath /home/sergio/MATLABCODE/COLORMAP
figure(1); clf
  pcolor(latx,plays,wv_ret*10); shading flat
  caxis([-0.025 +0.025]); colorbar; title('WV(lat,z)/decade')
  caxis([-0.01  +0.01]); colorbar; title('WV(lat,z)/decade')  
  caxis([-0.10  +0.10]); colorbar; title('WV(lat,z)/decade')  

figure(2); clf
  pcolor(latx,plays,temp_ret*10); shading flat;
  caxis([-0.15 +0.15]); colorbar; title('T(lat,z)/decade')
  caxis([-0.50 +0.50]); colorbar; title('T(lat,z)/decade')

figure(3); clf
  pcolor(latx,plays,ozone_ret*10); shading flat
  caxis([-0.025 +0.025]); colorbar; title('O3(lat,z)/decade')
  caxis([-0.02  +0.02]);  colorbar; title('O3(lat,z)/decade')  
  caxis([-0.2  +0.2]);  colorbar; title('O3(lat,z)/decade')  

figure(4); clf
  pcolor(latx,plays,wv_ret_sigs*10); shading flat
  caxis([0.0  +4e-3]); colorbar; title('WV(lat,z) sig/decade')  
  caxis([0.0  +4e-2]); colorbar; title('WV(lat,z) sig/decade')  

figure(5); clf
  pcolor(latx,plays,temp_ret_sigs*10); shading flat;
  caxis([0 0.15]); colorbar; title('T(lat,z) sig/decade')
  caxis([0 0.50]); colorbar; title('T(lat,z) sig/decade')

figure(6); clf
  pcolor(latx,plays,ozone_ret_sigs*10); shading flat
  caxis([0.0  +0.01]);  colorbar; title('O3(lat,z) sig/decade')  
  caxis([0.0  +0.10]);  colorbar; title('O3(lat,z) sig/decade')  

jett = jet; jett(1,:) = 1;
for ii = 1 : 6
  figure(ii)
  colormap(jett)
  if ii <= 6
    colormap usa2  
    llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); colormap(llsmap4.llsmap4);
  end
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
end

figure(7); clf
  plot(latx,traceNstemp(1:5,:),latx,traceNstemp(6,:)*100,'g','linewidth',2); grid
  hl = legend('CO2','N2O','CH4','CFC11','CFC12','stemp*100'); set(hl,'fontsize',10);
  title('tracegas and stemp rates/yr')

%{
figure(4); clf
  plot(latx,traceNstemp(1:5,:),latx,traceNstemp(6,:)*1,'g','linewidth',2); grid
  hl = legend('CO2','N2O','CH4','CFC11','CFC12','stemp*1'); set(hl,'fontsize',10);

%% these are absolute
figure(5); clf
  plot(latx,traceNstemp(7,:),'b',latx,traceNstemp(8,:),'r',...
       latx,traceNstemp(9,:),'c',latx,traceNstemp(10,:),'m','linewidth',2); grid
  hl = legend('cng1','cng2','csz1','csz2'); set(hl,'fontsize',10);

%% these are nominal fractional
figure(5); clf
  plot(latx,traceNstemp(7,:)/100,'b',latx,traceNstemp(8,:)/100,'r',...
       latx,traceNstemp(9,:)/50,'c',latx,traceNstemp(10,:)/20,'m','linewidth',2); grid
  hl = legend('frac cng1','frac cng2','frac csz1','frac csz2'); set(hl,'fontsize',10);

avgcloud = load('the_cloud_params.mat');
figure(5); clf
  plot(latx,traceNstemp(7,:)./avgcloud.cngwat,'b',latx,traceNstemp(8,:)./avgcloud.cngwat2,'r',...
       latx,traceNstemp(9,:)./avgcloud.cpsize,'c',latx,traceNstemp(10,:)./avgcloud.cpsize2,'m',...
       'linewidth',2); grid
  hl = legend('frac cng1','frac cng2','frac csz1','frac csz2'); set(hl,'fontsize',10);
axis([-100 +100 -4e-2 +2e-2])
%}

g = a.jacobian.chanset;
fairs = instr_chans;
hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
fairs = vchan2834;
load sarta_chans_for_l1c.mat
fairs = fairs(ichan);

iaLat = 1:40
iaLat = 2:39;
figure(8); clf
  plot(fairs(g),nanmean(input_rates(g,iaLat)'-fit_to_rates(g,iaLat)'),'b',...
       fairs(g),nanstd(input_rates(g,iaLat)'-fit_to_rates(g,iaLat)'),'r',...
       fairs(g),nanmean(input_rates(g,iaLat)'),'k','linewidth',2); grid
%  title('(b) bias (r) std (k) signal over all lats')
  hl = legend('bias','std dev','mean signal','location','best'); set(hl,'fontsize',10)
axis([min(fairs(g)) max(fairs(g)) -0.1 +0.1])

for ii = 1 : 3
  figure(ii); shading interp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iYes = input('show ERA rates (with AK?) : ');
if iYes > 0
  compute_geo_rates_ak
  figure(1); caxis([-0.1 +0.1]); colorbar; figure(9);  caxis([-0.1 +0.1]); colorbar; 
  figure(2); caxis([-0.5 +0.5]); colorbar; figure(10); caxis([-0.5 +0.5]); colorbar; 
  figure(7); figure(12);
end

iYes = input('show AIRS L3 rates? : ');
if iYes > 0
  %show_AIRSL3
  show_AIRSL3_rates_ak
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
m = mean(chanset');
s = std(chanset');
woo = find(m == 1 & s == 0); 

comment = 'see Plotutils/plot_all_latbins.m';
save good_chans.mat comment woo
%}
