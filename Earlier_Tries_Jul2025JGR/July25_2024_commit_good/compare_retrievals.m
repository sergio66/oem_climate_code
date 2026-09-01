function [] = compare_retrievals(mat1,mat2);

%{
compare_retrievals('results_nucal_nopert.mat','results_nucal_pert.mat');
compare_retrievals('results_nucal_LWMW2.mat','results_nucal_LWMW_era_ab0_v2.mat');
%}

a1 = load(mat1);
a2 = load(mat2);

for ii = 1 : 5
  figure(1); clf;
  if ii == 1
    str = 'co2';
  elseif ii == 2
    str = 'n2o';
  elseif ii == 3
    str = 'ch4';
  elseif ii == 4
    str = 'cfc';
  elseif ii == 5
    str = 'stemp';
  end
  subplot(211); plot(a1.latx,a1.traceNstemp(ii,:),'b',a1.latx,a2.traceNstemp(ii,:),'r'); title(num2str(ii)); title([str ' : 1-2']);
  subplot(212); plot(a1.latx,a1.traceNstemp(ii,:)-a2.traceNstemp(ii,:),'k');
  disp('ret to continue'); pause
end

figure(1); clf
  pcolor(a1.latx,a1.plays,a1.wv_ret*10); shading flat
  caxis([-0.025 +0.025]); colorbar; title('WV(lat,z)/decade')
  caxis([-0.01  +0.01]); colorbar; title('WV(lat,z)/decade')  
  caxis([-0.10  +0.10]); colorbar; title('WV(lat,z)/decade')  

figure(2); clf
  pcolor(a1.latx,a1.plays,a1.temp_ret*10); shading flat;
  caxis([-0.15 +0.15]); colorbar; title('T(lat,z)/decade')
  caxis([-0.50 +0.50]); colorbar; title('T(lat,z)/decade')

figure(3); clf
  pcolor(a1.latx,a1.plays,a1.ozone_ret*10); shading flat
  caxis([-0.025 +0.025]); colorbar; title('O3(lat,z) frac/decade')
  caxis([-0.02  +0.02]);  colorbar; title('O3(lat,z) frac/decade')  
  caxis([-0.2  +0.2]);  colorbar; title('O3(lat,z) frac/decade')  

llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); 
for ii = 1:3
  figure(ii)
  colormap usa2  
  colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp
end

disp('ret to continue'); pause
figure(4); clf
  pcolor(a1.latx,a1.plays,(a1.wv_ret./a2.wv_ret)-1); shading flat
  caxis([-0.025 +0.025]); colorbar; title('WV1(lat,z)/WV2(lat,z)-1')
  caxis([-0.01  +0.01]); colorbar; title('WV1(lat,z)/WV2(lat,z)-1')  
  caxis([-0.10  +0.10]); colorbar; title('WV1(lat,z)/WV2(lat,z)-1')  
  caxis([-10 +10]); colorbar; title('WV1(lat,z)/WV2(lat,z)-1')  

figure(5); clf
  pcolor(a1.latx,a1.plays,100*(a1.temp_ret-a2.temp_ret)./a1.temp_ret); shading flat;
  caxis([-0.15 +0.15]); colorbar; title('percent T1(lat,z) - T2(lat,z)/decade')
  caxis([-200.0 +200.0]); colorbar; title('percent T1(lat,z) - T2(lat,z)/decade')

figure(6); clf
  pcolor(a1.latx,a1.plays,(a1.ozone_ret./a2.ozone_ret)-1); shading flat
  caxis([-0.025 +0.025]); colorbar; title('O31(lat,z)/O32(lat,z)-1')
  caxis([-0.02  +0.02]);  colorbar; title('O31(lat,z)/O32(lat,z)-1')  
  caxis([-0.2  +0.2]);  colorbar; title('O31(lat,z)/O32(lat,z)-1')  
  caxis([-10 +10]);  colorbar; title('O31(lat,z)/O32(lat,z)-1')  

llsmap4 = load('//home/sergio/MATLABCODE/COLORMAP/llsmap4.mat'); 
for ii = 4:6
  figure(ii)
  colormap usa2  
  colormap(llsmap4.llsmap4);
  set(gca,'ydir','reverse')
  set(gca,'yscale','log')
  axis([-90 +90 9 1000])
  set(gca,'YTick',[10,100,1000])
  set(gca,'YTickLabel',{'10','100','1000'})  
  shading interp
end

