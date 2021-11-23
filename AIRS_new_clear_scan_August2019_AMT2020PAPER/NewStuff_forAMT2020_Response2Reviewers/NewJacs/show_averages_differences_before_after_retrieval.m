boo20 = find(pbefore.latbin == 20);
yyx = 2002:2017;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); subplot(211); plot(yyx,pafter.stemp(boo20)-pbefore.stemp(boo20)); title('\delta STEMP'); grid
           subplot(212); plot(yyx,pafter.colWV(boo20)-pbefore.colWV(boo20)); title('\delta COLWV'); grid

figure(2); pcolor(yyx,1:101,pafter.ptemp(:,boo20)-pbefore.ptemp(:,boo20)); title('\delta PTEMP');
  set(gca,'ydir','reverse'); shading interp; colormap(usa2); colorbar; caxis([-0.1 +0.1]); 

figure(3); pcolor(yyx,1:101,pafter.gas_1(:,boo20)./pbefore.gas_1(:,boo20)-1); title('\delta GAS1');
  set(gca,'ydir','reverse'); shading interp; colormap(usa2); colorbar; caxis([-0.25 +0.25]); 
figure(4); pcolor(yyx,1:101,pafter.RH(:,boo20) - pbefore.RH(:,boo20)-1); title('\delta RH');
  set(gca,'ydir','reverse'); shading interp; colormap(usa2); colorbar; caxis([-20 +20]); 

figure(5); pcolor(yyx,1:101,pafter.gas_2(:,boo20)./pbefore.gas_2(:,boo20)-1); title('\delta GAS2');
  set(gca,'ydir','reverse'); shading interp; colormap(usa2); colorbar; caxis([-0.25 +0.25]); 
figure(5); pcolor(yyx,1:101,pafter.gas_2(:,boo20)./(pbefore.gas_2(:,boo20(1)) * ones(1,length(yyx))) - 1); title('\delta GAS2 against 2002');
  set(gca,'ydir','reverse'); shading interp; colormap(usa2); colorbar; caxis([-0.25 +0.25]); 

figure(6); pcolor(yyx,1:101,pafter.gas_3(:,boo20)./pbefore.gas_3(:,boo20)-1); title('\delta GAS3');
  set(gca,'ydir','reverse'); shading interp; colormap(usa2); colorbar; caxis([-0.25 +0.25]); 
