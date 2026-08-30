function [] = sliceplot_vtk_paraview(rlon,rlat,z,palts,iType)

% sliceplot_vtk_paraview(rlon,rlat,compute_deltaRH.umbc,1);
% sliceplot_vtk_paraview(rlon,rlat,compute_deltaRH.umbc-compute_deltaRH.orig,p.palts(:,3000),2);

i2p5  = find(palts(:,1) <= 2500,1);
i5    = find(palts(:,1) <= 5000,1);
i7p5  = find(palts(:,1) <= 7500,1);
i10   = find(palts(:,1) <= 10000,1);
i12p5 = find(palts(:,1) <= 12500,1);
i15   = find(palts(:,1) <= 15000,1);
i17p5 = find(palts(:,1) <= 17500,1);
i20   = find(palts(:,1) <= 20000,1);
i30 = find(palts(:,1) <= 30000,1);
i40 = find(palts(:,1) <= 40000,1);
i50 = find(palts(:,1) <= 50000,1);
i60 = find(palts(:,1) <= 60000,1);

iaZ  = [2.5  5  7.5  10  12.5  15  17.5  20];
iAlt = [i2p5 i5 i7p5 i10 i12p5 i15 i17p5 i20];

iaZ  = [2.5  5  7.5  10];
iAlt = [i2p5 i5 i7p5 i10];

[X,Y,Z] = ndgrid(double(rlon),double(rlat),double(palts(iAlt,1)/1000));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mean(rlon) < 0
  azangle = +60;
  Yplot = min(Y(:));
else
  azangle = -120;
  Yplot = max(Y(:));
end

addpath Matlab3DPlots

if iType == 1
  figure(1); clf
  S = z(iAlt,:);
  S = reshape(S,length(iAlt),72,64);
  S = permute(S,[2 3 1]);
  vtkwrite('humidity.vtk', 'structured_grid', X, Y, Z, 'scalars', 'RH', S);

  h = slice(Y,X,Z,double(S), [Yplot], [], iaZ);
  shading interp; 
  set(h,'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
  % set transparency to correlate to the data values.
  %alpha('color');
  alpha(0.5);

  %h = contourslice(Y,X,Z,double(S), [Yplot], [], iaZ, 20);

  colormap jet; colorbar; caxis([0 120]); 
  ylabel('Longitude'); xlabel('Latitude'); zlabel('Hgt(km)'); title('RH(x,y,z)')
  zlim([0 max(iaZ)]);
  
  figure(2); clf
  isosurface(Y,X,Z,double(S),25); hold on
  isosurface(Y,X,Z,double(S),50); hold on
  isosurface(Y,X,Z,double(S),75); hold on
  isosurface(Y,X,Z,double(S),100); hold off
  %slice(Y,X,Z,double(S), [Yplot], [], iaZ);
  hold off
  colormap jet; colorbar; caxis([0 120]); 
  ylabel('Longitude'); xlabel('Latitude'); zlabel('Hgt(km)'); title('RH(x,y,z)')
  zlim([0 max(iaZ)]);

  figure(3); clf
  scatter3(X(:), Y(:), Z(:), 20, double(S(:)), 'filled')
  colorbar; caxis([0 120]); colormap jet
else
  figure(1); clf
  S = z(iAlt,:);
  S = reshape(S,length(iAlt),72,64);
  S = permute(S,[2 3 1]);
  vtkwrite('humidity_rate.vtk', 'structured_grid', X, Y, Z, 'scalars', 'RH', S);

  h = slice(Y,X,Z,double(S), [Yplot], [], iaZ);
  shading interp; 
  set(h,'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
  % set transparency to correlate to the data values.
  %alpha('color');
  alpha(0.5);

  %h = contourslice(Y,X,Z,double(S), [Yplot], [], iaZ, 40);

  colormap(usa2); colorbar; caxis([-1 +1]); 
  ylabel('Longitude'); xlabel('Latitude'); zlabel('Hgt(km)'); title('dRH/dt(x,y,z)')
  zlim([0 max(iaZ)]);
  
  figure(2); clf
  isosurface(Y,X,Z,double(S),-1.0); hold on
  isosurface(Y,X,Z,double(S),-0.5); hold on
  isosurface(Y,X,Z,double(S), 0.0); hold on
  isosurface(Y,X,Z,double(S),+0.5); hold on
  isosurface(Y,X,Z,double(S),+1.0); hold off
  %slice(Y,X,Z,double(S), [Yplot], [], iaZ);
  hold off
  colormap(usa2); colorbar; caxis([-1 +1]); 
  ylabel('Longitude'); xlabel('Latitude'); zlabel('Hgt(km)'); title('dRH/dt(x,y,z)')
  zlim([0 max(iaZ)]);

  figure(3); clf
  scatter3(X(:), Y(:), Z(:), 20, double(S(:)), 'filled')
  colorbar; caxis([-1 +1]); colormap(usa2)
end

figure(1); axis([-90 +90 -180 +180 0 max(iaZ)])
figure(2); axis([-90 +90 -180 +180 0 max(iaZ)])
figure(3); axis([-180 +180 -90 +90 0 max(iaZ)])
