figure(6); aslmap(6,rlat65,rlon73,maskLFmatr.*smoothn((reshape(results(:,6),72,64)'),1),[-90 +90],[-180 +180]); colormap(llsmap5);  
strXYZ = [];
if iNorD > 0
  strXYZ = 'Night ';
else
  strXYZ = 'Day ';
end

if dataset == 1
  strXYZ = [strXYZ 'Strow 02/09-19/08 Qtile' num2str(iQuantile)];
elseif dataset == -1
  strXYZ = [strXYZ 'Sergio 02/09-19/08 Qtile' num2str(iQuantile)];
elseif dataset == 2
  strXYZ = [strXYZ 'Sergio 02/09-21/08 Qtile' num2str(iQuantile)];
elseif dataset == 3
  strXYZ = [strXYZ 'Sergio 02/09-21/08 Extreme'];
end
colormap(llsmap5);
title([strXYZ ' d/dt ST']);    caxis([-0.15 +0.15])
