nan_bottom = ones(101,4608);
for junk = 1 : 4608
  nlevs = p.nlevs(junk);
  nan_bottom(nlevs:101,junk) = nan;
  deltaT(nlevs:101,junk) = 0.0;
  fracWV(nlevs:101,junk) = 0.0;
  fracO3(nlevs:101,junk) = 0.0;
end

nan_bottom100 = nan_bottom(1:100,:);
nan_bottom097 = nan_bottom(1:097,:);

figure(27); clf; aslmap(27,rlat65,rlon73,reshape(p.nlevs,72,64)', [-90 +90],[-180 +180]);
  title('Nlevs')

figure(28); 
junk = squeeze(nanmean(reshape(nan_bottom,101,72,64),2));
imagesc(junk); colorbar; title('NAN or not for SURFACE')
