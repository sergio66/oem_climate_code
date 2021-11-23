figure(6); clf
  junk = smoothn(reshape(((boo_lpsUMBC-boo_lps0)./results(:,6)'),72,64)',1);  aslmap(6,rlat65,rlon73,maskLFmatr.*junk,[-90 +90],[-180 +180]);   
  colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/K UMBC');  
figure(7); clf
  junk = smoothn(reshape(((boo_lpsAIRSL3-boo_lps0)./reshape(airsL3.thestats64x72.stemprate,1,72*64)),72,64)',1);  aslmap(7,rlat65,rlon73,maskLFmatr.*junk,[-90 +90],[-180 +180]);   
  colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/K AIRSL3');  
figure(8); clf
  junk = smoothn(reshape(((boo_lpsCMIP6-boo_lps0)./cmip6.trend_stemp),72,64)',1);  aslmap(8,rlat65,rlon73,maskLFmatr.*junk,[-90 +90],[-180 +180]);   
  colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/K CMIP6');  
figure(9); clf
  junk = smoothn(reshape(((boo_lpsERA5-boo_lps0)./era5.trend_stemp),72,64)',1);  aslmap(9,rlat65,rlon73,maskLFmatr.*junk,[-90 +90],[-180 +180]);   
  colormap(llsmap5); caxis([-1 +1]);  title('\delta lapse rate K/km/K ERA5');  
