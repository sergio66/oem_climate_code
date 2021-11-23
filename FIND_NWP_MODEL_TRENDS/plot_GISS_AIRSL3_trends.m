airsL3_2 = load('~/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ChrisHTrends/airsL3_trends.mat');
giss = load('~/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/ChrisHTrends/giss_trends.mat');

figure(51); aslmap(51,rlat65,rlon73,smoothn(airsL3_2.airsL3_trend4608',1), [-90 +90],[-180 +180]);  colormap(llsmap5); caxis([-0.15 +0.15]); title([strNorD ' RHsurf d/dt AIRS L3v2 /yr']);   
figure(52); aslmap(52,rlat65,rlon73,smoothn(giss.giss_trend4608',1), [-90 +90],[-180 +180]);  colormap(llsmap5); caxis([-0.15 +0.15]); title([strNorD ' stemp d/dt GISS K/yr']);   
