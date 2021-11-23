addpath /asl/matlib/maps
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

load llsmap5
llsmap5nan = llsmap5; llsmap5nan(1,:) = [0.6 0.6 0.6];

load hadISDH_trends.mat

moo = reshape(hadISDH_results.trend_t,36,72); oo1 = find(isnan(moo));
boo = smoothn(moo,1);
zoo = boo; zoo(oo1) = -9999;;
figure(1); aslmap(1,-90:5:+90,-180:5:+180,smoothn((reshape(hadISDH_results.trend_t,36,72)),1),[-90 +90],[-180 +180]); 
figure(1); aslmap(1,-90:5:+90,-180:5:+180,moo,[-90 +90],[-180 +180]); 
figure(1); aslmap(1,-90:5:+90,-180:5:+180,boo,[-90 +90],[-180 +180]); 
figure(1); aslmap(1,-90:5:+90,-180:5:+180,zoo,[-90 +90],[-180 +180]); 
  colormap(llsmap5nan); title('HadISDH dST/dt K/yr'); ccaxis([-1.5 +1.5]/10)

moo = reshape(hadISDH_results.trend_twet,36,72); oo2 = find(isnan(moo));
boo = smoothn(moo,1);
zoo = boo; zoo(oo2) = -9999;;
figure(2); aslmap(2,-90:5:+90,-180:5:+180,smoothn((reshape(hadISDH_results.trend_twet,36,72)),1),[-90 +90],[-180 +180]); 
figure(2); aslmap(2,-90:5:+90,-180:5:+180,moo,[-90 +90],[-180 +180]); 
figure(2); aslmap(2,-90:5:+90,-180:5:+180,boo,[-90 +90],[-180 +180]); 
figure(2); aslmap(2,-90:5:+90,-180:5:+180,zoo,[-90 +90],[-180 +180]); 
  colormap(llsmap5nan); title('HadISDH dSTwet/dt K/yr'); ccaxis([-1.5 +1.5]/10)

moo = reshape(hadISDH_results.trend_fracQ,36,72); oo3 = find(isnan(moo));
boo = smoothn(moo,1);
zoo = boo; zoo(oo3) = -9999;;
figure(3); aslmap(3,-90:5:+90,-180:5:+180,smoothn((reshape(hadISDH_results.trend_fracQ,36,72)),1),[-90 +90],[-180 +180]); 
figure(3); aslmap(3,-90:5:+90,-180:5:+180,moo,[-90 +90],[-180 +180]); 
figure(3); aslmap(3,-90:5:+90,-180:5:+180,boo,[-90 +90],[-180 +180]); 
figure(3); aslmap(3,-90:5:+90,-180:5:+180,zoo,[-90 +90],[-180 +180]); 
  colormap(llsmap5nan); title('HadISDH d fracQ/dt /yr'); ccaxis([-1 +1])

moo = reshape(hadISDH_results.trend_Q,36,72); oo3 = find(isnan(moo));
boo = smoothn(moo,1);
zoo = boo; zoo(oo3) = -9999;;
figure(3); aslmap(3,-90:5:+90,-180:5:+180,smoothn((reshape(hadISDH_results.trend_Q,36,72)),1),[-90 +90],[-180 +180]); 
figure(3); aslmap(3,-90:5:+90,-180:5:+180,moo,[-90 +90],[-180 +180]); 
figure(3); aslmap(3,-90:5:+90,-180:5:+180,boo,[-90 +90],[-180 +180]); 
figure(3); aslmap(3,-90:5:+90,-180:5:+180,zoo,[-90 +90],[-180 +180]); 
  colormap(llsmap5nan); title('HadISDH d Q/dt /yr'); ccaxis([-0.1 +0.1])

moo = reshape(hadISDH_results.trend_rh,36,72); oo4 = find(isnan(moo));
boo = smoothn(moo,1);
zoo = boo; zoo(oo4) = -9999;;
figure(4); aslmap(4,-90:5:+90,-180:5:+180,smoothn((reshape(hadISDH_results.trend_rh,36,72)),1),[-90 +90],[-180 +180]); 
figure(4); aslmap(4,-90:5:+90,-180:5:+180,moo,[-90 +90],[-180 +180]); 
figure(4); aslmap(4,-90:5:+90,-180:5:+180,boo,[-90 +90],[-180 +180]); 
figure(4); aslmap(4,-90:5:+90,-180:5:+180,zoo,[-90 +90],[-180 +180]); 
  colormap(llsmap5nan); title('HadISDH dRH/dt percent/yr'); ccaxis([-0.5 +0.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('hadCRU4_trends.mat')
  load hadCRU4_trends.mat
  
  moo = reshape(hadCRU4_results.trend_t,36,72); oo5 = find(isnan(moo));
  boo = smoothn(moo,1);
  zoo = boo; zoo(oo5) = -9999;;
  figure(5); aslmap(5,-90:5:+90,-180:5:+180,smoothn((reshape(hadCRU4_results.trend_t,36,72)),1),[-90 +90],[-180 +180]); 
  figure(5); aslmap(5,-90:5:+90,-180:5:+180,moo,[-90 +90],[-180 +180]); 
  figure(5); aslmap(5,-90:5:+90,-180:5:+180,boo,[-90 +90],[-180 +180]); 
  figure(5); aslmap(5,-90:5:+90,-180:5:+180,zoo,[-90 +90],[-180 +180]); 
    colormap(llsmap5nan); title('<HadCRU4> dST/dt K/yr'); ccaxis([-1.5 +1.5]/10)
end
