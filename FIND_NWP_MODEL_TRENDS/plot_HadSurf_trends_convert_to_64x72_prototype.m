addpath /asl/matlib/maps
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

load llsmap5
llsmap5nan = llsmap5; llsmap5nan(1,:) = [0.6 0.6 0.6];

load hadISDH_trends.mat

%% try to interp into 4608
[raaHadY,raaHadX] = meshgrid(-87.5:5:+87.5,-177.5:5:+177.5);
[raaLat64,raaLon72] = meshgrid(rlat,rlon);

[raaHadX,raaHadY] = meshgrid(-177.5:5:+177.5,-87.5:5:+87.5);
[raaLon72,raaLat64] = meshgrid(rlon,rlat);

%   fh = equal_area_map(fn, glat, glon, gval, latlim, lonlim, opts)
%
% INPUTS
%   fn    - matlab figure number
%   glat  - n+1 vector of latitude boundaries
%   glon  - m+1 vector of longitude boundaries
%   gval  - n x m array of map data values
%

for ii=1:6; figure(ii); clf; end

%% this works!!!! 36 x 72
moo = reshape(hadISDH_results.trend_t,36,72); oo2 = find(isnan(moo));
boo = smoothn(moo,1);
zoo = boo; zoo(oo2) = -9999;;
figure(2); aslmap(2,-90:5:+90,-180:5:+180,smoothn((reshape(hadISDH_results.trend_twet,36,72)),1),[-90 +90],[-180 +180]); 
figure(2); aslmap(2,-90:5:+90,-180:5:+180,moo,[-90 +90],[-180 +180]); 
figure(2); aslmap(2,-90:5:+90,-180:5:+180,boo,[-90 +90],[-180 +180]); 
figure(2); aslmap(2,-90:5:+90,-180:5:+180,zoo,[-90 +90],[-180 +180]); 
  colormap(llsmap5nan); title('HadISDH dST/dt K/yr'); ccaxis([-1.5 +1.5]/10)

moo = reshape(hadISDH_results.trend_t,36,72);
hadISDH_results.trend_t_64x72 = (interp2(raaHadX,raaHadY,moo,raaLon72,raaLat64));
moo = hadISDH_results.trend_t_64x72; oo1 = find(isnan(moo));
%boo = smoothn(moo,1);
boo = moo;
zoo = boo; zoo(oo1) = -9999;;
%zoo = moo;
figure(3); aslmap(3,rlat65,rlon73,zoo,[-90 +90],[-180 +180]); 
  colormap(llsmap5nan); title('HadISDH dST/dt K/yr'); ccaxis([-1.5 +1.5]/10)

moo = reshape(hadISDH_results.trend_t,36,72);
hadISDH_results.trend_t_64x72 = (interp2(raaHadX,raaHadY,moo,raaLon72,raaLat64));
moo = hadISDH_results.trend_t_64x72; oo1 = find(isnan(moo));
boo = smoothn(moo,1);
zoo = boo; zoo(oo1) = -9999;;
figure(4); aslmap(4,rlat65,rlon73,zoo,[-90 +90],[-180 +180]); 
  colormap(llsmap5nan); title('HadISDH dST/dt K/yr'); ccaxis([-1.5 +1.5]/10)

whos moo boo zoo poo raaHadX

poo = smoothn((reshape(results(:,6),64,72)),1);  %% completely wrong
poo = smoothn((reshape(results(:,6),72,64)'),1); %% correct  64x72

figure(6); aslmap(6,rlat65,rlon73,poo,[-90 +90],[-180 +180]); colormap(llsmap5);  
  caxis([-0.15 +0.15]); title('UMBC')

koo = poo; koo(oo1) = -9999;
figure(5); aslmap(5,rlat65,rlon73,koo,[-90 +90],[-180 +180]); colormap(llsmap5nan);  
  caxis([-0.15 +0.15]); title('UMBC')

error(';ksjgfkjs')
