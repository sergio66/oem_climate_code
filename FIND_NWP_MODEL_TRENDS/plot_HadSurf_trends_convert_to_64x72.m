addpath /asl/matlib/maps
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS

%% see plot_HadSurf_trends_convert_to_64x72_prototype.m
%   fh = equal_area_map(fn, glat, glon, gval, latlim, lonlim, opts)
%
% INPUTS
%   fn    - matlab figure number
%   glat  - n+1 vector of latitude boundaries
%   glon  - m+1 vector of longitude boundaries
%   gval  - n x m array of map data values
%

load llsmap5
llsmap5nan = llsmap5; llsmap5nan(1,:) = [0.6 0.6 0.6];

load hadISDH_trends.mat

%% try to interp into 4608
[raaHadX,raaHadY] = meshgrid(-177.5:5:+177.5,-87.5:5:+87.5);
[raaLon72,raaLat64] = meshgrid(rlon,rlat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

moo = reshape(hadISDH_results.trend_t,36,72);
hadISDH_results.trend_t_64x72 = (interp2(raaHadX,raaHadY,moo,raaLon72,raaLat64));
moo = hadISDH_results.trend_t_64x72; oo1 = find(isnan(moo));
boo = smoothn(moo,1);
zoo = boo; zoo(oo1) = -9999;;
figure(1); aslmap(1,rlat65,rlon73,zoo,[-90 +90],[-180 +180]); 
  colormap(llsmap5nan); title('HadISDH dST/dt K/yr'); ccaxis([-1.5 +1.5]/10)

poo = smoothn((reshape(results(:,6),72,64)'),1); %% correct  64x72
koo = poo; koo(oo1) = -9999;
figure(2); aslmap(2,rlat65,rlon73,koo,[-90 +90],[-180 +180]); colormap(llsmap5nan);  
  ccaxis([-0.15 +0.15]); title('UMBC dST/dt K/yr')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

moo = reshape(hadISDH_results.trend_rh,36,72);
hadISDH_results.trend_rh_64x72 = (interp2(raaHadX,raaHadY,moo,raaLon72,raaLat64));
moo = hadISDH_results.trend_rh_64x72; oo4 = find(isnan(moo));
boo = smoothn(moo,1);
zoo = boo; zoo(oo4) = -9999;;
figure(3); aslmap(3,rlat65,rlon73,zoo,[-90 +90],[-180 +180]); 
  colormap(llsmap5nan); title('HadISDH dRH/dt percent/yr'); ccaxis([-0.5 +0.5])

poo = smoothn((reshape(RHSurfpert-RHSurf0,72,64)'),1);
koo = poo; koo(oo1) = -9999;
figure(4); aslmap(4,rlat65,rlon73,koo,[-90 +90],[-180 +180]); colormap(llsmap5nan);  
  ccaxis([-0.5 +0.5]); title('UMBC dRH/dt percent/yr')

