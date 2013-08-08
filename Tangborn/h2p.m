function pin = p2h(hin)
%this function takes height in m, anc converts to pressure in mb

%load /asl/matlab/kcarta/airsheights.dat
%load /asl/matlab/kcarta/airslevels.dat
load airsheights.dat
load airslevels.dat

h=airsheights;
p=airslevels;
for ii=1:100
  pavg(ii)=(p(ii+1)-p(ii))/log(p(ii+1)/p(ii));
  end

pin=interp1(h,pavg,hin);
