addpath /home/sergio/MATLABCODE/oem_pkg

x = 0 : 1: 100;

k = 4*pi/100; y = cos(k*x);
k = 2*pi/100; y = sin(k*x);
k = 6*pi/100; y = sin(k*x);
%k = 2*pi/100; y = sin(k*x*k.*x);

dx = mean(diff(x));

dy1 = diff(y);
dy2 = diff(dy1);

L1 = get_l(length(x),1);
tik = L1'*L1;

knew = (k/dx)^2;
plot(x,y/400,'ko-',x,k*tik*y','rx-',x(2:end-1),k*dy2,'b'); grid
plot(x,k*tik*y','rx-',x(2:end-1),-k*dy2,'b'); grid


