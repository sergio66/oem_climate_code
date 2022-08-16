%% https://machinelearningmastery.com/arima-for-time-series-forecasting-with-python/

figure(7)
shampoo = load('shampoo.txt')
plot(1:length(shampoo),shampoo(:,3))

data = shampoo(:,3);

lx = xcorr(data,1,'coeff');
[lx,lagsx] = xcorr(data,18,'coeff');
[lx,lagsx] = xcorr(data,[],'coeff');
plot(lagsx(19:end),lx(19:end))

[la,lagsa] = autocorr(data,length(data)-1);
plot(lagsa,l); plotaxis2;
plot(lagsa(2:end),la(2:end),lagsx(19:end),lx(19:end)); plotaxis2;

plot(lagsa(2:length(la)),la(2:length(la)),lagsx(19:length(lagsx)),lx(19:length(lagsx))); plotaxis2;
plot(lagsa(2:length(la)),la(2:length(la)),lagsx(19:length(lagsx)),sqrt(lx(19:length(lagsx)))); plotaxis2;
plot(lagsa(2:length(la)),la(2:length(la)),lagsx(19:length(lagsx)),lx(19:length(lagsx)).^2); plotaxis2;

%% A random variable that is a time series is stationary if its
%% statistical properties are all constant over time.  A stationary
%% series has no trend, its variations around its mean have a constant
%% amplitude, and it wiggles in a consistent fashion,
