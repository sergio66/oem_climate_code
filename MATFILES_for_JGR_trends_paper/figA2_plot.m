%load figa2_data;  % mtime resid_670wn
load figA2_data;   % mtime resid_670wn

%figure;
figure(1); clf

plot(mtime,resid_670wn,'+-')
grid
ylabel('Fit BT Residual (K)')
