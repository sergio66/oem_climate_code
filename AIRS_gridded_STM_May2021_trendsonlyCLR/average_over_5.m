function y = average_over_5(x,N,iNlays)

for ii = 1:iNlays
  ind = (1:N) + (ii-1)*N;
  y(ii) = mean(x(ind));
end
