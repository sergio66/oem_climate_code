function y = average_over_5(x,N,iNlays)

for ii = 1:iNlays
  ind = (1:N) + (ii-1)*N;
  moo = find(ind <= length(x));
  ind = ind(moo);
  if length(ind) > 1
    y(ii) = mean(x(ind));
  else
    y(ii) = x(ind);
  end
end
