function [lag,ratio,deltaTuncx] = find_72_lag(deltaTunc);

[ll0 mm0 nn0] = size(deltaTunc);

for ll = 1 : ll0
  for latlat = 1 : 64
    data = squeeze(deltaTunc(ll,:,latlat));
    k = remove_nan(data);
    if length(k) > 10
      l = xcorr(data(k),1,'coeff');
      lag(ll,latlat) = l(1);
    else
      lag(ll,latlat) = NaN;
    end
  end
end
ratio = sqrt((1+lag)./(1-lag));

deltaTuncx = squeeze(nanmean(deltaTunc,2))/sqrt(72);
deltaTuncx = sqrt(squeeze(nansum(deltaTunc.^2,2)))/(72);

