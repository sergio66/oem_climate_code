function [zscore,pvalueT] = find_z_score(deltaTlat,xTunc,ratio);

%normcdf([-1.96 0 +1.96])

zscore = (abs(deltaTlat) - 0)./xTunc;
pvalueT = 1-normcdf(zscore);
pvalueT = 2*(1-tcdf(zscore,72-1));
pvalueT = 2*(1-tcdf(zscore,(72-2) * ratio'));
