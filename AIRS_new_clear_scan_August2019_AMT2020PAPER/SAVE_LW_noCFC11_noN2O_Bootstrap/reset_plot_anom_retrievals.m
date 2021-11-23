for ii = 1 : 6
  figure(ii); clf
end

iaaFound = 0 * iaaFound;
co2 = nan * co2;
n2o = nan * n2o;
ch4 = nan * ch4;
cfc11 = nan * cfc11;
stemp = nan * stemp;
dofs = nan * dofs;
bestloop = nan * bestloop;

if exist('aco2')
  aco2 = aco2 * nan;
  bco2 = bco2 * nan;
  cco2 = cco2 * nan;
end
