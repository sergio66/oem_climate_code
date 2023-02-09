pert     = p;
pert_unc = p; 

nanBelowSurf101 = ones(101,4608);

for ii = 1 : length(p.stemp)
  nlays = p.nlevs(ii)-1;
  playsjunk = p.plays(1:nlays,ii);

  nanBelowSurf101(p.nlevs(ii):101,ii) = nan;

  roo = interp1(log(pavg),resultsT(ii,:),log(playsjunk),[],'extrap');
  goo = find(isfinite(roo)); boo = find(~isfinite(roo));
  if length(boo) > 0
    roo(boo) = interp1(log(playsjunk(goo)),roo(goo),log(playsjunk(boo)),[],'extrap');
  end
    pert.ptemp(1:nlays,ii)     = pert.ptemp(1:nlays,ii) + roo;
    pert.ptemp_unc(:,ii)       = 0 * pert.ptemp(:,ii);
    roounc = interp1(log(pavg),resultsTunc(ii,:),log(playsjunk),[],'extrap');  
    pert.ptemp_unc(1:nlays,ii) = roounc;
    pert_unc.ptemp(1:nlays,ii) = pert_unc.ptemp(1:nlays,ii) + roo + roounc;

  roo = interp1(log(pavg),resultsWV(ii,:),log(playsjunk),[],'extrap');
  goo = find(isfinite(roo)); boo = find(~isfinite(roo));
  if length(boo) > 0
    roo(boo) = interp1(log(playsjunk(goo)),roo(goo),log(playsjunk(boo)),[],'extrap');
  end
    pert.gas_1(1:nlays,ii)     = pert.gas_1(1:nlays,ii) .* (1+roo);
    pert.gas_1_unc(:,ii)       = 0 * pert.gas_1(:,ii);
    roounc = interp1(log(pavg),resultsWVunc(ii,:),log(playsjunk),[],'extrap');
    pert.gas_1_unc(1:nlays,ii) = pert.gas_1(1:nlays,ii) .* (0+roounc);
    pert_unc.gas_1(1:nlays,ii) = pert_unc.gas_1(1:nlays,ii) .* (1+roo+roounc);

  roo = interp1(log(pavg),resultsO3(ii,:),log(playsjunk),[],'extrap');
  goo = find(isfinite(roo)); boo = find(~isfinite(roo));
  if length(boo) > 0
    roo(boo) = interp1(log(playsjunk(goo)),roo(goo),log(playsjunk(boo)),[],'extrap');
  end
    pert.gas_3(1:nlays,ii)     = pert.gas_3(1:nlays,ii) .* (1+roo);
    pert.gas_3_unc(:,ii)       = 0 * pert.gas_3(:,ii);
    roounc = interp1(log(pavg),resultsO3unc(ii,:),log(playsjunk),[],'extrap');
    pert.gas_3_unc(1:nlays,ii) = pert.gas_3(1:nlays,ii) .* (0+roounc);
    pert_unc.gas_3(1:nlays,ii) =  pert_unc.gas_3(1:nlays,ii) .* (1+roo+roounc);

  pert.stemp(ii) = pert.stemp(ii) + results(ii,6);
  pert.stemp_unc(ii) = results(ii,6);
  pert_unc.stemp(ii) = pert_unc.stemp(ii) + resultsunc(ii,6);

  pert.gas_2(1:nlays,ii) =  pert.gas_2(1:nlays,ii) .* (1+2.2/385);
  pert.gas_4(1:nlays,ii) =  pert.gas_4(1:nlays,ii) .* (1+0.8/300);
  pert.gas_6(1:nlays,ii) =  pert.gas_6(1:nlays,ii) .* (1+4.5/1700);
  pert_unc.gas_2(1:nlays,ii) =  pert_unc.gas_2(1:nlays,ii) .* (1+2.2/385);
  pert_unc.gas_4(1:nlays,ii) =  pert_unc.gas_4(1:nlays,ii) .* (1+0.8/300);
  pert_unc.gas_6(1:nlays,ii) =  pert_unc.gas_6(1:nlays,ii) .* (1+4.5/1700);
end
