for ii = 1 : length(p.stemp)
  nlays = p.nlevs(ii)-1;
  playsjunk = p.plays(1:nlays,ii);

  ptempjunk = p.ptemp(1:nlays,ii);
  roo0 = interp1(log(playsjunk),ptempjunk,log(p.spres(ii)),[],'extrap');
  roo  = interp1(log(pavg),resultsT(ii,:),log(p.spres(ii)),[],'extrap');
  pert.a2mtemp_orig(ii) =  roo0;
  pert.a2mtemp(ii) =  roo0 + roo;

  gas1junk = p.gas_1(1:nlays,ii);
  roo0 = interp1(log(playsjunk),log(gas1junk),log(p.spres(ii)),[],'extrap');
  roo0 = exp(roo0);
  pert.a2mgas_1_orig(ii) =  roo0;
  roo  = interp1(log(pavg),resultsT(ii,:),log(p.spres(ii)),[],'extrap');
  pert.a2mgas_1(ii) =  roo0 * (1+roo);
end
