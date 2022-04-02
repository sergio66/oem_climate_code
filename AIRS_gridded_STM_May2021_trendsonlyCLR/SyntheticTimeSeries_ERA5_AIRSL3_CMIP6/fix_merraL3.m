function   p72 = fix_merraL3(p72_0);

p72 = p72_0;

for ii = 1 : length(p72.stemp)
  plevs = p72.plevs(:,ii);

  tlevs = p72.ptemp(:,ii);
  gas_1 = p72.gas_1(:,ii);
  gas_3 = p72.gas_3(:,ii);
  spres = p72.spres(ii);
  stemp = p72.stemp(ii);

  %% just fix everything using extrapolations
  booBad = find(tlevs == -9999);
  booGood = find(tlevs > 150);
  tlevs(booBad) = interp1(log(plevs(booGood)),tlevs(booGood),log(plevs(booBad)),[],'extrap');
  gas_1(booBad) = exp(interp1(log(plevs(booGood)),log(gas_1(booGood)),log(plevs(booBad)),[],'extrap'));
  gas_3(booBad) = exp(interp1(log(plevs(booGood)),log(gas_3(booGood)),log(plevs(booBad)),[],'extrap'));
  p72.ptemp(:,ii) = tlevs;
  p72.gas_1(:,ii) = gas_1;
  p72.gas_3(:,ii) = gas_3;
end
