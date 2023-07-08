pertxb     = p;
pertxb_unc = p; 

xbunc   = 0 * xb;
xbTunc  = 0 * xbT;
xbWVunc = 0 * xbT;
xbO3unc = 0 * xbT;

nanBelowSurf101 = ones(101,4608);

for ii = 1 : length(p.stemp)
  nlays = p.nlevs(ii)-1;
  playsjunk = p.plays(1:nlays,ii);

  nanBelowSurf101(p.nlevs(ii):101,ii) = nan;

  roo = interp1(log(pavg),xbT(ii,:),log(playsjunk),[],'extrap');
  goo = find(isfinite(roo)); boo = find(~isfinite(roo));
  if length(boo) > 0
    roo(boo) = interp1(log(playsjunk(goo)),roo(goo),log(playsjunk(boo)),[],'extrap');
  end
    pertxb.ptemp(1:nlays,ii)     = pertxb.ptemp(1:nlays,ii) + roo;
    pertxb.ptemp_unc(:,ii)       = 0 * pertxb.ptemp(:,ii);
    roounc = interp1(log(pavg),xbTunc(ii,:),log(playsjunk),[],'extrap');  
    pertxb.ptemp_unc(1:nlays,ii) = roounc;
    pertxb_unc.ptemp(1:nlays,ii) = pertxb_unc.ptemp(1:nlays,ii) + roo + roounc;

  roo = interp1(log(pavg),xbWV(ii,:),log(playsjunk),[],'extrap');
  goo = find(isfinite(roo)); boo = find(~isfinite(roo));
  if length(boo) > 0
    roo(boo) = interp1(log(playsjunk(goo)),roo(goo),log(playsjunk(boo)),[],'extrap');
  end
    pertxb.gas_1(1:nlays,ii)     = pertxb.gas_1(1:nlays,ii) .* (1+roo);
    pertxb.gas_1_unc(:,ii)       = 0 * pertxb.gas_1(:,ii);
    roounc = interp1(log(pavg),xbWVunc(ii,:),log(playsjunk),[],'extrap');
    pertxb.gas_1_unc(1:nlays,ii) = pertxb.gas_1(1:nlays,ii) .* (0+roounc);
    pertxb_unc.gas_1(1:nlays,ii) = pertxb_unc.gas_1(1:nlays,ii) .* (1+roo+roounc);

  roo = interp1(log(pavg),xbO3(ii,:),log(playsjunk),[],'extrap');
  goo = find(isfinite(roo)); boo = find(~isfinite(roo));
  if length(boo) > 0
    roo(boo) = interp1(log(playsjunk(goo)),roo(goo),log(playsjunk(boo)),[],'extrap');
  end
    pertxb.gas_3(1:nlays,ii)     = pertxb.gas_3(1:nlays,ii) .* (1+roo);
    pertxb.gas_3_unc(:,ii)       = 0 * pertxb.gas_3(:,ii);
    roounc = interp1(log(pavg),xbO3unc(ii,:),log(playsjunk),[],'extrap');
    pertxb.gas_3_unc(1:nlays,ii) = pertxb.gas_3(1:nlays,ii) .* (0+roounc);
    pertxb_unc.gas_3(1:nlays,ii) =  pertxb_unc.gas_3(1:nlays,ii) .* (1+roo+roounc);

  pertxb.stemp(ii) = pertxb.stemp(ii) + xb(ii,6);
  pertxb.stemp_unc(ii) = xb(ii,6);
  pertxb_unc.stemp(ii) = pertxb_unc.stemp(ii) + xbunc(ii,6);

  pertxb.gas_2(1:nlays,ii) =  pertxb.gas_2(1:nlays,ii) .* (1+2.2/385);
  pertxb.gas_4(1:nlays,ii) =  pertxb.gas_4(1:nlays,ii) .* (1+0.8/300);
  pertxb.gas_6(1:nlays,ii) =  pertxb.gas_6(1:nlays,ii) .* (1+4.5/1700);
  pertxb_unc.gas_2(1:nlays,ii) =  pertxb_unc.gas_2(1:nlays,ii) .* (1+2.2/385);
  pertxb_unc.gas_4(1:nlays,ii) =  pertxb_unc.gas_4(1:nlays,ii) .* (1+0.8/300);
  pertxb_unc.gas_6(1:nlays,ii) =  pertxb_unc.gas_6(1:nlays,ii) .* (1+4.5/1700);
end
