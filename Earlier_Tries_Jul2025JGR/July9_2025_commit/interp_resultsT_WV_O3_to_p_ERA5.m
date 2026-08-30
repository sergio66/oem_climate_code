pertERA5     = p;
pertERA5_unc = p; 

nanBelowSurf101 = ones(101,4608);

pertERA5.stemp = pertERA5.stemp + era5.trend_stemp;
pertERA5.ptemp(1:100,:) = pertERA5.ptemp(1:100,:) + era5.trend_ptemp;
pertERA5.gas_1(1:100,:) = pertERA5.gas_1(1:100,:).*(1 + era5.trend_gas_1);
pertERA5.gas_3(1:100,:) = pertERA5.gas_3(1:100,:).*(1 + era5.trend_gas_3);

nlays = 101;
for ii = 1 : length(p.stemp)
  pertERA5.gas_2(1:nlays,ii) =  pertERA5.gas_2(1:nlays,ii) .* (1+2.2/385);
  pertERA5.gas_4(1:nlays,ii) =  pertERA5.gas_4(1:nlays,ii) .* (1+0.8/300);
  pertERA5.gas_6(1:nlays,ii) =  pertERA5.gas_6(1:nlays,ii) .* (1+4.5/1700);
  pertERA5_unc.gas_2(1:nlays,ii) =  pertERA5_unc.gas_2(1:nlays,ii) .* (1+2.2/385);
  pertERA5_unc.gas_4(1:nlays,ii) =  pertERA5_unc.gas_4(1:nlays,ii) .* (1+0.8/300);
  pertERA5_unc.gas_6(1:nlays,ii) =  pertERA5_unc.gas_6(1:nlays,ii) .* (1+4.5/1700);
end
