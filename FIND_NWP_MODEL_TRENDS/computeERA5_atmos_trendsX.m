%% this is only 100 layer trends ptemp,gas_1,gas_3,RH
%% in case computeERA5_atmos_trends.m forgot to do nanmean for gas_1,gas_3

disp(' ')
disp('<<< commencing atmos trends +=1000,x=100,.=10 >>> ')

warning off
for ii = 1 : 4608
  if mod(ii,1000) == 0
    fprintf(1,'+ \n')
  elseif mod(ii,100) == 0
    fprintf(1,'x')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end

  %fprintf(1,'doing 100 layers OP ptemp,rh,gas_1,gas_3 trends ii = %4i of 4608 \n',ii)
  for ll = 1 : 100    
    data = squeeze(pall.ptemp(:,ll,ii));  
    boo = find(isfinite(data));
    boo = intersect(boo,thetimeSeason);
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),iNumCycles); trend_ptemp(ll,ii) = B(2);  trend_ptemp_err(ll,ii) = stats.se(2);
    else
      trend_ptemp(ll,ii) = NaN; 
      trend_ptemp_err(ll,ii) = NaN;
    end

    data = squeeze(pall.gas_1(:,ll,ii));  data = data/nanmean(data); 
    boo = find(isfinite(data));
    boo = intersect(boo,thetimeSeason);
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),iNumCycles); trend_gas_1(ll,ii) = B(2);  trend_gas_1_err(ll,ii) = stats.se(2);
    else
      trend_gas_1(ll,ii) = NaN; 
      trend_gas_1_err(ll,ii) = NaN;
    end

    data = squeeze(pall.gas_3(:,ll,ii));  data = data/nanmean(data); 
    boo = find(isfinite(data));
    boo = intersect(boo,thetimeSeason);
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),iNumCycles); trend_gas_3(ll,ii) = B(2);  trend_gas_3_err(ll,ii) = stats.se(2);
    else
      trend_gas_3(ll,ii) = NaN; 
      trend_gas_3_err(ll,ii) = NaN;
    end

    data = squeeze(pall.RH(:,ll,ii));                             
    boo = find(isfinite(data));
    boo = intersect(boo,thetimeSeason);
    if length(boo) > 20
      [B, stats] = Math_tsfit_lin_robust(dayOFtime(boo),data(boo),iNumCycles); trend_RH(ll,ii) = B(2);     trend_RH_err(ll,ii) = stats.se(2);
    else
      trend_RH(ll,ii) = NaN; 
      trend_RH_err(ll,ii) = NaN;
    end
  end
end
fprintf(1,'\n')
warning on
