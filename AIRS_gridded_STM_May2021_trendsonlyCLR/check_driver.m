miaow.debug          = 'false';     %% no debuging
miaow.debug_dir      = '../Debug/'; %% no debuging
miaow.iQuantile         = +16;      %% 1-16 for coldest-warmest, 00 for mean
miaow.NorD = +1;                    %% [+1] night/desc   -1 : day/asc
miaow.iDebugRatesUseNWP = -1;       %% [-1] use AIRS observed trends
                                       %% 31,32 use reconstructed spectral trends from night/day AIRS L3
                                       %% 51,52 use reconstructed spectral trends from night/day ERA5
                                       %% 61,62 use reconstructed spectral trends from night/day CMIP6 (technically they are the same)

driver_allowedparams = [{'debug'},{'debug_dir'},{'iDebugRatesUseNWP'},{'NorD'},{'iQuantile'}];

%disp('miaow before')
%miaow

optvarD = fieldnames(driver);
optvarA = driver_allowedparams;
for i = 1 : length(optvarA)
 iFound = -1;
 for j = 1 : length(optvarD)
   if (length(intersect(optvarD{j},optvarA{i})) == 1)
     iFound = +1;
    end
  end
  if iFound == -1
    eval(sprintf('driver.%s = miaow.%s;', optvarA{i}, optvarA{i}));
    fprintf(1,'  check_driver.m : setting driver.%s \n',optvarA{i});
  end
end

%disp('driver after')
%driver

