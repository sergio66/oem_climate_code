iVers = 0;  %% orig
iVers = 1;  %% new

if iVers == 0
  trpi = 45 * ones(1,64);
  trpi = 9 * ones(1,64);       %% have reduced from 100 to 20 layers default, during AIRS STM Oct 2020, 20 layers ==> 45 --> 45/5 = 9
  
  %% trying these
  trpi = 5 * ones(1,64);    %% have reduced from 100 to 20 layers default, during AIRS STM Oct 2020
  trpi = 18 * ones(1,64);   %% have reduced from 100 to 20 layers default, during AIRS STM Oct 2020
  trpi = 1  * ones(1,64);   %% have reduced from 100 to 20 layers default, during AIRS STM Oct 2020, used till Aug 2021

elseif iVers == 1  
  junk = load('lps_save.mat');
  trpi = mean(reshape(junk.lps_index.trp_ind,72,64),1);
  
  clear junk
end
