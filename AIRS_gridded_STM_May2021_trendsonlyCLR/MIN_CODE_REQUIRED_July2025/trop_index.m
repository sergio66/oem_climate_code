iVers = 0;  %% orig
iVers = 1;  %% new

if iVers == 0
  trpi = [59 58 56 43 44 44 44 44 44 44 ...
          44 44 44 44 44 44 44 44 44 44 ...
          44 44 44 44 44 44 44 44 44 44 ...
          44 44 40 39 58 5959 59 59 59];
  trpi = 45 * ones(1,72*64);
  trpi = 9 * ones(1,72*64);   %% have reduced from 100 to 20 layers default, during AIRS STM Oct 2020, 20 layers ==> 45 --> 45/5 = 9
  
  %% trying these
  trpi = 5 * ones(1,72*64);    %% have reduced from 100 to 20 layers default, during AIRS STM Oct 2020
  trpi = 18 * ones(1,72*64);   %% have reduced from 100 to 20 layers default, during AIRS STM Oct 2020
  trpi = 1  * ones(1,72*64);   %% have reduced from 100 to 20 layers default, during AIRS STM Oct 2020, used till Aug 2021

elseif iVers == 1  
  %{
  %% after doing [h,ha,p,pa] = rtpread('avgprof');
  %% plot_profile_trends2.m:105:      [mmw0,lps0] = mmwater_rtp_pstop_lapse(h,p);
  lps_index.rlat = p.rlat;
  lps_index.rlon = p.rlon;
  lps_index.trp_ind = lps0.trp_ind;
  save lps_save.mat lps_index
  %}
  
  junk = load('lps_save.mat');
  %trpi = reshape(junk.lps_index.trp_ind,72,64);
  trpi = junk.lps_index.trp_ind;
  
  clear junk
end
