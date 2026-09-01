load /asl/s1/sergio/AIRS_L3/airsL3_v7_64x72_rates_stats_Sept2002_Aug2021_19yr_desc.mat

%% see eg driver_olr_fluxchanges_UMBC.m
cs = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX100_50fatlayers_AIRSL3_ERA5_CMIP6_feedback_olr_ceres_UMBC.mat');
cs = load('/asl/s1/sergio/JUNK/gather_tileCLRnight_Q16_newERA5_2021jacs_startwith0_uncX3_50fatlayers_AIRSL3_ERA5_CMIP6_feedback_TongaChans_olr_ceres_UMBC.mat');

figure(1); clf; 
plot(meanvaluebin(rlat),mean(thestats64x72_other.clrolrrate,1)*20,...
     cs.outflux.ceres_trends.lat,cs.outflux.ceres_trends.trend_ceres_lw_clr*20,'rx-',...
     cs.outflux.umbc_trends.rlat,cs.outflux.umbc_trends.umbc_all*20,'k','linewidth',2); plotaxis2;
hl = legend('AIRS L3','CERES','UMBC','location','best');

r64 = meanvaluebin(rlat);
junkceres = interp1(cs.outflux.ceres_trends.lat,cs.outflux.ceres_trends.trend_ceres_lw_clr,r64);

figure(1); clf; 
plot(r64,mean(thestats64x72_other.clrolrrate,1)*20,...
     r64,junkceres*20,'rx-',...
     r64,cs.outflux.umbc_trends.umbc_all*20,'k','linewidth',2); plotaxis2;
hl = legend('AIRS L3','CERES','UMBC','location','best');
title('OLR change in 20 years');

figure(2); clf
plot(r64,(junkceres' - mean(thestats64x72_other.clrolrrate,1))*20,...
     r64,(junkceres' - cs.outflux.umbc_trends.umbc_all)*20,'k','linewidth',2); plotaxis2;
hl = legend('AIRS L3','UMBC','location','best');
title('OLR change CERES-X')
