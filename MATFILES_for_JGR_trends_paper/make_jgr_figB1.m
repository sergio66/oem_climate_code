%% see  /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_compare_trends_Day_vs_Night.m  --> compare_SKT_trends_Day_vs_Night.m
load figB1.mat

load llsmap5

figure(1); clf;
    iFig = 1;
    aslmap_3x4tiledlayout(figB1_sktDNtrends.umbc_D,figB1_sktDNtrends.airsL3_D,figB1_sktDNtrends.climcaps_D,figB1_sktDNtrends.era5_D,...
                          figB1_sktDNtrends.umbc_N,figB1_sktDNtrends.airsL3_N,figB1_sktDNtrends.climcaps_N,figB1_sktDNtrends.era5_N,...
                          figB1_sktDNtrends.umbc_X,figB1_sktDNtrends.airsL3_X,figB1_sktDNtrends.climcaps_X,figB1_sktDNtrends.era5_X,...
                          iFig,figB1_sktDNtrends.plotoptions);
