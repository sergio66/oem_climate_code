[m,s,m0,s0] = weighted_mean_stddev(umbc.stemprate,cos(p.rlat*pi/180));               fprintf(1,'cosine weighted UMBC SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                 summary_olr_mmw_stats.umbc.stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate,cos(p.rlat*pi/180));                 fprintf(1,'cosine weighted UMBC mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                            summary_olr_mmw_stats.umbc.mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*umbc_mmwrate./umbc_mmw0,cos(p.rlat*pi/180));  fprintf(1,'cosine weighted UMBC mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                          summary_olr_mmw_stats.umbc.mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate+umbc_mmwrate_err,cos(p.rlat*pi/180));                     fprintf(1,'cosine weighted UMBC mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((umbc_mmwrate+umbc_mmwrate_err)./umbc_mmw0),cos(p.rlat*pi/180));  fprintf(1,'cosine weighted UMBC mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_stemp,cos(p.rlat*pi/180));                       fprintf(1,'cosine weighted ERA5 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                  summary_olr_mmw_stats.era5.stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw,cos(p.rlat*pi/180));                         fprintf(1,'cosine weighted ERA5 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                             summary_olr_mmw_stats.era5.mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*olr.trend_mmw./nanmean(olr.xmmw),cos(p.rlat*pi/180));  fprintf(1,'cosine weighted ERA5 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                           summary_olr_mmw_stats.era5.mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw + olr.trend_mmw_err,cos(p.rlat*pi/180));                    fprintf(1,'cosine weighted ERA5 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((olr.trend_mmw + olr.trend_mmw_err)./nanmean(olr.xmmw)),cos(p.rlat*pi/180));  fprintf(1,'cosine weighted ERA5 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(merra2_stemptrend,cos(p.rlat*pi/180));                       fprintf(1,'cosine weighted MERRA2 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                               summary_olr_mmw_stats.merra2.stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend,cos(p.rlat*pi/180));                         fprintf(1,'cosine weighted MERRA2 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                          summary_olr_mmw_stats.merra2.mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*merra2_mmwtrend./nanmean(olr.xmmw),cos(p.rlat*pi/180));  fprintf(1,'cosine weighted MERRA2 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                        summary_olr_mmw_stats.merra2.mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend + merra2_mmwtrend_err,cos(p.rlat*pi/180));                    fprintf(1,'cosine weighted MERRA2 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((merra2_mmwtrend + merra2_mmwtrend_err)./nanmean(olr.xmmw)),cos(p.rlat*pi/180));  fprintf(1,'cosine weighted MERRA2 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%
omi = load('omi_tcwv_trends.mat');
[omiY,omiX] = meshgrid(omi.omi_colwv.lat,omi.omi_colwv.lon);
[m,s,m0,s0] = weighted_mean_stddev(omi.omi_colwv.trend,cos(omiY*pi/180));                         fprintf(1,'cosine weighted OMI      mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(climcaps_mmwtrend,cos(p.rlat*pi/180));                         fprintf(1,'cosine weighted CLIMCAPS mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(airsL3_mmwtrend,cos(p.rlat*pi/180));                           fprintf(1,'cosine weighted AIRSL3   mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                         

boo = find(abs(omiY) <= 30);
[m,s,m0,s0] = weighted_mean_stddev(omi.omi_colwv.trend(boo),ones(size(boo)));                     fprintf(1,'tropical        OMI mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
boo = find(abs(p.rlat) <= 30);
[m,s,m0,s0] = weighted_mean_stddev(climcaps_mmwtrend(boo),ones(size(boo)));                       fprintf(1,'tropical        CLIMCAPS mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(airsL3_mmwtrend(boo),ones(size(boo)));                         fprintf(1,'tropical        AIRSL3   mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                         

%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
boo = find(abs(p.rlat) <= 30);
[m,s,m0,s0] = weighted_mean_stddev(umbc.stemprate(boo),ones(size(boo)));               fprintf(1,'tropical UMBC SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                 summary_olr_mmw_stats.umbc.tropical_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate(boo),ones(size(boo)));                 fprintf(1,'tropical UMBC mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                            summary_olr_mmw_stats.umbc.tropical_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(umbc_mmwrate(boo)./umbc_mmw0(boo)),ones(size(boo)));  fprintf(1,'tropical UMBC mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                          summary_olr_mmw_stats.umbc.tropical_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate(boo)+umbc_mmwrate_err(boo),ones(size(boo)));                     fprintf(1,'tropical UMBC mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((umbc_mmwrate(boo)+umbc_mmwrate_err(boo))./umbc_mmw0(boo)),ones(size(boo)));  fprintf(1,'tropical UMBC mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_stemp(boo),ones(size(boo)));                       fprintf(1,'tropical ERA5 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                  summary_olr_mmw_stats.era5.tropical_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw(boo),ones(size(boo)));                         fprintf(1,'tropical ERA5 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                             summary_olr_mmw_stats.era5.tropical_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(olr.trend_mmw(boo)./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'tropical ERA5 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                           summary_olr_mmw_stats.era5.tropical_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw(boo)+olr.trend_mmw_err(boo),ones(size(boo)));                            fprintf(1,'tropical ERA5 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((olr.trend_mmw(boo)+olr.trend_mmw_err(boo))./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'tropical ERA5 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(merra2_stemptrend(boo),ones(size(boo)));                       fprintf(1,'tropical MERRA2 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                               summary_olr_mmw_stats.merra2.tropical_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend(boo),ones(size(boo)));                         fprintf(1,'tropical MERRA2 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                          summary_olr_mmw_stats.merra2.tropical_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(merra2_mmwtrend(boo)./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'tropical MERRA2 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                        summary_olr_mmw_stats.merra2.tropical_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend(boo) + merra2_mmwtrend_err(boo),ones(size(boo)));                    fprintf(1,'tropical MERRA2 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((merra2_mmwtrend(boo) + merra2_mmwtrend_err(boo))./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'tropical MERRA2 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);

%%%%%
disp(' ')
boo = find(abs(p.rlat) > 30 & abs(p.rlat) <= 60);
[m,s,m0,s0] = weighted_mean_stddev(umbc.stemprate(boo),ones(size(boo)));               fprintf(1,'midlat UMBC SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                 summary_olr_mmw_stats.umbc.midlat_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate(boo),ones(size(boo)));                 fprintf(1,'midlat UMBC mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                            summary_olr_mmw_stats.umbc.midlat_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(umbc_mmwrate(boo)./umbc_mmw0(boo)),ones(size(boo)));  fprintf(1,'midlat UMBC mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                          summary_olr_mmw_stats.umbc.midlat_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate(boo)+umbc_mmwrate_err(boo),ones(size(boo)));                     fprintf(1,'midlat UMBC mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((umbc_mmwrate(boo)+umbc_mmwrate_err(boo))./umbc_mmw0(boo)),ones(size(boo)));  fprintf(1,'midlat UMBC mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_stemp(boo),ones(size(boo)));                       fprintf(1,'midlat ERA5 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                  summary_olr_mmw_stats.era5.midlat_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw(boo),ones(size(boo)));                         fprintf(1,'midlat ERA5 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                             summary_olr_mmw_stats.era5.midlat_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(olr.trend_mmw(boo)./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'midlat ERA5 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                           summary_olr_mmw_stats.era5.midlat_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw(boo)+olr.trend_mmw_err(boo),ones(size(boo)));                            fprintf(1,'midlat ERA5 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((olr.trend_mmw(boo)+olr.trend_mmw_err(boo))./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'midlat ERA5 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(merra2_stemptrend(boo),ones(size(boo)));                       fprintf(1,'midlat MERRA2 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                               summary_olr_mmw_stats.merra2.midlat_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend(boo),ones(size(boo)));                         fprintf(1,'midlat MERRA2 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                          summary_olr_mmw_stats.merra2.midlat_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(merra2_mmwtrend(boo)./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'midlat MERRA2 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                        summary_olr_mmw_stats.merra2.midlat_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend(boo) + merra2_mmwtrend_err(boo),ones(size(boo)));                    fprintf(1,'midlat MERRA2 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((merra2_mmwtrend(boo) + merra2_mmwtrend_err(boo))./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'midlat MERRA2 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);

%%%%%
disp(' ')
boo = find(abs(p.rlat) > 60);
[m,s,m0,s0] = weighted_mean_stddev(umbc.stemprate(boo),ones(size(boo)));               fprintf(1,'polar UMBC SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                 summary_olr_mmw_stats.umbc.polar_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate(boo),ones(size(boo)));                 fprintf(1,'polar UMBC mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                            summary_olr_mmw_stats.umbc.polar_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(umbc_mmwrate(boo)./umbc_mmw0(boo)),ones(size(boo)));  fprintf(1,'polar UMBC mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                          summary_olr_mmw_stats.umbc.polar_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate(boo)+umbc_mmwrate_err(boo),ones(size(boo)));                     fprintf(1,'polar UMBC mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((umbc_mmwrate(boo)+umbc_mmwrate_err(boo))./umbc_mmw0(boo)),ones(size(boo)));  fprintf(1,'polar UMBC mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_stemp(boo),ones(size(boo)));                       fprintf(1,'polar ERA5 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                  summary_olr_mmw_stats.era5.polar_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw(boo),ones(size(boo)));                         fprintf(1,'polar ERA5 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                             summary_olr_mmw_stats.era5.polar_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(olr.trend_mmw(boo)./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'polar ERA5 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                           summary_olr_mmw_stats.era5.polar_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw(boo)+olr.trend_mmw_err(boo),ones(size(boo)));                            fprintf(1,'polar ERA5 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((olr.trend_mmw(boo)+olr.trend_mmw_err(boo))./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'polar ERA5 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(merra2_stemptrend(boo),ones(size(boo)));                       fprintf(1,'polar MERRA2 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                               summary_olr_mmw_stats.merra2.polar_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend(boo),ones(size(boo)));                         fprintf(1,'polar MERRA2 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                          summary_olr_mmw_stats.merra2.polar_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(merra2_mmwtrend(boo)./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'polar MERRA2 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                        summary_olr_mmw_stats.merra2.polar_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend(boo) + merra2_mmwtrend_err(boo),ones(size(boo)));                    fprintf(1,'polar MERRA2 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((merra2_mmwtrend(boo) + merra2_mmwtrend_err(boo))./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'polar MERRA2 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);

%%%%%
disp(' ')
boo = find(p.landfrac == 0);
[m,s,m0,s0] = weighted_mean_stddev(umbc.stemprate(boo),ones(size(boo)));               fprintf(1,'ocean UMBC SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                 summary_olr_mmw_stats.umbc.ocean_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate(boo),ones(size(boo)));                 fprintf(1,'ocean UMBC mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                            summary_olr_mmw_stats.umbc.ocean_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(umbc_mmwrate(boo)./umbc_mmw0(boo)),ones(size(boo)));  fprintf(1,'ocean UMBC mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                          summary_olr_mmw_stats.umbc.ocean_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate(boo)+umbc_mmwrate_err(boo),ones(size(boo)));                     fprintf(1,'ocean UMBC mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((umbc_mmwrate(boo)+umbc_mmwrate_err(boo))./umbc_mmw0(boo)),ones(size(boo)));  fprintf(1,'ocean UMBC mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_stemp(boo),ones(size(boo)));                       fprintf(1,'ocean ERA5 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                  summary_olr_mmw_stats.era5.ocean_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw(boo),ones(size(boo)));                         fprintf(1,'ocean ERA5 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                             summary_olr_mmw_stats.era5.ocean_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(olr.trend_mmw(boo)./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'ocean ERA5 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                           summary_olr_mmw_stats.era5.ocean_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw(boo)+olr.trend_mmw_err(boo),ones(size(boo)));                            fprintf(1,'ocean ERA5 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((olr.trend_mmw(boo)+olr.trend_mmw_err(boo))./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'ocean ERA5 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(merra2_stemptrend(boo),ones(size(boo)));                       fprintf(1,'ocean MERRA2 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                               summary_olr_mmw_stats.merra2.ocean_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend(boo),ones(size(boo)));                         fprintf(1,'ocean MERRA2 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                          summary_olr_mmw_stats.merra2.ocean_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(merra2_mmwtrend(boo)./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'ocean MERRA2 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                        summary_olr_mmw_stats.merra2.ocean_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend(boo) + merra2_mmwtrend_err(boo),ones(size(boo)));                    fprintf(1,'ocean MERRA2 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((merra2_mmwtrend(boo) + merra2_mmwtrend_err(boo))./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'ocean MERRA2 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);

%%%%%
disp(' ')
boo = find(p.landfrac >= 0.99);
[m,s,m0,s0] = weighted_mean_stddev(umbc.stemprate(boo),ones(size(boo)));               fprintf(1,'land UMBC SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                 summary_olr_mmw_stats.umbc.land_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate(boo),ones(size(boo)));                 fprintf(1,'land UMBC mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                            summary_olr_mmw_stats.umbc.land_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(umbc_mmwrate(boo)./umbc_mmw0(boo)),ones(size(boo)));  fprintf(1,'land UMBC mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                          summary_olr_mmw_stats.umbc.land_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(umbc_mmwrate(boo)+umbc_mmwrate_err(boo),ones(size(boo)));                     fprintf(1,'land UMBC mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((umbc_mmwrate(boo)+umbc_mmwrate_err(boo))./umbc_mmw0(boo)),ones(size(boo)));  fprintf(1,'land UMBC mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_stemp(boo),ones(size(boo)));                       fprintf(1,'land ERA5 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                                  summary_olr_mmw_stats.era5.land_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw(boo),ones(size(boo)));                         fprintf(1,'land ERA5 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                             summary_olr_mmw_stats.era5.land_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(olr.trend_mmw(boo)./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'land ERA5 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                           summary_olr_mmw_stats.era5.land_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(olr.trend_mmw(boo)+olr.trend_mmw_err(boo),ones(size(boo)));                            fprintf(1,'land ERA5 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((olr.trend_mmw(boo)+olr.trend_mmw_err(boo))./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'land ERA5 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
disp(' ')
[m,s,m0,s0] = weighted_mean_stddev(merra2_stemptrend(boo),ones(size(boo)));                       fprintf(1,'land MERRA2 SKT absolute trends = %8.4f +/- %8.4f Kyr \n',m,s);                               summary_olr_mmw_stats.merra2.land_stemp = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend(boo),ones(size(boo)));                         fprintf(1,'land MERRA2 mmw absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);                          summary_olr_mmw_stats.merra2.land_mmw   = [m s];
[m,s,m0,s0] = weighted_mean_stddev(100*(merra2_mmwtrend(boo)./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'land MERRA2 mmw percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);                        summary_olr_mmw_stats.merra2.land_mmwpc = [m s];
[m,s,m0,s0] = weighted_mean_stddev(merra2_mmwtrend(boo) + merra2_mmwtrend_err(boo),ones(size(boo)));                    fprintf(1,'land MERRA2 mmw + err absolute trends = %8.4f +/- %8.4f mmH2O/yr \n',m,s);
[m,s,m0,s0] = weighted_mean_stddev(100*((merra2_mmwtrend(boo) + merra2_mmwtrend_err(boo))./nanmean(olr.xmmw(:,boo),1)),ones(size(boo)));  fprintf(1,'land MERRA2 mmw + err percent  trends = %8.4f +/- %8.4f percent/yr \n',m,s);
