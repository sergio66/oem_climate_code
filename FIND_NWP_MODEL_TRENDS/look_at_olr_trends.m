olr = load('OLR_ecRad/ERA5/all_era5_olr.mat');

junk = load('MERRA2_atm_data_2002_09_to_2022_08_trends_desc.mat','trend_mmw');     merra2_mmwtrend = junk.trend_mmw;
junk = load('MERRA2_atm_data_2002_09_to_2022_08_trends_desc.mat','trend_mmw_err'); merra2_mmwtrend_err = junk.trend_mmw_err;
junk = load('MERRA2_atm_data_2002_09_to_2022_08_trends_desc.mat','trend_stemp');   merra2_stemptrend = junk.trend_stemp;

junk = load('/asl/s1/sergio/AIRS_CLIMCAPS/airsclimcaps_64x72_rates_fastgrib_stats_Sept2002_Aug2022_20yr_desc.mat');
junk.thestats64x72.stemprate(isnan(junk.thestats64x72.stemprate)) = 0;
junk.thestats64x72.stempratestd(isnan(junk.thestats64x72.stempratestd)) = 0;
junk.thestats64x72.ptemprate(isnan(junk.thestats64x72.ptemprate)) = 0;
junk.thestats64x72.ptempratestd(isnan(junk.thestats64x72.ptempratestd)) = 0;
junk.thestats64x72.waterrate(isnan(junk.thestats64x72.waterrate)) = 0;
junk.thestats64x72.waterratestd(isnan(junk.thestats64x72.waterratestd)) = 0;
mmw0_junk = mmwater_rtp(h,p);

addpath /asl/matlib/science/
[salti,landfrac] =  usgs_deg10_dem(p.rlat,p.rlon);
p.landfrac = landfrac;

pjunk = p; 
  pjunk.stemp = pjunk.stemp + reshape(junk.thestats64x72.stemprate,1,4608);
  pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + (reshape(junk.thestats64x72.ptemprate,4608,100)');
  pjunk.gas_1((1:066)+34,:) = pjunk.gas_1((1:066)+34,:) .* (1 + (reshape(junk.thestats64x72.waterrate,4608,066)'));
mmwX = mmwater_rtp(h,pjunk);  climcaps_mmwtrend = mmwX - mmw0_junk;
%plot(rlat,nanmean(reshape(climcaps_mmwtrend,72,64),1))
pjunk = p; 
  pjunk.stemp = pjunk.stemp + reshape(junk.thestats64x72.stemprate,1,4608) + reshape(junk.thestats64x72.stempratestd,1,4608);
  pjunk.ptemp(1:100,:) = pjunk.ptemp(1:100,:) + (reshape(junk.thestats64x72.ptemprate,4608,100)' + reshape(junk.thestats64x72.ptempratestd,4608,100)');
  pjunk.gas_1((1:066)+34,:) = pjunk.gas_1((1:066)+34,:) .* (1 + (reshape(junk.thestats64x72.waterrate,4608,066)') + (reshape(junk.thestats64x72.waterratestd,4608,066)'));
mmwX = mmwater_rtp(h,pjunk);  climcaps_mmwtrend_err = mmwX - mmw0_junk;
clear mmw0_junk mmX

airsL3_mmwtrend = airsL3.trend_mmw;

umbc_mmw0 = mmwater_rtp(h,p);
pUMBC = p;
pUMBC.stemp = pUMBC.stemp + umbc.stemprate;
pUMBC.ptemp(1:100,:) = pUMBC.ptemp(1:100,:) + umbc.ptemprate;
pUMBC.gas_1(1:100,:) = pUMBC.gas_1(1:100,:) .*( 1 + umbc.waterrate);
umbc_mmwrate = mmwater_rtp(h,pUMBC) - umbc_mmw0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp('in OLR stuff')
  clear junk
  junkx.resultsunc = zeros(4608,6);
  junkx.resultsunc(:,6) = 0.2 * umbc.stemprate;
  try
    junk = load(strUMBC,'resultsunc');
  catch
    disp('no resultsunc in the saved retrievals, use constant value of 20 percent');
    junk.resultsunc(:,6) = junkx.resultsunc(:,6);
  end
  if size(junk) == [1 1]
    junk.resultsunc(:,6) = junkx.resultsunc(:,6);
  end

xpUMBC = p;
xpUMBC.stemp = xpUMBC.stemp + umbc.stemprate + junk.resultsunc(:,6)';
junk = load(strUMBC,'deltaTunc');
  xpUMBC.ptemp(1:100,:) = xpUMBC.ptemp(1:100,:) + umbc.ptemprate + junk.deltaTunc(1:100,:);
junk = load(strUMBC,'fracWVunc');
  xpUMBC.gas_1(1:100,:) = xpUMBC.gas_1(1:100,:) .*( 1 + umbc.waterrate(1:100,:) + junk.fracWVunc(1:100,:));
umbc_mmwrate_err = mmwater_rtp(h,xpUMBC) - umbc_mmw0;
umbc_mmwrate_err = abs(umbc_mmwrate_err - umbc_mmwrate);

disp('Recall density water = 1000 kg/m3 ==> 1 kg of water speread over area 1 m2 uses 1 mm ==> 1 mmH2O 1 kg/m2 of water')
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

%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFig = 50;
iFig = iFig + 1; figure(iFig); clf; aslmap(51,rlat65,rlon73,smoothn(reshape(nanmean(olr.xmmw,1),72,64)',1),  [-90 +90],[-180 +180]); colormap(jet); title('ERA5 mmw');
iFig = iFig + 1; figure(iFig); clf; aslmap(52,rlat65,rlon73,smoothn(reshape(nanmean(olr.xstemp,1),72,64)',1),[-90 +90],[-180 +180]); colormap(jet); title('ERA5 stemp');
iFig = iFig + 1; figure(iFig); clf; aslmap(53,rlat65,rlon73,smoothn(reshape(nanmean(olr.xolr,1),72,64)',1),  [-90 +90],[-180 +180]); colormap(jet); title('ERA5 olr clr');

iFig = iFig + 1; figure(iFig); clf; aslmap(54,rlat65,rlon73,smoothn(reshape(olr.trend_mmw,72,64)',1),  [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('ERA5 trend mmw');
iFig = iFig + 1; figure(iFig); clf; aslmap(55,rlat65,rlon73,smoothn(reshape(olr.trend_stemp,72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('ERA5 trend stemp');

iFig = iFig + 1; figure(iFig); clf; plot(1:4608,abs(umbc_mmwrate),'b',1:4608,umbc_mmwrate_err,'c',1:4608,abs(olr.trend_mmw),'r',1:4608,olr.trend_mmw_err,'m')

iFig = iFig + 1; figure(iFig); clf; aslmap(56,rlat65,rlon73,smoothn(reshape(umbc_mmwrate,72,64)',1),  [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('UMBC trend mmw');
iFig = iFig + 1; figure(iFig); clf; aslmap(57,rlat65,rlon73,smoothn(reshape(umbc.stemprate,72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('UMBC trend stemp');

iFig = iFig + 1; figure(iFig); clf; aslmap(58,rlat65,rlon73,smoothn(reshape(merra2_mmwtrend,72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('MERRA2 trend mmw');
iFig = iFig + 1; figure(iFig); clf; aslmap(59,rlat65,rlon73,smoothn(reshape(merra2_stemptrend,72,64)',1),[-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('MERRA2 trend stemp');

iFig = iFig + 1; figure(iFig); clf; aslmap(60,rlat65,rlon73,smoothn(reshape(olr.trend_olr,72,64)',1),   [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('ERA5 trend olr clr');

iFig = iFig + 1; figure(iFig); clf; 
plot(trend_rlat64,nanmean(reshape(olr.trend_mmw,72,64),1),'b',trend_rlat64,nanmean(reshape(merra2_mmwtrend,72,64),1),'g',trend_rlat64,nanmean(reshape(umbc_mmwrate,72,64),1),'r',...
     trend_rlat64,nanmean(reshape(airsL3_mmwtrend,72,64),1),'c',trend_rlat64,nanmean(reshape(climcaps_mmwtrend,72,64),1),'m','linewidth',2);
hold on
boo1 = nanmean(reshape(olr.trend_mmw,72,64),1); 
  boo2 = reshape(olr.trend_mmw_err,72,64) .* reshape(olr.trend_mmw_err,72,64); boo2 = sqrt(sum(boo2,1)); boo2 = boo2/72;
  errorbar(trend_rlat64,nanmean(reshape(olr.trend_mmw,72,64),1),boo2,'b','linewidth',2);
boo1 = nanmean(reshape(merra2_mmwtrend,72,64),1);
  boo2 = reshape(merra2_mmwtrend_err,72,64) .* reshape(merra2_mmwtrend_err,72,64); boo2 = sqrt(sum(boo2,1)); boo2 = boo2/72;
  errorbar(trend_rlat64,nanmean(reshape(merra2_mmwtrend,72,64),1),boo2,'g','linewidth',2);
boo1 = nanmean(reshape(umbc_mmwrate,72,64),1);
  boo2 = reshape(umbc_mmwrate_err,72,64) .* reshape(umbc_mmwrate_err,72,64); boo2 = sqrt(nansum(boo2,1)); boo2 = boo2/72;
  errorbar(trend_rlat64,nanmean(reshape(umbc_mmwrate,72,64),1),boo2,'r','linewidth',2);
boo1 = nanmean(reshape(airsL3_mmwtrend,72,64),1);
  %boo2 = reshape(airsL3_mmwtrend_err,72,64) .* reshape(airsL3_mmwtrend_err,72,64); boo2 = sqrt(nansum(boo2,1)); boo2 = boo2/72;
  boo2 = reshape(climcaps_mmwtrend_err,72,64) .* reshape(climcaps_mmwtrend_err,72,64); boo2 = sqrt(nansum(boo2,1)); boo2 = boo2/72;
  errorbar(trend_rlat64,nanmean(reshape(airsL3_mmwtrend,72,64),1),boo2,'c','linewidth',2);
boo1 = nanmean(reshape(climcaps_mmwtrend,72,64),1);
  boo2 = reshape(climcaps_mmwtrend_err,72,64) .* reshape(climcaps_mmwtrend_err,72,64); boo2 = sqrt(nansum(boo2,1)); boo2 = boo2/72;
  errorbar(trend_rlat64,nanmean(reshape(climcaps_mmwtrend,72,64),1),boo2,'m','linewidth',2);
hold off
plotaxis2; hl = legend('ERA5','MERRA2','UMBC','AIRS L3','CLIMCAPS L3','location','best','fontsize',10); 
xlabel('Latitude [deg]'); ylabel('d(mmw)/dt [mm/yr]'); xlim([-1 +1]*90) 

iFig = iFig + 1; figure(iFig); clf; plot(trend_rlat64,nanmean(reshape(olr.trend_stemp,72,64),1),'b',trend_rlat64,nanmean(reshape(merra2_stemptrend,72,64),1),'g',trend_rlat64,nanmean(reshape(umbc.stemprate,72,64),1),'r','linewidth',2);
  plotaxis2; hl = legend('ERA5','MERRA2','UMBC','location','best'); title('stemp trends')
iFig = iFig + 1; figure(iFig); clf; plot(nanmean(reshape(olr.trend_stemp,72,64),1),nanmean(reshape(olr.trend_mmw,72,64),1),'bs',nanmean(reshape(merra2_stemptrend,72,64),1),nanmean(reshape(merra2_mmwtrend,72,64),1),'gd',...
                      nanmean(reshape(umbc.stemprate,72,64),1),nanmean(reshape(umbc_mmwrate,72,64),1),'rx','markersize',5,'linewidth',2)
  plotaxis2; hl = legend('ERA5','MERRA2','UMBC','location','best'); xlabel('stemp trends'); ylabel('mmw trends')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
do_eofs

iFig = iFig + 1; figure(iFig); clf
  ymax = 4;
  imagesc(thecorrEOF.BT1231_ST); colorbar; caxis([-1 +1]); ylabel('BT1231 EOF'); ylim([0.5 0.5+ymax]); 
  title('dSKT/dt Correlations with ERA5 SKT EOF')
  set(gca,'ytick',[1:ymax],'yticklabel',num2str((1:ymax)'),'fontsize',10);
  set(gca,'xtick',[1:6],'xticklabel',{'ERA5','MERRA2','AIRS L3','CLIMCAPSL3','GISS','THIS WORK'},'fontsize',8);
  colormap(llsmap5)

iFig = iFig + 1; figure(iFig); clf
  ymax = 4;
  imagesc(thecorrEOF.era5skt_ST); colorbar; caxis([-1 +1]); ylabel('ERA5 SKT EOF'); ylim([0.5 0.5+ymax]); 
  title('dSKT/dt Correlations with ERA5 SKT EOF')
  set(gca,'ytick',[1:ymax],'yticklabel',num2str((1:ymax)'),'fontsize',10);
  set(gca,'xtick',[1:6],'xticklabel',{'ERA5','MERRA2','AIRS L3','CLIMCAPSL3','GISS','THIS WORK'},'fontsize',8);
  colormap(llsmap5)

iFig = iFig + 1; figure(iFig); clf
  ymax = 4;
  imagesc(thecorrEOF.BT1519_MMW); colorbar; caxis([-1 +1]); ylabel('BT 1519 EOF'); ylim([0.5 0.5+ymax]); 
  title('dMMW/dt Correlations with ERA5 SKT EOF')
  set(gca,'ytick',[1:ymax],'yticklabel',num2str((1:ymax)'),'fontsize',10);
  set(gca,'xtick',[1:6],'xticklabel',{'ERA5','MERRA2','AIRS L3','CLIMCAPSL3','GISS','THIS WORK'},'fontsize',8);
  colormap(llsmap5)

iFig = iFig + 1; figure(iFig); clf
  ymax = 4;
  imagesc(thecorrEOF.era5mmw_MMW); colorbar; caxis([-1 +1]); ylabel('ERA5 MMW EOF'); ylim([0.5 0.5+ymax]); 
  title('dMMW/dt Correlations with ERA5 MW EOF')
  set(gca,'ytick',[1:ymax],'yticklabel',num2str((1:ymax)'),'fontsize',10);
  set(gca,'xtick',[1:6],'xticklabel',{'ERA5','MERRA2','AIRS L3','CLIMCAPSL3','GISS','THIS WORK'},'fontsize',8);
  colormap(llsmap5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

umbcX = load(strUMBC,'rates');
fairs = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/h2645structure.mat','h');
fairs = fairs.h;
i900 = find(h.vchan >= 900,1);
i1227 = find(h.vchan >= 1227,1);

  %% see get_spectral_get_the_model_trends.m
  iEorMstrSPECTRA = 'ERA5';
  iEorMFstr = 'simulate64binsERA5_';
  iEorMFstr = 'reconstruct_era5_spectra_geo_rlat';
  for ii = 1 : 64
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iEorMstrSPECTRA '/' iEorMFstr num2str(ii,'%02i') '.mat'];
    fname = ['/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/' iEorMstrSPECTRA '/' iEorMFstr num2str(ii,'%02i') '_2002_09_2022_08.mat'];
    junk = load(fname,'fchanx');  fchanx   = junk.fchanx;
    junk = load(fname,'thesave'); thesave = junk.thesave.xtrendSpectral;
    ind = (ii-1)*72 + (1:72);
    era5.era5_spectral_rates(:,ind) = thesave;
  end

iFig = iFig + 1; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(umbcX.rates(1520,:),72,64)',1),   [-90 +90],[-180 +180]); caxis([-1 +1]*0.15); colormap(llsmap5); title('BT1231 obs');
iFig = iFig + 1; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(umbcX.rates(i900,:)-umbcX.rates(1520,:),72,64)',1),   [-90 +90],[-180 +180]); caxis([-1 +1]*0.01); colormap(llsmap5); title('BT900-BT1231 obs');
iFig = iFig + 1; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(umbcX.rates(1520,:)-umbcX.rates(i1227,:),72,64)',1),   [-90 +90],[-180 +180]); caxis([-1 +1]*0.015); colormap(llsmap5); title('BT1231-BT1227 obs');
iFig = iFig + 1; figure(iFig); clf; aslmap(iFig,rlat65,rlon73,smoothn(reshape(era5.era5_spectral_rates(1520,:)-era5.era5_spectral_rates(i1227,:),72,64)',1),   [-90 +90],[-180 +180]); caxis([-1 +1]*0.015); colormap(llsmap5); title('BT1231-BT1227 obs');
corr_mmw(1) = nanlinearcorrelation(olr.trend_mmw,olr.trend_stemp);
corr_mmw(2) = nanlinearcorrelation(olr.trend_mmw,umbc_mmwrate);
corr_mmw(3) = nanlinearcorrelation(umbc.stemprate,umbc_mmwrate);
corr_mmw(4) = nanlinearcorrelation(umbcX.rates(1520,:)-umbcX.rates(i1227,:),umbc_mmwrate);
corr_mmw(5) = nanlinearcorrelation(era5.era5_spectral_rates(1520,:)-era5.era5_spectral_rates(i1227,:),olr.trend_mmw);
corr_mmw(6) = nanlinearcorrelation(era5.era5_spectral_rates(1520,:)-era5.era5_spectral_rates(i1227,:),umbcX.rates(1520,:)-umbcX.rates(i1227,:));
fprintf(1,'correlation between : ERA5 mmw trend & ERA5 skt trend = %8.4f \n',corr_mmw(1))
fprintf(1,'                    : ERA5 mmw trend & UMBC mmw trend = %8.4f \n',corr_mmw(2))
fprintf(1,'                    : UMBC mmw trend & UMBC skt trend = %8.4f \n',corr_mmw(3))
fprintf(1,'                    : UMBC mmw trend & UMBC (AIRS OBS) BT1231-BT1227 = %8.4f \n',corr_mmw(4))
fprintf(1,'                    : ERA5 mmw trend & ERA5 (AIRS OBS) BT1231-BT1227 = %8.4f \n',corr_mmw(5))
fprintf(1,'                    : ERA5     trend & UMBC trend      BT1231-BT1227 = %8.4f \n',corr_mmw(6))
  summary_olr_mmw_stats.corr_mmw = corr_mmw;

disp(' ')
disp('global comparisons of ERA5 with UMBC : ')
disp(' ')
disp('global comparisons of ERA5 with UMBC : ')
disp('RH/WVfrac/T :        200 mb              |         500 mb             |     800 mb');
disp('-------------------------------------------------------------------------------------------')
fprintf(1,'CORRELATION : %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f  \n',squeeze(allXchi(:,1,4))');             summary_olr_mmw_stats.umbc.allXchi  = squeeze(allXchi(:,1,4))';
fprintf(1,'BIAS        : %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f  \n',squeeze(allXmean(:,1,4))');            summary_olr_mmw_stats.umbc.allXmean = squeeze(allXmean(:,1,4))';
fprintf(1,'STD         : %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f  \n',squeeze(allXstd(:,1,4))');             summary_olr_mmw_stats.umbc.allXstd  = squeeze(allXstd(:,1,4))';
fprintf(1,'FRAC +0-    : %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f  \n',squeeze(allX_frac_neg0pos(:,1,4))');   summary_olr_mmw_stats.umbc.allX_frac_neg0pos = squeeze(allX_frac_neg0pos(:,1,4))';
disp('-------------------------------------------------------------------------------------------')

disp(' ')
if length(comment) == 0
  fprintf(1,'<<<<<< these results used >>>>>>> %s \n',strUMBC)
  fprintf(1,'<<<<<< these results used >>>>>>> %s \n',strUMBC)
  fprintf(1,'<<<<<< these results used >>>>>>> %s \n',strUMBC)
else
  fprintf(1,'<<<<<< these results used >>>>>>> %s \n with comment %s \n',strUMBC,comment)
  fprintf(1,'<<<<<< these results used >>>>>>> %s \n with comment %s \n',strUMBC,comment)
  fprintf(1,'<<<<<< these results used >>>>>>> %s \n with comment %s \n',strUMBC,comment)
end

disp(' ')
