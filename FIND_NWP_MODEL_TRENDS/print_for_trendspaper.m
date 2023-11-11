addpath /asl/matlib/plotutils
dir0 = '//home/sergio/PAPERS/SUBMITPAPERS/trends/Figs/';

figure(37); aslprint([dir0 '/tz_trends_zonal_p_5panels20_years.pdf']);
figure(35); aslprint([dir0 '/fracWV_trends_zonal_p_5panels20_years.pdf']);
figure(33); aslprint([dir0 '/rh_trends_zonal_p_5panels20_years.pdf']);

figure(38); aslprint([dir0 '/tz_correl_slope_bias_stddev_20_years.pdf']);
figure(36); aslprint([dir0 '/fracWV_correl_slope_bias_stddev_20_years.pdf']);
figure(34); aslprint([dir0 '/rh_correl_slope_bias_stddev_20_years.pdf']);

figure(21); aslprint([dir0 '/200mb_rhtrends_lat_lon_5panels20_years.pdf']);
figure(24); aslprint([dir0 '/500mb_rhtrends_lat_lon_5panels20_years.pdf']);
figure(27); aslprint([dir0 '/800mb_rhtrends_lat_lon_5panels20_years.pdf']);

figure(22); aslprint([dir0 '/200mb_fracWVtrends_lat_lon_5panels20_years.pdf']);
figure(25); aslprint([dir0 '/500mb_fracWVtrends_lat_lon_5panels20_years.pdf']);
figure(28); aslprint([dir0 '/800mb_fracWVtrends_lat_lon_5panels20_years.pdf']);

figure(23); aslprint([dir0 '/200mb_tztrends_lat_lon_5panels20_years.pdf']);
figure(26); aslprint([dir0 '/500mb_tztrends_lat_lon_5panels20_years.pdf']);
figure(29); aslprint([dir0 '/800mb_tztrends_lat_lon_5panels20_years.pdf']);


