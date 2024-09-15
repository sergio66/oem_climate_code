addpath /asl/matlib/h4tools
addpath /asl/matlib/plotutils
addpath /asl/matlib/maps
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/NANROUTINES
addpath /asl/matlib/plotutils

clear all
for ii = 1 : 12
  figure(ii); clf
end

load llsmap5
cmap = usa2;
cmap = llsmap5;

do_XX_YY_from_X_Y

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');

%%%%%%%%%%%%%%%%%%%%%%%%%

%% section 5, Test 1 of trends paper_2024_07_26.tex
do_usualavg_VS_8timestepavg
disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

%% section 5, Test 2 of trends paper_2024_07_26.tex
do_Q90_VS_mean_VS_usual

disp('ret to continue'); pause
error('this is Section 5 of trends paper, Sept 2024')

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

asc_grid      = load('ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc.mat');
asc_surf_grid = load('ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_asc_surf.mat');

desc_grid      = load('ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc.mat');
desc_surf_grid = load('ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc_surf.mat');
desc_surf_rand = load('ERA5_atm_N_cld_data_2002_09_to_2022_08_trends_desc_randompt_surf.mat');

disp('Comparing Descending : usual grid center VS random point VS 273 center VS 273 Q90')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_surf_grid.trend_stemp,'desc usual');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_surf_rand.trend_stemp,'desc random');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_surf_grid.trend_stemp - desc_surf_rand.trend_stemp,       'desc : usual - rand');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_surf_grid.trend_stemp - desc_era5_273files.skt_trend_Q90, 'desc : usual - 273 center');
disp(' ')
stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_surf_grid.trend_stemp - desc_era5_273files.skt_trend_cntr,'desc : usual - 273 Q90');
disp(' ')

stats = print_simple_cosweighted_stats(p.rlat,p.landfrac,desc_surf_grid.trend_stemp - desc_era5_273files.skt_trend_mean,'desc : usual - 273 mean');
disp('see : global thediffs are typically 0 +/- 0.02 K for all .. put this into trends paper!!!!')
disp('  ')

oo = find(p.landfrac == 0 & abs(p.rlat) < 60);
junk = desc_surf_grid.trend_stemp(oo) - desc_era5_273files.skt_trend_mean(oo);
stats = print_simple_cosweighted_stats(p.rlat(oo),p.landfrac(oo),junk,'desc ocean trp/midlat : usual - 273 mean');
disp('see : midlat/ocean thediffs are typically 0 +/- 0.02 K for all .. put this into trends paper!!!!')

%%%%%%%%%%%%%%%%%%%%%%%%%

avg_asc_desc = 0.5*(asc_surf_grid.trend_stemp + desc_surf_grid.trend_stemp);
figure(1); clf; scatter_coast(p.rlon,p.rlat,50,avg_asc_desc);             title('U : Usual Grid Center (D+N)/2'); caxis([-1 +1]*0.150); colormap(cmap)
figure(2); clf; scatter_coast(p.rlon,p.rlat,50,dn_surf_grid.trend_stemp); title('A : 273 avg Grid Center DN');    caxis([-1 +1]*0.150); colormap(cmap)
figure(3); clf; scatter_coast(p.rlon,p.rlat,50,dn_surf_rand.trend_stemp); title('R : Random Grid Center  DN');    caxis([-1 +1]*0.150); colormap(cmap)
figure(4); clf; scatter_coast(p.rlon,p.rlat,50,dn_surf_grid.trend_stemp-dn_surf_rand.trend_stemp); title('A-R');  caxis([-1 +1]*0.025); colormap(cmap)
figure(5); clf; scatter_coast(p.rlon,p.rlat,50,dn_surf_rand.trend_stemp-avg_asc_desc);             title('A-U');  caxis([-1 +1]*0.025); colormap(cmap)

disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                              >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
disp('ret to continue !!!! look at these results !!!!!! '); pause
disp('ret to continue !!!! look at these results again !!!!!! '); pause
disp('ret to continue !!!! look at these results again again !!!!!! '); pause
disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                              >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

%%%%%%%%%%%%%%%%%%%%%%%%%

disp('fit look only at ERA5 SKT trends, descending')
disp('  take tile center desc_surf_grid.trend_stemp as truth')

aslmap(1,rlat65,rlon73,smoothn(reshape(desc_era5_273files.skt_trend_cntr,72,64)',1), [-90 +90],[-180 +180]);                            title('1.30 am Desc, tile mid pt of 273 points');
aslmap(2,rlat65,rlon73,smoothn(reshape(desc_era5_273files.skt_trend_Q90,72,64)',1), [-90 +90],[-180 +180]);                             title('1.30 am Desc, Q90  of 273 points');
aslmap(3,rlat65,rlon73,smoothn(reshape(desc_era5_273files.skt_trend_mean,72,64)',1), [-90 +90],[-180 +180]);                            title('1.30 am Desc, mean  of 273 points');
aslmap(4,rlat65,rlon73,smoothn(reshape(desc_surf_grid.trend_stemp,72,64)',1), [-90 +90],[-180 +180]);                                   title('1.30 am Desc, grid');
aslmap(5,rlat65,rlon73,smoothn(reshape(desc_surf_rand.trend_stemp,72,64)',1), [-90 +90],[-180 +180]);                                   title('1.30 am Desc, rand');
aslmap(6,rlat65,rlon73,smoothn(reshape(desc_surf_grid.trend_stemp-desc_surf_rand.trend_stemp,72,64)',1), [-90 +90],[-180 +180]);        title('1.30 am Desc, grid-rand');
aslmap(7,rlat65,rlon73,smoothn(reshape(desc_surf_grid.trend_stemp-desc_era5_273files.skt_trend_cntr,72,64)',1), [-90 +90],[-180 +180]); title('1.30 am Desc, grid-273mid');
aslmap(8,rlat65,rlon73,smoothn(reshape(desc_surf_grid.trend_stemp-desc_era5_273files.skt_trend_Q90,72,64)',1), [-90 +90],[-180 +180]);  title('1.30 am Desc, grid-273 Q90');
for ii = 1 : 5
  figure(ii); caxis([-1 +1]*0.15); colormap(cmap)
end
ii = 6;   figure(ii); caxis([-1 +1]*0.025); colormap(cmap)
ii = 7;   figure(ii); caxis([-1 +1]*0.025); colormap(cmap)
ii = 8;   figure(ii); caxis([-1 +1]*0.025); colormap(cmap)
figure(9); wah = desc_surf_grid.trend_stemp-desc_surf_rand.trend_stemp; dst = (-1:0.05:+1)*0.025; 
  oo = find(p.landfrac == 0 & abs(p.rlat) < 30);
  oo = find(p.landfrac == 0 & abs(p.rlat) < 60);
  fprintf(1,'mean ocean dSKT/dt in this abs(p.rat) < 60  range is %8.6f +/- %8.6f K/yr \n',mean(desc_surf_grid.trend_stemp(oo)),std(desc_surf_grid.trend_stemp(oo)))
  gah = smooth(histc(wah(oo),dst),5); gah = gah/length(oo); gah = gah/max(gah); plot(dst,gah); grid;
figure(10); plot(desc_surf_grid.trend_stemp,desc_surf_grid.trend_stemp,desc_surf_grid.trend_stemp,desc_surf_rand.trend_stemp,'r.')
figure(11); plot(desc_surf_grid.trend_stemp(oo),desc_surf_grid.trend_stemp(oo),desc_surf_grid.trend_stemp(oo),desc_surf_rand.trend_stemp(oo),'r.')
[r,chisqr,P,sigP,numpts] = nanlinearcorrelation(desc_surf_grid.trend_stemp(oo),desc_surf_rand.trend_stemp(oo));

oo = 1 : 4608;
oo = find(p.landfrac == 0);
oo = find(p.landfrac == 0 & abs(p.rlat) < 60);

x0 = desc_surf_grid.trend_stemp(oo);
x1 = desc_surf_grid.trend_stemp(oo)-desc_surf_rand.trend_stemp(oo);
x2 = desc_surf_grid.trend_stemp(oo)-desc_era5_273files.skt_trend_cntr(oo);
printarray([mean(x1) std(x1) mean(x2) std(x2)],'mean/std for usualpaper-paperrandompt and usualpaper-273pts center')

disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                              >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aslmap(1,rlat65,rlon73,smoothn(reshape(avg_asc_desc,72,64)',1), [-90 +90],[-180 +180]);               title('1.30 am/pm (Asc+Desc)/2');
aslmap(2,rlat65,rlon73,smoothn(reshape(dn_surf_grid.trend_stemp,72,64)',1), [-90 +90],[-180 +180]);   title('(24 hour avg, grid center');
aslmap(3,rlat65,rlon73,smoothn(reshape(dn_surf_rand.trend_stemp,72,64)',1), [-90 +90],[-180 +180]);   title('(24 hour avg, random pt');
for ii = 1 : 3
  figure(ii); caxis([-1 +1]*0.15); colormap(cmap)
end

aslmap(4,rlat65,rlon73,smoothn(reshape(dn_surf_grid.trend_stemp-dn_surf_rand.trend_stemp,72,64)',1), [-90 +90],[-180 +180]); title('24 hour avg : grid center - random pt');
aslmap(5,rlat65,rlon73,smoothn(reshape(dn_surf_rand.trend_stemp-avg_asc_desc,72,64)',1), [-90 +90],[-180 +180]);             title('24 hour avg grid center - \newline 1.30 am/pm (Asc+Desc)/2')
for ii = 4 : 5
  figure(ii); caxis([-1 +1]*0.025); colormap(cmap)
end

st_thresh = 1e-4;
st_thresh = 1e-5;

wah = (dn_surf_grid.trend_stemp-dn_surf_rand.trend_stemp)./dn_surf_grid.trend_stemp * 100;
  bah = find(abs(dn_surf_grid.trend_stemp) < st_thresh); wah(bah) = NaN;
  aslmap(6,rlat65,rlon73,smoothn(reshape(wah,72,64)',1), [-90 +90],[-180 +180]); title('24 hour avg : grid center - random pt : percent diff');
wah = (dn_surf_rand.trend_stemp-avg_asc_desc)./dn_surf_rand.trend_stemp * 100;
  bah = find(abs(dn_surf_grid.trend_stemp) < st_thresh); wah(bah) = NaN;
  aslmap(7,rlat65,rlon73,smoothn(reshape(wah,72,64)',1), [-90 +90],[-180 +180]); title('24 hour avg grid center - \newline 1.30 am/pm (Asc+Desc)/2 : percent diff')
for ii = 6 : 7
  figure(ii); caxis([-1 +1]*25); colormap(cmap)
end

dn = -100 : 5 : +100;
wah = (dn_surf_grid.trend_stemp-dn_surf_rand.trend_stemp)./dn_surf_grid.trend_stemp * 100;
  bah = find(abs(dn_surf_grid.trend_stemp) < st_thresh); wah(bah) = NaN;
  figure(6); clf; plot(dn,histc(wah,dn)); title('24 hour avg : grid center - random pt : percent diff'); grid on
wah = (dn_surf_rand.trend_stemp-avg_asc_desc)./dn_surf_rand.trend_stemp * 100;
  bah = find(abs(dn_surf_grid.trend_stemp) < st_thresh); wah(bah) = NaN;
  figure(7); clf; plot(dn,histc(wah,dn)); title('24 hour avg grid center - \newline 1.30 am/pm (Asc+Desc)/2 : percent diff'); grid on

moo = find(p.landfrac < eps & abs(p.rlat) <= 60);
moo = find(p.landfrac < eps & abs(p.rlat) <= 30);
dn = -100 : 5 : +100;
wah = (dn_surf_grid.trend_stemp-dn_surf_rand.trend_stemp)./dn_surf_grid.trend_stemp * 100;
  bah = find(abs(dn_surf_grid.trend_stemp) < st_thresh); wah(bah) = NaN;
  figure(6); clf; plot(dn,histc(wah(moo),dn)); title('24 hour avg : grid center - random pt : percent diff'); grid on
wah = (dn_surf_rand.trend_stemp-avg_asc_desc)./dn_surf_rand.trend_stemp * 100;
  bah = find(abs(dn_surf_grid.trend_stemp) < st_thresh); wah(bah) = NaN;
  figure(7); clf; plot(dn,histc(wah(moo),dn)); title('24 hour avg grid center - \newline 1.30 am/pm (Asc+Desc)/2 : percent diff'); grid on
