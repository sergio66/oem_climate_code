addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/TIME
addpath /asl/matlib/h4tools

%% see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/load_ceres_data.m

dir0 = '/asl/s1/sergio/MODIS_MONTHLY_L3/';
modis_fname = 'AQUA_MODIS.20230701_20230731.L3m.MO.SST.sst.4km.nc';
yymmdd = modis_fname(12:17);
modis = read_netcdf_lls([dir0 modis_fname]);

[h,ha,p,pa] = rtpread('../FIND_NWP_MODEL_TRENDS/summary_atm_N_cld_20years_all_lat_all_lon_2002_2022_monthlyERA5.ip.rtp');

hindex = 1 : 8;
hindex = 1;

F.s_longitude = ncread([dir0 '/' modis_fname],'lon');
F.s_latitude  = ncread([dir0 '/' modis_fname],'lat');

% [X,Y] = ndgrid(F.s_latitude,F.s_longitude);
% iX = flipud(X); iY = flipud(Y);
% figure(1); clf; simplemap(X,simplemap(wrapTo180(Y),squeeze(a.toa_lw_all_mon(:,:,1))))

[Y,X] = ndgrid(F.s_latitude,F.s_longitude);
Y = Y'; X = X';
iX = X; iY = Y;
figure(1); clf; simplemap(Y,X,a.sst);

moo = double(ncread([dir0 '/' modis_fname],'sst'));
F.junk.ig   = griddedInterpolant(iX,iY,moo,'linear');
miaow = F.junk.ig(p.rlon,p.rlat);;
figure(2); clf; simplemap(p.rlat,p.rlon,miaow,5); colorbar; shading flat; caxis([200 300]); colormap jet;
figure(2); clf; scatter_coast(p.rlon,p.rlat,100,miaow); colorbar; shading flat; caxis([200 300]); colormap jet;
figure(1); cx = caxis; figure(2); caxis(cx)

boo = dir([dir0 '/AQUA_MODIS*.nc']);
for hindex = 1 : length(boo)
  modis_fname = boo(hindex).name;
  yymmdd = modis_fname(12:17);
  yy = str2num(yymmdd(1:2));
  mm = str2num(yymmdd(3:4));
  dd = str2num(yymmdd(5:6));
  moo = double(ncread([dir0 '/' modis_fname],'sst'));
  F.junk.ig   = griddedInterpolant(iX,iY,moo,'linear');
  miaow = F.junk.ig(wrapTo360(p.rlon),p.rlat);;
  modis.sst_all_4608(hindex,:) = miaow;
  modis.daysSince2002(hindex) = change2days(yy,mm,dd,2002);
end


