function ceres = load_ceres_data(ceres_fname)

addpath /home/sergio/MATLABCODE
%% see /asl/s1/sergio/CERES_OLR_15year/Readme

%https://ceres.larc.nasa.gov/order_data.php
%
%Get monthly data set from 2002/09 to XXXX/08 : regional, monthly
%  CERES_EBAF-TOA_Ed4.1 Order

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = read_netcdf_lls(ceres_fname)
index = 1 : 228;

ceres.lon = a.lon;
ceres.lat = a.lat;

  data = a.toa_lw_all_mon(:,:,index);
  data = squeeze(mean(data,1));
  ceres.lwdata = data;

  data = a.toa_sw_all_mon(:,:,index);
  data = squeeze(mean(data,1));
  ceres.swdata = data;

  data = a.toa_net_all_mon(:,:,index);
  data = squeeze(mean(data,1));
  ceres.netdata = data;

  data = a.toa_lw_clr_c_mon(:,:,index);
  data = squeeze(mean(data,1));
  ceres.lwdata_clr = data;

  data = a.toa_sw_clr_c_mon(:,:,index);
  data = squeeze(mean(data,1));
  ceres.swdata_clr = data;

  data = a.toa_net_clr_c_mon(:,:,index);
  data = squeeze(mean(data,1));
  ceres.netdata_clr = data;

%  data = a.solar_mon(:,:,index);
%  data = squeeze(mean(data,1));
%  ceres.solardata = data;

%  data = a.cldarea_total_daynight_mon(:,:,index);
%  data = squeeze(mean(data,1));
%  ceres.cldareadata = data;

%  data = a.cldpress_total_daynight_mon(:,:,index);
%  data = squeeze(mean(data,1));
%  ceres.cldpressdata = data;

%  data = a.cldtemp_total_daynight_mon(:,:,index);
%  data = squeeze(mean(data,1));
%  ceres.cldtempdata = data;

%  data = a.cldtau_total_day_mon(:,:,index);
%  data = squeeze(mean(data,1));
%  ceres.cldtaudata = data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

EBAF 4.1

Parameters : 
I clicked on ALL params for
TOA fluxes
TOA CRE fluxes
Solar Flux
Cloud Params
Surface Fluxes
Surface CRE fluxes

Temmporal Res : monthly mean

Spatial : 1x1 grid

Time range : 09/2002 - 08/2018
%}

