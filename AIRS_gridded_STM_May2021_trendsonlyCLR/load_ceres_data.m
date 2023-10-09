function ceres = load_ceres_data(ceres_fname,iCorT)

addpath /home/sergio/MATLABCODE

%% see /asl/s1/sergio/CERES_OLR_15year/Readme

%https://ceres.larc.nasa.gov/order_data.php
%
%Get monthly data set from 2002/09 to XXXX/08 : regional, monthly
%  CERES_EBAF-TOA_Ed4.1 Order

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ryan suggested this
 
%% Interesting. I agree it often looks closer to the all-sky CERES but I
%% can’t think of a reason why. CERES does have biases of ~ 5W/m2 but
%% that’s for the raw values, not anomalies, where the errors are
%% believed to be much much smaller (and in line with ocean heat content
%% measurements) owing to the good stability of the instrument.  I wonder
%% if this is just a coincidence that your calculations match up with the
%% all-sky. Another possibility is it has to do with the clear-sky CERES
%% fluxes you are using. As of CERES EBAF version 4.1 they have a
%% traditional clear-sky product that uses clear-sky scenes, but they
%% also have a new clear-filled product where the total region is
%% clear-sky. I’m not sure how they do it, but it’s more lkike how AIRS
%% OLR and climate models define clear-sky, where there are no clouds at
%% all.  They must be doing some sort of correction using a radiative
%% transfer model or something.
%% 
%% Anyways, it may be worth comparing with that clear-filled version if
%% you aren’t already.  I’m not sure if it will improve things, but worth
%% a shot I think. You can download it from the CERES EBAF product rather
%% than the CERES EBAF TOA product on the CERES website. Here is a direct
%% link. Just click the little down arrow next to “TOA Fluxes” to see the
%% two clear-sky options.
%% https://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAF41Selection.jsp

%% https://ceres-tool.larc.nasa.gov/ord-tool/jsp/EBAF41Selection.jsp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
  error('need one argin')
  ceres_fname = '/asl/s1/sergio/CERES_OLR_15year/CERES_EBAF-TOA_Ed4.1_Subset_200209-202108.nc'; 
    Y0 = 2002; M0 = 09; YE = 2021; ME = 09; iNumY = 19; index = 1 : 228;
  ceres_fname = '/asl/s1/sergio/CERES_OLR_15year/CERES_EBAF-TOA_Ed4.2_Subset_200003-202307.nc'; 
    Y0 = 2000; M0 = 03; YE = 2023; ME = 07; iNumY = 20; index = 31 : 31-1+240;
  iCorT = 1;  %% just get down the clear filled region in a pixel .. rather than em[irically filed "Total" clear sky  
elseif nargin == 1
  iCorT = 1;  %% just get down the clear filled region in a pixel .. rather than empirically filed "Total" clear sky
end

a = read_netcdf_lls(ceres_fname)

[mmm] = length(a.time);

if ~isfield(a,'toa_lw_clr_t_mon')
  disp('toa_lw_clr_t_mon DNW so not getting the empirically filled total clear sky')
  iCorT = 1;
end

index = 1 : mmm;

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

  if iCorT == 1
    disp('getting _c_ == partially clear observed pixels') 
    data = a.toa_lw_clr_c_mon(:,:,index);
    data = squeeze(mean(data,1));
    ceres.lwdata_clr = data;
  
    data = a.toa_sw_clr_c_mon(:,:,index);
    data = squeeze(mean(data,1));
    ceres.swdata_clr = data;
  
    data = a.toa_net_clr_c_mon(:,:,index);
    data = squeeze(mean(data,1));
    ceres.netdata_clr = data;
   else
    disp('getting _t_ == empirically/astonishingly totally clear observed pixels') 
    data = a.toa_lw_clr_t_mon(:,:,index);
    data = squeeze(mean(data,1));
    ceres.lwdata_clr = data;
  
    data = a.toa_sw_clr_t_mon(:,:,index);
    data = squeeze(mean(data,1));
    ceres.swdata_clr = data;
  
    data = a.toa_net_clr_t_mon(:,:,index);
    data = squeeze(mean(data,1));
    ceres.netdata_clr = data;
   end
  
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

