disp('make sure you have run       module load netCDF-Fortran/4.4.4-intel-2018b      before starting Matlab')
disp('make sure you have run       module load netCDF-Fortran/4.4.4-intel-2018b      before starting Matlab')
disp('make sure you have run       module load netCDF-Fortran/4.4.4-intel-2018b      before starting Matlab')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RobinHoganECMWF/ECRAD_ECMWF_version_of_flux/ecRad/create_ecrad_inputSergio/
addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/DRIVER_CODE_RRTM_Band17/
addpath /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE/CRODGERS_FAST_CLOUD/
addpath /home/sergio/MATLABCODE/PLOTTER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('doing allfov')
[hh,hha,pp,ppa] = rtpread('/asl/s1/sergio/rtp/rtp_airicrad_v6/2019/04/25/cloudy_airs_l1c_ecm_sarta_baum_ice.2019.04.25.222_cumsum_-1.rtp');
ix = 1;
ix = 1 : 1000 : 12150;
ix = 1 : 10 : 12150;

iXStep = 25;  %% THIS WORKS
iXStep = 256; %% see nblocksize in /home/sergio/IR_NIR_VIS_UV_RTcodes/RobinHoganECMWF/ECRAD_ECMWF_version_of_flux/ecRad/create_ecrad_inputSergio/config_sergio.nam

if length(ix) <= 100
  [~,pjunk] = subset_rtp_allcloudfields(hh,pp,[],[],ix);
  junkall = driver_run_ecRad_rtp_loop_over_profiles(hh,pjunk,+1);
else
  ixlen = ceil(length(ix)/iXStep);
  for ixx = 1 : ixlen
    ind = (1:iXStep)  + (ixx-1)*iXStep;
    jind = find(ind <= length(ix));
    ind = ind(jind);
    [~,pjunk] = subset_rtp_allcloudfields(hh,pp,[],[],ix(ind));
    xjunkall = driver_run_ecRad_rtp_loop_over_profiles(hh,pjunk,+1);
    junkall.cld(ind) = xjunkall.cld;
    junkall.clr(ind) = xjunkall.clr;
  end
end
scatter_coast(pp.rlon(ix),pp.rlat(ix),50,junkall.clr-junkall.cld); title('Clear-CLoud OLR');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(' ')
disp(' ')
disp(' ')
disp('doing 4608 tiles')
load test_ecRad
px = p;
%umbc_spectral_olr.olr0 = compute_olr(h,px);
[~,pjunk] = subset_rtp_allcloudfields(h,px,[],[],1);
  junk = driver_run_ecRad_rtp_loop_over_profiles(h,pjunk,-1);
olr0_ecRad = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);               %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<
scatter_coast(px.rlon,px.rlat,50,olr0_ecRad.clr)
%umbc_spectral_olr.olr0_rrtm  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,px,0); %%% <<<<<<<<<<<<<<<<<<<<<<<<<<<
figure(1); scatter_coast(px.rlon,px.rlat,50,olr0_ecRad.clr); title('ecRad')

rrtm = load('olr0_rrtm.mat');
figure(2); scatter_coast(px.rlon,px.rlat,50,rrtm.olr0_rrtm); title('rrtm')
figure(3); hist(olr0_ecRad.clr./rrtm.olr0_rrtm,100); title('ecrad/rrtm')
figure(3); scatter_coast(px.rlon,px.rlat,50,olr0_ecRad.clr./rrtm.olr0_rrtm); title('ecrad/rrtm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KCARTA CHECK of RRTM vs ECRAD

bobo = find(rrtm.olr0_rrtm == max(rrtm.olr0_rrtm),1)
[~,pjunk] = subset_rtp_allcloudfields(h,px,[],[],bobo);
junkecRad = driver_run_ecRad_rtp_loop_over_profiles(h,pjunk,-1);
junkRRTM  = driver_rrtm_no_xsec_nocloud_twoslab_band17only_loop(h,pjunk,0);
%{
rtpwrite('testRRTM_ECRAD.rtp',h,[],pjunk,[]);
GO TO /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES, 
  edit   set_rtp to use this, 
  edit   set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm.m to have iFlux = 5    
  run    cluster_loop_allkcartabands
Move output to oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/TESTOLR_KCARTA
  load   allkc5bands_prof1.mat

trapz(wall,dall)/1000*pi = 366.81                        %% quick radiances --> OLR
trapz(wall,fluxall)/1000 = 329.0361  355.0080  353.1940  %%  ILR, trop(OLR),OLR

junkecRad =
    cld: 297.0296
    clr: 297.0296
junkRRTM = 
  359.6853

%}
