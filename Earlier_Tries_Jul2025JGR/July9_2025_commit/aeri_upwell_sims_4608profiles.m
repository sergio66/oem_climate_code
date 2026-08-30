%% see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/aeri_upwell_sims_rtp.m
%% see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/aeri_upwell_sims_4608profiles.m

addpath /asl/matlib/aslutil
addpath /asl/matlib/rtptools
addpath /asl/matlib/science/
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD
addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RobinHoganECMWF/ECRAD_ECMWF_version_of_flux/ecRad/create_ecrad_inputSergio/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% make the rtp files -- the raw + jacobian, and the one with perturbations
%%% make the rtp files -- the raw + jacobian, and the one with perturbations
%%% make the rtp files -- the raw + jacobian, and the one with perturbations

iMonthSoFar = 120;
loader = ['load /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/era5_tile_center_monthly_' num2str(iMonthSoFar) '.mat'];
if ~exist('pnew_op')
  eval(loader);
end

h5 = hnew_op;
p5 = pnew_op;

p5.upwell = 2 * ones(size(p5.stemp));

[p5.salti,p5.landfrac] = usgs_deg10_dem(p5.rlat,wrapTo180(p5.rlon));

[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2ppmv(h5,p5,1:length(p5.stemp),2);
fprintf(1,'mean CO2 ppmvAVG = %8.4f\n',mean(ppmvAVG));
ppmvAVG0 = ppmvAVG;
rCO2PPM = input('enter CO2 ppm (-1 for default 400 ppm): ');
if length(rCO2PPM) == 0
  rCO2PPM = -1;
end
if rCO2PPM > 0
  p5.gas_2 = p5.gas_2 * rCO2PPM/mean(ppmvAVG0);
  [ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2ppmv(h5,p5,1:length(p5.stemp),2);
end

%%% make the rtp files -- the raw + jacobian, and the one with perturbations
%%% make the rtp files -- the raw + jacobian, and the one with perturbations
%%% make the rtp files -- the raw + jacobian, and the one with perturbations

p5pertWV   = p5;   p5pertWV.gas_1   = p5pertWV.gas_1 * (1 + 0.001);
p5pertCO2  = p5;   p5pertCO2.gas_2  = p5pertCO2.gas_2 * (1 + 2.2/385);
p5pertT    = p5;   p5pertT.ptemp    = p5pertT.ptemp + 0.02;
p5pertALL  = p5;   p5pertALL.gas_1 = p5pertALL.gas_1 * (1 + 0.001);
                   p5pertALL.gas_2 = p5pertALL.gas_2 * (1 + 2.2/385);
                   p5pertALL.ptemp = p5pertALL.ptemp + 0.02;

p5pertX = p5pertWV;
[h5,p5pertX] = cat_rtp(h5,p5pertX,h5,p5pertCO2);
[h5,p5pertX] = cat_rtp(h5,p5pertX,h5,p5pertT);
[h5,p5pertX] = cat_rtp(h5,p5pertX,h5,p5pertALL);

p5pertT_TS = p5pertT;     p5pertT_TS.stemp = p5pertT_TS.stemp + 0.02;
p5pertALL_TS = p5pertALL; p5pertALL_TS.stemp = p5pertALL_TS.stemp + 0.02;
p5pertX_TS = p5pertWV;
[h5,p5pertX_TS] = cat_rtp(h5,p5pertX_TS,h5,p5pertCO2);
[h5,p5pertX_TS] = cat_rtp(h5,p5pertX_TS,h5,p5pertT_TS);
[h5,p5pertX_TS] = cat_rtp(h5,p5pertX_TS,h5,p5pertALL_TS);

%rtpwrite('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin4608.op_400ppm_uplook.rtp',h5,ha,p5,pa);
%rtpwrite('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin4608.op_400ppm_uplook_pert.rtp',h5,ha,p5pertX,pa);

%% for finite diff jacs, dQ = 0.001, dT = 0.01;
dQ = 0.001;
dT = 0.01;

%% see /umbc/xfs2/strow/asl/s1/sergio/home/IR_NIR_VIS_UV_RTcodes/RobinHoganECMWF/ECRAD_ECMWF_version_of_flux/ecRad/create_ecrad_inputSergio/wrapper_run_ecRad_rtp_loop_over_profiles.m
if ~exist('ecRad')
  ecRad0     = superdriver_run_ecRad_rtp_loop_over_profiles(h5,p5,-1,-1);
  %ecRad      = superdriver_run_ecRad_rtp_loop_over_profiles(h5,p5pertX,-1,-1);
  ecRad      = superdriver_run_ecRad_rtp_loop_over_profiles(h5,p5pertX_TS,-1,-1);
end

lat = p5.rlat;
ecRadWV     = ecRad.clr((1:4608)+0*4608)-ecRad0.clr;
ecRadCO2    = ecRad.clr((1:4608)+1*4608)-ecRad0.clr;
ecRadT      = ecRad.clr((1:4608)+2*4608)-ecRad0.clr;
ecRadtotal  = ecRad.clr((1:4608)+3*4608)-ecRad0.clr;
ecRadtotal2 = ecRadWV + ecRadCO2 + ecRadT;

figure(6);
plot(lat,ecRadWV,'b',lat,ecRadCO2,'g',lat,ecRadT,'r',lat,ecRadtotal,'kx-','linewidth',2)
ylim([0 0.2])
title('Actual ecRad Flux Changes!!!')
plotaxis2; hl = legend('WV','CO2','T','total','location','best');
ylabel('Flux change W/m2'); xlabel('Latitude'); xlim([-max(abs(lat)) +max(abs(lat))]);

figure(7)
plot(lat,ecRadWV./ecRadtotal,'b',lat,ecRadCO2./ecRadtotal,'g',lat,ecRadT./ecRadtotal,'r','linewidth',2)
plotaxis2;
hl = legend('WV','CO2','T','location','best','fontsize',10);; xlim([-max(abs(lat)) +max(abs(lat))]);
ylabel('Fractional Contribution \newline to Flux change'); title('for realistic annual change : ecRad'); xlabel('Latitude');

comment = 'see aeri_upwell_sims_4608profiles.m';
if rCO2PPM == -1
  saver = ['save ilr_4608_ecRad.mat ecRad* comment']; 
else
  saver = ['save ilr_4608_ecRad_co2ppm_' num2str(rCO2PPM) '.mat ecRad* comment']; 
end
% eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iPrint = -1;
if iPrint > 0
  print_results_co2ppm_aeri_brutsaert
end
