addpath /asl/matlib/plotutils
addpath /asl/matlib/rtptools
addpath /asl/matlib/maps
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/PLOTTER/TILEDPLOTS
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /asl/matlib/science/
addpath /home/sergio/MATLABCODE/SHOWSTATS
addpath /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/
addpath /home/sergio/MATLABCODE/NANROUTINES
addpath /asl/matlib/aslutil

load llsmap5

ceres_olr = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/ceres_trends_20year_T.mat');
ceres_ilr = load('../AIRS_gridded_STM_May2021_trendsonlyCLR/ceres_trends_20year_T_ilr.mat');

[h,ha,p,pa] = rtpread('summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.op.rtp');
[salti,landfrac] =  usgs_deg10_dem(p.rlat,p.rlon);
p.landfrac = landfrac;

mmw0 = mmwater_rtp(h,p);
RH0  = layeramt2RH(h,p);

for ii = 1 : length(p.stemp)
  nlay = p.nlevs(ii) - 1;
  RHSurf(ii) = RH0(nlay,ii);
end

% umbc_night_file  = '/asl/s1/sergio/JUNK/smallgather_tileCLRnight_SEQN_dataset09_Q03_newERA5_2021jacs_startwith0_50fatlayers_CarbonTrackerCO2_test2D.mat';
% junk = load(umbc_night_file,'results');
% junkSKTtrend = junk.results(:,6);
% 
% dRH_Ts = jeevanjee_PBL_deltaRH_Ts(p.stemp,RHSurf/100,precipitation_vs_skt_changes(p.rlat,p.landfrac),junkSKTtrend,p);
% figure(1); clf; scatter_coast(p.rlon,p.rlat,100,dRH_Ts); title('Jevanjee dRH per K'); colormap(llsmap5)

ocean = find(p.landfrac == 0);
land  = find(p.landfrac == 1);

do_XX_YY_from_X_Y
mu = cos(YY'*pi/180);
load llsmap5

tropics = find(abs(YY) <= 30);          
midlatsNtropics = find(abs(YY) <= 60); 
midlats = setdiff(midlatsNtropics,tropics); 
poles = find(abs(YY) > 60);             

tropics = nan(size(YY));         junk = find(abs(YY) <= 30);                tropics(junk) = 1;
midlatsNtropics = nan(size(YY)); junk = find(abs(YY) <= 60);                midlatsNtropics(junk) = 1;
midlats = nan(size(YY));         junk = find(abs(YY) > 30 & abs(YY) <= 60); midlats(junk) = 1;
poles = nan(size(YY));           junk = find(abs(YY) > 60);                 poles(junk) = 1;

plays100 = load('/home/sergio/MATLABCODE/airslevels.dat');
plays100 = flipud(plevs2plays(plays100));
plays100 = plays100(2:end);

usa2x = usa2; usa2x = usa2x(61:120,:);
jett = jet; jett(1,:) = 1;

% https://uw.pressbooks.pub/fundamentalsofclimatechange/chapter/moving-energy/
%
%  The atmospheric energy transport can be decomposed into several
%  components. The flux of dry static energy, s = cp T + g z represents
%  the internal and gravitational potential energy that is transported,
%  with cp the specific heat of dry air at constant pressure, T the
%  temperature, g=9.8 m/s2 the gravitational acceleration, and z the
%  height. Latent energy also contributes to the poleward transport of
%  heat, and comes into the moist static energy m = cp T + g z + Lv q ,
%  with Lv = 2.5 x 106 J/kg the latent heat of vaporization and q the
%  specific humidity. Latent energy is an important component of the
%  poleward transport of energy. When water vapor is evaporated in the
%  subtropics, it requires energy. If the moisture is transported
%  polewards in a baroclinic eddy, before condensing out at higher
%  latitudes, it releases its latent energy where the condensation
%  occurs. Thus energy is transported poleward by the movement of water
%  vapor, just like when warmer air moves poleward and is
%  mixed. Additional latent energy is released when freezing occurs, but
%  this is smaller than the condensational heating. The transport of
%  kinetic energy is much smaller compared with the other components.

cp = 1005;       %% J/kg/K
g  = 9.81;       %% m/s2
Lv = 2.50084e6;  %% J/kg
%% so Q = cp T + g z + Lv q = [J/kg/K][K] + g z + [J/jg][kg/kg] = J/kg + g z     
%%   g z = phi = geopotential = work done to raise 1 kg from sea level to that point (m=1) === m g z = (1) g z = J for 1 kg = J/kg
Q0 = 1000 * layers2gg(h,p,1:4608,1);
[mm,nn] = size(Q0);
Q0 = Lv * Q0 + cp * p.ptemp(1:mm,:) + g * p.palts(1:mm,:);

clear f*
prep_colWV_T_WV_trends_Day_vs_Night

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('ret to continue to <<< simulated simulated simulated >>> SPECTRAL trends'); pause
compare_spectral_trends_Day_or_Night
disp('run "driver_compare_Day_vs_Night_spectral_trends.m" for the real spectral comparison')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('ret to continue to SKT trends'); pause
compare_SKT_trends_Day_vs_Night
compare_SKT_trends_Day_vs_Night_stats

disp('ret to continue to colWV trends'); pause
compare_colWV_trends_Day_vs_Night
compare_colWV_trends_Day_vs_Night_stats

disp('ret to continue to T trends'); pause
compare_T_trends_Day_vs_Night

disp('ret to continue to fracWV trends'); pause
compare_WV_trends_Day_vs_Night

disp('ret to continue to RH trends'); pause
compare_RH_trends_Day_vs_Night

disp('ret to continue to RHsurf trends'); pause
compare_RHsurf_trends_Day_vs_Night

disp('ret to continue to misc surface stuff such as MSE trends, ILR trends, compared to dSKT/dt'); pause
compare_misc_stuff_trends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
print_Day_vs_Night_trendspaper
%}
