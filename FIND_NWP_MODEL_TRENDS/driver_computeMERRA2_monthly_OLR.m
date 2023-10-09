%% monthly, 18 years x 12 months/year = 216
%% monthly, 19 years x 12 months/year = 228
%% monthly, 12 years x 12 months/year = 144

addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/StrowCodeforTrendsAndAnomalies/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/Strow_humidity/convert_humidity/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/IDL_WV_ROUTINES/atmos_phys/MATLAB/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/create_ecrad_inputSergio/
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD

load('llsmap5.mat');

%% note this code only handles 2002_09 to YYYY_08 sp
%%      this code does not currently handle eg OCO2 MERRA2_atm_N_cld_data_2012_05_to_2019_04_trends_desc.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('assumes you have run driver_computeMERRA2_monthly_trends_desc_or_asc.m');

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 240 for months

iNumYears = 20;
for JOB = 1 : iNumYears*12

  fnameOUT = ['OLR_ecRad/MERRA2/merra2_olr_' num2str(JOB,'%03d') '.mat'];
  if ~exist(fnameOUT)    
    fnameIN = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/MERRA2/Tile_Center/DESC/merra2_tile_center_monthly_' num2str(JOB,'%03d') '.mat'];
    loader = ['load ' fnameIN];
    eval(loader);
    clear fnameIN
    
    h = hnew_op;
    px = pnew_op;
    olr = superdriver_run_ecRad_rtp_loop_over_profiles(h,px,-1);
    
    stemp = px.stemp;
    saver = ['save OLR_ecRad/MERRA2/merra2_olr_' num2str(JOB,'%03d') '.mat stemp olr'];
    eval(saver)
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computeMERRA2_OLR_trend
