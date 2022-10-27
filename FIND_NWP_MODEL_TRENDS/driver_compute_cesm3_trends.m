%% copied from /driver_compute_AIRSL3_trends_desc_or_asc.m

addpath /asl/matlib/aslutil/
addpath /asl/matlib/science
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/FIND_TRENDS/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

clear all

for ii = 1 : 6
  figure(ii); colormap jet
end

timespan = 12;
timespan = 07;
timespan = 18;
timespan = 19;
fprintf(1,'timespan = %2i years \n',timespan)

if timespan == 07
  savestr_version = 'May2012_07yr';
  StartY = 2012; StartYM = 5;   %% start 05/2012
  StopY  = 2019; StopYM  = 4;   %% stop  04/2019  
elseif timespan == 16
  savestr_version = 'Sept2017_15yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2018; StopYM  = 8;   %% stop  08/2017  
elseif timespan == 18
  savestr_version = 'Sept2002_Aug2020_18yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2020; StopYM  = 8;   %% stop  08/2020
elseif timespan == 12
  savestr_version = 'Sept2002_Aug2014_12yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2014; StopYM  = 8;   %% stop  08/2014
elseif timespan == 19
  savestr_version = 'Sept2002_Jul2021_19yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2021; StopYM  = 7;   %% stop  08/2021  
  savestr_version = 'Sept2002_Aug2021_19yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2021; StopYM  = 8;   %% stop  08/2021  
else
  error('huh check timespan')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
maindir = '/asl/models/cesm3/';
cesm_files = findfiles([maindir '/*reduced_cam*.nc']);

cesm_Date_Start = datenum(1998,01+(1:length(cesm_files))-1,01);   %%%% we start reading data potentially from 2002
%%%%%%%%%%%%%%%%%%%%%%%%%

cesm_Date = cesm_Date_Start;

for i = 1:length(cesm_files)  
  fin = cesm_files{i};
  yy(i) = str2num(fin(length(fin)-9:length(fin)-6));
  mm(i) = str2num(fin(length(fin)-4:length(fin)-3));

  read_this_file(i) = -1;
  if yy(i) == StartY & mm(i) >= StartYM
    read_this_file(i) = +1;
  elseif yy(i) == StopY & mm(i) <= StopYM
    read_this_file(i) = +1;
  elseif yy(i) > StartY & yy(i) < StopY
    read_this_file(i) = +1;
  end

  fprintf(1,'%3i : %4i / %2i Y/N = %2i \n',i,yy(i),mm(i),read_this_file(i))
end

iOK = find(read_this_file > 0);
if i >= abs(timespan)*12 & length(iOK) == abs(timespan)*12
  fprintf(1,'%4i/%2i to %4i/%2i %s \n',[StartY StartYM StopY StopYM],savestr_version)  
  disp('ok, found the files needed for the YEAR timespan you chose')
else
  iStop = input('WARNING : oops need to get in more CESM L3 data, +1 to stop???????? ')
  if iStop == 1
    error('quitting')
  else
    disp('hmm, proceed at your OWN RISK!')
  end
end

disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

%% this is suppose we save 2002/09 to 2021/08 but only want eg AMIP/CMIP time = 2002/09 to 2014/08
iSkipTo_64x72_trend = input('Skip directly to trends by reading in earlier files????? (-1 default /+1) : ');
if length(iSkipTo_64x72_trend) == 0
  iSkipTo_64x72_trend = -1;
end
if iSkipTo_64x72_trend == +1

  load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
  drlon = 5; 
  rlon = -180 : drlon : +180;      %% 73 edges, so 72 bins
  rlat = latB2;                    %% 65 edges, so 64 bins
  %save_lat64x72 = 0.5*(rlat(1:end-1)+rlat(2:end));
  %save_lon64x72 = 0.5*(rlon(1:end-1)+rlon(2:end));

  %iDorA = input('Enter Asc(-1) or Desc (+1) : ');
  iDorA = 1;
  if length(iDorA) == 0
    iDorA = +1;
  end
  iL3orCLIMCAPS = +1;

  %savestr_version_big = 'Sept2002_Aug2021_19yr';
  %loader = ['load /asl/s1/sergio/CESM3/cesm3_' savestr_version_big '.mat'];

  savestr_version_big = 'Sept2002_Aug2021_19yr';
  loader = ['load /asl/s1/sergio/CESM3/cesm_64x72_rates_' savestr_version_big '.mat'];

  eval(loader)
  do_the_fits_cesmL3_rates_tiles  %% which is part of do_the_CESM3_trends.m, but assumes you have successfully run "plot_cesm3_data_tiles" to completion
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

woo = find(read_this_file > 0);  %% typically 1 -- 12*timespan

%% from standard
cesm_plev  = zeros(length(woo),58,192,288,'single');
cesm_ptemp  = zeros(length(woo),58,192,288,'single');
cesm_gas_1   = zeros(length(woo),58,192,288,'single');  %% really MMR in g/kg dry air
cesm_gas_3 = zeros(length(woo),58,192,288,'single');
cesm_rh    = zeros(length(woo),58,192,288,'single');

cesm_CO2   = zeros(length(woo),58,192,288,'single');  %% vmr
cesm_CH4   = zeros(length(woo),58,192,288,'single');  %% vmr

cesm_stemp = zeros(length(woo),192,288,'single');
cesm_spres = zeros(length(woo),192,288,'single');
cesm_wspd  = zeros(length(woo),192,288,'single');

cesm_ciwwc = zeros(length(woo),58,192,288,'single');
cesm_clwc = zeros(length(woo),58,192,288,'single');
cesm_cc   = zeros(length(woo),58,192,288,'single');
cesm_tcc  = zeros(length(woo),192,288,'single');

clear yy mm

iDo = input('(+1) read in saved mat file or (-1) read in CESM3 files, one at a time  ? ');
if iDo < 0
  clear yy mm
  for ix = 1:length(woo)
    i = woo(ix);
    fprintf(1,'reading in %3i of %3i required files -- file %3i in overall list of length %3i \n',ix,length(woo),i,length(read_this_file))
    
    fin = cesm_files{i};
    disp(fin)
    yy(i) = str2num(fin(length(fin)-9:length(fin)-6));
    mm(i) = str2num(fin(length(fin)-4:length(fin)-3));

    a = read_netcdf_lls(fin);

    if isfield(a,'lat')
      cesm_Lat = a.lat;
      cesm_Lon = a.lon;
    end

    cesm_clwc(ix,:,:,:)  = permute(a.CLDLIQ,[3 2 1]);
    cesm_ciwcC(ix,:,:,:)  = permute(a.CLDICE,[3 2 1]);
    cesm_cc(ix,:,:,:)    = permute(a.CLOUD,[3 2 1]);
    cesm_tcc(ix,:,:)     = permute(a.CLDTOT,[2 1]);

    cesm_plev(ix,:,:,:)  = permute(a.PMID,[3 2 1]);
    cesm_ptemp(ix,:,:,:) = permute(a.T,[3 2 1]);
    cesm_gas_1(ix,:,:,:) = permute(a.Q,[3 2 1]);
    cesm_gas_3(ix,:,:,:) = permute(a.O3,[3 2 1]);
    cesm_rh(ix,:,:,:)    = permute(a.RELHUM,[3 2 1]);

    cesm_spres(ix,:,:)   = permute(a.PS,[2 1]);
    cesm_stemp(ix,:,:)   = permute(a.TS,[2 1]);
    cesm_wspd(ix,:,:)    = permute(a.U10,[2 1]);

    cesm_CO2(ix,:,:)   = a.co2vmr*1e6;
    cesm_CH4(ix,:,:)   = a.ch4vmr*1e6;

    %%%%%%%%%%%%%%%%%%%%%%%%%
  end

  cesm_landfrac = a.LANDFRAC;

  iSave = input('save these mega huge files (+1) or store in memory??? (-1) : ? ');
  if iSave > 0
    saver = ['save -v7.3 /asl/s1/sergio/CESM3/cesm3_' savestr_version '.mat  cesm* yy mm'];
    eval(saver);
  end
  
else

  disp('loading in data A')
  %% should be same as /cesmv7_March2014.mat

  loader = ['load /asl/s1/sergio/CESM3/cesm3_' savestr_version '.mat'];
  eval(loader);

end

error('ooo')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_the_CESM3_trends %%  basically same as do_the_AIRSL3_trends
