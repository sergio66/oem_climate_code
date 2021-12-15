%% copied from /home/sergio/MATLABCODE/AIRS_L3/driver_airsL3_ratesv6_Sept2017_equalarea_latbins.m

addpath /asl/matlib/aslutil/
addpath /home/sergio/MATLABCODE/TIME
addpath /asl/matlib/science
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/FIND_TRENDS/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

clear all

for ii = 1 : 6
  figure(ii); colormap jet
end

timespan = 18;
timespan = 19;
fprintf(1,'timespan = %2i years \n',timespan)

if timespan == 16
  savestr_version = 'Sept2017_15yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2018; StopYM  = 8;   %% stop  08/2017  
elseif timespan == 18
  savestr_version = 'Sept2002_Aug2020_18yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2020; StopYM  = 8;   %% stop  08/2020
elseif timespan == 19
  savestr_version = 'Sept2002_Aug2021_19yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2021; StopYM  = 8;   %% stop  08/2021  
  savestr_version = 'Sept2002_Jul2021_19yr';
  StopY  = 2021; StopYM  = 7;   %% stop  08/2021  
else
  error('huh check timespan')
end

%%% TEMPERATURE

Airs_PT = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, ...
70, 50, 30, 20, 15, 10, 7, 5, 3, 2, 1.5,1];
Airs_PQ = [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100];

%%%%%%%%%%%%%%%%%%%%%%%%%
maindir = '/asl/data/airs/AIRX3STM/';
maindir = '/asl/data/airs/Aqua_AIRS_Level3/AIRS3STM.006/';
maindir = '/asl/s1/sergio/AIRS_L3/acdisc.gsfc.nasa.gov/ftp/data/s4pa/Aqua_AIRS_Level3/AIRX3STM.006/';
maindir = '/asl/airs/AIRS3STM/v7/'; %% https://acdisc.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_Level3/AIRS3STM.006/2018/
Airs_files = findfiles([maindir '/*/*.hdf']);

Airs_Date_Start = datenum(StartY,09+(1:length(Airs_files))-1,01);
Airs_Date_Start = datenum(2002,09+(1:length(Airs_files))-1,01);   %%%% we start reading data potentially from 2002
%Airs_Date_End   = datenum(StopY,08+(1:length(Airs_files))-1,01);

%maindirS = '/asl/s1/sergio/AIRS_L3/acdisc.gsfc.nasa.gov/ftp/data/s4pa/Aqua_AIRS_Level3/AIRX3SPM.006/';
%%%%%%%%%%%%%%%%%%%%%%%%%

Airs_Date = Airs_Date_Start;

for i = 1:length(Airs_files)  
  fin = Airs_files{i};
  oops = findstr(fin,'.L3.RetStd');
  yy(i) = str2num(fin(oops-10+0:oops-10+3));
  mm(i) = str2num(fin(oops-10+5:oops-10+6));

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
  disp('ok, found the files needed for the YEAR timespan you chose')
else
  iStop = input('WARNING : oops need to get in more AIRS L3 data, +1 to stop???????? ')
  if iStop == 1
    error('quitting')
  else
    disp('hmm, proceed at your OWN RISK!')
  end
end

disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

woo = find(read_this_file > 0);  %% typically 1 -- 12*timespan

%% from standard
Airs_Temp_A = zeros(length(woo),24,180,360,'single');
Airs_Temp_D = zeros(length(woo),24,180,360,'single');
Airs_STemp_A = zeros(length(woo),180,360,'single');
Airs_STemp_D = zeros(length(woo),180,360,'single');
Airs_H2OVap_A = zeros(length(woo),12,180,360,'single');  %% really MMR in g/kg dry air
Airs_H2OVap_D = zeros(length(woo),12,180,360,'single');  %% really MMR in g/kg dry air
Airs_RH_A     = zeros(length(woo),12,180,360,'single');  
Airs_RH_D     = zeros(length(woo),12,180,360,'single');  
Airs_RHSurf_A = zeros(length(woo),180,360,'single');  
Airs_RHSurf_D = zeros(length(woo),180,360,'single');  

Airs_Ozone_A = zeros(length(woo),24,180,360,'single');
Airs_Ozone_D = zeros(length(woo),24,180,360,'single');
Airs_CO_A = zeros(length(woo),24,180,360,'single');
Airs_CO_D = zeros(length(woo),24,180,360,'single');
Airs_CH4_A = zeros(length(woo),24,180,360,'single');
Airs_CH4_D = zeros(length(woo),24,180,360,'single');
Airs_SPres_A = zeros(length(woo),180,360,'single');
Airs_SPres_D = zeros(length(woo),180,360,'single');
Airs_OLR_A = zeros(length(woo),180,360,'single');
Airs_OLR_D = zeros(length(woo),180,360,'single');
Airs_ClrOLR_A = zeros(length(woo),180,360,'single');
Airs_ClrOLR_D = zeros(length(woo),180,360,'single');
Airs_CldPres_A = zeros(length(woo),3,180,360,'single');
Airs_CldPres_D = zeros(length(woo),3,180,360,'single');
Airs_CldFrac_A = zeros(length(woo),3,180,360,'single');
Airs_CldFrac_D = zeros(length(woo),3,180,360,'single');

%% from MW
Airs_LiqWater_A = zeros(length(woo),180,360,'single');
Airs_LiqWater_D = zeros(length(woo),180,360,'single');

%% from support
Airs_IceOD_A = zeros(length(woo),180,360,'single');  %% ice_cld_opt_dpth
Airs_IceOD_D = zeros(length(woo),180,360,'single');
Airs_IceSze_A = zeros(length(woo),180,360,'single');  %% ice_cld_eff_diam
Airs_IceSze_D = zeros(length(woo),180,360,'single');
Airs_IceT_A = zeros(length(woo),180,360,'single');  %% ice_cld_temp_eff
Airs_IceT_D = zeros(length(woo),180,360,'single');

clear yy mm

iDo = input('(+1) read in saved mat file or (-1) read in L3 files, one at a time  ? ');
if iDo < 0
  clear yy mm
  for ix = 1:length(woo)
    i = woo(ix);
    fprintf(1,'reading in %3i of %3i required files -- file %3i in overall list of length %3i \n',ix,length(woo),i,length(read_this_file))
    
    fin = Airs_files{i};
    disp(fin)
    oops = findstr(fin,'.L3.RetStd');
    yy(ix) = str2num(fin(oops-10+0:oops-10+3));
    mm(ix) = str2num(fin(oops-10+5:oops-10+6));

    Airs_Temp_A(ix,:,:,:) = hdfread(fin, 'ascending', 'Fields','Temperature_A');
    Airs_Temp_D(ix,:,:,:) = hdfread(fin, 'descending','Fields','Temperature_D');
    if ~ exist('Airs_Lat','var')
      Airs_Lat = hdfread(fin, 'location', 'Fields', 'Latitude');
      Airs_Lon = hdfread(fin, 'location', 'Fields', 'Longitude');
    end

    Airs_STemp_A(ix,:,:) = hdfread(fin, 'ascending', 'Fields','SurfSkinTemp_A');
    Airs_STemp_D(ix,:,:) = hdfread(fin, 'descending', 'Fields','SurfSkinTemp_D');

    %Airs_H2OVap_A(ix,:,:,:) = hdfread(fin, 'ascending', 'Fields','H2OVapMMR_A');
    %Airs_H2OVap_D(ix,:,:,:) = hdfread(fin, 'descending','Fields','H2OVapMMR_D');

    Airs_H2OVap_A(ix,:,:,:) = hdfread(fin, 'ascending', 'Fields','H2O_MMR_A');
    Airs_H2OVap_D(ix,:,:,:) = hdfread(fin, 'descending','Fields','H2O_MMR_D');

    Airs_RH_A(ix,:,:,:) = hdfread(fin, 'ascending', 'Fields','RelHum_A');
    Airs_RH_D(ix,:,:,:) = hdfread(fin, 'descending','Fields','RelHum_D');
    Airs_RHSurf_A(ix,:,:) = hdfread(fin, 'ascending', 'Fields','RelHumSurf_A');
    Airs_RHSurf_D(ix,:,:) = hdfread(fin, 'descending','Fields','RelHumSurf_D');

    Airs_Ozone_A(ix,:,:,:) = hdfread(fin, 'ascending', 'Fields','O3_VMR_A');
    Airs_Ozone_D(ix,:,:,:) = hdfread(fin, 'descending','Fields','O3_VMR_D');

    Airs_CH4_A(ix,:,:,:) = hdfread(fin, 'ascending', 'Fields','CH4_VMR_A');
    Airs_CH4_D(ix,:,:,:) = hdfread(fin, 'descending','Fields','CH4_VMR_D');

    Airs_CO_A(ix,:,:,:) = hdfread(fin, 'ascending', 'Fields','CO_VMR_A');
    Airs_CO_D(ix,:,:,:) = hdfread(fin, 'descending','Fields','CO_VMR_D');

    %if ~exist('Airs_Lat','var')
    %  Airs_Lat = hdfread(fin, 'location', 'Fields', 'Latitude');
    %  Airs_Lon = hdfread(fin, 'location', 'Fields', 'Longitude');
    %end

    Airs_SPres_A(ix,:,:) = hdfread(fin, 'ascending', 'Fields','SurfPres_Forecast_A');
    Airs_SPres_D(ix,:,:) = hdfread(fin, 'descending', 'Fields','SurfPres_Forecast_D');

    Airs_OLR_A(ix,:,:) = hdfread(fin, 'ascending',  'Fields','OLR_A');
    Airs_OLR_D(ix,:,:) = hdfread(fin, 'descending', 'Fields','OLR_D');

    Airs_ClrOLR_A(ix,:,:) = hdfread(fin, 'ascending',  'Fields','ClrOLR_A');
    Airs_ClrOLR_D(ix,:,:) = hdfread(fin, 'descending', 'Fields','ClrOLR_D');

    Airs_CldPres_A(ix,:,:,:) = hdfread(fin, 'ascending',  'Fields','CoarseCloudPres_A');
    Airs_CldPres_D(ix,:,:,:) = hdfread(fin, 'descending', 'Fields','CoarseCloudPres_D');

    Airs_CldFrac_A(ix,:,:,:) = hdfread(fin, 'ascending',  'Fields','CoarseCloudFrc_A');
    Airs_CldFrac_D(ix,:,:,:) = hdfread(fin, 'descending', 'Fields','CoarseCloudFrc_D');

    Airs_LiqWater_A(ix,:,:,:) = hdfread(fin, 'ascending_MW_only',  'Fields','TotCldLiqH2O_MW_A');
    Airs_LiqWater_D(ix,:,:,:) = hdfread(fin, 'descending_MW_only', 'Fields','TotCldLiqH2O_MW_D');

    %%%%%%%%%%%%%%%%%%%%%%%%%
%    finS = Airs_filesS{i};
%    disp(finS)

%    Airs_IceOD_A(ix,:,:,:) = hdfread(finS, 'ascending',  'Fields','ice_cld_opt_dpth_A');
%    Airs_IceOD_D(ix,:,:,:) = hdfread(finS, 'descending', 'Fields','ice_cld_opt_dpth_D');
%    Airs_IceSze_A(ix,:,:,:) = hdfread(finS, 'ascending',  'Fields','ice_cld_eff_diam_A');
%    Airs_IceSze_D(ix,:,:,:) = hdfread(finS, 'descending', 'Fields','ice_cld_eff_diam_D');
%    Airs_IceT_A(ix,:,:,:) = hdfread(finS, 'ascending',  'Fields','ice_cld_temp_eff_A');
%    Airs_IceT_D(ix,:,:,:) = hdfread(finS, 'descending', 'Fields','ice_cld_temp_eff_D');

  end
  
  %% save /asl/s1/sergio/AIRS_L3/airs_L3v6_March2014.mat Airs_Date* Airs_Temp* Airs_STemp* Airs_H2OVap* Airs_Lat Airs_Lon yy mm  
  %% should be same as /airs_L3v6_March2014.mat

  figure(1); pcolor(flipud(squeeze(mean(double(Airs_RHSurf_D),1)))); shading interp; colorbar; caxis([0 120]); title('Surf RH')

  iSave = input('save these mega huge files (+1) or store in memory??? (-1) : ? ');
  if iSave > 0
    saver = ['save /asl/s1/sergio/AIRS_L3/airs_L3v7_' savestr_version '.mat  Airs_Date* Airs_Temp* Airs_STemp* Airs_H2OVap* Airs_RH* Airs_Lat Airs_Lon Airs_CO* Airs_CH4* yy mm'];
    eval(saver);
    saver = ['save /asl/s1/sergio/AIRS_L3/airs_L3v7_extra_' savestr_version '.mat Airs_Oz* Airs_PQ Airs_PT Airs_SPres* Airs_Date* Airs_Cl* Airs_Ice* Airs_Lat Airs_Lon Airs_Liq*  yy mm'];
    eval(saver);
  end
  
else

  disp('loading in data A')
  %% should be same as /airs_L3v7_March2014.mat
  loader = ['load /asl/s1/sergio/AIRS_L3/airs_L3v7_' savestr_version '.mat'];
  eval(loader);

  disp('loading in data B')
  loader = ['load /asl/s1/sergio/AIRS_L3/airs_L3v7_extra_' savestr_version '.mat'];
  eval(loader);

end

error('ooo')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_the_AIRSL3_trends
