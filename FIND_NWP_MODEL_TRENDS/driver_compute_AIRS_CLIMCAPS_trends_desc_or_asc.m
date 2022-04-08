%% copied from /home/sergio/MATLABCODE/AIRS_L3/driver_airsL3_ratesv6_Sept2017_equalarea_latbins.m

addpath /asl/matlib/aslutil/
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/TIME
addpath /asl/matlib/science
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/FIND_TRENDS/
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

clear all

for ii = 1 : 6
  figure(ii); colormap jet
end

timespan = 19;
fprintf(1,'timespan = %2i years \n',timespan)

if timespan == 19
  savestr_version = 'Sept2002_Aug2021_19yr';
  StartY = 2002; StartYM = 9;   %% start 09/2002
  StopY  = 2021; StopYM  = 8;   %% stop  08/2021  
else
  error('huh check timespan')
end

%%% TEMPERATURE

a = read_netcdf_lls('/asl/models/CLIMCAPS_SNDR_AIRS_L3/2002/SNDR.AQUA.AIRS.20021001.M01.L3_CLIMCAPS_QCC.std.v02_38.G.210413175740.nc');

Airs_PT = a.air_pres;
Airs_PQ = a.air_pres_h2o;

%%%%%%%%%%%%%%%%%%%%%%%%%
maindir = '/asl/airs/AIRS3STM/v7/';
maindir = '/asl/models/CLIMCAPS_SNDR_AIRS_L3/'; %% https://acdisc.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_Level3/AIRS3STM.006/2018/
Airs_files = findfiles([maindir '/*/*.nc']);

Airs_Date_Start = datenum(StartY,09+(1:length(Airs_files))-1,01);
Airs_Date_Start = datenum(2002,09+(1:length(Airs_files))-1,01);   %%%% we start reading data potentially from 2002
%Airs_Date_End   = datenum(StopY,08+(1:length(Airs_files))-1,01);

%maindirS = '/asl/s1/sergio/AIRS_L3/acdisc.gsfc.nasa.gov/ftp/data/s4pa/Aqua_AIRS_Level3/AIRX3SPM.006/';
%%%%%%%%%%%%%%%%%%%%%%%%%

Airs_Date = Airs_Date_Start;

for i = 1:length(Airs_files)  
  fin = Airs_files{i};
  oops = findstr(fin,'.M01.L3_CLIMCAPS');
  yy(i) = str2num(fin(oops-8+0:oops-8+3));
  mm(i) = str2num(fin(oops-8+4:oops-8+5));

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
Airs_Temp_A = zeros(length(woo),length(Airs_PT),180,360,'single');
Airs_Temp_D = zeros(length(woo),length(Airs_PT),180,360,'single');
Airs_STemp_A = zeros(length(woo),180,360,'single');
Airs_STemp_D = zeros(length(woo),180,360,'single');
Airs_H2OVap_A = zeros(length(woo),length(Airs_PQ),180,360,'single');  %% really MMR in g/kg dry air
Airs_H2OVap_D = zeros(length(woo),length(Airs_PQ),180,360,'single');  %% really MMR in g/kg dry air
Airs_RH_A     = zeros(length(woo),length(Airs_PQ),180,360,'single');  
Airs_RH_D     = zeros(length(woo),length(Airs_PQ),180,360,'single');  
%Airs_RHSurf_A = zeros(length(woo),180,360,'single');  
%Airs_RHSurf_D = zeros(length(woo),180,360,'single');  

%% Airs_Ozone_A = zeros(length(woo),length(Airs_PT),180,360,'single');
%% Airs_Ozone_D = zeros(length(woo),length(Airs_PT),180,360,'single');
%% Airs_CO_A = zeros(length(woo),length(Airs_PT),180,360,'single');
%% Airs_CO_D = zeros(length(woo),length(Airs_PT),180,360,'single');
%% Airs_CH4_A = zeros(length(woo),length(Airs_PT),180,360,'single');
%% Airs_CH4_D = zeros(length(woo),length(Airs_PT),180,360,'single');
Airs_Ozone_A = zeros(length(woo),180,360,'single');
Airs_Ozone_D = zeros(length(woo),180,360,'single');
Airs_CO_A = zeros(length(woo),180,360,'single');
Airs_CO_D = zeros(length(woo),180,360,'single');
Airs_CH4_A = zeros(length(woo),180,360,'single');
Airs_CH4_D = zeros(length(woo),180,360,'single');
Airs_SPres_A = zeros(length(woo),180,360,'single');
Airs_SPres_D = zeros(length(woo),180,360,'single');
%Airs_OLR_A = zeros(length(woo),180,360,'single');
%Airs_OLR_D = zeros(length(woo),180,360,'single');
%Airs_ClrOLR_A = zeros(length(woo),180,360,'single');
%Airs_ClrOLR_D = zeros(length(woo),180,360,'single');
Airs_CldPres_A = zeros(length(woo),1,180,360,'single');
Airs_CldPres_D = zeros(length(woo),1,180,360,'single');
%Airs_CldFrac_A = zeros(length(woo),3,180,360,'single');
%Airs_CldFrac_D = zeros(length(woo),3,180,360,'single');

%% from MW
%Airs_LiqWater_A = zeros(length(woo),180,360,'single');
%Airs_LiqWater_D = zeros(length(woo),180,360,'single');

%% from support
%Airs_IceOD_A = zeros(length(woo),180,360,'single');  %% ice_cld_opt_dpth
%Airs_IceOD_D = zeros(length(woo),180,360,'single');
%Airs_IceSze_A = zeros(length(woo),180,360,'single');  %% ice_cld_eff_diam
%Airs_IceSze_D = zeros(length(woo),180,360,'single');
%Airs_IceT_A = zeros(length(woo),180,360,'single');  %% ice_cld_temp_eff
%Airs_IceT_D = zeros(length(woo),180,360,'single');

clear yy mm

iDo = input('(+1) read in saved mat file or (-1) read in L3 NUCAPS files, one at a time  ? ');
if iDo < 0
  clear yy mm
  for ix = 1:length(woo)
    i = woo(ix);
    fprintf(1,'reading in %3i of %3i required files -- file %3i in overall list of length %3i \n',ix,length(woo),i,length(read_this_file))
    
    fin = Airs_files{i};
    disp(fin)

    oops = findstr(fin,'.M01.L3_CLIMCAPS');
    yy(i) = str2num(fin(oops-8+0:oops-8+3));
    mm(i) = str2num(fin(oops-8+4:oops-8+5));

    a = read_netcdf_lls(fin);

    Airs_Temp_A(ix,:,:,:) = permute(squeeze(a.air_temp(:,:,:,1)),[3 2 1]);
    Airs_Temp_D(ix,:,:,:) = permute(squeeze(a.air_temp(:,:,:,2)),[3 2 1]);
    if ~ exist('Airs_Lat','var')
      Airs_Lat = a.lat;
      Airs_Lon = a.lon;
    end

    Airs_STemp_A(ix,:,:) = squeeze(a.surf_temp(:,:,1))';
    Airs_STemp_D(ix,:,:) = squeeze(a.surf_temp(:,:,2))';

    Airs_H2OVap_A(ix,:,:,:) = permute(squeeze(a.spec_hum(:,:,:,1)),[3 2 1]);
    Airs_H2OVap_D(ix,:,:,:) = permute(squeeze(a.spec_hum(:,:,:,2)),[3 2 1]);

    Airs_RH_A(ix,:,:,:) = permute(squeeze(a.rel_hum(:,:,:,1)),[3 2 1]);
    Airs_RH_D(ix,:,:,:) = permute(squeeze(a.rel_hum(:,:,:,2)),[3 2 1]);
    %Airs_RHSurf_A(ix,:,:) = 
    %Airs_RHSurf_D(ix,:,:) = 

    Airs_Ozone_A(ix,:,:) = squeeze(a.o3_tot(:,:,1))';
    Airs_Ozone_D(ix,:,:) = squeeze(a.o3_tot(:,:,2))';

    Airs_CH4_A(ix,:,:) = squeeze(a.ch4_mmr_midtrop(:,:,1))';
    Airs_CH4_D(ix,:,:) = squeeze(a.ch4_mmr_midtrop(:,:,2))';

    Airs_CO_A(ix,:,:,:) = squeeze(a.co_mmr_midtrop(:,:,1))';
    Airs_CO_D(ix,:,:,:) = squeeze(a.co_mmr_midtrop(:,:,2))';

    Airs_CO2_A(ix,:,:,:) = squeeze(a.co2_vmr_uppertrop(:,:,1))';
    Airs_CO2_D(ix,:,:,:) = squeeze(a.co2_vmr_uppertrop(:,:,2))';

    Airs_SPres_A(ix,:,:) = squeeze(a.prior_surf_pres(:,:,1))';
    Airs_SPres_D(ix,:,:) = squeeze(a.prior_surf_pres(:,:,2))';

    %Airs_OLR_A(ix,:,:) = 
    %Airs_OLR_D(ix,:,:) = 

    %Airs_ClrOLR_A(ix,:,:) = 
    %Airs_ClrOLR_D(ix,:,:) = 

    Airs_CldPres_A(ix,:,:,:) = squeeze(a.cld_top_pres(:,:,1))';
    Airs_CldPres_D(ix,:,:,:) = squeeze(a.cld_top_pres(:,:,2))';

    %Airs_CldFrac_A(ix,:,:,:) = 
    %Airs_CldFrac_D(ix,:,:,:) = 

    %Airs_LiqWater_A(ix,:,:,:) = 
    %Airs_LiqWater_D(ix,:,:,:) = 

  end
  
  %% save /asl/s1/sergio/AIRS_L3/airs_L3v6_March2014.mat Airs_Date* Airs_Temp* Airs_STemp* Airs_H2OVap* Airs_Lat Airs_Lon yy mm  
  %% should be same as /airs_L3v6_March2014.mat

  figure(1); pcolor(flipud(squeeze(mean(double(Airs_STemp_D),1)))); shading interp; colorbar; caxis([220 320]); title('Surf Temp')
  figure(1); pcolor((squeeze(mean(double(Airs_STemp_D),1)))); shading interp; colorbar; caxis([220 320]); title('Surf Temp')

  iSave = input('save these mega huge files (+1) or store in memory??? (-1) : ? ');
  if iSave > 0
    saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airs_climcaps_' savestr_version '.mat  Airs_Date* Airs_Temp* Airs_STemp* Airs_H2OVap* Airs_RH* Airs_Lat Airs_Lon Airs_CO* Airs_CH4* yy mm'];
    eval(saver);
    saver = ['save /asl/s1/sergio/AIRS_CLIMCAPS/airs_climcaps_extra_' savestr_version '.mat Airs_Oz* Airs_PQ Airs_PT Airs_SPres* Airs_Date* Airs_Cl* Airs_Lat Airs_Lon yy mm'];
    eval(saver);
  end
  
else

  disp('loading in data A')
  %% should be same as /airs_L3v7_March2014.mat
  loader = ['load /asl/s1/sergio/AIRS_CLIMCAPS/airs_climcaps_' savestr_version '.mat'];
  eval(loader);

  disp('loading in data B')
  loader = ['load /asl/s1/sergio/AIRS_CLIMCAPS/airs_climcaps_extra_' savestr_version '.mat'];
  eval(loader);

end

error('ooo')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); pcolor((squeeze(mean(double(Airs_STemp_D),1)))); shading interp; colorbar; caxis([220 320]); title('Surf Temp'); colormap(jet);
for ii = 1 : 180
  if mod(ii,100) == 0
    fprintf(1,'+')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end
  for jj = 1 : 360
    junk = squeeze(Airs_STemp_D(:,ii,jj));
    boo = Math_tsfit_lin_robust((1:228)*30,junk,4);
    quickSTrate(ii,jj) = boo(2);
  end
end
fprintf(1,'\n');
figure(2); pcolor(quickSTrate); shading interp; colorbar; caxis([-1 +1]*0.15); title('Surf Temp Rate'); colormap(usa2);
addpath /umbc/xfs2/strow/asl/matlib/maps/aslmap.m
addpath /home/sergio/MATLABCODE/COLORMAP/LLSMAPS
load llsmap5
figure(2); aslmap(2,-90:1:+90,-180:1:+180,quickSTrate,[-90 +90],[-180 +180]);  colormap(llsmap5); caxis([-0.15 +0.15]);

iL3orCLIMCAPS = -1;
do_the_AIRSL3_trends