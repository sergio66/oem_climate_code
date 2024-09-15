%% see driver_computeERA5_monthly_trends_desc_or_asc.m
%% see driver_computeERA5_monthly_trends_NIGHT_or_DAY_or_BOTH.m

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
addpath /asl/matlib/plotutils
addpath /asl/matlib/maps

load('llsmap5.mat');

%% note this code only handles 2002_09 to YYYY_08 sp
%%      this code does not currently handle eg OCO2 ERA5_atm_N_cld_data_2012_05_to_2019_04_trends_desc.mat

can start afresh from here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('takes about 5 hours for trends')

iEnsureSave = -1;   %% only save if file DNE (big fat file, or trends)
iEnsureSave = +1;   %% this makes sure big fat file is saved, as are trend files, whether or not files have already been made

if iEnsureSave < 0
  disp('iEnsureSave = -1 so files only saved (big fat file, trends etc) if they DNE')
elseif iEnsureSave > 0
  disp('iEnsureSave = +1 so files always saved (big fat file, trends etc) even if they already exist')
end
disp('if you have ALREADY saved BIG FAT FILE also look for lines that say     can start afresh from here       ')
disp('if you have ALREADY saved BIG FAT FILE also look for lines that say     can start afresh from here       ')
disp('<RET> to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% 1 : 20 for the iNumYears 
if length(JOB) == 0
  JOB = 20;
  JOB = 22;
end

if ~exist('iTrendsOrAnoms')
  iTrendsOrAnoms = -1;  %% surface and atmospheric anomalies
  iTrendsOrAnoms = -10; %% surface                 anomalies

  iTrendsOrAnoms = +10; %% surface                 trends
  iTrendsOrAnoms = +1;  %% surface and atmospheric trends
end

if ~exist('iDorA')
  iDorA = -10; %% asc,  random pert about tile center
  iDorA = +10; %% desc, random pert about tile center
  iDorA = -1;  %% asc,  tile center
  iDorA = +1;  %% desc, tile center
end
fprintf(1,'iDorA = %2i \n',iDorA)

clear iaFound
%iaMax = 18*12; %% 18 year
%iaMax = 19*12; %% 19 year

iNumYears = 18; %% 2002/09-2020/08
iNumYears = 19; %% 2002/09-2021/08
iNumYears = 12; %% 2002/09-2014/08
iNumYears = 070; %% 2012/05-2019/04 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iNumYears = 20; %% 2002/09-2022/08
iNumYears = 22; %% 2002/09-2024/06

iCnt = 0;
for yy = 2002 : 2024
  mmS = 1; mmE = 12;
  if yy == 2002
    mmS = 09;
  elseif yy == 2024
    mmE = 06;
    mmE = 08;
  end
  for mm = mmS : mmE
    iCnt = iCnt + 1;
    yysave(iCnt) = yy;
    mmsave(iCnt) = mm;
  end
end

yymmS = [2020 07]; yymmE = [2024 06];
yymmS = [2002 09]; yymmE = [2024 06];
yymmS = [2002 09]; yymmE = [2024 08];

iaMin = find(yymmS(1) == yysave & yymmS(2) == mmsave);
iaMax = find(yymmE(1) == yysave & yymmE(2) == mmsave);

iAllorSeasonal = -4;  %% SON
iAllorSeasonal = -3;  %% JJA
iAllorSeasonal = -2;  %% MAM
iAllorSeasonal = -1;  %% DJF
iAllorSeasonal = +1;  %% all

if ~exist('iCldORClr')
  iCldORClr = -1; %% clear only
  iCldORClr = +1; %% include cloud fields
end

find_computeERA5_monthly_trends_foutname

fout_trendjunk_surface = [fout_trendjunk(1:end-4) '_surf.mat'];
if exist(fout_trendjunk) & exist(fout_trendjunk_surface)
  fprintf(1,'  fout_trendjunk = %s \n',fout_trendjunk)
  fprintf(1,'  fout_trendjunk_surface = %s \n',fout_trendjunk_surface)
  error('fout_trendjunk,fout_trendjunk_surface exist')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% see /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/clust_loop_make_monthly_tile_center_asc_or_desc.m
fprintf(1,'seeing how many of the expected %3i files exist \n',iaMax);
disp('  "+" are the 100, "x" are tens .. ')

iOLR = +1;
for ii = iaMin : iaMax
  % iForwardFrom2002_or_BackwardFrom2022 = +1;
  %if mod(ii,100) == 0 & iForwardFrom2002_or_BackwardFrom2022 > 0
  %  fprintf(1,'+ \n')
  %elseif mod(ii,10) == 0 & iForwardFrom2002_or_BackwardFrom2022 > 0
  %  fprintf(1,'x')
  %elseif iForwardFrom2002_or_BackwardFrom2022 > 0
  %  fprintf(1,'.')
  %end

  if iOLR < 0
    if iDorA == 1
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -1
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == 10
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -10
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    end
  else
    if iDorA == +1
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -1
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC_WithOLR/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == +10
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -10
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC_WithOLR/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    end
  end

  if exist(fin)
    fprintf(1,'%3i : %s exists \n',ii,fin)
    iaFound(ii) = 1;
    moo = load(fin,'pnew_op');
    iaNumProf(ii) = length(moo.pnew_op.stemp);
    if length(moo.pnew_op.stemp) ~= 4608
      fprintf(1,'OH NO %s %4i has %4i profiles instead of 4608 \n',fin,ii,iaNumProf(ii))
    end
  else
    iaFound(ii) = 0;
    iaNumProf(ii) = -1;
    fprintf(1,'%3i : %s DNE \n',ii,fin)
  end
end

fprintf(1,'iNumYears = %3i : found %3i of expected %3i files \n',[iNumYears sum(iaFound) length(iaFound)])
bad = find(iaFound > 0 & iaNumProf ~= 4608);
if length(bad) > 0
  fprintf(1,'found the following %3i bad files ... please redo \n',length(bad))
  bad
  error('bad files ... redo')
end

%iForwardFrom2002_or_BackwardFrom2022 = +1;
plot(1:iaMax,iaFound,'+-')
if sum(iaFound) < length(iaFound(iaMin:iaMax)) 
  printarray(find(iaFound == 0))
  error('not enough')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('reading in monthly ERA5 data in /asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ ')
disp('  made by /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_ERA5/clust_loop_make_monthly_tile_center_asc_or_desc.m')
fprintf(1,' "+" are the 100, "x" are tens .. need to read in %3i \n',iaMax)
for ii = iaMin : iaMax
  if iOLR < 0
    if iDorA == 1
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -1
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == 10
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -10
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    end
  else
    if iDorA == +1
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -1
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC_WithOLR/era5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == +10
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/DESC_WithOLR/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    elseif iDorA == -10
      fin = ['/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center/ASC_WithOLR/randomptera5_tile_center_monthly_' num2str(ii,'%03d') '.mat'];
    end
  end

  iii = ii - iaMin + 1;
  if exist(fin)
    fprintf(1,'%3i : reading %s \n',ii,fin)
    % if mod(ii,100) == 0
    %   fprintf(1,'+ \n')
    % elseif mod(ii,10) == 0
    %   fprintf(1,'x')
    % else
    %   fprintf(1,'.')
    % end
    iaFound(ii) = 1;

    a = load(fin);
    %% already done this in /home/sergio/MATLABCODE/RTPMAKE/CLUST_RTPMAKE/CLUSTMAKE_CMIP6/clust_compute_cmip6_profile_rtpfiles.m
    %a.pnew_ip.rh = convert_humidity(a.pnew_ip.plevs*100,a.pnew_ip.ptemp,a.pnew_ip.gas_1,'mixing ratio','relative humidity');
    if ~isfield(a.pnew_ip,'rh')
      a.pnew_ip.rh = convert_humidity(a.pnew_ip.plevs*100,a.pnew_ip.ptemp,a.pnew_ip.gas_1,'specific humidity','relative humidity');
    end

    all.yy(iii) = a.thedateS(1);
    all.mm(iii) = a.thedateS(2);
    all.dd(iii) = a.thedateS(3);

    all.nwp_ptemp(iii,:,:) = a.pnew_ip.ptemp;
    all.nwp_gas_1(iii,:,:) = a.pnew_ip.gas_1;
    all.nwp_gas_3(iii,:,:) = a.pnew_ip.gas_3;
    all.nwp_rh(iii,:,:)    = a.pnew_ip.rh;
    all.nwp_plevs(iii,:,:) = a.pnew_ip.plevs;

    all.gas_1(iii,:,:) = a.pnew_op.gas_1;
    all.gas_3(iii,:,:) = a.pnew_op.gas_3;
    all.ptemp(iii,:,:) = a.pnew_op.ptemp;
    all.stemp(iii,:)   = a.pnew_op.stemp;
    all.mmw(iii,:)     = a.pnew_op.mmw;
    all.nlays(iii,:)   = a.pnew_op.nlevs-1;
    all.RH(iii,:,:)    = a.pnew_op.RH;
    all.TwSurf(iii,:)  = a.pnew_op.TwSurf;
    all.RHSurf(iii,:)  = a.pnew_op.RHSurf;
    
    if iOLR > 0
      if isfield(a.pnew_op,'d2m')
        iYes2m = +1;
        all.d2m(iii,:)     = a.pnew_op.d2m;
        all.t2m(iii,:)     = a.pnew_op.t2m;

        %% Improved Magnus Form Approximation of Saturation Vapor Pressure
        %% Oleg A. Alduchov  and Robert E. Eskridge
        %% JAMS 1996, v35 DOI: https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2  Page(s): 601–609
        all.RH2m(iii,:) = 100 * exp((17.625*(a.pnew_op.d2m-273.13)./(243.04+(a.pnew_op.d2m-273.13))));
        all.RH2m(iii,:) = all.RH2m(iii,:) ./ exp((17.625*(a.pnew_op.t2m-273.13)./(243.04+(a.pnew_op.t2m-273.13))));;

        %% from "Understanding variations in downwelling longwave radiation using Brutsaert’s equation"
        %% Earth Syst. Dynam., 14, 1363–1374, 2023 https://doi.org/10.5194/esd-14-1363-2023, eqn 6
        all.e2a(iii,:) = 6.1079*exp(17.269 * (a.pnew_op.d2m-273.13)./(237.3 + (a.pnew_op.d2m-273.13)));            
        all.ecs(iii,:) = 1.24*(all.e2a(iii,:)./a.pnew_op.t2m).^(1/7);
        all.Rld(iii,:) = 5.67e-8 * all.ecs(iii,:) .* a.pnew_op.t2m.^(4); 

      else
        iYes2m = -1;
        all.d2m(iii,:)     = zeros(size(a.pnew_op.stemp));
        all.t2m(iii,:)     = zeros(size(a.pnew_op.stemp));

        all.RH2m(iii,:)     = zeros(size(a.pnew_op.stemp));
        all.e2a(iii,:)     = zeros(size(a.pnew_op.stemp));
        all.ecs(iii,:)     = zeros(size(a.pnew_op.stemp));
        all.Rld(iii,:)     = zeros(size(a.pnew_op.stemp));
        
      end

      if ~isfield(a.pnew_op,'olr') & iOLR > 0
        fprintf(1,'%s \n',fin)
        error('no olr field but iOLR > 0')
      end
      if isfield(a.pnew_op,'olr')
        all.olr(iii,:)     = a.pnew_op.olr;
        all.olr_clr(iii,:) = a.pnew_op.olr_clr;
      else
        all.olr(iii,:)     = nan;
        all.olr_clr(iii,:) = nan;
      end
      
      if ~isfield(a.pnew_op,'ilr') & iOLR > 0
        disp('WARNING : no ilr field but iOLR > 0')
      end
      if isfield(a.pnew_op,'ilr')
        %% now do adjustments, since I think ILR is a net flux instead of actual downwelling flux
        all.ilr_adj(iii,:) = -a.pnew_op.ilr_clr + 5.67e-8 * a.pnew_op.stemp.^(4);
        all.ilr(iii,:)     = a.pnew_op.ilr;
        all.ilr_clr(iii,:) = a.pnew_op.ilr_clr;
      else
        all.ilr_adj(iii,:) = nan;
        all.ilr(iii,:)     = nan;
        all.ilr_clr(iii,:) = nan;
      end
    end

    if iCldORClr == +1
      hunk = a.hnew_op;
      junk = a.pnew_op;
      %junk = find_average_rtp(hunk,junk);
      all.ctype(iii,:)  = a.pnew_op.ctype;
      all.cfrac(iii,:)  = a.pnew_op.cfrac;
      all.cpsize(iii,:) = a.pnew_op.cpsize;
      all.cngwat(iii,:) = a.pnew_op.cngwat;
      all.cprtop(iii,:) = a.pnew_op.cprtop;
      all.cprbot(iii,:) = a.pnew_op.cprbot;

      all.ctype2(iii,:)  = a.pnew_op.ctype2;
      all.cfrac2(iii,:)  = a.pnew_op.cfrac2;
      all.cpsize2(iii,:) = a.pnew_op.cpsize2;
      all.cngwat2(iii,:) = a.pnew_op.cngwat2;
      all.cprtop2(iii,:) = a.pnew_op.cprtop2;
      all.cprbot2(iii,:) = a.pnew_op.cprbot2;
      
      all.cfrac12(iii,:)  = a.pnew_op.cfrac12;
    end
  else
    iaFound(iii) = 0;
  end
end
fprintf(1,'\n');
all.rlon = a.pnew_op.rlon;
all.rlat = a.pnew_op.rlat;

monitor_memory_whos

comment = 'see computeERA5_trends.m';
comment = 'see driver_computeERA5_monthly_trends_desc_or_asc.m';

find_computeERA5_monthly_foutname

iSave = -1;  %%% unless now you do a whole new set of files Steve brings down from ERA5
             %%% be careful though, first time you run this, you need it DUDE and iSave should be +1   
             %%%     eg when I did the daytime   "asc"  2002-2022 I forgot to save, so could not generate spectral trends boo hoo
             %%%     eg when I did the nighttime "desc" 2002-2024 I forgot to save, so could not generate spectral trends boo hoo

if iNumYears <= 100
  iSave = -1;
else
  eeee = exist(foutjunk);
  if eeee > 0
    fprintf(1,'foutjunk = %s a HUGE FILE WITH ALL THE 20 years of PROFILE data already exists!!! \n',foutjunk)
  else
    fprintf(1,'foutjunk = %s a HUGE FILE WITH ALL THE 20 years of PROFILE data DNE DNE DNE ! \n',foutjunk)
  end
  iSave = input('Save huge file (-1/+1) : ');
  %%iSave = +1;
end

if ~exist(foutjunk)
  iSave = +1;   %% this is doing it first time
  fprintf(1,'BIG FAT FILE with all data %s DNE and so will NOT be able to do spectral trends using driver_check_WV_T_RH_ERA5_geo_and_spectral_rates2.m \n',foutjunk)
  iSave = input('Save BIG FAT FILES???? (-1/+1) : ');
end

fprintf(1,'iSave   iEnsureSave    foutjunk  = %2i %s  \n',iSave,iEnsureSave,foutjunk);
disp('RET to continue'); pause

if iSave > 0 | iEnsureSave > 0
  %%% foutjunk = ['ERA5_atm_data_2002_09_to_*.mat'];
  fprintf(1,'saving huge file : can type in a separate window         watch "ls -lt %s " \n',foutjunk)
  saver = ['save -v7.3 ' foutjunk ' comment all'];
  eval(saver);
end

figure(1); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(all.stemp,1)); colormap(jet); title('ERA5 mean stemp')
figure(2); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(all.RHSurf,1)); colormap(jet); title('ERA5 mean RHsurf DO NOT BELIEVE')
figure(3); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(all.TwSurf,1)); colormap(jet); title('ERA5 mean TWSurf DO BELIEVE')
figure(4); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(all.mmw,1)); colormap(jet); title('ERA5 mean mmw')

figure(5); clf; scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.stemp); colormap(jet); title('ERA5 mean stemp')
figure(6); clf; scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.RHSurf); colormap(jet); title('ERA5 mean RHsurf DO NOT BELIEVE')
figure(7); clf; scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.TwSurf); colormap(jet); title('ERA5 mean TWSurf DO BELIEVE')
figure(8); clf; scatter_coast(a.pnew_op.rlon,a.pnew_op.rlat,40,a.pnew_op.mmw); colormap(jet); title('ERA5 mean mmw')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('can start afresh from here eg    clear all; load ERA5_atm_N_cld_data_2002_09_to_2024_08_desc.mat');
disp('can start afresh from here eg    clear all; load ERA5_atm_N_cld_data_2002_09_to_2024_08_desc.mat');
disp('can start afresh from here eg    clear all; load ERA5_atm_N_cld_data_2002_09_to_2024_08_desc.mat');
%% check these :::: iOLR = +1; iCldORClr = +1; iDorA = +1; yymmS = [2002 09]; yymmE = [2024 08];  iTrendsOrAnoms = +1; find_computeERA5_monthly_trends_foutname

plevs = load('/home/sergio/MATLABCODE/airslevels.dat');
pN = plevs(1:end-1)-plevs(2:end);
pD = log(plevs(1:end-1)./plevs(2:end));
plays = flipud(pN./pD);

load /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_For_HowardObs_TimeSeries/latB64.mat
rlat65 = latB2; rlon73 = -180 : 5 : +180;
rlon = -180 : 5 : +180;  rlat = latB2;
rlon = 0.5*(rlon(1:end-1)+rlon(2:end));
rlat = 0.5*(rlat(1:end-1)+rlat(2:end));

figure(9); clf; junk = reshape(a.pnew_op.ptemp,101,72,64); junk = squeeze(nanmean(junk,2)); junk = junk(1:100,:); pcolor(rlat,plays,junk);
  caxis([200 300]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
  title('Mean T')
figure(10); clf; junk = reshape(a.pnew_op.RH,100,72,64); junk = squeeze(nanmean(junk,2)); junk = junk(1:100,:); pcolor(rlat,plays,junk);
  caxis([00 100]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar
  title('Mean RH')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dayOFtime = change2days(all.yy,all.mm,all.dd,2002);

if iTrendsOrAnoms > 0
  computeERA5_surface_trends

  figure(1); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_stemp,1)); title('ERA5 trend  stemp K/yr');        caxis([-0.2 +0.2]); colormap(usa2);
  figure(2); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_RHSurf,1)); title('ERA5 trend  RHsurf UGH pc/yr'); caxis([-0.4 +0.4]); colormap(usa2);
  figure(3); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_TwSurf,1)); title('ERA5 trend  TWSurf UGH K/yr');  caxis([-0.2 +0.2]); colormap(usa2);
  figure(4); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_mmw,1)); title('ERA5 trend  colwater mm/yr');      caxis([-0.2 +0.2]); colormap(usa2);

  figure(1); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_stemp,1)); title('ERA5 trend  stemp K/yr');        caxis([-1 +1]*0.15); colormap(llsmap5);
  figure(2); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_RHSurf,1)); title('ERA5 trend  RHsurf UGH pc/yr'); caxis([-1 +1]*0.4);  colormap(llsmap5);
  figure(3); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_TwSurf,1)); title('ERA5 trend  TWSurf UGH K/yr');  caxis([-1 +1]*0.1);  colormap(llsmap5);
  figure(4); clf; scatter_coast(all.rlon,all.rlat,40,nanmean(trend_mmw,1)); title('ERA5 trend  colwater mm/yr');      caxis([-1 +1]*0.2);  colormap(llsmap5);

  figure(1); clf; aslmap(1,rlat65,rlon73,smoothn((reshape(trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]);, title('ERA5 trend  stemp K/yr');        caxis([-1 +1]*0.15); colormap(llsmap5);
  figure(2); clf; aslmap(2,rlat65,rlon73,smoothn((reshape(trend_RHSurf,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 trend  RHsurf UGH pc/yr');  caxis([-1 +1]*0.4);  colormap(llsmap5);
  figure(3); clf; aslmap(3,rlat65,rlon73,smoothn((reshape(trend_TwSurf,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 trend  TWSurf UGH K/yr');   caxis([-1 +1]*0.1);  colormap(llsmap5);
  figure(4); clf; aslmap(4,rlat65,rlon73,smoothn((reshape(trend_mmw,72,64)') ,1), [-90 +90],[-180 +180]);    title('ERA5 trend  colwater mm/yr');    caxis([-1 +1]*0.2);  colormap(llsmap5);
  pause(0.1)

  if ~exist(fout_trendjunk_surface) | iEnsureSave > 0
    fprintf(1,'saving surface/scalar trend file : can type in a separate window         watch "ls -lt %s " \n',fout_trendjunk_surface)
    saver = ['save ' fout_trendjunk_surface ' comment trend*'];
    eval(saver);
  end

elseif iTrendsOrAnoms < 0
  computeERA5_surface_anoms

%  fout_trendjunk_surface_anom = ['anom_' fout_trendjunk_surface_anom];

  if ~exist(fout_trendjunk_surface_anom | iEnsureSave > 0
    fprintf(1,'saving surface/scalar anomaly file : can type in a separate window         watch "ls -lt %s " \n',fout_trendjunk_surface_anom)
    saver = ['save ' fout_trendjunk_surface_anom ' comment anom_*'];
    eval(saver);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iTrendsOrAnoms == 1 
  computeERA5_atmos_trends

  if isfield(all,'nwp_plevs')
    trend_nwp_plevs_mean = mean(squeeze(nanmean(all.nwp_plevs,1)),2);
  end
  
  trend_plays = flipud(pN./pD);
  
  trend_rlat = all.rlat;
  trend_rlon = all.rlon;
  trend_rlat64 = rlat; trend_rlon72 = rlon;
  %trend_plevs37 = permute(all.nwp_plevs,[2 1 3]); trend_plevs37 = reshape(trend_plevs37,37,227*4608); trend_plevs37 = mean(trend_plevs37,2);
  
  %find_computeERA5_monthly_trends_foutname
  if ~exist(fout_trendjunk) | iEnsureSave > 0
    fprintf(1,'saving trend file : can type in a separate window         watch "ls -lt %s " \n',fout_trendjunk)
    saver = ['save ' fout_trendjunk ' comment trend*'];
    eval(saver);
  end
  
  figure(5); clf; pcolor_sin(trend_rlat64,trend_plays,squeeze(nanmean(reshape(trend_ptemp,100,72,64),2))); caxis([-1 +1]*0.150); set(gca,'ydir','reverse'); colormap(llsmap5); title('ERA5 dT/dt'); shading interp
    set(gca,'yscale','log'); ylim([10 1000]); colorbar('horizontal')
  figure(6); clf; pcolor_sin(trend_rlat64,trend_plays,squeeze(nanmean(reshape(trend_gas_1,100,72,64),2))); caxis([-1 +1]*0.015); set(gca,'ydir','reverse'); colormap(llsmap5); title('ERA5 dWVfrac/dt'); shading interp
    set(gca,'yscale','linear'); ylim([100 1000]); colorbar('horizontal')

elseif iTrendsOrAnoms == -1 
  error('too long and painful, forget this')
  computeERA5_atmos_anoms

  if isfield(all,'nwp_plevs')
    anom_nwp_plevs_mean = mean(squeeze(nanmean(all.nwp_plevs,1)),2);
  end
  
  anom_plays = flipud(pN./pD);
  
  anom_rlat = all.rlat;
  anom_rlon = all.rlon;
  anom_rlat64 = rlat; anom_rlon72 = rlon;
  %anom_plevs37 = permute(all.nwp_plevs,[2 1 3]); anom_plevs37 = reshape(anom_plevs37,37,227*4608); anom_plevs37 = mean(anom_plevs37,2);
  
  %find_computeERA5_monthly_anoms_foutname
  fout_anomjunk = ['anom_' fout_trendjunk];
  if ~exist(fout_anomjunk) | iEnsureSave > 0
    fprintf(1,'saving anom file : can type in a separate window         watch "ls -lt %s " \n',fout_anomjunk)
    saver = ['save ' fout_anomjunk ' comment anom*'];
    eval(saver);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PARTIALLY BELIEVE RHSURF stuff since I use t2m and d2m for RHsurf in the CLUSTMAKE_RTP ... so did silly conversions, interpolating RH(p) to surf .. need to go back to orig rtps and fix')
disp('PARTIALLY BELIEVE RHSURF stuff since I use t2m and d2m for RHsurf in the CLUSTMAKE_RTP ... so did silly conversions, interpolating RH(p) to surf .. need to go back to orig rtps and fix')
disp('PARTIALLY BELIEVE RHSURF stuff since I use t2m and d2m for RHsurf in the CLUSTMAKE_RTP ... so did silly conversions, interpolating RH(p) to surf .. need to go back to orig rtps and fix')

figure(1); clf; scatter_coast(all.rlon,all.rlat,40,trend_stemp); title('ERA5 trend  stemp K/yr');    caxis([-1 +1]*0.15); colormap(llsmap5);
figure(2); clf; scatter_coast(all.rlon,all.rlat,40,trend_RHSurf); title('ERA5 trend  UGH RHsurf pc/yr'); caxis([-1 +1]*0.4); colormap(llsmap5);
addpath /asl/matlib/maps/
aslmap(1,rlat65,rlon73,smoothn((reshape(trend_stemp,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 dST/dt');      caxis([-1 +1]*0.15); colormap(llsmap5)
aslmap(2,rlat65,rlon73,smoothn((reshape(trend_RHSurf,72,64)') ,1), [-90 +90],[-180 +180]); title('ERA5 UGH dRHSurf/dt'); caxis([-1 +1]*0.25); colormap(llsmap5)

figure(3); junk = reshape(trend_ptemp,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer trend ptemp K/yr');  caxis([-1 +1]*0.15); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(4); junk = reshape(trend_RH,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer trend RH percent/yr');  caxis([-1 +1]*0.15); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar
figure(5); junk = reshape(trend_gas_1,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer trend WVfrac /yr');  caxis([-1 +1]*0.01); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

figure(6); junk = squeeze(nanmean(all.ptemp,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer mean ptemp K');  caxis([200 300]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(7); junk = squeeze(nanmean(all.RH,1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer mean RH percent');  caxis([0 100]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

%{
figure(8); junk = squeeze(nanstd(all.ptemp,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanstd(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev ptemp K');  caxis([00 20]; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(9); junk = squeeze(nanstd(all.RH,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanstd(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev RH percent');  caxis([0 20]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar

figure(10); junk = squeeze(max(all.ptemp,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev ptemp K');  caxis([00 20]; colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([10 1000]); colorbar
figure(11); junk = squeeze(max(all.RH,[],1)); junk = junk(1:100,:); junk = reshape(junk,100,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_plays,junk); title('ERA5 100 layer stddev RH percent');  caxis([0 20]); colormap(jet); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading flat; ylim([100 1000]); colorbar
%}

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(8); junk = reshape(trend_nwp_ptemp,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend ptemp K/yr');  caxis([-1 +1]*0.15); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([10 1000]); colorbar
figure(9); junk = reshape(trend_nwp_rh,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend RH percent/yr');  caxis([-1 +1]*0.15); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

figure(10); junk = reshape(trend_nwp_gg,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend SH g/g/yr');  caxis([0 +2.5]*1e-5); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar
figure(11); junk = reshape(trend_nwp_ppmv,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  trend PPMV /yr');  caxis([0 40]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar
figure(12); junk = reshape(trend_nwp_frac,37,72,64); junk = squeeze(nanmean(junk,2)); 
  pcolor(trend_rlat64,trend_nwp_plevs_mean,junk); title('ERA5 37 lvl  frac /yr');  caxis([-10 +10]*1e-3); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

pause(0.1);

junk1 = squeeze(nanmean(all.nwp_gas_1,1)); junk1 = reshape(junk1,37,72,64); junk1 = squeeze(nanmean(junk1,2)); 
junk2 = reshape(trend_nwp_frac,37,72,64); junk2 = squeeze(nanmean(junk2,2)); 
figure(13); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1.*junk2); title('ERA5 37 lvl  SH g/g/yr VERS2');  caxis([0 +2.5]*1e-5); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar

junk1 = toppmv(all.nwp_plevs,all.nwp_ptemp,all.nwp_gas_1,18,21); junk1 = squeeze(nanmean(junk1,1));
junk1 = reshape(junk1,37,72,64); junk1 = squeeze(nanmean(junk1,2));
junk2 = reshape(trend_nwp_frac,37,72,64); junk2 = squeeze(nanmean(junk2,2)); 
figure(14); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1);
figure(14); loglog(nanmean(junk1,2),trend_nwp_plevs_mean);  set(gca,'ydir','reverse'); xlim([1 1e4]); grid
figure(14); pcolor(trend_rlat64,trend_nwp_plevs_mean,junk1.*junk2); title('ERA5 37 lvl  PPMV/yr VERS2');  caxis([0 40]); colormap(llsmap5); set(gca,'ydir','reverse'); set(gca,'yscale','log'); shading interp; ylim([100 1000]); colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% this is MLS comparisons, Frank Werner JPL
figure(11); ylim([1 300]); colormap(llsmap5)
figure(07); ylim([1 300]);
figure(09); ylim([1 300]); caxis([-1 1]*1e-7); colormap(llsmap5)
figure(12); ylim([1 300]); caxis([-1 1]*1e-7); colormap(llsmap5)
figure(10); ylim([1 300]); caxis([-1 1]*1e-1); colormap(llsmap5)
figure(13); ylim([1 300]); caxis([-1 1]*1e-1); colormap(llsmap5)
%}
