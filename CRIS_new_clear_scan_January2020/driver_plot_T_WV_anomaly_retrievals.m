%% now plot all results
addpath /home/sergio/MATLABCODE/TIME

clear all

i1231 = 731; % from Larrabee

iOBSorCAL = 0;  %% OBS
iOBSorCAL = 1;  %% CAL

change_important_topts_settings;
iOBSorCAL = topts.ocb_set;
clear topts

iMaxTimeSteps = 157;

iAvgNumDays = 180;
iAvgNumDays = 000;

dataset = 2; %% IASI
dataset = 1; %% AIRS

co2lays = 3; %% lower trop, upper trop, strat
co2lays = 1; %% coljac

if iOBSorCAL == 0 & dataset == 1
  disp('iOBSorCAL = 0;  %% AIRS OBS')
elseif iOBSorCAL == 1 & dataset == 1
  disp('iOBSorCAL = 1;  %% AIRS CAL')
elseif iOBSorCAL == 0 & dataset == 2
  disp('iOBSorCAL = 0;  %% IASI OBS')
elseif iOBSorCAL == 1 & dataset == 2
  disp('iOBSorCAL = 1;  %% IASI CAL')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iLatbin = 20;
if iOBSorCAL == 0 & dataset == 1
  fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
elseif iOBSorCAL == 1 & dataset == 1
  fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
elseif iOBSorCAL == 0 & dataset == 2
  fin = ['IASI_ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
elseif iOBSorCAL == 1 & dataset == 2
  fin = ['IASI_ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
end

load(fin)
save_days = avg_doy_since2012;
save_rtime = avg16_rtime;
save_dat_1231 = zeros(length(save_days),40);
save_days_map = zeros(length(save_days),40);

for iLatbin = 1 : 40
  if iOBSorCAL == 0 & dataset == 1
    fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
  elseif iOBSorCAL == 1 & dataset == 1
    fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
  elseif iOBSorCAL == 0 & dataset == 2
    fin = ['IASI_ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
  elseif iOBSorCAL == 1 & dataset == 2
    fin = ['IASI_ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
  end
  load(fin)

  figure(2); clf
  pcolor(2012+(avg_doy_since2012)/365,1:1305,avg16_btanom'); colormap jet; colorbar; shading interp;
  caxis([-2 +2]);
  xlabel('days since 1/1/2012'); ylabel('Channel ID')
  title(num2str(iLatbin))

  save_dat_day(:,iLatbin) = [avg_doy_since2012(1) avg_doy_since2012(end) length(avg_doy_since2012)];

  if length(avg_doy_since2012) == length(save_days)
    save_dat_1231(:,iLatbin) = avg16_btanom(:,i1231);
    save_days_map(:,iLatbin) = 1 : iMaxTimeSteps;
  else
    for dd = 1 : length(avg_doy_since2012)
      moo = abs(avg_doy_since2012(dd)-save_days);
      moo = find(moo == min(moo));
      save_days_map(moo,iLatbin) = moo;
      save_dat_1231(moo,iLatbin) = avg16_btanom(dd,i1231);
    end
  end
  pause(0.1);
end

figure(1); clf
pcolor(2012+save_days/365,1:40,save_dat_1231'); caxis([-1 +1]); colorbar; title('ANOM BT1231');
shading flat

% save save_365datemaps.mat save*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/MATLABCODE/PLOTTER

iaaFound = zeros(size(save_dat_1231))';

disp('now keep running read_the_anom_retrievals  (and/or read_the_anom_retrievals_spectra)')
read_the_anom_T_WV_retrievals
%read_the_anom_retrievals_spectra
