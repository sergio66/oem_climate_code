%% now plot all results
addpath /home/sergio/MATLABCODE/TIME

clear all

i1231 = 731; % from Larrabee

iOBSorCAL = 0;  %% OBS
iOBSorCAL = 1;  %% CAL

change_important_topts_settings;
iOBSorCAL = topts.ocb_set;
clear topts

iAvgNumDays = 180;
iAvgNumDays = 000;

dataset = 1; %% CRIS

co2lays = 3; %% lower trop, upper trop, strat
co2lays = 1; %% coljac

if iOBSorCAL == 0 & dataset == 1
  disp('iOBSorCAL = 0;  %% CRIS OBS')
elseif iOBSorCAL == 1 & dataset == 1
  disp('iOBSorCAL = 1;  %% CRIS CAL')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iLatbin = 20;
if iOBSorCAL == 0 & dataset == 1
  fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
elseif iOBSorCAL == 1 & dataset == 1
  fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
end

load(fin)
save_days = avg_doy_since2012;
save_rtime = avg16_rtime;
save_dat_1231 = zeros(length(save_days),40);
save_days_map = zeros(length(save_days),40);

iReadLatbin = 30;
iReadLatbin = 20;
for iLatbin = iReadLatbin
  if iOBSorCAL == 0 & dataset == 1
    fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '.mat'];
  elseif iOBSorCAL == 1 & dataset == 1
    fin = ['ANOM_16dayavg/latbin_' num2str(iAvgNumDays) 'dayavg_' num2str(iLatbin) '_cal.mat'];
  end
  load(fin)

  figure(2); clf
  pcolor(2012+(avg_doy_since2012)/365,1:1305,avg16_btanom'); colormap jet; colorbar; shading interp;
  caxis([-2 +2]);
  xlabel('days since 1/1/2012'); ylabel('Channel ID')
  title(num2str(iLatbin))

  save_dat_day(:,iLatbin) = [avg_doy_since2012(1) avg_doy_since2012(end) length(avg_doy_since2012)];

  if length(avg_doy_since2012) == length(save_days)
    save_dat_1231 = avg16_btanom(:,i1231);
    save_days_map = 1 : 365;
  else
    for dd = 1 : length(avg_doy_since2012)
      moo = abs(avg_doy_since2012(dd)-save_days);
      moo = find(moo == min(moo));
      save_days_map(moo) = moo;
      save_dat_1231(moo) = avg16_btanom(dd,i1231);
    end
  end
  pause(0.1);
end

figure(1)
plot(2012+save_days/365,save_dat_1231'); title('ANOM BT1231');
shading flat

% save save_365datemaps.mat save*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/MATLABCODE/PLOTTER

iaFound = zeros(size(save_dat_1231))';

disp('now keep running read_the_anom_retrievals  (and/or read_the_anom_retrievals_spectra)')
read_the_anom_T_WV_retrievals_onelatbin
%read_the_anom_retrievals_spectra
