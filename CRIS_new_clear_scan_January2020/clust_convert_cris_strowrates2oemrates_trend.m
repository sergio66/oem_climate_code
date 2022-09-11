addpath  /asl/matlib/time
addpath /home/sergio/MATLABCODE/TIME

iPlot = -1;
iDebug = -1;

%i_obsORcal = -1; %% cal
%i_obsORcal = +1; %% obs

%iAvgNumDays = 180;
%iAvgNumDays = 000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i1231 = 731;

iLS = 1; iLE = 40;

iLS = 1; iLE = 40;

for iLatbin = iLS : iLE
  clear avg16_rtime avg_doy_since2012 avg16_btanom
  fnameO = ['/home/strow/Work/Cris/Stability/Data/Desc_fits/fit_robs_lat' num2str(iLatbin) '.mat'];
  fnameC = ['/home/strow/Work/Cris/Stability/Data/Desc_fits/fit_rclr_lat' num2str(iLatbin) '.mat'];
  aO = load(fnameO);
  aC = load(fnameC);

  %% figure out "bad" data
  fprintf(1,'latbin %2i fname %s \n',iLatbin,fnameO);

  b_cal(iLatbin,:,2) = aC.dbt * 1000;
  b_err_cal(iLatbin,:,2) = aC.dbt_err  * 1000;
  lagcor_cal_anom(iLatbin,:) = aC.lag;

  b_obs(iLatbin,:,2) = aO.dbt  * 1000;
  b_err_obs(iLatbin,:,2) = aO.dbt_err  * 1000;
  lagcor_obs_anom(iLatbin,:) = aO.lag;
  
end         %% loop over latbins  iLatbin

fsave = ['cris_npp_2012_o5_to_2019_04_trends_fron_data.mat'];

comment = 'see clust_convert_cris_strowrates2oemrates_trend.m';
load f1305
saver = ['save ' fsave ' comment f1305 b_* lagcor*'];

%if ~exist(fsave)
  eval(saver)
%else
%  fprintf(1,'%s already exists \n',fsave)
%end
