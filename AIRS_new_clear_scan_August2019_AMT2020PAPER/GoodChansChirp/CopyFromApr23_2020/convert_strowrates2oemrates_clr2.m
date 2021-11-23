addpath /home/sergio/MATLABCODE

%fairs = instr_chans;    %% 2378

load sarta_chans_for_l1c.mat

hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
fairs = vchan2834;     %% 2834

fairs = fairs(ichan);  %% 2645  

for ix = 1 : 40
  fprintf(1,'ix = %2i \n',ix);

%% for Mar 23, 2019
  fname1 = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/Small/fit_nucal_robs_lat' num2str(ix) '.mat'];   %% exists
  fname2 = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/Small/fit_nucal_rclr_lat' num2str(ix) '.mat'];   %% DNE
  fname3 = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/Small/fit_nucal_bias_lat' num2str(ix) '.mat'];   %% DNE

%% for May 31, 2019
  fname1 = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/Small/fit_nucal_robs_lat' num2str(ix) '.mat'];  %% exists
  fname2 = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/Small/fit_rclr_lat' num2str(ix) '.mat'];        %% exists
  fname3 = ['/home/strow/Work/Airs/Stability/Data/Desc_fits/Small/fit_bias_lat' num2str(ix) '.mat'];        %% exists

  clear a
  if exist(fname1)
    a = load(fname1);
    b_obs(ix,:,2)     = a.dbt;
    b_err_obs(ix,:,2) = a.dbt_err;    
    lagcor_obs_anom(ix,:) = a.lag;
  end

  clear a
  if exist(fname2)
    a = load(fname2);
    b_cal(ix,:,2)     = a.dbt;
    b_err_cal(ix,:,2) = a.dbt_err;
    lagcor_cal_anom(ix,:) = a.lag;
  end

  %% for Mar 23, 2019
%{
  clear a
  if exist(fname3)
    a = load(fname3);
    b_bias(ix,:,2)     = a.dbt;
    b_err_bias(ix,:,2) = a.dbt_err;
    lagcor_bias_anom(ix,:) = a.lag;
  end
%}

  %% for May 31, 2019
  clear a
  if exist(fname1) & exist(fname2)
    a1 = load(fname1);
    a2 = load(fname2);
    fprintf(1,'doing bias for May 31, 2019 version, latbin %2i \n',ix)
    b_bias(ix,:,2)     = a1.dbt-a2.dbt;
    b_err_bias(ix,:,2) = a1.dbt_err;
    lagcor_bias_anom(ix,:) = a1.lag;
  end
end

comment = 'from /home/strow/Work/Airs/Stability/Data/Desc_fits/Small Feb 26, 2019 and convert_strowrates2oemrates_clr2.m';
%% for Mar 23, 2019
%% save convert_strowrates2oemrates_random_16_year_v32_clear_nucal.mat  b_* comment lagcor*

%% for May 29. 2019
save convert_strowrates2oemrates_random_16_year_v32_clear_nucal_obs_cal_bias.mat  b_* comment lagcor*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); plot(fairs,squeeze(b_obs(:,:,2))); ax = axis; ax(3) = -0.25; ax(4) = +0.25; axis(ax);
figure(2); plot(fairs,lagcor_obs_anom); ax = axis; ax(3) = -1; ax(4) = +1; axis(ax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
if ye wanna have dem 2378 chans

cp /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/map2645to2378.mat
mapp = load('map2645to2378.mat');
plot(mapp.closest2645to2378_freq,squeeze(b_obs(:,mapp.closest2645to2378_ind,2)));

b_bias     = b_bias(:,mapp.closest2645to2378_ind,:);
b_err_bias = b_err_bias(:,mapp.closest2645to2378_ind,:);
b_cal     = b_cal(:,mapp.closest2645to2378_ind,:);
b_err_cal = b_err_cal(:,mapp.closest2645to2378_ind,:);
b_obs     = b_obs(:,mapp.closest2645to2378_ind,:);
b_err_obs = b_err_obs(:,mapp.closest2645to2378_ind,:);
save convert_strowrates2oemrates_random_16_year_v32.mat  b_* comment
%}
