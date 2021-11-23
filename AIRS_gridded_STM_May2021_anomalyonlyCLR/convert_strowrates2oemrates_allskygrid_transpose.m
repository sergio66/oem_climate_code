addpath /home/sergio/MATLABCODE

%fairs = instr_chans;    %% 2378

load sarta_chans_for_l1c.mat

hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
fairs = vchan2834;     %% 2834

fairs = fairs(ichan);  %% 2645  

%a = load('/home/strow/Work/Airs/Tiles/Data/grid_fit_var_trends_all.mat');
x = load('/home/strow/Work/Airs/Tiles/Data/grid_fit_trends_all.mat');
a.dbt     = permute(x.dbt,[3 1 2]);
  a.dbt     = reshape(a.dbt,2645,64*72);
a.dbt_err = permute(x.dbt_err,[3 1 2]);
  a.dbt_err = reshape(a.dbt_err,2645,64*72);
a.lag     = permute(x.lag,[3 1 2]);
  a.lag     = reshape(a.lag,2645,64*72);

b_obs(:,:,2)     = a.dbt';
b_err_obs(:,:,2) = a.dbt_err';    
lagcor_obs_anom(:,:) = a.lag';
b_obs = real(b_obs);
b_err_obs = real(b_err_obs);

%{
b_cal(:,:,2)     = b.dbt';
b_err_cal(:,:,2) = b.dbt_err';
lagcor_cal_anom(:,:) = b.lag';
b_cal = real(b_cal);
b_err_cal = real(b_err_cal);

b_bias(:,:,2)     = c.dbt';
b_err_bias(:,:,2) = c.dbt_err';
lagcor_bias_anom(:,:) = c.lag';
b_bias = real(b_bias);
b_err_bias = real(b_err_bias);
%}

comment = 'from /home/strow/Work/Airs/Stability/Data/Desc_fits_iasi_times/Small/fit_* and convert_strowrates2oemrates_allskygrid_transpose.m';
save convert_strowrates2oemrates_allskygrid_obsonly_transpose  b_* comment lagcor*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); plot(fairs,squeeze(b_obs(:,:,2))); ax = axis; ax(3) = -0.25; ax(4) = +0.25; axis(ax);
figure(2); plot(fairs,lagcor_obs_anom); ax = axis; ax(3) = -1; ax(4) = +1; axis(ax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

