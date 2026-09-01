addpath /home/sergio/MATLABCODE

%fairs = instr_chans;    %% 2378

load sarta_chans_for_l1c.mat

hdffile = '/home/sergio/MATLABCODE/airs_l1c_srf_tables_lls_20181205.hdf';   % what he gave in Dec 2018
vchan2834 = hdfread(hdffile,'freq');
fairs = vchan2834;     %% 2834

fairs = fairs(ichan);  %% 2645  

iasi = load('/asl/s1/strow/iasi_desc_airs_srf.mat');
for ix = 1 : 40
  fprintf(1,'ix = %2i \n',ix);

  clear a
    a = iasi.o2a; 
    b_obs(ix,:,2)     = real(a.dbt(ix,:));
    b_err_obs(ix,:,2) = real(a.dbt_err(ix,:));    
    lagcor_obs_anom(ix,:) = real(a.lag(ix,:));
end

comment = 'from /asl/s1/strow/iasi_desc_airs_srf.mat and convert_strowrates2oemrates_iasi2airs_iasitimespan.m'
%% for June 04, 2019
save convert_strowrates2oemrates_clear_11_year_iasi2airs_obs.mat  b_* comment lagcor*

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
