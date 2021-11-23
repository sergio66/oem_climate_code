%% simple code to produce the list of L1c channels used in the AMT 2020 stability paper

load ../SAVE_LW_noCFC11_noN2O_covx10/OutputAnomaly_OBS/20/anomtest_timestep265.mat
chanset = jacobian.chanset;
load ../f2645.mat

chris = load('chanset_chris.mat');

[Y,iSergio,iChris] = intersect(chanset,chris.chanset);

comment = 'see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/GoodChansChirp/amt2020_climateQA.m';
save climateQA_LW_noCFC.mat chanset comment
