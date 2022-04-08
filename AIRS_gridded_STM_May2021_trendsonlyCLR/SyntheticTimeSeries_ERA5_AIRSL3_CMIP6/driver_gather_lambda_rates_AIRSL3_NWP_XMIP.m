addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/COLORMAP/LLS
addpath /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS

%% see driver_gather_spectralrates_AIRSL3_NWP_XMIP6.m

load('llsmap5');
if length(llsmap5) == 64
  %% need to center the white 1.0 1.0 1.0 .. right now it is at position 33, so need 65 points, or remove first ... choose that
  llsmap5 = llsmap5(2:64,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

junk = load('../../FIND_NWP_MODEL_TRENDS/MLS_atm_data_2004_09_to_2020_08_trends.mat');
rlat = junk.trend_rlat64;

set_gather_savename_rates_AIRSL3_NWP_XMIP

savename1 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/olr_feedbacks_AIRSL3_ERA5_CMIP6.mat';     %% this is base lambda for AIRSL3_ERA5_CMIP6
savename2 = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/olr_feedbacks_CLIMCAPS_MERRA2_AMIP6.mat'; %% this is base lambda for CLIMCAPS_MERRA2_AMIP6

%% see ~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/do_avg_feedback2cos.m
junk = load(savename, 'umbc_spectral_olr');    umbc_spectral_olr = junk.umbc_spectral_olr;
junk = load(savename1,'cmip6_spectral_olr');   cmip6_spectral_olr = junk.cmip6_spectral_olr;
junk = load(savename1,'era5_spectral_olr');    era5_spectral_olr = junk.era5_spectral_olr;
junk = load(savename1,'airsL3_spectral_olr');  airsL3_spectral_olr = junk.airsL3_spectral_olr;

junk = load(savename, 'umbc_spectral_olr');   umbc2_spectral_olr = junk.umbc_spectral_olr;
junk = load(savename2,'cmip6_spectral_olr');  amip6_spectral_olr = junk.cmip6_spectral_olr;
junk = load(savename2,'era5_spectral_olr');   merra2_spectral_olr = junk.era5_spectral_olr;
junk = load(savename2,'airsL3_spectral_olr'); climcaps_spectral_olr = junk.airsL3_spectral_olr;

junk = load(savename,'rlat'); rlat = junk.rlat;
cosrlat = cos(rlat'*pi/180);
nmcos = 1/nanmean(cosrlat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = 1;
  wonk11 = era5_spectral_olr.feedback.planck_ecRad; wonk11(wonk11 < -10) = NaN; wonk11(wonk11 > +00) = NaN; wonk11 = reshape(wonk11,72,64); wonk11 = nanmean(wonk11,1); plot(rlat,wonk11); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk11);
  wonk12 = era5_spectral_olr.feedback.lapse_ecRad;  wonk12(wonk12 < -05) = NaN; wonk12(wonk12 > +05) = NaN; wonk12 = reshape(wonk12,72,64); wonk12 = nanmean(wonk12,1); plot(rlat,wonk12); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk12);
  wonk13 = era5_spectral_olr.feedback.wv_ecRad;     wonk13(wonk13 < -05) = NaN; wonk13(wonk13 > +05) = NaN; wonk13 = reshape(wonk13,72,64); wonk13 = nanmean(wonk13,1); plot(rlat,wonk13); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk13);
  wonk14 = era5_spectral_olr.feedback.skt_ecRad;    wonk14(wonk14 < -05) = NaN; wonk14(wonk14 > +00) = NaN; wonk14 = reshape(wonk14,72,64); wonk14 = nanmean(wonk14,1); plot(rlat,wonk14); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk14);
model = 2;
  wonk21 = merra2_spectral_olr.feedback.planck_ecRad; wonk21(wonk21 < -10) = NaN; wonk21(wonk21 > +00) = NaN; wonk21 = reshape(wonk21,72,64); wonk21 = nanmean(wonk21,1); plot(rlat,wonk21); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk21);
  wonk22 = merra2_spectral_olr.feedback.lapse_ecRad;  wonk22(wonk22 < -05) = NaN; wonk22(wonk22 > +05) = NaN; wonk22 = reshape(wonk22,72,64); wonk22 = nanmean(wonk22,1); plot(rlat,wonk22); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk22);
  wonk23 = merra2_spectral_olr.feedback.wv_ecRad;     wonk23(wonk23 < -05) = NaN; wonk23(wonk23 > +05) = NaN; wonk23 = reshape(wonk23,72,64); wonk23 = nanmean(wonk23,1); plot(rlat,wonk23); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk23);
  wonk24 = merra2_spectral_olr.feedback.skt_ecRad;    wonk24(wonk24 < -05) = NaN; wonk24(wonk24 > +00) = NaN; wonk24 = reshape(wonk24,72,64); wonk24 = nanmean(wonk24,1); plot(rlat,wonk24); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk24);

model = 3;
  wonk31 = airsL3_spectral_olr.feedback.planck_ecRad; wonk31(wonk31 < -10) = NaN; wonk31(wonk31 > +00) = NaN; wonk31 = reshape(wonk31,72,64); wonk31 = nanmean(wonk31,1); plot(rlat,wonk31); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk31);
  wonk32 = airsL3_spectral_olr.feedback.lapse_ecRad;  wonk32(wonk32 < -05) = NaN; wonk32(wonk32 > +05) = NaN; wonk32 = reshape(wonk32,72,64); wonk32 = nanmean(wonk32,1); plot(rlat,wonk32); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk32);
  wonk33 = airsL3_spectral_olr.feedback.wv_ecRad;     wonk33(wonk33 < -05) = NaN; wonk33(wonk33 > +05) = NaN; wonk33 = reshape(wonk33,72,64); wonk33 = nanmean(wonk33,1); plot(rlat,wonk33); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk33);
  wonk34 = airsL3_spectral_olr.feedback.skt_ecRad;    wonk34(wonk34 < -05) = NaN; wonk34(wonk34 > +00) = NaN; wonk34 = reshape(wonk34,72,64); wonk34 = nanmean(wonk34,1); plot(rlat,wonk34); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk34);
model = 4;
  wonk41 = climcaps_spectral_olr.feedback.planck_ecRad; wonk41(wonk41 < -10) = NaN; wonk41(wonk41 > +00) = NaN; wonk41 = reshape(wonk41,72,64); wonk41 = nanmean(wonk41,1); plot(rlat,wonk41); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk41);
  wonk42 = climcaps_spectral_olr.feedback.lapse_ecRad;  wonk42(wonk42 < -05) = NaN; wonk42(wonk42 > +05) = NaN; wonk42 = reshape(wonk42,72,64); wonk42 = nanmean(wonk42,1); plot(rlat,wonk42); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk42);
  wonk43 = climcaps_spectral_olr.feedback.wv_ecRad;     wonk43(wonk43 < -05) = NaN; wonk43(wonk43 > +05) = NaN; wonk43 = reshape(wonk43,72,64); wonk43 = nanmean(wonk43,1); plot(rlat,wonk43); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk43);
  wonk44 = climcaps_spectral_olr.feedback.skt_ecRad;    wonk44(wonk44 < -05) = NaN; wonk44(wonk44 > +00) = NaN; wonk44 = reshape(wonk44,72,64); wonk44 = nanmean(wonk44,1); plot(rlat,wonk44); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk44);

model = 5;
  wonk51 = cmip6_spectral_olr.feedback.planck_ecRad; wonk51(wonk51 < -10) = NaN; wonk51(wonk51 > +00) = NaN; wonk51 = reshape(wonk51,72,64); wonk51 = nanmean(wonk51,1); plot(rlat,wonk51); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk51);
  wonk52 = cmip6_spectral_olr.feedback.lapse_ecRad;  wonk52(wonk52 < -05) = NaN; wonk52(wonk52 > +05) = NaN; wonk52 = reshape(wonk52,72,64); wonk52 = nanmean(wonk52,1); plot(rlat,wonk52); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk52);
  wonk53 = cmip6_spectral_olr.feedback.wv_ecRad;     wonk53(wonk53 < -05) = NaN; wonk53(wonk53 > +05) = NaN; wonk53 = reshape(wonk53,72,64); wonk53 = nanmean(wonk53,1); plot(rlat,wonk53); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk53);
  wonk54 = cmip6_spectral_olr.feedback.skt_ecRad;    wonk54(wonk54 < -05) = NaN; wonk54(wonk54 > +00) = NaN; wonk54 = reshape(wonk54,72,64); wonk54 = nanmean(wonk54,1); plot(rlat,wonk54); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk54);
model = 6;
  wonk61 = amip6_spectral_olr.feedback.planck_ecRad; wonk61(wonk61 < -10) = NaN; wonk61(wonk61 > +00) = NaN; wonk61 = reshape(wonk61,72,64); wonk61 = nanmean(wonk61,1); plot(rlat,wonk61); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk61);
  wonk62 = amip6_spectral_olr.feedback.lapse_ecRad;  wonk62(wonk62 < -05) = NaN; wonk62(wonk62 > +05) = NaN; wonk62 = reshape(wonk62,72,64); wonk62 = nanmean(wonk62,1); plot(rlat,wonk62); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk62);
  wonk63 = amip6_spectral_olr.feedback.wv_ecRad;     wonk63(wonk63 < -05) = NaN; wonk63(wonk63 > +05) = NaN; wonk63 = reshape(wonk63,72,64); wonk63 = nanmean(wonk63,1); plot(rlat,wonk63); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk63);
  wonk64 = amip6_spectral_olr.feedback.skt_ecRad;    wonk64(wonk64 < -05) = NaN; wonk64(wonk64 > +00) = NaN; wonk64 = reshape(wonk64,72,64); wonk64 = nanmean(wonk64,1); plot(rlat,wonk64); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk64);

model = 7;
  wonk71 = umbc_spectral_olr.feedback.planck_ecRad; wonk71(wonk71 < -10) = NaN; wonk71(wonk71 > +00) = NaN; wonk71 = reshape(wonk71,72,64); wonk71 = nanmean(wonk71,1); plot(rlat,wonk71); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk71);
  wonk72 = umbc_spectral_olr.feedback.lapse_ecRad;  wonk72(wonk72 < -05) = NaN; wonk72(wonk72 > +05) = NaN; wonk72 = reshape(wonk72,72,64); wonk72 = nanmean(wonk72,1); plot(rlat,wonk72); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk72);
  wonk73 = umbc_spectral_olr.feedback.wv_ecRad;     wonk73(wonk73 < -05) = NaN; wonk73(wonk73 > +05) = NaN; wonk73 = reshape(wonk73,72,64); wonk73 = nanmean(wonk73,1); plot(rlat,wonk73); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk73);
  wonk74 = umbc_spectral_olr.feedback.skt_ecRad;    wonk74(wonk74 < -05) = NaN; wonk74(wonk74 > +00) = NaN; wonk74 = reshape(wonk74,72,64); wonk74 = nanmean(wonk74,1); plot(rlat,wonk74); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk74);
%{
model = 8;
  wonk81 = mls_spectral_olr.feedback.planck_ecRad; wonk81(wonk81 < -10) = NaN; wonk81(wonk81 > +00) = NaN; wonk81 = reshape(wonk81,72,64); wonk81 = nanmean(wonk81,1); plot(rlat,wonk81); mean_feedback(model,1) = nmcos * nanmean(cosrlat.*wonk81);
  wonk82 = mls_spectral_olr.feedback.lapse_ecRad;  wonk82(wonk82 < -05) = NaN; wonk82(wonk82 > +05) = NaN; wonk82 = reshape(wonk82,72,64); wonk82 = nanmean(wonk82,1); plot(rlat,wonk82); mean_feedback(model,2) = nmcos * nanmean(cosrlat.*wonk82);
  wonk83 = mls_spectral_olr.feedback.wv_ecRad;     wonk83(wonk83 < -05) = NaN; wonk83(wonk83 > +05) = NaN; wonk83 = reshape(wonk83,72,64); wonk83 = nanmean(wonk83,1); plot(rlat,wonk83); mean_feedback(model,3) = nmcos * nanmean(cosrlat.*wonk83);
  wonk84 = mls_spectral_olr.feedback.skt_ecRad;    wonk84(wonk84 < -05) = NaN; wonk84(wonk84 > +00) = NaN; wonk84 = reshape(wonk84,72,64); wonk84 = nanmean(wonk84,1); plot(rlat,wonk84); mean_feedback(model,4) = nmcos * nanmean(cosrlat.*wonk84);
%}

for ii = 1 : 7
  for jj = 1 : 4
    str = ['junk = wonk' num2str(ii) num2str(jj) ';'];
    eval(str);
    junk = smooth(junk,4)';
    %junk = smooth(junk,8)';
    str = ['swonk' num2str(ii) num2str(jj) ' = junk;'];
    eval(str);
  end
end

disp('<<<<<<<<<<<<< ---- do_avg_feedback2cos.m ----- >>>>>>>>>>>>>>')
disp('           Planck          Lapse           WV            SKT')
fprintf(1,'ERA5     %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(1,:))
fprintf(1,'MERRA2   %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(2,:))
fprintf(1,'AIRS L3  %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(3,:))
fprintf(1,'CLIMCAPS %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(4,:))
fprintf(1,'CMIP6    %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(5,:))
fprintf(1,'AMIP6    %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(6,:))
fprintf(1,'UMBC     %9.2f      %9.2f     %9.2f    %9.2f \n',mean_feedback(7,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear plotoptions;
iFig = 21; figure(iFig); clf; plot(rlat,swonk11,'b',rlat,swonk21,'c',rlat,swonk31,'r',rlat,swonk41,'m',rlat,swonk51,'g',rlat,swonk61,'ys',rlat,swonk71,'k','linewidth',2); title('Planck Feedback W/m2/K');
 plotaxis2; hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS','CMIP6','AMIP6','UMBC','location','best','fontsize',8); xlabel('Latitude'); ylabel('\lambda');
iFig = 22; figure(iFig); clf; plot(rlat,swonk12,'b',rlat,swonk22,'c',rlat,swonk32,'r',rlat,swonk42,'m',rlat,swonk52,'g',rlat,swonk62,'ys',rlat,swonk72,'k','linewidth',2); title('Lapse Feedback W/m2/K');
 plotaxis2; hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS','CMIP6','AMIP6','UMBC','location','best','fontsize',8); xlabel('Latitude'); ylabel('\lambda');
iFig = 23; figure(iFig); clf; plot(rlat,swonk13,'b',rlat,swonk23,'c',rlat,swonk33,'r',rlat,swonk43,'m',rlat,swonk53,'g',rlat,swonk63,'ys',rlat,swonk73,'k','linewidth',2); title('WV Feedback W/m2/K');
 plotaxis2; hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS','CMIP6','AMIP6','UMBC','location','best','fontsize',8); xlabel('Latitude'); ylabel('\lambda');
iFig = 24; figure(iFig); clf; plot(rlat,swonk14,'b',rlat,swonk24,'c',rlat,swonk34,'r',rlat,swonk44,'m',rlat,swonk54,'g',rlat,swonk64,'ys',rlat,swonk74,'k','linewidth',2); title('SKT Feedback W/m2/K');
 plotaxis2; hl = legend('ERA5','MERRA2','AIRS L3','CLIMCAPS','CMIP6','AMIP6','UMBC','location','best','fontsize',8); xlabel('Latitude'); ylabel('\lambda');

if iSave > 0
  comment1 = 'variables are wonkXY, swonkXY where  wonk* = raw lambda, swonk* = smoothed lambda';
  comment2 = ' X=1--7   1=ERA5 2=MERRA2    3=AIRSL3 4=CLIMCAPS    5=CMIP6 6=CMIP6    7=UMBC'; 
  comment3 = ' Y=1--4   1=planck, 2=lapse  3=wv     4=skt';
  saver = ['save FIGS/Figs_JPL_Apr2022/strow_jpl_Apr2022_lambda' savestr '.mat wonk* swonk* rlat comment*'];
  eval(saver)
end

