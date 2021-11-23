addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/COLORMAP
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

latbins = equal_area_spherical_bands(20);

iShowStats = +1; %% do     run plot_wv_co2_o3_co_ch4_bias_qa
iShowStats = -1; %% do not run plot_wv_co2_o3_co_ch4_bias_qa

for yy = 2002 : 2017
  fprintf(1,'year = %4i \n',yy)

  clear poemNew poem btobs btcal*

  if mod(yy,4) ~= 0
    fname = ['/asl/s1/sergio/rtp/singlefootprintretrievals_airicrad_v6/' num2str(yy) '/12/12/era_airicrad_day_clearonly_346_4_iDET_4_iStemp_ColWV_5_iCenterFov_-1_retrclr_no_clouds.mat'];
  elseif mod(yy,4) == 0
    fname = ['/asl/s1/sergio/rtp/singlefootprintretrievals_airicrad_v6/' num2str(yy) '/12/12/era_airicrad_day_clearonly_347_4_iDET_4_iStemp_ColWV_5_iCenterFov_-1_retrclr_no_clouds.mat'];
  end

  loader = ['load ' fname];
  eval(loader)

  if iShowStats > 0
    plot_wv_co2_o3_co_ch4_bias_qa
    disp('ret to continue'); pause
  else
    if ~exist('poem')
      disp('warning .. no poem so just setting poem to be orig stuff');
      poem = poemNew;
      poem.gas_1 = poemNew.gas_1_orig;
      poem.gas_2 = poemNew.gas_2_orig;
      poem.gas_3 = poemNew.gas_3_orig;
      %  poem.gas_4 = poemNew.gas_4_orig;
      poem.gas_5 = poemNew.gas_5_orig;
      poem.gas_6 = poemNew.gas_6_orig;
      poem.ptemp = poemNew.ptemp_orig;
      poem.stemp = poemNew.stemp_orig;
    end
  end

  [RH0,RH1km0,colwater0] = layeramt2RH(hoemNew,poem);
  [RHF,RH1kmF,colwaterF] = layeramt2RH(hoemNew,poemNew);

  RH0(101,:) = 0;
  RHF(101,:) = 0;
  do_averages_before_after_retrieval

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(pafter.diff_mean_gas_1,1:101,'c',nanmean(pafter.diff_mean_gas_1'),1:101,'b'); set(gca,'ydir','reverse'); axis([0.6 1.4 1 97])
title('mean(gas1after/gas1before)')

figure(2)
plot(pafter.diff_std_gas_1,1:101,'c',nanmean(pafter.diff_std_gas_1'),1:101,'b'); set(gca,'ydir','reverse'); axis([0 0.6 1 97])
title('std(gas1after/gas1before)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(pafter.diff_mean_ptemp,1:101,'c',nanmean(pafter.diff_mean_ptemp'),1:101,'b'); set(gca,'ydir','reverse'); axis([-0.1 +0.1 1 97])
title('mean(Tafter-Tbefore)')

figure(2)
plot(pafter.diff_std_ptemp,1:101,'c',nanmean(pafter.diff_std_ptemp'),1:101,'b'); set(gca,'ydir','reverse'); axis([0 0.1 1 97])
title('std(Tafter-Tbefore)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

show_averages_differences_before_after_retrieval

%% fill_in_the_fields

hx = hoemNew;
hx.pfields = 1;  %% only has profile
hx.ptype = 1;

%{
save_the_rtp
%}
